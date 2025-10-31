/*
 *   This file is part of the OpenPhase (R) software library.
 *  
 *  Copyright (c) 2009-2025 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2025 OpenPhase Solutions GmbH,
 *                Universitaetsstrasse 136, D-44799 Bochum, Germany.
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *     
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo
 *
 */

#include "UserDrivingForce.h"
#include "Settings.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Temperature.h"
#include "Composition.h"

namespace openphase
{
using namespace std;

UserDrivingForce::UserDrivingForce(Settings& locSettings,
                                               const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void UserDrivingForce::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "UserDrivingForce";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;

    Mode.Allocate(Nphases,Nphases);
    Value.Allocate(Nphases,Nphases);
    Teq.Allocate(Nphases,Nphases);
    Slope.Allocate(Nphases,Nphases);
    Component.Allocate(Nphases,Nphases);
    LatentHeat.Allocate(Nphases,Nphases);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void UserDrivingForce::ReadInput(std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("UserDrivingForce input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void UserDrivingForce::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    // Reading driving force modes and parameters for each phase pair
    for (size_t n = 0; n < Nphases; n++)
    for (size_t m = n; m < Nphases; m++)
    {
        stringstream idx;
        idx << "_" << n << "_" << m;
        if(n == m)
        {
            Mode(n, m) = UserDrivingForceModes::None;
        }
        else
        {
            string locMode = FileInterface::ReadParameterK(inp, moduleLocation, "UDF_Mode" + idx.str(), false, "NONE");

            if(locMode == "NONE")    Mode(n, m) = Mode(m, n) = UserDrivingForceModes::None;
            if(locMode == "VALUE")   Mode(n, m) = Mode(m, n) = UserDrivingForceModes::Value;
            if(locMode == "FORMULA") Mode(n, m) = Mode(m, n) = UserDrivingForceModes::Formula;
        }

        switch(Mode(n, m))
        {
            case UserDrivingForceModes::Value:
            {
                Value(n,m) = FileInterface::ReadParameterD(inp, moduleLocation, "UDF_Value" + idx.str(), true, 0.0);
                Value(m,n) = -Value(n,m);

                LatentHeat(n,m) = LatentHeat(m,n) = 0.0;
                Teq(n,m)        = Teq(m,n)        = 1.0;                        /// Set to 1 to avoid accidental division by zero
                Slope(n,m)      = Slope(m,n)      = 0.0;

                break;
            }
            case UserDrivingForceModes::Formula:
            {
                Value(n,m) = Value(m,n) = 0.0;

                LatentHeat(n,m) = FileInterface::ReadParameterD(inp, moduleLocation, "UDF_LatentHeat" + idx.str(), true, 0.0);
                LatentHeat(m,n) = -LatentHeat(n,m);

                Teq(n,m) = FileInterface::ReadParameterD(inp, moduleLocation, "UDF_Teq" + idx.str(), true, 1.0);
                Teq(m,n) = Teq(n,m);

                Slope(n,m) = FileInterface::ReadParameterD(inp, moduleLocation, "UDF_Slope" + idx.str(), false, 0.0);
                Slope(m,n) = Slope(n,m);

                Component(n,m) = FileInterface::ReadParameterI(inp, moduleLocation, "UDF_Component" + idx.str(), false, 0);
                Component(m,n) = Component(n,m);

                break;
            }
            case UserDrivingForceModes::None:
            default:
            {
                Value(n,m)       = Value(m,n)      = 0.0;
                LatentHeat(n,m)  = LatentHeat(m,n) = 0.0;
                Teq(n,m)         = Teq(m,n)        = 1.0;                       /// Set to 1 to avoid accidental division by zero
                break;
            }
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase,
                                       DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                int pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                if(Mode(pIndexA,pIndexB) == UserDrivingForceModes::Value)
                {
                    dGab.Force(i,j,k).add_raw(alpha->index, beta->index, Value(pIndexA,pIndexB));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase,
                                       DrivingForce& dGab,
                                       Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                switch(Mode(pIndexA,pIndexB))
                {
                    case UserDrivingForceModes::None:
                    {
                        break;
                    }
                    case UserDrivingForceModes::Formula:
                    {
                        double dG_AB = LatentHeat(pIndexA,pIndexB)*
                                       (Tx(i,j,k) - Teq(pIndexA,pIndexB))/
                                        Teq(pIndexA,pIndexB);

                        dGab.Force(i,j,k).add_raw(alpha->index, beta->index, dG_AB);
                        break;
                    }
                    case UserDrivingForceModes::Value:
                    {
                        dGab.Force(i,j,k).add_raw(alpha->index, beta->index, Value(pIndexA,pIndexB));
                        break;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase,
                                       DrivingForce& dGab,
                                       Temperature& Tx,
                                       Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                switch(Mode(pIndexA,pIndexB))
                {
                    case UserDrivingForceModes::None:
                    {
                        break;
                    }
                    case UserDrivingForceModes::Formula:
                    {
                        double locTeq = Teq(pIndexA,pIndexB)
                                      + Slope(pIndexA,pIndexB)*Cx.MoleFractions(i,j,k,{pIndexA,Component(pIndexA,pIndexB)});
                        double dG_AB = LatentHeat(pIndexA,pIndexB)*
                                       (Tx(i,j,k) - locTeq)/locTeq;

                        dGab.Force(i,j,k).add_raw(alpha->index, beta->index, dG_AB);
                        break;
                    }
                    case UserDrivingForceModes::Value:
                    {
                        dGab.Force(i,j,k).add_raw(alpha->index, beta->index, Value(pIndexA,pIndexB));
                        break;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetDrivingForce(PhaseField& Phase, DrivingForce& dGab,
                                       int indexA, int indexB, double dGvalue)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            if (Phase.Fields(i,j,k).present(indexA) and
                Phase.Fields(i,j,k).present(indexB))
            {
                dGab.Force(i,j,k).add_raw(indexA, indexB, dGvalue);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

}// namespace openphase
