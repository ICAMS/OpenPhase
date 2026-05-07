/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2018
 *  Main contributors :   Matthias Stratmann; Oleg Shchyglo
 *
 */

#include "ThermodynamicInterfaces/ThermokineticInterface.h"
#include "Diffusion/DiffusionProperties.h"
#include "Settings.h"
#include "Composition.h"
#include "Temperature.h"
#include "PhysicalConstants.h"

namespace openphase
{

using namespace std;

void ThermokineticInterface::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ThermokineticInterface";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;

    PhaseNames = locSettings.PhaseNames;
    ElementNames = locSettings.ElementNames;

    DiffusionCoefficients.Allocate({Nphases,Ncomp,Ncomp});
    ActivationEnergies.Allocate({Nphases,Ncomp,Ncomp});
    
    BaseDir = locSettings.BaseDir;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ThermokineticInterface::ReadInput(const string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    ReadInput(data);

    ConsoleOutput::WriteLine();
}

void ThermokineticInterface::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    DBfile = FileInterface::ReadParameterFW(inp, moduleLocation, "DBfile", true, "");

	std::filesystem::path DB(DBfile);
	
	if (DB.is_relative())
        DB = BaseDir / DB;

    DB = std::filesystem::weakly_canonical(DB);

	DBfile = DB.string();
    ReadTDB();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void ThermokineticInterface::ReadTDB()
{
    fstream inp(DBfile.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + DBfile + " could not be opened", thisclassname, "ReadTDB()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();

    ReadTDB(data);

    inp.close();
}

void ThermokineticInterface::ReadTDB(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Thermo-kinetic Properties");
    ConsoleOutput::WriteStandard("Source", DBfile);

    // Read components available in the database
    int idx = 0;
    bool continue_reading = true;
    std::vector<std::string> DB_elements;
    while(continue_reading)
    {
        string comp = string("Comp_") + to_string(idx);

        if(FileInterface::FindParameter(inp, moduleLocation, comp) != -1)
        {
            string tmp_element = FileInterface::ReadParameterK(inp, moduleLocation, comp, true, "NN");
            DB_elements.push_back(tmp_element);
            idx++;
        }
        else
        {
            continue_reading = false;
        }
    }

    // Compare available components with system settings (assuming alphabetic order)
    if (DB_elements != ElementNames)
    {
        ConsoleOutput::WriteExit("Elements in the database file " + DBfile + " do not correspond to system settings", thisclassname, "ReadTDB()");
        OP_Exit(EXIT_FAILURE);
    }

    // Read phases available in the database
    idx = 0;
    continue_reading = true;
    std::vector<std::string> DB_phases;
    while(continue_reading)
    {
        string phase = string("Phase_") + to_string(idx);

        if(FileInterface::FindParameter(inp, moduleLocation, phase) != -1)
        {
            string tmp_phase = FileInterface::ReadParameterK(inp, moduleLocation, phase, true, "NN");
            DB_phases.push_back(tmp_phase);
            idx++;
        }
        else
        {
            continue_reading = false;
        }
    }

    // Compare available phases with system settings
    size_t success_counter = 0;
    for(size_t n = 0; n < PhaseNames.size(); n++)
    for(size_t m = 0; m < DB_phases.size(); m++)
    if (PhaseNames[n] == DB_phases[m])
    {
        success_counter++;
    }
    if(success_counter != PhaseNames.size())
    {
        ConsoleOutput::WriteExit("Not all required phases are present in the database file " + DBfile, thisclassname, "ReadTDB()");
        OP_Exit(EXIT_FAILURE);
    }

    //Reading diffusion coefficients and activation energies for phase and component
    for(size_t n = 0; n < Nphases; n++)
    for(size_t i = 0; i < Ncomp; i++)
    for(size_t j = 0; j < Ncomp; j++)
    {
        string counter = to_string(n) + "_" + ElementNames[i] + "_" + ElementNames[j];
        DiffusionCoefficients({n,i,j}) = FileInterface::ReadParameterD(inp, moduleLocation, string("DC_") + counter, false, 0.0);
        ActivationEnergies({n, i, j})  = FileInterface::ReadParameterD(inp, moduleLocation, string("AE_") + counter, false, 0.0);
    }

    ConsoleOutput::WriteLine();
}

void ThermokineticInterface::SetDiffusionData(PhaseField& Phase,
                                              Composition& Cx,
                                              Temperature& Tx,
                                              DiffusionProperties& DP)
{
    /** This function will populate diffusion properties based on the given
    parameters for temperature and composition for all active phases. */
    const double R = PhysicalConstants::R;
    double maxDC = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,DP.DiffusionCoefficients,DP.DiffusionCoefficients.Bcells(),reduction(max:maxDC))
    {
        for(size_t m = 0; m < Nphases; m++)
        if(Phase.Fractions(x,y,z,{m}) != 0.0)
        {
            for(size_t p = 0; p < Ncomp; p++)
            for(size_t q = 0; q < Ncomp; q++)
            if(!Cx.Phase[m].Component[p].isReference and
               !Cx.Phase[m].Component[p].isFastDiffusor and
               !Cx.Phase[m].Component[q].isReference and
               !Cx.Phase[m].Component[q].isFastDiffusor)
            {
                double locDC = DP.ScalingFactor({m,p,q})*DiffusionCoefficients({m,p,q})*exp(-ActivationEnergies({m,p,q})/(R*Tx(x,y,z)));
                DP.DiffusionCoefficients(x,y,z,{m,p,q}) = locDC;
                maxDC = std::max(maxDC,locDC);
            }
        }

        if(Phase.Fields(x,y,z).interface())
        {
            for(size_t m = 0; m < Nphases; m++)
            for(size_t n = 0; n < Nphases; n++)
            if(Phase.Fractions(x,y,z,{m}) != 0.0 and
               Phase.Fractions(x,y,z,{n}) != 0.0)
            {
                for(size_t p = 0; p < Ncomp; p++)
                for(size_t q = 0; q < Ncomp; q++)
                if(!Cx.Phase[m].Component[p].isReference and
                   !Cx.Phase[m].Component[p].isFastDiffusor and
                   !Cx.Phase[m].Component[q].isReference and
                   !Cx.Phase[m].Component[q].isFastDiffusor)
                {
                    double locDC = DP.DiffusionCoefficients(x,y,z,{m,p,q}) * DP.GrainBoundaryFactor({m,n});
                    DP.DiffusionCoefficients(x,y,z,{m,p,q}) = locDC;
                    maxDC = std::max(maxDC,locDC);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &maxDC, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
    DP.MaxDiffusionCoefficient = maxDC;
}

}//namespace openphase

