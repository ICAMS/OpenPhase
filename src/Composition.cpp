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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Matthias Stratmann
 *
 */

#include "Composition.h"
#include "ElasticProperties.h"
#include "Settings.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include "Velocities.h"
#include "AdvectionHR.h"
#include "Thermodynamics/PeriodicTable.h"
#include "H5Interface.h"

namespace openphase
{
using namespace std;
/*************************************************************************/

void Composition::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Composition";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Ncomp   = locSettings.Ncomp;
    Nphases = locSettings.Nphases;

    for(size_t i = 0; i < Ncomp; i++)
    {
        Element tempComp;
        tempComp.Index = i;
        tempComp.Value = 0.0;
        tempComp.Min = 0.0;
        tempComp.Max = 1.0;
        tempComp.Name = locSettings.ElementNames[i];
        tempComp.isVacancy = false;
        tempComp.isStoichiometric = false;
        tempComp.isReference = false;
        Component.push_back(tempComp);
    }
    Phase.resize(Nphases);
    for(size_t n = 0; n < Nphases; n++)
    {
        Phase[n].Name = locSettings.PhaseNames[n];
        Phase[n].Index = n;
        Phase[n].Number = n;
        Phase[n].Ncomp = Ncomp;
        Phase[n].State = locSettings.PhaseAggregateStates[n];
    }

    size_t Bcells = Grid.Bcells;

    MoleFractions.Allocate(Grid, {Nphases, Ncomp}, Bcells);

    NormIn.Allocate(Grid, {Nphases, Ncomp}, Bcells);
    NormOut.Allocate(Grid, {Nphases, Ncomp}, Bcells);

    MoleFractionsDotIn.Allocate(Grid, {Nphases, Ncomp}, Bcells);
    MoleFractionsDotOut.Allocate(Grid, {Nphases, Ncomp}, Bcells);

    NormTotal.Allocate(Grid, {Ncomp}, Bcells);

    MoleFractionsTotal.Allocate(Grid, {Ncomp}, Bcells);
    MoleFractionsTotalDot.Allocate(Grid, {Ncomp}, Bcells);

    MassFractionsTotal.Allocate(Grid, {Ncomp}, Bcells);
    MassFractionsTotalOld.Allocate(Grid, {Ncomp}, Bcells);
    MolecularWeight.Allocate({Ncomp});

    Limiting.Allocate(Grid, Bcells);

    Initial.Allocate({Nphases, Ncomp});

    MoleFractionsAverage.Allocate({Nphases, Ncomp});
    MoleFractionsTotalAverage.Allocate({Ncomp});
    MoleFractionsInterfaceAverage.Allocate({Nphases, Nphases, Ncomp});
    Interface.Allocate({Nphases, Nphases});

    TotInitial.Allocate({Ncomp});

    PT.Initialize();

    AtStart = true;

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Composition::ReadInput(const string InputFileName)
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

void Composition::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    //Read reference elements
    std::vector<std::string> RefElements;
    if(Component.size() != 0)
    for (size_t n = 0; n < Nphases; n++)
    {
        stringstream converter;
        converter << string("RefElement_") << n;
        if (FileInterface::FindParameter(inp, moduleLocation, converter.str()))
        {
            string tmp = FileInterface::ReadParameterK(inp, moduleLocation, converter.str(),false,Component[0].Name);
            tmp.erase(remove(tmp.begin(), tmp.end(), ' '), tmp.end());
            RefElements.push_back(tmp);
        }
    }
    // Read initial compositions
    for(size_t n = 0; n < Ncomp; n++)
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream converter;
        converter << alpha << string("_") << Component[n].Name;
        string counter = converter.str();

        Initial({alpha, n}) = FileInterface::ReadParameterD(inp, moduleLocation, string("C0_") + converter.str());
    }

    // Set reference element flags
    for(size_t i = 0; i < Nphases; i++)
    {
        Phase[i].Index = i;
        Phase[i].Component = Component;

        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            if(Phase[i].Component[comp].Name == RefElements[i])
            {
                Phase[i].Component[comp].isReference = true;
                Phase[i].ReferenceElement = comp;
            }
        }
    }

    // Set fast diffusors
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t  comp = 0;  comp < Ncomp;   comp++)
    {
        string identifier = "FastDiffusor_" + to_string(alpha) + "_"
                                            + Phase[alpha].Component[comp].Name;

        Phase[alpha].Component[comp].isFastDiffusor = FileInterface::ReadParameterI(
                                                    inp,moduleLocation,identifier,
                                                    false,false);
    }

    for(size_t n = 0; n < Ncomp; n++)
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream converter;
        converter << alpha << string("_") << Component[n].Name;
        string counter = converter.str();

        Phase[alpha].Component[n].isStoichiometric = FileInterface::ReadParameterB(inp, moduleLocation, string("Stoichiometric_") + converter.str(), false, false);
    }
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        size_t nStoich = 0;
        for(size_t n = 0; n < Ncomp; n++)
        {
            if (Phase[alpha].Component[n].isStoichiometric) nStoich++;
        }
        if (nStoich == Ncomp-1)
        {
            for(size_t n = 0; n < Ncomp; n++)
            {
                Phase[alpha].Component[n].isStoichiometric = true;
            }
        }
    }

    // Read composition limits
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        string identifier = "Cmin_" + to_string(alpha) + "_"
                                    + Phase[alpha].Component[comp].Name;

        Phase[alpha].Component[comp].Min = FileInterface::ReadParameterD(inp,
                                           moduleLocation,identifier,false,
                                           Phase[alpha].Component[comp].Min);

        identifier = "Cmax_" + to_string(alpha) + "_"
                             + Phase[alpha].Component[comp].Name;

        Phase[alpha].Component[comp].Max = FileInterface::ReadParameterD(inp,
                                           moduleLocation,identifier,false,
                                           Phase[alpha].Component[comp].Max);
    }

    // Detect stoichiometric components
    /*for (size_t alpha = 0; alpha < Nphases; alpha++)
    {
        Phase[alpha].StoichiometricFlag = false;

        size_t nStoich = 0;
        for (size_t comp = 0; comp < Ncomp; comp++)
        if (fabs(Phase[alpha].Component[comp].Max
            - Phase[alpha].Component[comp].Min) <= DBL_EPSILON)
        {
            Phase[alpha].Component[comp].isStoichiometric = true;
            Phase[alpha].StoichiometricFlag = true;
            nStoich++;

            string tmp = "Component " + Component[comp].Name
                + " in phase " + Phase[alpha].Name
                + " will be regarded as a stoichiometric component!";
            cout << tmp << endl;
        }

        if (nStoich == Ncomp - 1)
        for (size_t comp = 0; comp < Ncomp; comp++)
        if (!Phase[alpha].Component[comp].isStoichiometric)
        {
            Phase[alpha].Component[comp].isStoichiometric = true;

            string tmp = "Component " + Component[comp].Name
                + " in phase " + Phase[alpha].Name
                + " has no composition limits but seems to be a"
                " stoichiometric component!";
            ConsoleOutput::WriteWarning(tmp, thisclassname, "Initialize");
        }
    }*/

    // Enable 1D extensions
    if(Grid.dNx)
    {
        int X0_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_X0", false, 0);
        if(X0_size > 0)
        {
            int X0_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_X0_Phase", true, 0);
            ExtensionX0.Initialize(X0_size,X0_phase, Ncomp,(iVector3){-1, 0, 0});
            ExtensionsActive = true;
        };

        int XN_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_XN", false, 0);
        if(XN_size > 0)
        {
            int XN_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_XN_Phase", true, 0);
            ExtensionXN.Initialize(XN_size,XN_phase,Ncomp,(iVector3){ 1, 0, 0});
            ExtensionsActive = true;
        };
    }
    if(Grid.dNy)
    {
        int Y0_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Y0", false, 0);
        if(Y0_size > 0)
        {
            int Y0_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Y0_Phase", true, 0);
            ExtensionX0.Initialize(Y0_size,Y0_phase, Ncomp,(iVector3){0, -1, 0});
            ExtensionsActive = true;
        };

        int YN_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_YN", false, 0);
        if(YN_size > 0)
        {
            int YN_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_YN_Phase", true, 0);
            ExtensionYN.Initialize(YN_size,YN_phase,Ncomp,(iVector3){0, 1, 0});
            ExtensionsActive = true;
        };
    }
    if(Grid.dNz)
    {
        int Z0_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Z0", false, 0);
        if(Z0_size > 0)
        {
            int Z0_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Z0_Phase", true, 0);
            ExtensionX0.Initialize(Z0_size,Z0_phase, Ncomp,(iVector3){0, 0, -1});
            ExtensionsActive = true;
        };

        int ZN_size  = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_ZN", false, 0);
        if(ZN_size > 0)
        {
            int ZN_phase = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_ZN_Phase", true, 0);
            ExtensionXN.Initialize(ZN_size,ZN_phase,Ncomp,(iVector3){0, 0, 1});
            ExtensionsActive = true;
        };
    }

    SortedElementMatrix = CalculateSortedElementMatrix(); 

    PrintData();

    //CalculateMolefractionLimits();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Composition::SetInitiaMoleFractions1Dextension(Composition1Dextension& CxExt)
{
    for(size_t x = 0; x < CxExt.Data.total_size(); x++)
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        CxExt.Data[x]({comp}) = Initial({CxExt.PhaseIndex,comp});
    }
}
void Composition::SetInitialMoleFractions(PhaseField& Phi)
{
    /* Simple composition set up, uniform component fractions per phase. */

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal, MoleFractionsTotal.Bcells(),)
    {
        MoleFractions(i,j,k) = Initial;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    CalculateTotalMoleFractions(Phi);
    CalculateMoleFractionsTotalAverage();
    CalculateMoleFractionsAverage(Phi);
    CalculateTotalMolarVolume();

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    if(ExtensionX0.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionX0);}
#ifdef MPI_PARALLEL
    if(MPI_RANK == MPI_SIZE - 1)
#endif
    if(ExtensionXN.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionXN);}

    if(ExtensionY0.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionY0);}
    if(ExtensionYN.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionYN);}
    if(ExtensionZ0.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionZ0);}
    if(ExtensionZN.isActive()) {SetInitiaMoleFractions1Dextension(ExtensionZN);}

    ConsoleOutput::WriteLine("Total molar volume: " + to_string(TotalMolarVolume));
}

Tensor<double, 1> Composition::WeightFractionsTotal(int x, int y, int z) const
{
    /**This function calculates weight fractions of all elements in a given
     * grid cell.*/

    Tensor<double, 1> result({Ncomp});

    double div = 0.0;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        double locWeight = MoleFractionsTotal(x,y,z,{comp})*PT.GetData(Component[comp].Name).AtomicWeight;
        div += locWeight;
        result({comp}) = locWeight;
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        result({comp}) /= div;
    }
    return result;
}

void Composition::CalculateTotalMoleFractions(PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,MoleFractionsTotal.Bcells(),)
    {
        MoleFractionsTotal(i,j,k).set_to_zero();

        for(size_t comp = 0; comp < Ncomp; comp++)
        for(size_t alpha = 0; alpha < Nphases; alpha++)
        {
            MoleFractionsTotal(i,j,k,{comp}) += Phi.Fractions(i,j,k,{alpha})*
                                             MoleFractions(i,j,k,{alpha,comp});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Composition::CalculateTotalMolarVolume(void)
{
    TotalMolarVolume = 0.0;

    for(size_t i = 0; i < Ncomp; i++)
    {
        TotalMolarVolume += PT.GetData(Component[i].Name).MolarVolume
                            *MoleFractionsTotalAverage[i];
    }
}

void Composition::CalculateMoleFractionsTotalAverage(void)
{
    Tensor<double, 1> locMoleFractionsTotalAverage({Ncomp});
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0,reduction(TensorD1Sum: locMoleFractionsTotalAverage))
    {
        locMoleFractionsTotalAverage += MoleFractionsTotal(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    MoleFractionsTotalAverage = locMoleFractionsTotalAverage;

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(locMoleFractionsTotalAverage.data(), MoleFractionsTotalAverage.data(), Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif
    MoleFractionsTotalAverage /= double(Grid.TotalNumberOfCells());
}

void Composition::CalculateMoleFractionsAverage(PhaseField& Phase)
{
    Tensor<double, 3> locMoleFractionsInterfaceAverage({Nphases,Nphases,Ncomp});
    Tensor<double, 2> locMoleFractionsAverage({Nphases,Ncomp});

    Tensor<double, 2> NpointsAB({Nphases,Nphases});
    Tensor<double, 1> NpointsA({Nphases});

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,MoleFractions,0,reduction(TensorD2Sum: locMoleFractionsAverage) reduction(TensorD3Sum: locMoleFractionsInterfaceAverage) reduction(TensorD1Sum: NpointsA) reduction(TensorD2Sum: NpointsAB))
    {
        for(size_t alpha = 0; alpha < Nphases; alpha++)
        if(Phase.Fractions(i,j,k,{alpha}) != 0.0)
        {
            NpointsA({alpha}) += 1;
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                locMoleFractionsAverage({alpha, comp}) += MoleFractions(i,j,k,{alpha, comp});
            }
        }

        if(Phase.Fields(i,j,k).interface())
        {
            for(size_t alpha = 0; alpha < Nphases; alpha++)
            if(Phase.Fractions(i,j,k,{alpha}) != 0.0)
            for(size_t beta = 0; beta < Nphases; beta++)
            if(Phase.Fractions(i,j,k,{beta}) != 0.0)
            {
                NpointsAB({alpha, beta}) += 1;
                for(size_t comp = 0; comp < Ncomp; comp++)
                {
                    locMoleFractionsInterfaceAverage({alpha, beta, comp}) +=
                            Phase.Fractions(i,j,k,{alpha})*MoleFractions(i,j,k,{alpha, comp}) +
                            Phase.Fractions(i,j,k,{ beta})*MoleFractions(i,j,k,{ beta, comp});
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    MoleFractionsAverage = locMoleFractionsAverage;
    MoleFractionsInterfaceAverage = locMoleFractionsInterfaceAverage;

#ifdef MPI_PARALLEL
    Tensor<double, 1> tmpNpointsA = NpointsA;
    Tensor<double, 2> tmpNpointsAB = NpointsAB;
    OP_MPI_Allreduce(tmpNpointsA.data(), NpointsA.data(), Nphases, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(tmpNpointsAB.data(), NpointsAB.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);

    OP_MPI_Allreduce(locMoleFractionsAverage.data(), MoleFractionsAverage.data(), Nphases*Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(locMoleFractionsInterfaceAverage.data(), MoleFractionsInterfaceAverage.data(), Nphases*Nphases*Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    for(size_t alpha = 0; alpha < Nphases; alpha++)
    {
        if(NpointsA({alpha}) != 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsAverage({alpha, comp}) /= NpointsA({alpha});
            }
        }
        else
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsAverage({alpha, comp}) = MoleFractionsTotalAverage({comp});
            }
        }
    }

    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t  beta = 0;  beta < Nphases; beta++)
    {
        if(NpointsAB({alpha, beta}) != 0)
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsInterfaceAverage({alpha, beta, comp}) /= NpointsAB({alpha, beta});
            }
        }
        else
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                MoleFractionsInterfaceAverage({alpha, beta, comp}) = 0.5*(MoleFractionsAverage({alpha, comp}) + MoleFractionsAverage({beta, comp}));
            }
        }
    }
    
    for(size_t alpha = 0; alpha < Nphases; alpha++)
    for(size_t  beta = 0;  beta < Nphases; beta++)
    {
        if(NpointsAB({alpha, beta}) != 0)
        {
            Interface({alpha, beta}) = true;
        }
        else
        {
            Interface({alpha, beta}) = false;
        }
    }
}

void Composition::Remesh(const int newNx, const int newNy, const int newNz, const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy,newNz);

    MoleFractions.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    MoleFractionsTotal.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    MoleFractionsDotIn.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    MoleFractionsDotOut.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    MoleFractionsTotalDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    Limiting.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    NormIn.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    NormOut.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void Composition::MoveFrame(const int dx, const int dy, const int dz,
                            const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    int xEnd = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    int yEnd = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    int zEnd = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        MoleFractions(i,j,k) = MoleFractions(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
    }
    SetBoundaryConditions(BC);

    if(ExtensionX0.isActive()) ExtensionX0.moveFrame(-dx*Grid.dNx, BC.BC0X);
    if(ExtensionXN.isActive()) ExtensionXN.moveFrame( dx*Grid.dNx, BC.BCNX);
    if(ExtensionY0.isActive()) ExtensionY0.moveFrame(-dy*Grid.dNy, BC.BC0Y);
    if(ExtensionYN.isActive()) ExtensionYN.moveFrame( dy*Grid.dNy, BC.BCNY);
    if(ExtensionZ0.isActive()) ExtensionZ0.moveFrame(-dz*Grid.dNz, BC.BC0Z);
    if(ExtensionZN.isActive()) ExtensionZN.moveFrame( dz*Grid.dNz, BC.BCNZ);

    ConsoleOutput::WriteStandard(thisclassname, "Frame moved");
}

void Composition::ConsumePlane(const int dx, const int dy, const int dz,
                               const int  x, const int  y, const int  z,
                               const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    int xEnd = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    int yEnd = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    int zEnd = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz >= 0)
    {
        MoleFractions(i,j,k) = MoleFractions(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
    }
    xBeg = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    xEnd = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    xInc = 2*(dx < 0) - 1;

    yBeg = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    yEnd = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    yInc = 2*(dy < 0) - 1;

    zBeg = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    zEnd = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    zInc = 2*(dz < 0) - 1;

    for(int i = xBeg; ((dx >= 0) and (i >= xEnd)) or ((dx < 0) and (i <= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j >= yEnd)) or ((dy < 0) and (j <= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k >= zEnd)) or ((dz < 0) and (k <= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz < 0)
    {
        MoleFractions(i,j,k) = MoleFractions(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
        MoleFractionsTotal(i,j,k) = MoleFractionsTotal(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
    }
    SetBoundaryConditions(BC);

    if(ExtensionX0.isActive()) ExtensionX0.moveFrame( dx*Grid.dNx, BC.BC0X);
    if(ExtensionXN.isActive()) ExtensionXN.moveFrame( dx*Grid.dNx, BC.BCNX);
    if(ExtensionY0.isActive()) ExtensionY0.moveFrame( dy*Grid.dNy, BC.BC0Y);
    if(ExtensionYN.isActive()) ExtensionYN.moveFrame( dy*Grid.dNy, BC.BCNY);
    if(ExtensionZ0.isActive()) ExtensionZ0.moveFrame( dz*Grid.dNz, BC.BC0Z);
    if(ExtensionZN.isActive()) ExtensionZN.moveFrame( dz*Grid.dNz, BC.BCNZ);

    ConsoleOutput::WriteStandard(thisclassname, "Plane consumed");
}

void Composition::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(MoleFractions);
    BC.SetY(MoleFractions);
    BC.SetZ(MoleFractions);

    BC.SetX(MoleFractionsTotal);
    BC.SetY(MoleFractionsTotal);
    BC.SetZ(MoleFractionsTotal);
}

void Composition::SetLimitsBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(NormIn);
    BC.SetY(NormIn);
    BC.SetZ(NormIn);

    BC.SetX(NormOut);
    BC.SetY(NormOut);
    BC.SetZ(NormOut);

    BC.SetXFlags(Limiting);
    BC.SetYFlags(Limiting);
    BC.SetZFlags(Limiting);
}

void Composition::SetTotalLimitsBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(NormTotal);
    BC.SetY(NormTotal);
    BC.SetZ(NormTotal);

    BC.SetXFlags(Limiting);
    BC.SetYFlags(Limiting);
    BC.SetZFlags(Limiting);
}

bool Composition::Write(const Settings& locSettings, const int tStep) const
{

#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_", tStep, ".dat");
#endif

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be created";
        ConsoleOutput::WriteExit(message.str(), thisclassname, "Write()");
        OP_Exit(EXIT_FAILURE);
    };

    int tmp = Grid.Nx;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    tmp = Grid.Ny;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    tmp = Grid.Nz;
    out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    size_t tmp2 = Nphases;
    out.write(reinterpret_cast<char*>(&tmp2), sizeof(size_t));
    tmp2 = Ncomp;
    out.write(reinterpret_cast<char*>(&tmp2), sizeof(size_t));

    STORAGE_LOOP_BEGIN(i,j,k,MoleFractions,0)
        out.write(reinterpret_cast<const char*>(MoleFractions(i,j,k).data()), MoleFractions(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        out.write(reinterpret_cast<const char*>(MoleFractionsTotal(i,j,k).data()), MoleFractionsTotal(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    return true;
}

bool Composition::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{

#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir, thisclassname + "_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be opened";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    int locNx = Grid.Nx;
    int locNy = Grid.Ny;
    int locNz = Grid.Nz;
    size_t locNphases = Nphases;
    size_t locNcomp = Ncomp;
    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNphases), sizeof(size_t));
    inp.read(reinterpret_cast<char*>(&locNcomp), sizeof(size_t));
    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz or locNphases != Nphases or locNcomp != Ncomp)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Input data Nphases: " << locNphases << "\n"
                << "Input data Ncomp: " << locNcomp << "\n"
                << "Required data dimensions: (" << Grid.Nx << ", " << Grid.Ny << ", " << Grid.Nz << ") grid points.\n"
                << "Required data Nphases: " << Nphases << "\n"
                << "Required data Ncomp: " << Ncomp << "\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    STORAGE_LOOP_BEGIN(i,j,k,MoleFractions,0)
        inp.read(reinterpret_cast<char*>(MoleFractions(i,j,k).data()), MoleFractions(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        inp.read(reinterpret_cast<char*>(MoleFractionsTotal(i,j,k).data()), MoleFractionsTotal(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,MoleFractionsTotal,0)
        for(size_t n = 0; n < Ncomp; n++)
        {
            TotInitial({n}) += MoleFractionsTotal(i, j, k)({n});
        }
    STORAGE_LOOP_END
    for(size_t n = 0; n < Ncomp; n++)
    {
        TotInitial({n}) /= double(Grid.LocalNumberOfCells());
    }
    SetBoundaryConditions(BC);

    // Calculation of MoleFractionsTotalAverage is needed for CalculateTotalMolarVolume()
    CalculateMoleFractionsTotalAverage();
    CalculateTotalMolarVolume();
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

void Composition::WriteH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Grid.Nx);
    dbuffer.push_back(Grid.Ny);
    dbuffer.push_back(Grid.Nz);
    dbuffer.push_back(Nphases);
    dbuffer.push_back(Ncomp);
    H5.WriteCheckPoint(tStep, "CxDomain", dbuffer);
    dbuffer.clear();
    dbuffer = MoleFractions.pack();
    H5.WriteCheckPoint(tStep, "MoleFractions", dbuffer);
    dbuffer.clear();
    dbuffer = MoleFractionsTotal.pack();
    H5.WriteCheckPoint(tStep, "MoleFractionsTotal", dbuffer);
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "H5Interface()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

bool Composition::ReadH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "CxDomain", dbuffer);
    int locNx = dbuffer[0];
    int locNy = dbuffer[1];
    int locNz = dbuffer[2];
    size_t locNphases = dbuffer[3];
    size_t locNcomp = dbuffer[4];
    dbuffer.clear();
    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz or locNphases != Nphases or locNcomp != Ncomp)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Input data Nphases: " << locNphases << "\n"
                << "Input data Ncomp: " << locNcomp << "\n"
                << "Required data dimensions: (" << Grid.Nx << ", " << Grid.Ny << ", " << Grid.Nz << ") grid points.\n"
                << "Required data Nphases: " << Nphases << "\n"
                << "Required data Ncomp: " << Ncomp << "\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }
    H5.ReadCheckPoint(tStep, "MoleFractions", dbuffer);
    MoleFractions.unpack(dbuffer);
    dbuffer.clear();
    H5.ReadCheckPoint(tStep, "MoleFractionsTotal", dbuffer);
    MoleFractionsTotal.unpack(dbuffer);
    return true;
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "ReadH5()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
    return false;
}

void Composition::WriteDistortedVTK(Settings& locSettings, const ElasticProperties& EP, const int tStep)
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        const std::string nameComp1 = "WeightFractionsTotal_" + Component[comp].Name;
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return WeightFractionsTotal(i,j,k)({comp});}});
        const std::string nameComp2 = "MoleFractionsTotal_" + Component[comp].Name;
        ListOfFields.push_back((VTK::Field_t) {nameComp2, [this, comp](int i,int j,int k){return MoleFractionsTotal(i,j,k,{comp});}});
        for (size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            const std::string namePhase = "MoleFractionsPhase_" + Component[comp].Name + "(" + to_string(alpha) + ")";
            ListOfFields.push_back((VTK::Field_t) {namePhase, [this, comp, alpha](int i,int j,int k){return MoleFractions(i,j,k,{alpha,comp});}});
        }
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + "Distorted_", tStep, ".vts");

    VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with distorted composition data", "application/xml");
    #endif
}

void Composition::WriteVTK(Settings& locSettings, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        const std::string nameComp1 = "WeightFractionsTotal_" + Component[comp].Name;
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return WeightFractionsTotal(i,j,k)({comp});}});
        const std::string nameComp2 = "MoleFractionsTotal_" + Component[comp].Name;
        ListOfFields.push_back((VTK::Field_t) {nameComp2, [this, comp](int i,int j,int k){return MoleFractionsTotal(i,j,k,{comp});}});
        for (size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            const std::string namePhase = "MoleFractionsPhase_" + Component[comp].Name + "(" + to_string(alpha) + ")";
            ListOfFields.push_back((VTK::Field_t) {namePhase, [this, comp, alpha](int i,int j,int k){return MoleFractions(i,j,k,{alpha,comp});}});
        }
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + '_', tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with composition data", "application/xml");
    #endif
}

void Composition::PrintPointStatistics(int x, int y, int z)
{
    cout << "Components:\t";
    for(size_t comp = 0; comp < Ncomp; comp++) cout << Component[comp].Name << "\t";
    cout << endl;
    cout << "Total:\t\t";
    for(size_t comp = 0; comp < Ncomp; comp++) cout << MoleFractionsTotal(x, y, z,{comp}) << "\t";
    cout << endl;
    cout << "Phase:  " << endl;

    for (size_t n = 0; n < Nphases; n++)
    {
        cout << "Phase " << n << ":\t";
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            cout << MoleFractions(x,y,z,{n, comp}) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Composition::WriteStatistics(const Settings& locSettings, const int tStep, double sim_time)
{
    vector<double> total(Ncomp, 0);
    vector<double> deviation(Ncomp, 0);

    CalculateMoleFractionsTotalAverage();
    if (AtStart)
    {
        TotInitial = MoleFractionsTotalAverage;
    }
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        deviation[comp] = TotInitial({comp}) - MoleFractionsTotalAverage({comp});
    }

    AtStart = false;
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    ofstream output_file;
    if (tStep == 0)
    {
        output_file.open(locSettings.TextDir + "CompositionStatistics.opd", ios::out);
        output_file << "sim_time" << "   ";
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            output_file << "total_" << Component[comp].Name << "   "
                        << "deviation_" << Component[comp].Name << "   ";
        }
        output_file << endl;
        output_file.close();
    }

    output_file.open(locSettings.TextDir + "CompositionStatistics.opd", ios::app);
    output_file << sim_time << "   ";
    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        output_file << MoleFractionsTotalAverage({comp})  << "   " << deviation[comp] << "   ";
    }
    output_file << endl;
    output_file.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void Composition::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    Adv.AdvectField(MoleFractions, MoleFractionsDotIn, Vel, BC, dt);
    Adv.AdvectField(MoleFractionsTotal, MoleFractionsTotalDot, Vel, BC, dt);
    //CalculateTotalMoleFractions(Phi);
}

//======================== Moved from ChemicalProperties =======================

void Composition::PrintData(void)
{
    /**This function prints the names of all components stored in the local
    copy of the thermodynamic system. The indices of each component as stored in
    the TDB/GES5/manual input file, the names of all thermodynamic phases and
    its index, as well as the sublattice model. Vacancies are represented with
    a -1.*/

    ConsoleOutput::WriteLineInsert("Data stored in " + thisclassname +
                                               ". Position in TDB in brackets");

    string ctmp;

    for(size_t i = 0; i < Ncomp; i++)
    {
        ctmp += "[" + to_string(Component[i].Index) + "] " + Component[i].Name;
        if(i < Ncomp-1)
        ctmp += ", ";
    }

    ConsoleOutput::WriteStandard("Components in the system", ctmp);

    string ptmp;

    for(size_t i = 0; i < Nphases; i++)
    {
        ptmp += "[" + to_string(Phase[i].Index) + "] " + Phase[i].Name;
        if(i < Nphases-1)
        ptmp += ", ";
    }

    ConsoleOutput::WriteStandard("Phases in the system", ptmp);

    ConsoleOutput::WriteLine();
}

size_t Composition::PhaseNumber(string name)
{
    /**This function is used to address thermodynamic phases by their names. It
    returns the index, given the phase name as a string.*/
    transform(name.begin(), name.end(), name.begin(), ::toupper);
    long phasenr = -1;
    for(size_t i = 0; i < Nphases; i++)
    if(name == Phase[i].Name)
    {
        phasenr = i;
    }
    if(phasenr < 0)
    {
        cerr << "Error in Composition::PhaseNumber(), no phase available"
             << " with the name of " << name << endl;
        OP_Exit(EXIT_FAILURE);
    }
    return size_t(phasenr);
}

void Composition::CalculateMolefractionLimits(void)
{
    /**This function analyzes the sublattice models of each phase and
    calculates the minimum and maximum limits for the mole fraction of each
    component.*/
    size_t locNphases = Nphases;
    size_t locNcomps  = Ncomp;
    for(size_t alpha = 0; alpha < locNphases; alpha++)
    for(size_t comp = 0; comp < locNcomps; comp++)
    {
        Phase[alpha].Component[comp].Min = 0.0;
        Phase[alpha].Component[comp].Max = 1.0;

        size_t Nsubs = Phase[alpha].Nsubs;
        double totalSitesMAX = 0.0;
        double totalSitesMIN = 0.0;
        double sitesWithCompMAX = 0.0;
        double sitesWithCompMIN = 0.0;
        for(size_t sub = 0; sub < Nsubs; sub++)
        {
            //COUNT MINIMUM
            if(Phase[alpha].Sublattice[sub].isElementPresent(Component[comp].Index)
            and(Phase[alpha].Sublattice[sub].Ncons < 2))
            {
                sitesWithCompMIN += Phase[alpha].Sublattice[sub].Site;
            }
            if(!(Phase[alpha].Sublattice[sub].hasVacancies
            and(Phase[alpha].Sublattice[sub].Ncons < 2)))
            {
                totalSitesMIN += Phase[alpha].Sublattice[sub].Site;
            }
            //COUNT MAXIMUM
            if(Phase[alpha].Sublattice[sub].isElementPresent(Component[comp].Index))
            {
                sitesWithCompMAX += Phase[alpha].Sublattice[sub].Site;
            }
            if(!(Phase[alpha].Sublattice[sub].hasVacancies)
            and!(Phase[alpha].isInterstitial(comp)))
            {
                totalSitesMAX += Phase[alpha].Sublattice[sub].Site;
            }
            else if(Phase[alpha].isInterstitial(comp))
            {
                totalSitesMAX += Phase[alpha].Sublattice[sub].Site;
            }
        }
        if(totalSitesMIN > 0)
        Phase[alpha].Component[comp].Min = sitesWithCompMIN/totalSitesMIN;
        if(totalSitesMAX > 0)
        Phase[alpha].Component[comp].Max = sitesWithCompMAX/totalSitesMAX;
    }
}

vector<size_t> Composition::CalculateSortedElementMatrix(void)
{
    vector<size_t> result(Ncomp);
    vector<string> Names(Ncomp);
    vector<string> SortedNames;

    for(size_t n = 0; n < Ncomp; n++)
    Names[n] = Component[n].Name;

    SortedNames = Names;
    sort(SortedNames.begin(), SortedNames.end());

    for(size_t n = 0; n < Ncomp; n++)
    for(size_t m = 0; m < Ncomp; m++)
    if(Names[n] == SortedNames[m])
    result[m] = n;

    ElementsAreSorted = true;
    for(size_t n = 0; n < Ncomp; n++)
    if(Names[n] != SortedNames[n])
    ElementsAreSorted = false;

    return result;
}

}// namespace openphase
