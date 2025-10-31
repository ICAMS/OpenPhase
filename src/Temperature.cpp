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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#include "Temperature.h"
#include "AdvectionHR.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"
#include "RunTimeControl.h"
#include "Velocities.h"

namespace openphase
{
using namespace std;

inline pair<int, int> three_state_bounds_selector(int selector, int lower_bound, int upper_bound)
{
    pair<int, int> result;
    switch (selector)
    {
        case -1: // lower bound only
        {
            result = make_pair(lower_bound, lower_bound + 1);
            break;
        }
        case 0: // both bounds
        {
            result = make_pair(lower_bound, upper_bound);
            break;
        }
        case 1: // upper bound only
        {
            result = make_pair(upper_bound - 1, upper_bound);
            break;
        }
    }
    return result;
}

Temperature::Temperature(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Temperature::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Temperature";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    Tmin = 0.0;
    Tmax = 0.0;
    Tavg = 0.0;

    Tiavg.Allocate({Nphases,Nphases});

    ReadFromFile = false;
    ExtensionsActive = false;

    LatentHeatMode = LatentHeatModes::Off;

    size_t Bcells = Grid.Bcells;
    Tx   .Allocate(Grid, Bcells);
    TxOld.Allocate(Grid, Bcells);

    HeatCapacity.Allocate(Nphases);
    LatentHeat.Allocate({Nphases,Nphases});

    dT_dr.set_to_zero();
    r0.set_to_zero();

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Temperature::ReadInput(const string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());

    ReadInput(data);

    ConsoleOutput::WriteLine();
}

void Temperature::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    ReadFromFile = FileInterface::ReadParameterB(inp, moduleLocation,"ReadfromFile", false, false);

    r0[0] = FileInterface::ReadParameterD(inp, moduleLocation, string("R0X"), false, 0.0);
    r0[1] = FileInterface::ReadParameterD(inp, moduleLocation, string("R0Y"), false, 0.0);
    r0[2] = FileInterface::ReadParameterD(inp, moduleLocation, string("R0Z"), false, 0.0);

    dT_dr[0] = FileInterface::ReadParameterD(inp, moduleLocation, string("DT_DRX"), false, 0.0);
    dT_dr[1] = FileInterface::ReadParameterD(inp, moduleLocation, string("DT_DRY"), false, 0.0);
    dT_dr[2] = FileInterface::ReadParameterD(inp, moduleLocation, string("DT_DRZ"), false, 0.0);

    TBC0X = FileInterface::ReadParameterD(inp, moduleLocation, string("TBC0X"), false, 0.0);
    TBCNX = FileInterface::ReadParameterD(inp, moduleLocation, string("TBCNX"), false, 0.0);
    TBC0Y = FileInterface::ReadParameterD(inp, moduleLocation, string("TBC0Y"), false, 0.0);
    TBCNY = FileInterface::ReadParameterD(inp, moduleLocation, string("TBCNY"), false, 0.0);
    TBC0Z = FileInterface::ReadParameterD(inp, moduleLocation, string("TBC0Z"), false, 0.0);
    TBCNZ = FileInterface::ReadParameterD(inp, moduleLocation, string("TBCNZ"), false, 0.0);

    if(ReadFromFile)
    {
        string DataFile = FileInterface::ReadParameterF(inp, moduleLocation, "DataFile");

        fstream inp(DataFile.c_str(), ios::in);
        if (!inp)
        {
            ConsoleOutput::WriteExit("File \"" + DataFile + "\" could not be opened",
                            thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        };

        double activation_time;
        double value;
        while(inp >> activation_time >> value)
        {
            pair<double,double> tempProfile;
            tempProfile.first  = activation_time;
            tempProfile.second = value;
            TemperatureProfile.push_back(tempProfile);

            ConsoleOutput::WriteSimple("Temperature at " + to_string(tempProfile.first) + " s is set to: " + to_string(tempProfile.second));
        }

        if(TemperatureProfile.size() == 0)
        {
            ConsoleOutput::WriteExit("Could not extract temperature information from file \"" + DataFile + "\".",
                            thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }

        T0 = TemperatureProfile[0].second;

        inp.close();
    }
    else
    {
        T0 = FileInterface::ReadParameterD(inp, moduleLocation, string("T0"));

        int n = 0;
        while(FileInterface::FindParameter(inp, moduleLocation, "ControlMode_" + to_string(n)) > 0)
        {
            TemperatureParameters localParameters;

            std::string sControlMode = FileInterface::ReadParameterK(inp, moduleLocation,
                                       string("ControlMode_" + to_string(n)), false, "NONE");
            bool ControlModeSet = false;
            if (sControlMode == "LINEAR")
            {
                localParameters.Mode = ControlModes::Linear;
                localParameters.LinearRate = FileInterface::ReadParameterD(inp, moduleLocation,
                                                    string("LinearRate_" + to_string(n)), false, 0.0);
                ControlModeSet = true;
            }
            if (sControlMode == "NEWTON")
            {
                localParameters.Mode = ControlModes::Newton;
                localParameters.NewtonRate = FileInterface::ReadParameterD(inp, moduleLocation,
                                                    string("NewtonRate_" + to_string(n)), false, 0.0);
                ControlModeSet = true;
            }
            if (sControlMode == "NONE" or sControlMode == "NO" or sControlMode == "OFF")
            {
                localParameters.Mode = ControlModes::None;
                ControlModeSet = true;
            }

            if(!ControlModeSet)
            {
                std::cerr << "Unknown ControlMode selected." << std::endl;
                OP_Exit(EXIT_CONTROLMODE_ERROR);
            }

            if(localParameters.Mode != ControlModes::None)
            {
                localParameters.HoldingTemperature = FileInterface::ReadParameterD(inp, moduleLocation,
                                                     string("HoldingTemperature_" + to_string(n)), false, 0);
                localParameters.ActivationTime     = FileInterface::ReadParameterD(inp, moduleLocation,
                                                     string("ActivationTime_" + to_string(n)), false, 0.0);
            }

            ControlParameters.push_back(localParameters);
            n++;
        }

        std::string sLatentHeatMode = FileInterface::ReadParameterK(inp, moduleLocation,"LatentHeatMode", false, "OFF");
        bool LatentHeatModeSet = false;
        if (sLatentHeatMode == "OFF" or sLatentHeatMode == "NO")
        {
            LatentHeatMode = LatentHeatModes::Off;
            LatentHeatModeSet = true;
        }
        if (sLatentHeatMode == "LOCAL")
        {
            LatentHeatMode = LatentHeatModes::Local;
            LatentHeatModeSet = true;
        }
        if (sLatentHeatMode == "GLOBAL")
        {
            LatentHeatMode = LatentHeatModes::Global;
            LatentHeatModeSet = true;
        }
        if(!LatentHeatModeSet)
        {
            std::cerr << "Unknown LatentHeatMode selected. " << std::endl;
            OP_Exit(EXIT_LATENTHEATMODE_ERROR);
        }

        if(LatentHeatMode != LatentHeatModes::Off)
        {
            for (size_t n = 0; n < Nphases; n++)
            {
                stringstream converter;
                converter << n;
                string counter = converter.str();

                ConsoleOutput::WriteBlankLine();

                HeatCapacity[n] = FileInterface::ReadParameterD(inp, moduleLocation, {string("VolumetricHeatCapacity_") + counter, string("HeatCapacity_") + counter}, false, 0.0);
            }
            for (size_t n = 0; n < Nphases; n++)
            for (size_t m = n; m < Nphases; m++)
            {
                if(n != m)
                {
                    stringstream idx;
                    idx << "_" << n << "_" << m;
                    LatentHeat({n,m}) = FileInterface::ReadParameterD(inp, moduleLocation, "LatentHeat" + idx.str(), false, 0.0);
                    LatentHeat({m,n}) = -LatentHeat({n,m});
                }
            }
        }
    }

    if(Grid.dNx)
    {
        int X0_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_X0", false, 0);
        int XN_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_XN", false, 0);

        if(X0_size > 0) {ExtensionX0.Initialize(X0_size,(iVector3){-1, 0, 0}); ExtensionsActive = true;};
        if(XN_size > 0) {ExtensionXN.Initialize(XN_size,(iVector3){ 1, 0, 0}); ExtensionsActive = true;};
    }
    if(Grid.dNy)
    {
        int Y0_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Y0", false, 0);
        int YN_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_YN", false, 0);
        if(Y0_size > 0) {ExtensionY0.Initialize(Y0_size,(iVector3){ 0,-1, 0}); ExtensionsActive = true;};
        if(YN_size > 0) {ExtensionYN.Initialize(YN_size,(iVector3){ 0, 1, 0}); ExtensionsActive = true;};
    }
    if(Grid.dNz)
    {
        int Z0_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_Z0", false, 0);
        int ZN_size = FileInterface::ReadParameterI(inp, moduleLocation, "Extension_ZN", false, 0);
        if(Z0_size > 0) {ExtensionZ0.Initialize(Z0_size,(iVector3){ 0, 0,-1}); ExtensionsActive = true;};
        if(ZN_size > 0) {ExtensionZN.Initialize(ZN_size,(iVector3){ 0, 0, 1}); ExtensionsActive = true;};
    }
}
double Temperature::CalculateLatentHeatEffect(const PhaseField& Phase,
                                              const double dt)
{
    size_t Nthreads = 1;

#ifdef _OPENMP
    Nthreads = omp_get_max_threads();
#endif

    /* Calculate phase fractions increments. */

    std::vector<Matrix<double> > nFractionsDot(Nthreads);

    for (size_t i = 0; i < Nthreads; ++i)
    {
        nFractionsDot[i].Allocate(Phase.Nphases,Phase.Nphases);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
    {
        size_t thread = 0;

        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it  = Phase.FieldsDot(i,j,k).cbegin();
                 it != Phase.FieldsDot(i,j,k).cend(); ++it)
        {
            size_t PIdxA = Phase.FieldsProperties[it->indexA].Phase;
            size_t PIdxB = Phase.FieldsProperties[it->indexB].Phase;

            nFractionsDot[thread](PIdxA,PIdxB) += it->value1;
            nFractionsDot[thread](PIdxB,PIdxA) -= it->value1;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Matrix<double> FractionsDot(Phase.Nphases,Phase.Nphases);

    for (size_t t = 0; t < Nthreads; ++t)
    for (size_t n = 0; n < Phase.Nphases; ++n)
    for (size_t m = 0; m < Phase.Nphases; ++m)
    {
        FractionsDot(n,m) += nFractionsDot[t](n,m);
    }

    /* Calculate global effective latent heat. */
    double AverageLatentHeat = 0.0;                                             // average volumetric latent heat [J/(m^3)]

    for (size_t n = 0; n < Phase.Nphases; ++n)
    for (size_t m = n; m < Phase.Nphases; ++m)
    {
        AverageLatentHeat += LatentHeat({n,m})*FractionsDot(n,m)*dt;
    }

    /* Calculate global effective volumetric heat capacity. */
    double AverageHeatCapacity = 0.0;                                           // average volumetric heat capacity [J/(m^3 K)]

    for (size_t n = 0; n < Phase.Nphases; ++n)
    {
        AverageHeatCapacity += Phase.FractionsTotal[n]*HeatCapacity[n];
    }

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &AverageLatentHeat, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    AverageLatentHeat /= double(Grid.TotalNumberOfCells());

    return AverageLatentHeat/AverageHeatCapacity;
}

void Temperature::CalculateMinMaxAvg()
{
    double locTmin = 6000.0;
    double locTmax =    0.0;
    double locTavg =    0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,reduction(min:locTmin) reduction(max:locTmax) reduction(+:locTavg))
    {
        locTmin = min(locTmin,Tx(i,j,k));
        locTmax = max(locTmax,Tx(i,j,k));
        locTavg += Tx(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locTmin, 1, OP_MPI_DOUBLE, OP_MPI_MIN, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locTmax, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locTavg, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    Tmin = locTmin;
    Tmax = locTmax;
    Tavg = locTavg/double(Grid.TotalNumberOfCells());
    Tiavg.set_to_value(Tavg);
}

void Temperature::SetInitial1Dextension(Temperature1Dextension& TxExt)
{
    // Get average temperature at the boundary of interest
    std::pair<int, int> Xbounds = three_state_bounds_selector(TxExt.Direction[0], 0, Grid.Nx);
    std::pair<int, int> Ybounds = three_state_bounds_selector(TxExt.Direction[1], 0, Grid.Ny);
    std::pair<int, int> Zbounds = three_state_bounds_selector(TxExt.Direction[2], 0, Grid.Nz);

    double boundaryValue = 0.0;
    double area          = 0.0;

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        boundaryValue += Tx(i,j,k);
        area++;
    }

#ifdef MPI_PARALLEL
    if(TxExt.Direction[0] == 0)
    {
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &boundaryValue, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &area, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
#endif

    TxExt.Data(-1) = boundaryValue/area;

    // Set initial temperature in the extension
    int offset_x = 0.5*Grid.Nx*TxExt.Direction[1]*TxExt.Direction[2];
    int offset_y = 0.5*Grid.Ny*TxExt.Direction[0]*TxExt.Direction[2];
    int offset_z = 0.5*Grid.Nz*TxExt.Direction[0]*TxExt.Direction[1];

    bool limit = false;

    for (size_t i = 0; i <= TxExt.size(); ++i)
    {
        TxExt.Data(i) = TxExt.Data(-1) + (dT_dr[0]*(offset_x + i*TxExt.Direction[0] - r0[0]) +
                                          dT_dr[1]*(offset_y + i*TxExt.Direction[1] - r0[1]) +
                                          dT_dr[2]*(offset_z + i*TxExt.Direction[2] - r0[2]))*Grid.dx;
        if(TxExt.Data(i) < 0.0)
        {
            limit = true;
        }
    }

    if(limit)
    {
        ConsoleOutput::WriteExit("Negative temperature detected!\n", thisclassname, "SetInitial1Dextension()");
        OP_Exit(EXIT_FAILURE);
    }
}

void Temperature::SetInitial(const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,Tx.Bcells(),)
    {
        Tx(i,j,k) = T0 + (dT_dr[0]*(i - r0[0]) +
                          dT_dr[1]*(j - r0[1]) +
                          dT_dr[2]*(k - r0[2]))*Grid.dx;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    //NX0
    if (TBC0X != 0.0)
    for(long int i = -Tx.BcellsX(); i < 0; i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBC0X;
    }

    //NY0
    if (TBC0Y != 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < 0; j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBC0Y;
    }

    //NZ0
    if (TBC0Z != 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < 0; k++)
    {
        Tx(i,j,k) = TBC0Z;
    }

    //NXN
    if (TBCNX != 0.0)
    for(long int i =    Tx.sizeX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNX;
    }

    //NYN
    if (TBCNY != 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j =    Tx.sizeY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k = -Tx.BcellsZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNY;
    }

    //NZN
    if (TBCNZ != 0.0)
    for(long int i = -Tx.BcellsX(); i < Tx.sizeX() + Tx.BcellsX(); i++)
    for(long int j = -Tx.BcellsY(); j < Tx.sizeY() + Tx.BcellsY(); j++)
    for(long int k =    Tx.sizeZ(); k < Tx.sizeZ() + Tx.BcellsZ(); k++)
    {
        Tx(i,j,k) = TBCNZ;
    }

#ifdef MPI_PARALLEL
    if(MPI_3D_DECOMPOSITION)
    {
        BC.CommunicateX(Tx);
        BC.CommunicateY(Tx);
        BC.CommunicateZ(Tx);
    }
    else
    {
        BC.Communicate(Tx);
    }
#endif

    CalculateMinMaxAvg();
    SetBoundaryConditions(BC);

    if(Tmin < 0.0)
    {
        stringstream message;
        message << "Negative temperature detected!\n";
        ConsoleOutput::WriteExit(message.str(), thisclassname, "SetInitial()");
        OP_Exit(EXIT_FAILURE);
    }

#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    if(ExtensionX0.isActive()) {SetInitial1Dextension(ExtensionX0);}
#ifdef MPI_PARALLEL
    if(MPI_RANK == MPI_SIZE - 1)
#endif
    if(ExtensionXN.isActive()) {SetInitial1Dextension(ExtensionXN);}

    if(ExtensionY0.isActive()) {SetInitial1Dextension(ExtensionY0);}
    if(ExtensionYN.isActive()) {SetInitial1Dextension(ExtensionYN);}
    if(ExtensionZ0.isActive()) {SetInitial1Dextension(ExtensionZ0);}
    if(ExtensionZN.isActive()) {SetInitial1Dextension(ExtensionZN);}
}

void Temperature::Set(const BoundaryConditions& BC, const PhaseField& Phase,
                      const double simulation_time, const double dt)
{
    if(ReadFromFile)
    {
        /* Read temperature profile over time from file:
        In this case the simulation time is used to interpolate the temperature
        using the input data from the supplied external file.
        In the input the parameters $ReadfromFile and $DataFile should be
        specified as well as an external file with the only content of two
        double values separated by blank spaces. The first column should contain
        the time and the second should provide the corresponding temperature
        value, e.g.:
            0.0    1050
            1E3    900
            1E5    293
        Outside of the given time range, the last given temperature value is
        used. */

        double TemperatureInput = 293;
        size_t nEntries = TemperatureProfile.size();

        if(nEntries == 0)
        {
            ConsoleOutput::WriteExit("Could not extract temperature information from "
                            "internal Temperature Profile.",
                            thisclassname, "Set()");
            OP_Exit(EXIT_FAILURE);
        }

        if(nEntries == 1)
        {
            TemperatureInput = TemperatureProfile[0].second;
        }
        else
        {
            if(simulation_time <= TemperatureProfile[0].first)
            {
                TemperatureInput = TemperatureProfile[0].second;
            }
            else if(simulation_time >= TemperatureProfile[nEntries-1].first)
            {
                TemperatureInput = TemperatureProfile[nEntries-1].second;
            }
            else
            {
                for(size_t i = 1; i < nEntries; i++)
                if(simulation_time >= TemperatureProfile[i-1].first and
                   simulation_time <  TemperatureProfile[i].first)
                {
                    TemperatureInput = TemperatureProfile[i-1].second;
                    double locTimeRange  = TemperatureProfile[i].first - TemperatureProfile[i-1].first;

                    if(locTimeRange > DBL_EPSILON)
                    {
                        double locTimeOffset = simulation_time - TemperatureProfile[i-1].first;
                        double locTempRange  = TemperatureProfile[i].second - TemperatureProfile[i-1].second;

                        TemperatureInput += locTimeOffset*locTempRange/locTimeRange;
                    }

                    break;
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,Tx.Bcells(),)
        {
            Tx(i, j, k) = TemperatureInput + (dT_dr[0] * (i - r0[0]) + dT_dr[1] * (j - r0[1]) + dT_dr[2] * (k - r0[2]))*Grid.dx;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        /* Calculate latent heat effect */
        if(LatentHeatMode == LatentHeatModes::Global)
        {
            double LatentHeatEffect = CalculateLatentHeatEffect(Phase, dt);
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
            {
                Tx(i,j,k) += LatentHeatEffect;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        }

        long int ControlCycle = -1;
        for(size_t n = 0; n < ControlParameters.size(); n++)
        if(simulation_time >= ControlParameters[n].ActivationTime)
        {
            ControlCycle = n;
        }

        /* Apply temperature increment in all points. */
        if(ControlCycle != -1)
        switch(ControlParameters[ControlCycle].Mode)
        {
            case ControlModes::Linear:
            {
                if((ControlParameters[ControlCycle].HoldingTemperature > Tavg and
                    ControlParameters[ControlCycle].LinearRate > 0) or
                   (ControlParameters[ControlCycle].HoldingTemperature < Tavg and
                    ControlParameters[ControlCycle].LinearRate < 0))
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
                    {
                        Tx(i,j,k) += dt*ControlParameters[ControlCycle].LinearRate;
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                }
                break;
            }
            case ControlModes::Newton:
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
                {
                    Tx(i,j,k) += dt * ControlParameters[ControlCycle].NewtonRate * (ControlParameters[ControlCycle].HoldingTemperature - Tx(i,j,k));
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                break;
            }
            case ControlModes::None:
            default:
            {
                break;
            }
        }
    }

    CalculateMinMaxAvg();
    SetBoundaryConditions(BC);

    if(Tmin < 0.0)
    {
        stringstream message;
        message << "Negative temperature detected!\n";
        ConsoleOutput::WriteExit(message.str(), thisclassname, "Set()");
        OP_Exit(EXIT_FAILURE);
    }
}

void Temperature::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if (Grid.dNx > 0) BC.SetX(Tx);
    if (Grid.dNy > 0) BC.SetY(Tx);
    if (Grid.dNz > 0) BC.SetZ(Tx);
}

void Temperature::MoveFrame(const int dx, const int dy, const int dz,
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
         Tx(i,j,k) = Tx(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
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

void Temperature::ConsumePlane(const int dx, const int dy, const int dz,
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
        Tx(i,j,k) = Tx(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
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
        Tx(i,j,k) = Tx(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
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

bool Temperature::Write(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,thisclassname+"_", tStep, ".dat");
#endif
    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        return false;
    };

    out.write(reinterpret_cast<const char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    out.close();

    if(ExtensionsActive)
    {
        /* Write Extensions */
#ifdef MPI_PARALLEL
        string FileName2 = FileInterface::MakeFileName(locSettings.RawDataDir,"TemperatureExtensions_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
        string FileName2 = FileInterface::MakeFileName(locSettings.RawDataDir,"TemperatureExtensions_", tStep, ".dat");
#endif

        ofstream out2(FileName2.c_str(), ios::out | ios::binary);

        if (!out2)
        {
            ConsoleOutput::WriteWarning("File \"" + FileName2 + "\" could not be created", thisclassname, "Write()");
            return false;
        };

        if(ExtensionX0.isActive()) ExtensionX0.write(out2);
        if(ExtensionXN.isActive()) ExtensionXN.write(out2);
        if(ExtensionY0.isActive()) ExtensionY0.write(out2);
        if(ExtensionYN.isActive()) ExtensionYN.write(out2);
        if(ExtensionZ0.isActive()) ExtensionZ0.write(out2);
        if(ExtensionZN.isActive()) ExtensionZN.write(out2);
    }
    return true;
}

bool Temperature::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
        return false;
    };

    inp.read(reinterpret_cast<char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    inp.close();

    CalculateMinMaxAvg();
    SetBoundaryConditions(BC);

    if(ExtensionsActive)
    {
        /* Read Extensions */
#ifdef MPI_PARALLEL
        string FileName2 = FileInterface::MakeFileName(locSettings.InputRawDataDir,"TemperatureExtensions_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
        string FileName2 = FileInterface::MakeFileName(locSettings.InputRawDataDir,"TemperatureExtensions_", tStep, ".dat");
#endif

        ifstream inp2(FileName2.c_str(), ios::in | ios::binary);

        if (!inp2)
        {
            ConsoleOutput::WriteWarning("File \"" + FileName2 + "\" could not be opened", thisclassname, "Read()");
            return false;
        };

        if(ExtensionX0.isActive()) ExtensionX0.read(inp2);
        if(ExtensionXN.isActive()) ExtensionXN.read(inp2);
        if(ExtensionY0.isActive()) ExtensionY0.read(inp2);
        if(ExtensionYN.isActive()) ExtensionYN.read(inp2);
        if(ExtensionZ0.isActive()) ExtensionZ0.read(inp2);
        if(ExtensionZN.isActive()) ExtensionZN.read(inp2);
    }
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");

    return true;
}

void Temperature::WriteH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Nx);
    dbuffer.push_back(Ny);
    dbuffer.push_back(Nz);
    H5.WriteCheckPoint(tStep, "TxDomain", dbuffer);
    dbuffer.clear();

    dbuffer = Tx.pack();
    H5.WriteCheckPoint(tStep, "Tx", dbuffer);

    if(ExtensionX0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExX0", ExtensionX0.Data);
    }
    if(ExtensionXN.isActive()){
        H5.WriteCheckPoint(tStep, "TxExXN", ExtensionXN.Data);
    }
    if(ExtensionY0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExY0", ExtensionY0.Data);
    }
    if(ExtensionYN.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExYN", ExtensionYN.Data);
    }
    if(ExtensionZ0.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExZ0", ExtensionZ0.Data);
    }
    if(ExtensionZN.isActive()) {
        H5.WriteCheckPoint(tStep, "TxExZN", ExtensionZN.Data);
    }

    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteH5()");
    OP_Exit(EXIT_FAILURE);
    #endif
}

bool Temperature::ReadH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "PFDomain", dbuffer);
    int locNx = dbuffer[0];
    int locNy = dbuffer[1];
    int locNz = dbuffer[2];
    dbuffer.clear();
    if(locNx != Nx or locNy != Ny or locNz != Nz)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    H5.ReadCheckPoint(tStep, "Tx", dbuffer);
    Tx.unpack(dbuffer);

    if(ExtensionX0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExX0", ExtensionX0.Data);
    }
    if(ExtensionXN.isActive()){
        H5.ReadCheckPoint(tStep, "TxExXN", ExtensionXN.Data);
    }
    if(ExtensionY0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExY0", ExtensionY0.Data);
    }
    if(ExtensionYN.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExYN", ExtensionYN.Data);
    }
    if(ExtensionZ0.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExZ0", ExtensionZ0.Data);
    }
    if(ExtensionZN.isActive()) {
        H5.ReadCheckPoint(tStep, "TxExZN", ExtensionZN.Data);
    }

    return true;
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "ReadH5()");
    OP_Exit(EXIT_FAILURE);
    #endif
    return false;
}

void Temperature::Remesh(const int newNx, const int newNy, const int newNz,
                                                    const BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);

    Grid.SetDimensions(newNx, newNy, newNz);

    Tx.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    if(TxDot.IsAllocated())
    {
        TxDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    }

    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void Temperature::PrintPointStatistics(const int x, const int y, const int z) const
{
    ConsoleOutput::WriteStandard("Point", iVector3{x,y,z});
    ConsoleOutput::WriteStandard(thisclassname, Tx(x, y, z));
    ConsoleOutput::WriteBlankLine();
}
void Temperature::PrintStatistics() const
{
    std::string message  = "\n";
                message += ConsoleOutput::GetStandard("Tavg = ", Tavg);
                message += ConsoleOutput::GetStandard("Tmin = ", Tmin);
                message += ConsoleOutput::GetStandard("Tmax = ", Tmax);
    ConsoleOutput::WriteStandard(thisclassname, message);
}
void Temperature::WriteVTK(Settings& locSettings, const int tStep) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"T", [this](int i,int j,int k){return Tx(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with temperature data", "application/xml");
    #endif
}
void Temperature::WriteGradientVTK(Settings& locSettings, const int tStep) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"Gradient_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"dT_dx", [this](int i,int j,int k){return (dVector3){0.5*(Tx(i+1,j,k)-Tx(i-1,j,k))/Grid.dx,0.5*(Tx(i,j+1,k)-Tx(i,j-1,k))/Grid.dx,0.5*(Tx(i,j,k+1)-Tx(i,j,k-1))/Grid.dx};}});
    VTK::Write(Filename, locSettings, ListOfFields);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with temperature gradient data", "application/xml");
    #endif
}

void Temperature::WriteMinMaxAverage(int time_step, double time, string filename) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        string separator = " ";

        if(!std::filesystem::exists(filename) or time_step <= 0)
        {
            ofstream file(filename, ios::out);
            file << "TimeStep" << separator
                 << "Time"     << separator
                 << "Min"      << separator
                 << "Max"      << separator
                 << "Average"  << endl;
            file.close();
        }

        ofstream file(filename, ios::app);
        file << time_step << separator
             << time      << separator
             << Tmin      << separator
             << Tmax      << separator
             << Tavg      << endl;
        file.close();
    }
}

Temperature& Temperature::operator=(const Temperature& rhs)
{
    // protect against invalid self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "Temperature")
    {
        thisclassname = rhs.thisclassname;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;

        dT_dr = rhs.dT_dr;
        r0    = rhs.r0;
        T0    = rhs.T0;

        Tmin  = rhs.Tmin;
        Tmax  = rhs.Tmax;
        Tavg  = rhs.Tavg;
        Tiavg = rhs.Tiavg;

        Tx    = rhs.Tx;
        TxDot = rhs.TxDot;
        TxOld = rhs.TxOld;

        ExtensionX0 = rhs.ExtensionX0;
        ExtensionXN = rhs.ExtensionXN;
        ExtensionY0 = rhs.ExtensionY0;
        ExtensionYN = rhs.ExtensionYN;
        ExtensionZ0 = rhs.ExtensionZ0;
        ExtensionZN = rhs.ExtensionZN;
        ExtensionsActive = rhs.ExtensionsActive;

        LatentHeatMode   = rhs.LatentHeatMode;
        HeatCapacity     = rhs.HeatCapacity;
        LatentHeat       = rhs.LatentHeat;

        ReadFromFile       = rhs.ReadFromFile;
        TemperatureProfile = rhs.TemperatureProfile;

        ControlParameters  = rhs.ControlParameters;
    }
    return *this;
}

void Temperature::Advect(AdvectionHR& Adv, const Velocities& Vel,
                         PhaseField& Phi, const BoundaryConditions& BC,
                         const double dt, const double tStep)
{
    if(TxDot.IsNotAllocated()) TxDot.Allocate(Grid,Tx.Bcells());

    Adv.AdvectField(Tx, TxDot, Vel, BC, dt);
}

void Temperature::CalculateInterfaceAverage(PhaseField& Phi)
{
    Tensor<double, 2> locTiavg({Nphases,Nphases});
    Tensor<double, 2> NpointsAB({Nphases,Nphases});

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Tx,0,reduction(TensorD2Sum:locTiavg,NpointsAB))
    {
        if(Phi.Fields(x,y,z).interface())
        {
            for (size_t alpha =     0; alpha < Nphases; ++alpha)
            for (size_t  beta = alpha;  beta < Nphases; ++beta)
            if(Phi.Fractions(x,y,z,{alpha}) != 0.0 and
               Phi.Fractions(x,y,z,{ beta}) != 0.0)
            {
                locTiavg({alpha,beta}) += Tx(x,y,z);
                NpointsAB({alpha,beta})++;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Tiavg = locTiavg;

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, NpointsAB.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, Tiavg.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    for(size_t alpha =     0; alpha < Nphases; alpha++)
    for(size_t  beta = alpha;  beta < Nphases; beta++)
    {
        if(NpointsAB({alpha, beta}) != 0)
        {
            Tiavg({alpha, beta}) /= NpointsAB({alpha, beta});
        }
        else
        {
            Tiavg({alpha, beta}) = Tavg;
        }
        Tiavg({beta, alpha}) = Tiavg({alpha, beta});
    }
}

}// namespace openphase
