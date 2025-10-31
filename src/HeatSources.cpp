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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann
 *
 */

#include "HeatSources.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "HeatDiffusion.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;

void HeatSources::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "HeatSources";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    ConsiderHeatSources = false;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void HeatSources::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("HeatSources input");
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

    ConsoleOutput::WriteLine();
}

void HeatSources::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    bool endofsources = false;
    size_t idx = 0;
    while(!endofsources)
    {
        stringstream converter;
        converter << idx;
        string counter = converter.str();
        if(FileInterface::FindParameter(inp, moduleLocation, string("Source_") + counter) != -1)
        {
            ConsoleOutput::WriteLine("=");
            bool SourceEnabled = FileInterface::ReadParameterB(inp, moduleLocation, string("Source_") + counter, true, false);

            if(SourceEnabled)
            {
                HeatSourceStructure locHeatSource;
                locHeatSource.Active = false;
                locHeatSource.Value  = FileInterface::ReadParameterD(inp, moduleLocation, string("Value_") + counter, true, 0.0);

                bool SourceTypeSet = false;
                string SourceType = FileInterface::ReadParameterK(inp, moduleLocation, string("Type_") + counter, true, "NONE");
                if(SourceType == "BCX0"       ) {locHeatSource.Type = HeatSourceTypes::BCX0;        SourceTypeSet = true;}
                if(SourceType == "BCXN"       ) {locHeatSource.Type = HeatSourceTypes::BCXN;        SourceTypeSet = true;}
                if(SourceType == "BCY0"       ) {locHeatSource.Type = HeatSourceTypes::BCY0;        SourceTypeSet = true;}
                if(SourceType == "BCYN"       ) {locHeatSource.Type = HeatSourceTypes::BCYN;        SourceTypeSet = true;}
                if(SourceType == "BCZ0"       ) {locHeatSource.Type = HeatSourceTypes::BCZ0;        SourceTypeSet = true;}
                if(SourceType == "BCZN"       ) {locHeatSource.Type = HeatSourceTypes::BCZN;        SourceTypeSet = true;}
                if(SourceType == "PHASE"      ) {locHeatSource.Type = HeatSourceTypes::Phase;       SourceTypeSet = true;}
                if(SourceType == "ELLIPSOIDAL") {locHeatSource.Type = HeatSourceTypes::Ellipsoidal; SourceTypeSet = true;}

                if(not SourceTypeSet)
                {
                    string message = "No valid type is provided for the heat source " + counter + " --> " + SourceType;
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                switch(locHeatSource.Type)
                {
                    case HeatSourceTypes::Phase:
                    {
                        locHeatSource.PhaseIndex = FileInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndex_") + counter, true, 0);
                        break;
                    }
                    case HeatSourceTypes::Ellipsoidal:
                    {
                        locHeatSource.Position = FileInterface::ReadParameterV3(inp, moduleLocation, string("Position_") + counter, true, {0.0,0.0,0.0});
                        locHeatSource.Size     = FileInterface::ReadParameterV3(inp, moduleLocation, string("Radii_"   ) + counter, true, {0.0,0.0,0.0});
                        locHeatSource.Velocity = FileInterface::ReadParameterV3(inp, moduleLocation, string("Velocity_") + counter, false, {0.0,0.0,0.0});

                        string EnergyDistribution = FileInterface::ReadParameterK(inp, moduleLocation, string("EnergyDistribution_") + counter, false, "UNIFORM");
                        bool EnergyDistributionSet = false;
                        if(EnergyDistribution == "UNIFORM" ) {locHeatSource.Distribution = HeatSourceEnergyDistributions::Uniform;  EnergyDistributionSet = true;}
                        if(EnergyDistribution == "GAUSSIAN") {locHeatSource.Distribution = HeatSourceEnergyDistributions::Gaussian; EnergyDistributionSet = true;}

                        if(not EnergyDistributionSet)
                        {
                            string message = "No valid energy distribution is provided for the heat source " + counter + " --> " + EnergyDistribution;
                            ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                            OP_Exit(EXIT_FAILURE);
                        }
                        break;
                    }
                    default:
                    {
                        break;
                    }
                }

                bool SourceONtriggerSet = false;
                string SourceONtrigger = FileInterface::ReadParameterK(inp, moduleLocation, string("TriggerON_") + counter, true, "USER");
                if(SourceONtrigger == "USER"            ) {locHeatSource.TriggerON = EventTriggers::User;             SourceONtriggerSet = true;}
                if(SourceONtrigger == "TMAX"            ) {locHeatSource.TriggerON = EventTriggers::Tmax;             SourceONtriggerSet = true;}
                if(SourceONtrigger == "TMIN"            ) {locHeatSource.TriggerON = EventTriggers::Tmin;             SourceONtriggerSet = true;}
                if(SourceONtrigger == "TIME"            ) {locHeatSource.TriggerON = EventTriggers::Time;             SourceONtriggerSet = true;}
                if(SourceONtrigger == "TIMESTEP"        ) {locHeatSource.TriggerON = EventTriggers::TimeStep;         SourceONtriggerSet = true;}
                if(SourceONtrigger == "PHASEFRACTIONMAX") {locHeatSource.TriggerON = EventTriggers::PhaseFractionMax; SourceONtriggerSet = true;}
                if(SourceONtrigger == "PHASEFRACTIONMIN") {locHeatSource.TriggerON = EventTriggers::PhaseFractionMin; SourceONtriggerSet = true;}

                if(not SourceONtriggerSet)
                {
                    string message = "No valid ON trigger is provided for the heat source " + counter + " --> " + SourceONtrigger;
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                if(locHeatSource.TriggerON != EventTriggers::User)
                {
                    locHeatSource.ONtriggerValue = FileInterface::ReadParameterD(inp, moduleLocation, string("TriggerONvalue_") + counter, true, 0);
                }
                if(locHeatSource.TriggerON == EventTriggers::PhaseFractionMax or locHeatSource.TriggerON == EventTriggers::PhaseFractionMin)
                {
                    locHeatSource.PhaseIndexON = FileInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndexON_") + counter, true, 0);
                }

                bool SourceOFFtriggerSet = false;
                string SourceOFFtrigger = FileInterface::ReadParameterK(inp, moduleLocation, string("TriggerOFF_") + counter, false, SourceONtrigger);
                if(SourceOFFtrigger == "USER"            ) {locHeatSource.TriggerOFF = EventTriggers::User;             SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "TMAX"            ) {locHeatSource.TriggerOFF = EventTriggers::Tmax;             SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "TMIN"            ) {locHeatSource.TriggerOFF = EventTriggers::Tmin;             SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "TIME"            ) {locHeatSource.TriggerOFF = EventTriggers::Time;             SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "TIMESTEP"        ) {locHeatSource.TriggerOFF = EventTriggers::TimeStep;         SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "PHASEFRACTIONMAX") {locHeatSource.TriggerOFF = EventTriggers::PhaseFractionMax; SourceOFFtriggerSet = true;}
                if(SourceOFFtrigger == "PHASEFRACTIONMIN") {locHeatSource.TriggerOFF = EventTriggers::PhaseFractionMin; SourceOFFtriggerSet = true;}

                if(not SourceOFFtriggerSet)
                {
                    string message = "No valid OFF trigger is provided for the heat source " + counter + " --> " + SourceOFFtrigger;
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                if(locHeatSource.TriggerOFF != EventTriggers::User)
                {
                    locHeatSource.OFFtriggerValue = FileInterface::ReadParameterD(inp, moduleLocation, string("TriggerOFFvalue_") + counter, true, 0);
                }

                if(locHeatSource.TriggerOFF == EventTriggers::PhaseFractionMax or locHeatSource.TriggerOFF == EventTriggers::PhaseFractionMin)
                {
                    locHeatSource.PhaseIndexOFF = FileInterface::ReadParameterI(inp, moduleLocation, string("PhaseIndexOFF_") + counter, true, 0);
                }

                locHeatSource.Repeat = FileInterface::ReadParameterI(inp, moduleLocation, string("Repeat_") + counter, false, 1);
                if(locHeatSource.Repeat <= 0)
                {
                    string message = "Zero or negative \"Repeat\" parameter for heat source " + counter + " --> " + to_string(locHeatSource.Repeat);
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                switch(locHeatSource.Type)
                {
                    case HeatSourceTypes::Phase:
                    {
                        locHeatSource.Mode = HeatSourceModes::ThreeDimensional;
                        break;
                    }
                    case HeatSourceTypes::Ellipsoidal:
                    {
                        if(locHeatSource.Size[0] == 0 or
                           locHeatSource.Size[1] == 0 or
                           locHeatSource.Size[2] == 0)
                        {
                            locHeatSource.Mode = HeatSourceModes::TwoDimensional;
                        }
                        else
                        {
                            locHeatSource.Mode = HeatSourceModes::ThreeDimensional;
                        }
                        if(locHeatSource.Size[0] == 0)
                        {
                            locHeatSource.Size[0] = 0.5;
                        }
                        if(locHeatSource.Size[1] == 0)
                        {
                            locHeatSource.Size[1] = 0.5;
                        }
                        if(locHeatSource.Size[2] == 0)
                        {
                            locHeatSource.Size[2] = 0.5;
                        }
                        break;
                    }
                    default:
                    {
                        locHeatSource.Mode = HeatSourceModes::TwoDimensional;
                    }
                }

                Sources.push_back(locHeatSource);
            }
            idx++;
        }
        else
        {
            endofsources = true;
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

bool HeatSources::Read(const Settings& OPSettings, const int tStep)
{
    std::string FileName = FileInterface::MakeFileName(OPSettings.InputRawDataDir, thisclassname + "_", tStep,".dat");
    return Read(FileName);
}

bool HeatSources::Read(const std::string FileName)
{
    std::ifstream inp(FileName.c_str(), std::ios::in | std::ios::binary);
    if (!inp)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
        return false;
    }

    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));

    if(Sources.size() != size)
    {
        Sources.resize(size);
    }

    for(size_t n = 0; n != size; n++) Sources[n].Read(inp);

    inp.close();
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

bool HeatSources::Write(const Settings& OPSettings, const int tStep) const
{
    std::string FileName =
        FileInterface::MakeFileName(OPSettings.RawDataDir, thisclassname + "_", tStep,".dat");
    return Write(FileName);
}

bool HeatSources::Write(const std::string FileName) const
{
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
    {
#endif
    std::ofstream outp(FileName.c_str(), std::ios::out | std::ios::binary);
    if (!outp)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be created!\n", thisclassname, "Write");
        return false;
    };

    const size_t size = Sources.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(size_t n = 0; n != size; n++) Sources[n].Write(outp);
    outp.close();
#ifdef MPI_PARALLEL
    }
#endif
    return true;
}

void HeatSources::Apply(PhaseField& Phase, Temperature& Tx, HeatDiffusion& HD)
{
    for(auto source = Sources.begin(); source != Sources.end(); source++)
    if(source->Active)
    {
        double source_value = source->Value;
        if(source->Mode == HeatSourceModes::TwoDimensional)
        {
            source_value /= HD.Grid.dx;
        }

        switch(source->Type)
        {
            case HeatSourceTypes::BCX0:
            {
                if(Tx.ExtensionX0.isActive())
                {
                    Tx.ExtensionX0.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetX == 0)
                {
                    int x = 0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int y = 0; y < HD.Grid.Ny; y++)
                    for(int z = 0; z < HD.Grid.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCXN:
            {
                if(Tx.ExtensionXN.isActive())
                {
                    Tx.ExtensionXN.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetX + Tx.Grid.Nx == Tx.Grid.TotalNx)
                {
                    int x = HD.Grid.Nx-1;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int y = 0; y < HD.Grid.Ny; y++)
                    for(int z = 0; z < HD.Grid.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCY0:
            {
                if(Tx.ExtensionY0.isActive())
                {
                    Tx.ExtensionY0.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetY == 0)
                {
                    int y = 0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int x = 0; x < HD.Grid.Nx; x++)
                    for(int z = 0; z < HD.Grid.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCYN:
            {
                if(Tx.ExtensionYN.isActive())
                {
                    Tx.ExtensionYN.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetY + Tx.Grid.Ny == Tx.Grid.TotalNy)
                {
                    int y = HD.Grid.Ny-1;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int x = 0; x < HD.Grid.Nx; x++)
                    for(int z = 0; z < HD.Grid.Nz; z++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCZ0:
            {
                if(Tx.ExtensionZ0.isActive())
                {
                    Tx.ExtensionZ0.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetZ == 0)
                {
                    int z = 0;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int x = 0; x < HD.Grid.Nx; x++)
                    for(int y = 0; y < HD.Grid.Ny; y++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::BCZN:
            {
                if(Tx.ExtensionZN.isActive())
                {
                    Tx.ExtensionZN.Qdot = source_value;
                }
                else if(Tx.Grid.OffsetZ + Tx.Grid.Nz == Tx.Grid.TotalNz)
                {
                    int z = HD.Grid.Nz-1;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                    for(int x = 0; x < HD.Grid.Nx; x++)
                    for(int y = 0; y < HD.Grid.Ny; y++)
                    {
                        HD.Qdot(x,y,z) += source_value;
                    }
                }
                break;
            }
            case HeatSourceTypes::Phase:
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
                {
                    if(Phase.Fractions(i,j,k,{source->PhaseIndex}) > DBL_EPSILON)
                    {
                        HD.Qdot(i,j,k) += source_value * Phase.Fractions(i,j,k,{source->PhaseIndex});
                    }
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                break;
            }
            case HeatSourceTypes::Ellipsoidal:
            {
                double radiusX2 = pow(source->Size[0],2);
                double radiusY2 = pow(source->Size[1],2);
                double radiusZ2 = pow(source->Size[2],2);

                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
                {
                    double radX2 = pow((i+Phase.Grid.OffsetX - source->Position[0]), 2);
                    double radY2 = pow((j+Phase.Grid.OffsetY - source->Position[1]), 2);
                    double radZ2 = pow((k+Phase.Grid.OffsetZ - source->Position[2]), 2);

                    if ((radX2/radiusX2 + radY2/radiusY2 + radZ2/radiusZ2) <= 1.0)
                    {
                        switch(source->Distribution)
                        {
                            case HeatSourceEnergyDistributions::Uniform:
                            {
                                HD.Qdot(i,j,k) += source_value;
                                break;
                            }
                            case HeatSourceEnergyDistributions::Gaussian:
                            {
                                HD.Qdot(i,j,k) += source_value * exp(-0.5*9.0*(radX2/radiusX2 + radY2/radiusY2 + radZ2/radiusZ2));
                                break;
                            }
                        }
                    }
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                break;
            }

            default:
            {
                break;
            }
        }
    }
}

void HeatSources::Activate(PhaseField& Phase, Temperature& Tx, RunTimeControl& RTC)
{
    int counter = 0;
    for(auto source = Sources.begin(); source != Sources.end(); source++)
    {
        if(source->Active)
        {
            switch(source->TriggerOFF)
            {
                case EventTriggers::Tmin:
                {
                    if(Tx.Tmin < source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::Tmax:
                {
                    if(Tx.Tmax > source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }

                case EventTriggers::TimeStep:
                {
                    if(RTC.TimeStep == source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > source->OFFtriggerValue - RTC.dt and
                       RTC.SimulationTime < source->OFFtriggerValue + RTC.dt)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMax:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexOFF] > source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMin:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexOFF] < source->OFFtriggerValue)
                    {
                        source->Active = false;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(!source->Active)
            {
                std::stringstream message;
                message << " Activate Heat Sources: heat source " << counter << " is deactivated.";
                ConsoleOutput::WriteStandard(thisclassname, message.str());
            }
        }
        else if (source->Repeat != 0)
        {
            switch(source->TriggerON)
            {
                case EventTriggers::Tmin:
                {
                    if(Tx.Tmin < source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::Tmax:
                {
                    if(Tx.Tmax > source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }

                case EventTriggers::TimeStep:
                {
                    if(RTC.TimeStep == source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > source->ONtriggerValue - RTC.dt and
                       RTC.SimulationTime < source->ONtriggerValue + RTC.dt)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMax:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexON] > source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::PhaseFractionMin:
                {
                    if(Phase.FractionsTotal[source->PhaseIndexON] < source->ONtriggerValue)
                    {
                        source->Active = true;
                        source->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(source->Active)
            {
                std::stringstream message;
                message << " Activate Heat Sources: heat source " << counter << " is activated.";
                ConsoleOutput::WriteStandard(thisclassname, message.str());
            }
        }
        // Move heat source
        if(source->Active and source->Type == HeatSourceTypes::Ellipsoidal)
        {
            source->Position += source->Velocity*RTC.dt/Tx.Grid.dx;
        }

        counter++;
    }
}

}// namespace openphase
