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

 *   File created :   2023
 *   Main contributors :   Oleg Shchyglo
 *
 */

#include "MechanicalLoads.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "ElasticProperties.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;

void MechanicalLoads::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "MechanicalLoads";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    ConsiderLoads = false;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void MechanicalLoads::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("MechanicalLoads input");
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

void MechanicalLoads::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    bool endofloads = false;
    size_t idx = 0;
    while(!endofloads)
    {
        stringstream converter;
        converter << idx;
        string counter = converter.str();
        if(FileInterface::FindParameter(inp, moduleLocation, string("Load_") + counter) != -1)
        {
            ConsoleOutput::WriteLine("=");

            bool LoadEnabled = FileInterface::ReadParameterB(inp, moduleLocation, string("Load_") + counter, true, false);

            if(LoadEnabled)
            {
                MechanicalLoadStructure locLoad;
                locLoad.Active = false;

                bool LoadONtriggerSet = false;
                string LoadONtrigger = FileInterface::ReadParameterK(inp, moduleLocation, string("TriggerON_") + counter, true, "USER");
                if(LoadONtrigger == "USER")     {locLoad.TriggerON = EventTriggers::User;     LoadONtriggerSet = true;}
                if(LoadONtrigger == "TIME")     {locLoad.TriggerON = EventTriggers::Time;     LoadONtriggerSet = true;}
                if(LoadONtrigger == "TIMESTEP") {locLoad.TriggerON = EventTriggers::TimeStep; LoadONtriggerSet = true;}
                if(LoadONtrigger == "STRESS")   {locLoad.TriggerON = EventTriggers::Stress;   LoadONtriggerSet = true;}
                if(LoadONtrigger == "STRAIN")   {locLoad.TriggerON = EventTriggers::Strain;   LoadONtriggerSet = true;}

                if(not LoadONtriggerSet)
                {
                    string message = "No valid ON trigger is provided for the load " + counter + " --> " + LoadONtrigger;
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                if(locLoad.TriggerON != EventTriggers::User)
                {
                    locLoad.ONtriggerValue = FileInterface::ReadParameterD(inp, moduleLocation, string("TriggerONvalue_") + counter, true, 0);
                }
                if(locLoad.TriggerON == EventTriggers::Strain or locLoad.TriggerON == EventTriggers::Stress)
                {
                    string OnComp = FileInterface::ReadParameterK(inp, moduleLocation, string("ONindex_") + counter, true, "XX");
                    if (OnComp == "XX") locLoad.ONindex = 0;
                    if (OnComp == "YY") locLoad.ONindex = 1;
                    if (OnComp == "ZZ") locLoad.ONindex = 2;
                    if (OnComp == "YZ" or OnComp == "ZY") locLoad.ONindex = 3;
                    if (OnComp == "XZ" or OnComp == "ZX") locLoad.ONindex = 4;
                    if (OnComp == "XY" or OnComp == "YX") locLoad.ONindex = 5;
                }

                bool LoadOFFtriggerSet = false;
                string LoadOFFtrigger = FileInterface::ReadParameterK(inp, moduleLocation, string("TriggerOFF_") + counter, false, LoadONtrigger);
                if(LoadOFFtrigger == "USER")     {locLoad.TriggerOFF = EventTriggers::User;     LoadOFFtriggerSet = true;}
                if(LoadOFFtrigger == "TIME")     {locLoad.TriggerOFF = EventTriggers::Time;     LoadOFFtriggerSet = true;}
                if(LoadOFFtrigger == "TIMESTEP") {locLoad.TriggerOFF = EventTriggers::TimeStep; LoadOFFtriggerSet = true;}
                if(LoadOFFtrigger == "STRESS")   {locLoad.TriggerOFF = EventTriggers::Stress;   LoadOFFtriggerSet = true;}
                if(LoadOFFtrigger == "STRAIN")   {locLoad.TriggerOFF = EventTriggers::Strain;   LoadOFFtriggerSet = true;}

                if(not LoadOFFtriggerSet)
                {
                    string message = "No valid OFF trigger is provided for the load " + counter + " --> " + LoadOFFtrigger;
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                if(locLoad.TriggerOFF != EventTriggers::User)
                {
                    locLoad.OFFtriggerValue = FileInterface::ReadParameterD(inp, moduleLocation, string("TriggerOFFvalue_") + counter, true, 0);
                }

                if(locLoad.TriggerOFF == EventTriggers::Strain or locLoad.TriggerOFF == EventTriggers::Stress)
                {
                    string OffComp = FileInterface::ReadParameterK(inp, moduleLocation, string("OFFindex_") + counter, true, "XX");

                    if (OffComp == "XX") locLoad.OFFindex = 0;
                    if (OffComp == "YY") locLoad.OFFindex = 1;
                    if (OffComp == "ZZ") locLoad.OFFindex = 2;
                    if (OffComp == "YZ" or OffComp == "ZY") locLoad.OFFindex = 3;
                    if (OffComp == "XZ" or OffComp == "ZX") locLoad.OFFindex = 4;
                    if (OffComp == "XY" or OffComp == "YX") locLoad.OFFindex = 5;
                }

                locLoad.Repeat = FileInterface::ReadParameterI(inp, moduleLocation, string("Repeat_") + counter, false, 1);
                if(locLoad.Repeat <= 0)
                {
                    string message = "Zero or negative \"Repeat\" parameter for load " + counter + " --> " + to_string(locLoad.Repeat);
                    ConsoleOutput::WriteExit(message,thisclassname,"ReadInput()");
                    OP_Exit(EXIT_FAILURE);
                }

                /// Mechanical boundary conditions

                std::vector<std::string> BCnames {"BCX", "BCY","BCZ","BCYZ","BCXZ","BCXY"};
                std::vector<std::string> BCvalueNames {"BCValueX", "BCValueY","BCValueZ","BCValueYZ","BCValueXZ","BCValueXY"};

                for(int n = 0; n < 6; n++)
                {
                    string tmp2 = FileInterface::ReadParameterK(inp, moduleLocation, BCnames[n] + "_" + counter, false);

                    if (tmp2 == "STRAIN")
                    {
                        locLoad.AppliedStrain[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n] + "_" + counter);
                        locLoad.AppliedStrainMask[n] = 1.0;
                    }
                    if (tmp2 == "STRAINRATE")
                    {
                        locLoad.AppliedStrainRate[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n] + "_" + counter);
                        locLoad.AppliedStrainRateMask[n] = 1.0;
                    }
                    if (tmp2 == "STRESS")
                    {
                        locLoad.AppliedStressMask[n] = 1.0;
                        locLoad.AppliedStress[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n] + "_" + counter);
                    }
                }

                Loads.push_back(locLoad);
            }
            idx++;
        }
        else
        {
            endofloads = true;
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

bool MechanicalLoads::Read(const Settings& OPSettings, const int tStep)
{
    std::string FileName = FileInterface::MakeFileName(OPSettings.InputRawDataDir, thisclassname + "_", tStep,".dat");
    return Read(FileName);
}

bool MechanicalLoads::Read(const std::string FileName)
{
    std::ifstream inp(FileName.c_str(), std::ios::in | std::ios::binary);
    if (!inp)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
        return false;
    }

    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));

    if(Loads.size() != size)
    {
        Loads.resize(size);
    }

    for(size_t n = 0; n != size; n++) Loads[n].Read(inp);

    inp.close();
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

bool MechanicalLoads::Write(const Settings& OPSettings, const int tStep) const
{
    std::string FileName =
        FileInterface::MakeFileName(OPSettings.RawDataDir, thisclassname + "_", tStep,".dat");
    return Write(FileName);
}

bool MechanicalLoads::Write(const std::string FileName) const
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

    const size_t size = Loads.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(size_t n = 0; n != size; n++) Loads[n].Write(outp);
    outp.close();
#ifdef MPI_PARALLEL
    }
#endif
    return true;
}

void MechanicalLoads::Apply(PhaseField& Phase, ElasticProperties& EP, RunTimeControl& RTC)
{
    for(auto load = Loads.begin(); load != Loads.end(); load++)
    if(load->Active)
    {
        if(!load->Set)
        {
            for(int n = 0; n < 6; n++)
            {
                if(load->AppliedStrainMask[n])
                {
                    EP.AppliedStrainMask[n] = load->AppliedStrainMask[n];
                    EP.AppliedStrain[n]     = load->AppliedStrain[n];
                }
                if(load->AppliedStrainRateMask[n])
                {
                    EP.AppliedStrainRateMask[n] = load->AppliedStrainRateMask[n];
                    EP.AppliedStrainRate[n]     = load->AppliedStrainRate[n];
                }
                if(load->AppliedStressMask[n])
                {
                    EP.AppliedStressMask[n] = load->AppliedStressMask[n];
                    EP.AppliedStress[n]     = load->AppliedStress[n];
                }
            }
            load->Set = true;
        }
    }
    else
    {
        if(load->Set)
        {
            for(int n = 0; n < 6; n++)
            {
                if(load->AppliedStrainMask[n])
                {
                    EP.AppliedStrainMask[n] = 0.0;
                    EP.AppliedStrain[n]     = 0.0;
                }
                if(load->AppliedStrainRateMask[n])
                {
                    EP.AppliedStrainRateMask[n] = 0.0;
                    EP.AppliedStrainRate[n] = 0.0;
                }
                if(load->AppliedStressMask[n])
                {
                    EP.AppliedStressMask[n] = 0.0;
                    EP.AppliedStress[n]     = 0.0;
                }
            }
            load->Set = false;
        }
    }
}

void MechanicalLoads::Activate(PhaseField& Phase, ElasticProperties& EP, RunTimeControl& RTC)
{
    int counter = 0;
    for(auto load = Loads.begin(); load != Loads.end(); load++)
    {
        if(load->Active)
        {
            switch(load->TriggerOFF)
            {
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > load->OFFtriggerValue - RTC.dt and
                       RTC.SimulationTime < load->OFFtriggerValue + RTC.dt)
                    {
                        load->Active = false;
                    }
                    break;
                }
                case EventTriggers::TimeStep:
                {
                    if(RTC.TimeStep == load->OFFtriggerValue)
                    {
                        load->Active = false;
                    }
                    break;
                }
                case EventTriggers::Stress:
                {
                    if(( signbit(load->OFFtriggerValue) and EP.AverageStress[load->OFFindex] < load->OFFtriggerValue) or
                       (!signbit(load->OFFtriggerValue) and EP.AverageStress[load->OFFindex] > load->OFFtriggerValue))
                    {
                        load->Active = false;
                    }
                    break;
                }
                case EventTriggers::Strain:
                {
                    if(( signbit(load->OFFtriggerValue) and EP.AverageStrain[load->OFFindex] < load->OFFtriggerValue) or
                       (!signbit(load->OFFtriggerValue) and EP.AverageStrain[load->OFFindex] > load->OFFtriggerValue))
                    {
                        load->Active = false;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(!load->Active)
            {
                std::stringstream message;
                message << "load " << counter << " is deactivated.";
                ConsoleOutput::WriteStandard(thisclassname, message.str());
            }
        }
        else if (load->Repeat != 0)
        {
            switch(load->TriggerON)
            {
                case EventTriggers::Time:
                {
                    if(RTC.SimulationTime > load->ONtriggerValue - RTC.dt and
                       RTC.SimulationTime < load->ONtriggerValue + RTC.dt)
                    {
                        load->Active = true;
                        load->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::TimeStep:
                {
                    if(RTC.TimeStep == load->ONtriggerValue)
                    {
                        load->Active = true;
                        load->Repeat -= 1;
                    }
                    break;
                }

                case EventTriggers::Stress:
                {
                    if(( signbit(load->ONtriggerValue) and EP.AverageStress[load->ONindex] < load->ONtriggerValue) or
                       (!signbit(load->ONtriggerValue) and EP.AverageStress[load->ONindex] > load->ONtriggerValue))
                    {
                        load->Active = true;
                        load->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::Strain:
                {
                    if(( signbit(load->ONtriggerValue) and EP.AverageStrain[load->ONindex] < load->ONtriggerValue) or
                       (!signbit(load->ONtriggerValue) and EP.AverageStrain[load->ONindex] > load->ONtriggerValue))
                    {
                        load->Active = true;
                        load->Repeat -= 1;
                    }
                    break;
                }
                case EventTriggers::User:
                default:
                {
                    break;
                }
            }
            if(load->Active)
            {
                std::stringstream message;
                message << "load " << counter << " is activated.";
                ConsoleOutput::WriteStandard(thisclassname, message.str());
            }
        }
        counter++;
    }
}

}// namespace openphase
