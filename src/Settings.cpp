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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "Settings.h"

#include "GridParameters.h"
#include "RunTimeControl.h"
#include "BoundaryConditions.h"
#include "OPObject.h"
#include "PhaseField.h"

namespace openphase
{

using namespace std;

void Settings::Initialize(std::string ObjectNameSuffix)
{
    thisclassname = "Settings";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = 0;
    Ncomp = 0;

    if(VTKDir.size() == 0) VTKDir = DefaultVTKDir;                              //Only overwrite default directory if VTKDir is empty. This way, VTKDir can be defined earlier.
    if(RawDataDir.size() == 0) RawDataDir = DefaultRawDataDir;                  //Only overwrite default directory if RawDataDir is empty. This way, RawDataDir can be defined earlier.
    if(InputRawDataDir.size() == 0) InputRawDataDir = DefaultRawDataDir;        //Only overwrite default directory if InputRawDataDir is empty. This way, RawDataDir can be defined earlier.
    if(TextDir.size() == 0) TextDir = DefaultTextDir;                           //Only overwrite default directory if TextDir is empty. This way, TextDir can be defined earlier.

    initialized = true;
    ConsoleOutput::WriteStartScreen();
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Settings::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::string filetype = FileInterface::getFileExtension(InputFileName); 
    if (filetype == "opi")
    {
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
    else
    if (filetype == "json")
    {
        ReadJSON(InputFileName);
    }
    else
    {
        std::cerr << "Filetype " << filetype << " not recognized. Filetype must be opi or json." << std::endl;
        OP_Exit(1);
    }
}

void Settings::ReadInput(std::stringstream& inp)
{
    Grid.ReadInput(inp);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Settings");

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    ConsoleOutput::WriteLineInsert("Active phases");

    bool endofnames = false;
    size_t n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Phase_") << n;
        if(FileInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = FileInterface::ReadParameterK(inp, moduleLocation, converter.str());
            PhaseNames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }
    Nphases = PhaseNames.size();

    // Reading states of matter for all phases.
    PhaseAggregateStates.resize(Nphases);
    for(size_t m = 0; m < Nphases; m++)
    {
        stringstream converter;
        converter << string("State_") << m;
        string state_of_matter = FileInterface::ReadParameterK(inp, moduleLocation, converter.str(), false, "SOLID");
        bool state_set = false;
        if(state_of_matter == "SOLID")
        {
            PhaseAggregateStates[m] = AggregateStates::Solid;
            state_set = true;
        }
        if(state_of_matter == "LIQUID" )
        {
            PhaseAggregateStates[m] = AggregateStates::Liquid;
            state_set = true;
        }
        if(state_of_matter == "GAS")
        {
            PhaseAggregateStates[m] = AggregateStates::Gas;
            state_set = true;
        }
        if(!state_set)
        {
            ConsoleOutput::WriteExit("Wrong state of matter is selected for phase " + to_string(m) + " -> " + state_of_matter, "Settings", "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }
    }

//    // Reading phase activation modes.
//    PhaseActivationModes.resize(Nphases);
//    for(size_t m = 0; m < Nphases; m++)
//    {
//        stringstream converter;
//        converter << string("PhaseActivation_") << m;
//        string state_of_matter = FileInterface::ReadParameterK(inp, moduleLocation, converter.str(), false, "ENABLED");
//        bool state_set = false;
//        if(state_of_matter == "ENABLED")
//        {
//            PhaseActivationModes[m] = ActivationModes::Enabled;
//            state_set = true;
//        }
//        if(state_of_matter == "DISABLED")
//        {
//            PhaseActivationModes[m] = ActivationModes::Disabled;
//            state_set = true;
//        }
//        if(!state_set)
//        {
//            ConsoleOutput::WriteExit("Wrong phase interaction mode is selected for phase " + to_string(m) + " -> " + state_of_matter, "Settings", "ReadInput()");
//            OP_Exit(EXIT_FAILURE);
//        }
//    }

    // Reading interaction modes for phase pairs.
//    PhaseInteractions.Allocate(Nphases, Nphases);
//    for(size_t n = 0; n < Nphases; n++)
//    for(size_t m = n; m < Nphases; m++)
//    if(PhaseActivationModes[n] == ActivationModes::Enabled and
//       PhaseActivationModes[m] == ActivationModes::Enabled)
//    {
//        stringstream converter;
//        converter << string("PhaseInteractions_") << n << "_" << m;
//        PhaseInteractions(n,m) = FileInterface::ReadParameterB(inp, moduleLocation, converter.str(), false, true);
//        PhaseInteractions(m,n) = PhaseInteractions(n,m);
//    }
//    else
//    {
//        PhaseInteractions(n,m) = false;
//        PhaseInteractions(m,n) = false;
//    }

    // Reading equilibrium density for all phases (Rigid-Body Motion).
    PhaseEquilibriumDensities.resize(Nphases);
    for(size_t m = 0; m < Nphases; m++)
    {
        stringstream converter;
        converter << string("MassDensity_") << m;
        PhaseEquilibriumDensities[m] = FileInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
    }

    // Reading number of crystallographic variants of each phase
    Nvariants.resize(Nphases);
    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "Nvariants_" << pIndex;
        Nvariants[pIndex] = FileInterface::ReadParameterI(inp, moduleLocation, converter.str(),false,1);
        if (Nvariants[pIndex] < 1)
        {
            ConsoleOutput::WriteExit("Wrong number of crystallographic variants for phase " + std::to_string(pIndex) + ": " + std::to_string(Nvariants[pIndex]) + " (Minimum 1)", "Settings", "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }
    }

    endofnames = false;
    n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Comp_") << n;
        if(FileInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = FileInterface::ReadParameterK(inp, moduleLocation, converter.str());
            /*if(tmp == "VA")
            {
                cerr << "\"VA\" as a component is not allowed in "
                     << "ChemicalProperties.opi" << endl;
                OP_Exit(EXIT_SUCCESS);
            }*/
            ElementNames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }
    Ncomp = ElementNames.size();

    VTKDir     = FileInterface::ReadParameterF(inp, moduleLocation,"VTKDir",false,VTKDir);
    RawDataDir = FileInterface::ReadParameterF(inp, moduleLocation,"RAWDir",false,RawDataDir);
    InputRawDataDir = FileInterface::ReadParameterF(inp, moduleLocation,"InputRAWDir",false,RawDataDir);
    TextDir    = FileInterface::ReadParameterF(inp, moduleLocation,"DATADir",false,TextDir);

    string from = "\\";
    string to = "/";

#ifdef _WIN32
    from = "/";
    to = "\\";
#endif

    size_t start_pos = 0;
    while((start_pos = VTKDir.find(from, start_pos)) != std::string::npos)
    {
        VTKDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    start_pos = 0;
    while((start_pos = RawDataDir.find(from, start_pos)) != std::string::npos)
    {
        RawDataDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    start_pos = 0;
    while((start_pos = TextDir.find(from, start_pos)) != std::string::npos)
    {
        TextDir.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    if (!VTKDir.empty())
    {
        stringstream ss; ss << *VTKDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) VTKDir += dirSeparator;
    }

    if (!RawDataDir.empty())
    {
        stringstream ss; ss << *RawDataDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) RawDataDir += dirSeparator;
    }

    if (!TextDir.empty())
    {
        stringstream ss; ss << *TextDir.rbegin();
        string s; ss >> s;
        if(s != dirSeparator) TextDir += dirSeparator;
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Settings::ReadJSON(const string InputFileName)
{
    Grid.ReadJSON(InputFileName);
    std::ifstream f(InputFileName);
    
    json data = json::parse(f);
    if (data.contains(thisclassname))
    {
        json settings = data[thisclassname];

        ConsoleOutput::WriteLine();
        ConsoleOutput::WriteLineInsert("Settings");
        ConsoleOutput::WriteLineInsert("Active phases");

        size_t number_of_phases = FileInterface::ReadParameter<size_t>(settings, {"Phases","Number"});
        for (size_t i = 0; i < number_of_phases; ++i)
        {
            string tmp = FileInterface::ReadParameter<std::string>(settings, {"Phases",i,"Name"});
            PhaseNames.push_back(tmp);
        }
        Nphases = PhaseNames.size();

        // Reading states of matter for all phases.
        PhaseAggregateStates.resize(Nphases);
        for(size_t m = 0; m < Nphases; m++)
        {

            string state_of_matter = FileInterface::ReadParameter<std::string>(settings, {"Phases",m,"State"},"SOLID");
            bool state_set = false;
            if(state_of_matter == "SOLID")
            {
                PhaseAggregateStates[m] = AggregateStates::Solid;
                state_set = true;
            }
            if(state_of_matter == "LIQUID" )
            {
                PhaseAggregateStates[m] = AggregateStates::Liquid;
                state_set = true;
            }
            if(state_of_matter == "GAS")
            {
                PhaseAggregateStates[m] = AggregateStates::Gas;
                state_set = true;
            }
            if(!state_set)
            {
                ConsoleOutput::WriteExit("Wrong state of matter is selected for phase " + to_string(m) + " -> " + state_of_matter, "Settings", "ReadInput()");
                OP_Exit(EXIT_FAILURE);
            }
        }

        // Reading equilibrium density for all phases (Rigid-Body Motion).
        PhaseEquilibriumDensities.resize(Nphases);
        for(size_t m = 0; m < Nphases; m++)
        {
            PhaseEquilibriumDensities[m] = FileInterface::ReadParameter<double>(settings, {"Phases",m,"EquilibriumDensity"}, 0.0);
        }

        // Reading number of crystallographic variants of each phase
        Nvariants.resize(Nphases);
        for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            Nvariants[pIndex] = FileInterface::ReadParameter<int>(settings, {"Phases",pIndex,"Variants"}, 1);
            if (Nvariants[pIndex] < 1)
            {
                ConsoleOutput::WriteExit("Wrong number of crystallographic variants for phase " + std::to_string(pIndex) + ": " + std::to_string(Nvariants[pIndex]) + " (Minimum 1)", "Settings", "ReadInput()");
                OP_Exit(EXIT_FAILURE);
            }
        }

        size_t number_of_components = FileInterface::ReadParameter<size_t>(settings, {"Components","Number"},0);
        for (size_t i = 0; i < number_of_components; ++i)
        {
            string tmp = FileInterface::ReadParameter<std::string>(settings, {"Components",i,"Name"});
            PhaseNames.push_back(tmp);
        }
        Ncomp = ElementNames.size();

        VTKDir     = FileInterface::ReadParameter<std::string>(settings, {"VTKDir"}, VTKDir);
        RawDataDir = FileInterface::ReadParameter<std::string>(settings, {"RawDataDir"}, RawDataDir);
        InputRawDataDir = FileInterface::ReadParameter<std::string>(settings, {"InputRawDataDir"}, InputRawDataDir);
        TextDir    = FileInterface::ReadParameter<std::string>(settings, {"TextDir"}, TextDir);

        string from = "\\";
        string to = "/";

#ifdef _WIN32
        from = "/";
        to = "\\";
#endif

        size_t start_pos = 0;
        while((start_pos = VTKDir.find(from, start_pos)) != std::string::npos)
        {
            VTKDir.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }

        start_pos = 0;
        while((start_pos = RawDataDir.find(from, start_pos)) != std::string::npos)
        {
            RawDataDir.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }

        start_pos = 0;
        while((start_pos = TextDir.find(from, start_pos)) != std::string::npos)
        {
            TextDir.replace(start_pos, from.length(), to);
            start_pos += to.length();
        }

        if (!VTKDir.empty())
        {
            stringstream ss; ss << *VTKDir.rbegin();
            string s; ss >> s;
            if(s != dirSeparator) VTKDir += dirSeparator;
        }

        if (!RawDataDir.empty())
        {
            stringstream ss; ss << *RawDataDir.rbegin();
            string s; ss >> s;
            if(s != dirSeparator) RawDataDir += dirSeparator;
        }

        if (!TextDir.empty())
        {
            stringstream ss; ss << *TextDir.rbegin();
            string s; ss >> s;
            if(s != dirSeparator) TextDir += dirSeparator;
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Settings::RemeshAll(int newNx, int newNy, int newNz, int tStep, const BoundaryConditions& BC)
{
    if(RemeshingAllowed)
    {
        if(!Grid.dNx) newNx = 0;
        if(!Grid.dNy) newNy = 0;
        if(!Grid.dNz) newNz = 0;

        for(size_t n = 0; n < ObjectsToRemesh.size(); n++)
        {
            ObjectsToRemesh[n]->Remesh(newNx, newNy, newNz, BC);
        }

        Grid.SetDimensions(newNx, newNy, newNz);

        Grid.Idx = tStep;

        GridHistory.push_back(Grid);
    }
}

void Settings::AddForRemeshing(OPObject& Obj)
{
    ObjectsToRemesh.push_back(&Obj);
}

void Settings::AddForAdvection(OPObject& Obj)
{
    ObjectsToAdvect.push_back(&Obj);
}

void Settings::AdvectAll(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, double dt, int tStep)
{
    for(size_t n = 0; n < ObjectsToAdvect.size(); n++)
    if(ObjectsToAdvect[n]->thisclassname != "PhaseField")
    {
        ObjectsToAdvect[n]->Advect(Adv, Vel, Phi, BC, dt, tStep);
    }
    else
    {
        Phi.AdvectALE(Adv, Vel, BC, dt, tStep);
    }
}

void Settings::Advect(std::string ObjNameBase, AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, double dt, int tStep)
{
    for(size_t n = 0; n < ObjectsToAdvect.size(); n++)
    if(ObjectsToAdvect[n]->thisobjectname == ObjNameBase)
    {
        ObjectsToAdvect[n]->Advect(Adv, Vel, Phi, BC, dt, tStep);
    }
}

void Settings::AddForReading(OPObject& Obj)
{
    ObjectsToRead.push_back(&Obj);
}

bool Settings::ReadAll(const BoundaryConditions& BC, const int tStep)
{
    bool read_status = true;
    for(size_t n = 0; n < ObjectsToRead.size(); n++)
    {
        read_status = read_status && ObjectsToRead[n]->Read(*this,BC, tStep);
    }
    return read_status;
}

bool Settings::Read(std::string ObjNameBase, const BoundaryConditions& BC, const int tStep)
{
    bool read_status = true;
    for(size_t n = 0; n < ObjectsToRead.size(); n++)
    if(ObjectsToRead[n]->thisobjectname == ObjNameBase)
    {
        read_status = read_status && ObjectsToRead[n]->Read(*this,BC, tStep);
    }
    return read_status;
}

Settings& Settings::operator= (const Settings& rhs)
{
    if (this != &rhs) // protect against self-assignment
    {
        thisclassname  = rhs.thisclassname;
        thisobjectname = rhs.thisobjectname;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;
        Ncomp   = rhs.Ncomp;

        RemeshingAllowed = rhs.RemeshingAllowed;
        initialized = rhs.initialized;

        ObjectsToRemesh = rhs.ObjectsToRemesh;
        ObjectsToAdvect = rhs.ObjectsToAdvect;
        ObjectsToRead = rhs.ObjectsToRead;

        PhaseNames = rhs.PhaseNames;
        ElementNames = rhs.ElementNames;

        PhaseAggregateStates = rhs.PhaseAggregateStates;
        PhaseEquilibriumDensities = rhs.PhaseEquilibriumDensities;

        //PhaseActivationModes = rhs.PhaseActivationModes;

        VTKDir = rhs.VTKDir;
        RawDataDir = rhs.RawDataDir;
        InputRawDataDir = rhs.InputRawDataDir;
        TextDir = rhs.TextDir;

        GridHistory = rhs.GridHistory;
    }
    return *this;
}

GridParameters Settings::GridHistoryParameters(int time_step) const
{
    if(time_step < 0)
    {
        time_step = 0;
    }

    size_t history_record_idx = 0;
    for(size_t n = 0; n < GridHistory.size(); n++)
    if(n <= size_t(time_step))
    {
        history_record_idx = n;
    }
    else
    {
        break;
    }
    return GridHistory[history_record_idx];
}

} //namespace openphase
