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
 *  File created :   2025
 *  Main contributors :  Oleg Shchyglo;
 *
 */

#include "ThermodynamicInterfaces/ThermodynamicInterfaceLPD.h"
#include "Settings.h"
#include "Composition.h"
#include "Temperature.h"
#include "PhysicalConstants.h"
#include "Containers/StringManipulations.h"
#include "Thermodynamics/ThermodynamicPropertiesEQP.h"
#include <filesystem>

namespace openphase
{

using namespace std;

void ThermodynamicInterfaceLPD::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ThermodynamicInterface";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;
    Cnom.resize(Ncomp, 0.0);

    PhaseNames = locSettings.PhaseNames;
    ElementNames = locSettings.ElementNames;

    Data.Allocate(Nphases,Nphases);
    DataShortcuts.Allocate(Nphases,Nphases);
    
    BaseDir = locSettings.BaseDir;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ThermodynamicInterfaceLPD::ReadInput(const string InputFileName)
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

void ThermodynamicInterfaceLPD::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    DBfile = FileInterface::ReadParameterFW(inp, moduleLocation, {std::string("TQfile"),std::string("DBfile")}, true, "");

	std::filesystem::path DB(DBfile);
	
	if (DB.is_relative())
        DB = BaseDir / DB;

    DB = std::filesystem::weakly_canonical(DB);

	DBfile = DB.string();
    ReadTDB();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void ThermodynamicInterfaceLPD::ReadTDB()
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

void ThermodynamicInterfaceLPD::ReadTDB(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Thermodynamic Properties LPD");
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

    // Read nominal composition from the database
    for(size_t n = 0; n < ElementNames.size(); n++)
    {
        Cnom[n] = FileInterface::ReadParameterD(inp, moduleLocation, string("Cnom_") + ElementNames[n], true, 0.0);
    }

    // Read database entries for phase pairs
    for(size_t m =   0; m < PhaseNames.size()-1; m++)
    for(size_t n = m+1; n < PhaseNames.size(); n++)
    {
        int locationA = FileInterface::FindParameterLocation(inp, moduleLocation, "EQ_" + PhaseNames[m] + "_" + PhaseNames[n]);
        int locationB = FileInterface::FindParameterLocation(inp, moduleLocation, "EQ_" + PhaseNames[n] + "_" + PhaseNames[m]);

        if(locationA == -1 and locationB == -1)
        {
            ConsoleOutput::WriteExit("Equilibrium entries for phase pair " + PhaseNames[m] + " and " + PhaseNames[n] + " are not found in database file " + DBfile, thisclassname, "ReadTDB");
            OP_Exit(EXIT_FAILURE);
        }

        size_t idxA = m;
        size_t idxB = n;

        if(locationA != -1)
        {
            inp.seekg(locationA);
        }

        if(locationB != -1)
        {
            inp.seekg(locationB);
            idxA = n;
            idxB = m;
        }

        int  entry_counter = 0;
        bool valid_input = true;
        bool empty_line  = false;
        bool random_text = false;

        while(valid_input or empty_line or random_text)
        {
            valid_input = false;
            empty_line  = false;
            random_text = false;

            // Get entire line
            string tmp;
            //cout << "Raw string: " << tmp << endl;

            if (!getline(inp,tmp)) break;

            // remove leading/training whitespaces
            RemoveLeadingTrailingWhitespaces(tmp);

            // Exit if end of data set is reached
            if(StringIsSectionName(tmp) or StringIsInputParameter(tmp))
            {
                //cout << "Next section/input parameter: " << tmp << endl;
                break;
            }

            // Skip empty lines and arbitrary text (or comments)
            if(StringIsEmpty(tmp))
            {
                empty_line = true;
                //cout << "Empty line: " << tmp << endl;
                continue;
            }

            EquilibriumData tmpDataA(Ncomp);
            EquilibriumData tmpDataB(Ncomp);
            // Extract string between brackets, continue otherwise
            if(StringIsInBrackets(tmp, "</","/>"))
            {
                //cout << "Valid input: " << tmp << endl;

                tmp = ExtractStringBetweenBrackets(tmp, "</","/>");
                vector<string> tmp2 = ExtractSubstrings(tmp);
                //cout << "Substrings: " << endl;

                //for(size_t idx = 0; idx < tmp2.size(); idx++)
                //{
                //    cout << tmp2[idx] << endl;
                //}

                if(tmp2.size() == Ncomp*2 + 3)
                {
                    valid_input = true;

                    try
                    {
                        tmpDataA.eq_temperature = std::stod(tmp2[0]);
                        //cout << "Teq: " << tmpDataA.eq_temperature << endl;
                    }
                    catch (const std::invalid_argument&)
                    {
                       ConsoleOutput::WriteExit("Argument for temperature is invalid", thisclassname, "ReadTDB()");
                       OP_Exit(EXIT_FAILURE);
                    }
                    tmpDataB.eq_temperature = tmpDataA.eq_temperature;
                    //getchar();
                    try
                    {
                        tmpDataA.eq_entropy = std::stod(tmp2[1]);
                        //cout << "Ent 1: " << tmpDataA.eq_entropy << endl;
                    }
                    catch (const std::invalid_argument&)
                    {
                       ConsoleOutput::WriteExit("Argument for entropy of phase " + PhaseNames[m] + " is invalid", thisclassname, "ReadTDB()");
                       OP_Exit(EXIT_FAILURE);
                    }

                    try
                    {
                        tmpDataB.eq_entropy = std::stod(tmp2[2]);
                        //cout << "Ent 2: " << tmpDataB.eq_entropy << endl;
                    }
                    catch (const std::invalid_argument&)
                    {
                       ConsoleOutput::WriteExit("Argument for entropy of phase " + PhaseNames[n] + " is invalid", thisclassname, "ReadTDB()");
                       OP_Exit(EXIT_FAILURE);
                    }
                    //getchar();

                    for(size_t i = 3, ii = 0; i < Ncomp+3; i++, ii++)
                    {
                        //cout << tmp2[i] << endl;
                        try
                        {
                            tmpDataA.eq_composition({ii}) = std::stod(tmp2[i]);
                            //cout << "CeqA: " << ii << " " << tmpDataA.eq_composition({ii}) << endl;
                        }
                        catch (const std::invalid_argument&)
                        {
                           ConsoleOutput::WriteExit("Entry " + tmp2[i] + " for equilibrium composition is invalid", thisclassname, "ReadTDB()");
                           OP_Exit(EXIT_FAILURE);
                        }

                    }
                    //getchar();

                    for(size_t i = Ncomp+3, ii = 0; i < 2*Ncomp+3; i++, ii++)
                    {
                        //cout << tmp2[i] << endl;
                        try
                        {
                            tmpDataB.eq_composition({ii}) = std::stod(tmp2[i]);
                            //cout << "CeqB: " << ii << " " << tmpDataB.eq_composition({ii}) << endl;
                        }
                        catch (const std::invalid_argument&)
                        {
                           ConsoleOutput::WriteExit("Entry " + tmp2[i] + " for equilibrium composition is invalid", thisclassname, "ReadTDB()");
                           OP_Exit(EXIT_FAILURE);
                        }
                    }
                    //getchar();

                    Data(idxA,idxB).push_back(tmpDataA);
                    Data(idxB,idxA).push_back(tmpDataB);

                    entry_counter++;
                }
            }
            else
            {
                //cout << "Random text: " << tmp << endl;
                random_text = true;
                continue;
            }
        }
        if(entry_counter < 2)
        {
            ConsoleOutput::WriteExit("Insufficient number of database entries (less than two) for phase pair " + PhaseNames[m] + " and " + PhaseNames[n], thisclassname, "ReadTDB()");
            OP_Exit(EXIT_FAILURE);
        }
    }
    //getchar();
    SortDatabase();

    ConsoleOutput::WriteLine();
}

void ThermodynamicInterfaceLPD::SortDatabase(void)
{
    for(size_t m = 0; m < PhaseNames.size(); m++)
    for(size_t n = 0; n < PhaseNames.size(); n++)
    if(m != n)
    {
        // Sort entries for a given phase pair by ascending equilibrium temperature
        for(size_t i =   0; i < Data(m,n).size()-1; i++)
        for(size_t j = i+1; j < Data(m,n).size()  ; j++)
        {
            if(Data(m,n)[i].eq_temperature > Data(m,n)[j].eq_temperature)
            {
                EquilibriumData tmpEntry1 = Data(m,n)[j];
                Data(m,n)[j] = Data(m,n)[i];
                Data(m,n)[i] = tmpEntry1;
            }
        }
        // Set lowest and highest available temperatures
        double Tmin = Data(m,n).front().eq_temperature;
        double Tmax = Data(m,n).back().eq_temperature;

        // Create temperature shortcuts for maximum temperature range from 0K to 6000K
        DataShortcuts(m,n).resize(6001);

        // Use lower temperature of the range as the shortcut position.
        // The nearest range is always between the shortcut temperature
        // and the next higher temperature.
        for(size_t i = 0; i < floor(Tmin); i++)
        {
            DataShortcuts(m,n)[i] = 0;
        }

        for(size_t i = ceil(Tmax); i < 6001; i++)
        {
            DataShortcuts(m,n)[i] = Data(m,n).size()-2;
        }

        for(size_t i = floor(Tmin); i < ceil(Tmax); i++)
        for(size_t j = 1; j < Data(m,n).size() - 1;)
        {
            if(i < Data(m,n)[j].eq_temperature)
            {
                DataShortcuts(m,n)[i] = j-1;
                break;
            }
            else
            {
                DataShortcuts(m,n)[i] = j;
                j++;
            }
        }
        /*for(size_t i = 0; i < 6001; i++)
        {
            cout << i << ": " << DataShortcuts(m,n)[i] << endl;
            if(!(i%50)) getchar();
        }*/
    }
}

void ThermodynamicInterfaceLPD::SetEquilibriumData(PhaseField& Phase,
                                                       Composition& Cx,
                                                       Temperature& Tx,
                                                       ThermodynamicPropertiesEQP& TP)
{
    /** This function will populate thermodynamic equilibrium data from the
    given phase diagram with the given parameters for T, P and X for all active
    phase pairs. Output is the entropy difference, phase diagram
    slopes, equilibrium compositions for all active phase pairs and
    equilibrium volume fractions.*/

    switch(TP.ExtrapolationMode)
    {
        case ExtrapolationModes::Global:
        {
            for(size_t m = 0; m < Nphases; m++)
            for(size_t n = 0; n < Nphases; n++)
            if(m != n)
            if( TP.GlobalExtrapolationData({m,n}).out_of_T_range or
               !TP.GlobalExtrapolationData({m,n}).is_set)
            {
                size_t idx = DataShortcuts(m,n)[(int)Tx.Tavg];

                TP.GlobalExtrapolationData({m,n}) = Data(m,n)[idx];
                TP.GlobalExtrapolationData({m,n}).is_set = true;
                TP.GlobalExtrapolationData({m,n}).out_of_T_range = false;

                double dT = (Data(m,n)[idx+1].eq_temperature - Data(m,n)[idx].eq_temperature);
                for(size_t i = 0; i < Ncomp; i++)
                {
                    double dC = Data(m,n)[idx+1].eq_composition({i}) - Data(m,n)[idx].eq_composition({i});

                    if(dC != 0.0)
                    {
                        TP.GlobalExtrapolationData({m,n}).eq_slope({i}) = dT/dC;
                    }
                    else
                    {
                        TP.GlobalExtrapolationData({m,n}).eq_slope({i}) = 0.0;
                    }
                }
            }

            for(size_t m = 0; m < Nphases; m++)
            for(size_t n = 0; n < Nphases; n++)
            if(m != n)
            {
                for(size_t i = 0; i < Ncomp; i++)
                if(TP.GlobalExtrapolationData({m,n}).eq_slope({i}) != 0.0 and
                   TP.GlobalExtrapolationData({n,m}).eq_slope({i}) != 0.0)
                {
                    TP.GlobalExtrapolationData({m,n}).part_coefficient({i}) =
                            fabs(TP.GlobalExtrapolationData({m,n}).eq_slope({i})
                                /TP.GlobalExtrapolationData({n,m}).eq_slope({i}));
                }
                else
                {
                    TP.GlobalExtrapolationData({m,n}).part_coefficient({i}) = 1.0;
                }
            }

            break;
        }
        case ExtrapolationModes::Local:
        case ExtrapolationModes::None:
        default:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,TP.LocalExtrapolationData,TP.LocalExtrapolationData.Bcells(),)
            {
                if(Phase.Fields(x,y,z).wide_interface())
                {
                    for(size_t m = 0; m < Nphases; m++)
                    for(size_t n = 0; n < Nphases; n++)
                    if(m != n and (TP.LocalExtrapolationData(x,y,z,{m,n}).out_of_T_range or
                                  !TP.LocalExtrapolationData(x,y,z,{m,n}).is_set))
                    {
                        size_t idx = DataShortcuts(m,n)[(int)Tx.Tx(x,y,z)];

                        TP.LocalExtrapolationData(x,y,z,{m,n}) = Data(m,n)[idx];
                        TP.LocalExtrapolationData(x,y,z,{m,n}).is_set = true;
                        TP.LocalExtrapolationData(x,y,z,{m,n}).out_of_T_range = false;

                        double dT = (Data(m,n)[idx+1].eq_temperature - Data(m,n)[idx].eq_temperature);

                        for(size_t i = 0; i < Ncomp; i++)
                        {
                            double dC = Data(m,n)[idx+1].eq_composition({i}) - Data(m,n)[idx].eq_composition({i});

                            if(dC != 0.0)
                            {
                                TP.LocalExtrapolationData(x,y,z,{m,n}).eq_slope({i}) = dT/dC;
                            }
                            else
                            {
                                TP.LocalExtrapolationData(x,y,z,{m,n}).eq_slope({i}) = 0.0;
                            }
                        }
                    }

                    for(size_t m = 0; m < Nphases; m++)
                    for(size_t n = 0; n < Nphases; n++)
                    if(m != n)
                    {
                        for(size_t i = 0; i < Ncomp; i++)
                        if(TP.LocalExtrapolationData(x,y,z,{m,n}).eq_slope({i}) != 0.0 and
                           TP.LocalExtrapolationData(x,y,z,{n,m}).eq_slope({i}) != 0.0)
                        {
                            TP.LocalExtrapolationData(x,y,z,{m,n}).part_coefficient({i}) =
                                    fabs(TP.LocalExtrapolationData(x,y,z,{m,n}).eq_slope({i})
                                        /TP.LocalExtrapolationData(x,y,z,{n,m}).eq_slope({i}));
                        }
                        else
                        {
                            TP.LocalExtrapolationData(x,y,z,{m,n}).part_coefficient({i}) = 1.0;
                        }
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            break;
        }
    }
}

}//namespace openphase

