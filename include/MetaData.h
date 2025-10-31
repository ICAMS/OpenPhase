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

 *   File created :   2024
 *   Main contributors :   Marvin Tegeler
 *
 */

#ifndef METADATA_H
#define METADATA_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "Globals.h"
#include "Includes.h"

namespace openphase
{

struct AuthorInformation                                                    
{
    std::string Name;
    std::string Affiliation;
    std::string ORCID;
};

struct SimulationInformation                                                    
{
    std::string Title;
    std::string Description;
    std::string License;
    std::vector<std::string> Keywords;
    std::vector<std::string> RelatedArticles;
};

struct ProgramInformation                                                  
{
    std::string Title;
    std::string GitHashCore;
    std::string VersionCore;
    std::string GitHashStudio;
    std::string VersionStudio;
};

struct Part                                                    
{
    std::string Name;
    std::string Type;
    std::string Description;
    std::string Format;
};

class OP_EXPORTS MetaData 
{
  public:

    MetaData();
    //MetaData(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);
    //void Initialize(Settings& locSettings);                            ///< Initializes storages, sets internal variables.
    //void ReadInput(const std::string InputFileName);                   ///< Reads input values from file
    //void ReadInput(std::stringstream& inp);                            ///< Reads input values from file
    std::string to_iso8601(std::time_t t);
    std::string make_relative_path(const std::string& path_to_file,
                               const std::string& path_to_working_directory);
    void WriteJSON(std::string FileName);
    void AddPart(std::string FileName, std::string type, std::string description, std::string format);
    void ReadAuthorInformation(std::stringstream& inp);
    void ReadSimulationInformation(std::stringstream&);
    std::vector<AuthorInformation> Authors;
    SimulationInformation Simulation;
    ProgramInformation Program;
    std::vector<Part> Parts;
    std::string StartTime;
};

}// namespace openphase

#endif // MetaData_H
