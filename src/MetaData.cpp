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

#include "MetaData.h"
#include "json.hpp"
#include "BuildInfo.h"
#ifdef __linux__
#include <unistd.h>         // readlink
#include <linux/limits.h>   // PATH_MAX
#endif
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace openphase
{
    using json = nlohmann::json;
    
    std::string MetaData::to_iso8601(std::time_t t) {
        std::tm tm = *std::gmtime(&t); // use gmtime for UTC ("Z"), localtime for local offset
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
        return oss.str();
    }
  
std::string MetaData::make_relative_path(const std::string& path_to_file,
                               const std::string& path_to_working_directory)
{
    try {
        fs::path file_path(path_to_file);
        fs::path base_path(path_to_working_directory);

        // Ensure both paths are absolute for correct relative calculation
        file_path = fs::absolute(file_path);
        base_path = fs::absolute(base_path);

        fs::path relative = fs::relative(file_path, base_path);
        return "./"+relative.generic_string();  // Use forward slashes
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return path_to_file; // fallback to original
    }
}

    MetaData::MetaData()
    {
        time_t timestamp;
        time(&timestamp);
        Program.Title = "OPExe"; 
        Program.GitHashCore = GIT_COMMIT_SHA;               
        StartTime = to_iso8601(timestamp);
        #ifdef _WIN32
		    char exepath[MAX_PATH];
		    std::string s_exepath(exepath);
		    HMODULE hModule = GetModuleHandle(NULL);
		    if (hModule != NULL)
		    {
		        // Use GetModuleFileName() with module handle to get the path
		        GetModuleFileName(hModule, exepath, (sizeof(exepath)));
                s_exepath = exepath;
		    }
		    Program.Title = s_exepath; 
        #else
		    std::string result;
		    result.resize(PATH_MAX);
			ssize_t count = readlink("/proc/self/exe", result.data(), PATH_MAX);
			result.resize(count);
			Program.Title = result; 
        #endif                           
    }

    void MetaData::WriteJSON(std::string FileName)
    {
        time_t timestamp;
        time(&timestamp);
        fs::path p(FileName);
        fs::path dir = p.parent_path();
        json j;
        j["@context"][0] = "https://w3id.org/ro/crate/1.1/context";
        j["@context"][1] = "https://w3id.org/ro/terms/workflow-run/context";
        j["@graph"][0]["@id"] = "./ro-crate-metadata.json";
        j["@graph"][0]["@type"] = "CreativeWork";
        j["@graph"][0]["conformsTo"]["@id"] = "https://w3id.org/ro/crate/1.1";
        j["@graph"][0]["about"]["@id"] = "./";
        j["@graph"][1]["@id"] = "./";
        j["@graph"][1]["@type"] = "Dataset";
        j["@graph"][1]["name"] = Simulation.Title;
        j["@graph"][1]["description"] = Simulation.Description;
        j["@graph"][1]["license"] = Simulation.License;
        j["@graph"][1]["datePublished"] = to_iso8601(timestamp);
        j["@graph"][1]["author"] = json::array();
        j["@graph"][1]["conformsTo"]["@id"] = "https://w3id.org/ro/wfrun/process/0.4";
        j["@graph"][1]["mentions"]["@id"] = "Workflow";
        j["@graph"][2]["@id"] = "Workflow"; 
        j["@graph"][2]["@type"] = "WorkflowRun"; 
        j["@graph"][2]["name"] = Simulation.Title;
        j["@graph"][2]["description"] = Simulation.Description;
        j["@graph"][2]["startTime"] = StartTime;
        j["@graph"][2]["endTime"] = to_iso8601(timestamp);
        j["@graph"][2]["instrument"]["@id"]  = Program.Title; 
        json jApp;
        jApp["@id"] = Program.Title;
        jApp["@type"] = "Application";
        jApp["softwareVersion"] = Program.VersionStudio;
        jApp["identifier"] = Program.GitHashStudio;
        jApp["dateCreated"] = BUILD_TIME;
        j["@graph"].push_back(jApp);
        j["@graph"][0]["datePublished"] = to_iso8601(timestamp);
        j["@graph"][0]["description"] = Simulation.Description;
        j["@graph"][0]["hasPart"];
        for (size_t i = 0; i < Parts.size(); ++i)
        {
            j["@graph"][0]["hasPart"][i]["@id"] = make_relative_path(Parts[i].Name,dir.string());
            j["@graph"][1]["result"][i]["@id"] = make_relative_path(Parts[i].Name,dir.string());
            json j1;
            j1["@id"] = make_relative_path(Parts[i].Name,dir.string());
            j1["@type"] = Parts[i].Type;
            j1["description"] = Parts[i].Description;
            j1["encodingFormat"] = Parts[i].Format;
            j["@graph"].push_back(j1);
        }
        for (size_t i = 0; i < Authors.size(); ++i)
        {
            j["@graph"][0]["author"][i]["@id"] = Authors[i].Name;
            json j1;
            j1["@id"] = Authors[i].Name;
            j1["@type"] = "Person";
            j1["name"] = Authors[i].Name;
            j1["affiliation"] = Authors[i].Affiliation;
            j1["identifier"] = Authors[i].ORCID;
            j["@graph"].push_back(j1);
        }
        for (size_t i = 0; i < Simulation.Keywords.size(); ++i)
        {
            j["@graph"][0]["keywords"][i] = Simulation.Keywords[i];
        }
        std::cout << "FileName: " << FileName << std::endl;
        std::ofstream o(FileName);
        o << std::setw(4) << j << std::endl;
        o.close();
    }
    
    void MetaData::AddPart(std::string FileName, std::string type, std::string description, std::string format)
    {
        Part p;
        p.Name = FileName;
        p.Type = type;
        p.Description = description;
        p.Format = format;
        Parts.push_back(p);
    }
    
    void MetaData::ReadAuthorInformation(std::stringstream& inp)
    {
        std::string module = " ";
        int moduleLocation = FileInterface::FindModuleLocation(inp, "AuthorInformation");
        if (moduleLocation != 0)
        {
            int i = 0;
            do {
                module = "Author_" + std::to_string(i);
                std::string tmpMd = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"NOAUTHOR");          
                if (tmpMd != "NOAUTHOR")
                {
                    AuthorInformation NewAuthor;
                    NewAuthor.Name  = tmpMd;
                    module = "Affiliation_" + std::to_string(i);
                    tmpMd = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"");
                    NewAuthor.Affiliation = tmpMd;
                    module = "ORCID_" + std::to_string(i);
                    tmpMd = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"");
                    NewAuthor.ORCID = tmpMd;
                    Authors.push_back(NewAuthor);
                    ++i;
                }
                else
                {
                    break;
                }
            }            
            while(true);
        }
    }
    
    void MetaData::ReadSimulationInformation(std::stringstream& inp)
    {
        std::string module = " ";
        int moduleLocation = FileInterface::FindModuleLocation(inp, "SimulationInformation");
        if (moduleLocation != 0)
        {
            module = "Title";
            Simulation.Title = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"NOTITLE");      
            module = "Description";
            Simulation.Description = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"NODESCRIPTION");     
            module = "GUIGitHash";
            Program.GitHashStudio = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"0000");      
            module = "GUIVersion";
            Program.VersionStudio = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"0000");   
            
            int i = 0;
            do {
                module = "Keyword_" + std::to_string(i);
                std::string tmpMd  = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"NOKEYWORD");          
                if (tmpMd != "NOKEYWORD")
                {
                    Simulation.Keywords.push_back(tmpMd);
                    ++i;
                }
                else
                {
                    break;
                }
            }            
            while(true);
            
            i = 0;
            do {
                module = "Article_" + std::to_string(i);
                std::string tmpMd  = FileInterface::ReadParameterS(inp,moduleLocation,module,false,"NOARTICLE");          
                if (tmpMd != "NOARTICLE")
                {
                    Simulation.RelatedArticles[i] = tmpMd;
                    ++i;
                }
                else
                {
                    break;
                }
            }            
            while(true);
        }
        
    }

}// namespace openphase

