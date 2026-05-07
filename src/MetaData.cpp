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
 *  File created :   2024
 *  Main contributors :   Marvin Tegeler
 *
 */

#include "MetaData.h"
#include "Settings.h" 
#include "RunTimeControl.h" 
#include "json.hpp"
#include "BuildInfo.h"
#ifdef __linux__
#include <unistd.h>         // readlink
#include <linux/limits.h>   // PATH_MAX
#endif
#include <filesystem>
#include <string>
#include <unordered_set>
#ifdef OP_HAS_OPENSSL
#include <openssl/sha.h>
#include <openssl/evp.h>
#include <iomanip>
#endif

namespace fs = std::filesystem;

namespace openphase
{
    using json = nlohmann::json;
    
    inline std::string MetaData::trim(const std::string& str)
	{
		size_t first = str.find_first_not_of(" \t\n\r");
		if (first == std::string::npos) return ""; // all spaces
		size_t last = str.find_last_not_of(" \t\n\r");
		return str.substr(first, last - first + 1);
	}
    
    std::string MetaData::to_iso8601(std::time_t t) {
        std::tm tm = *std::gmtime(&t); // use gmtime for UTC ("Z"), localtime for local offset
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
        return oss.str();
    }
    
	void MetaData::IndexOutputFiles()
	{
		// Keep track of absolute paths we've already added
		std::unordered_set<std::string> seen_files;

		// Lambda to determine type and MIME format
		auto determine_type_and_format = [](const std::string& path) -> std::pair<std::string,std::string> {
		    std::string ext = fs::path(path).extension().string();
		    if(ext == ".vts" || ext == ".vtk") return {"File","jApp/xml"};
		    else if(ext == ".opi" || ext == ".txt" || ext == ".csv") return {"File","text/plain"};
		    else return {"File","jApp/octet-stream"};
		};

		// Lambda to generate a human-readable description
		auto generate_description = [](const std::string& path) -> std::string {
		    std::string name = fs::path(path).filename().string();
		    size_t pos = name.find_last_of('_');
		    size_t dot = name.find_last_of('.');
		    if(pos != std::string::npos && dot != std::string::npos && pos < dot)
		    {
		        std::string field = name.substr(0,pos);
		        std::string step  = name.substr(pos+1, dot-pos-1);
		        return field + " field output at step " + step;
		    }
		    return "Automatically indexed file: " + name;
		};

		// Directories to scan
		std::vector<std::string> dirs = {OutPutDir, VTKDir, RawDataDir, TextDir};

		for(const auto& dirPath : dirs)
		{
		    if(dirPath.empty()) continue;

		    fs::path dir(dirPath);

		    if(fs::exists(dir) && fs::is_directory(dir))
		    {
		        for(const auto& entry : fs::recursive_directory_iterator(dir))
		        {
		            if(!fs::is_regular_file(entry)) continue;

		            std::string absPath = fs::absolute(entry.path()).string();

		            // Skip duplicates
		            if(seen_files.count(absPath)) continue;
		            seen_files.insert(absPath);

		            auto [type, format] = determine_type_and_format(absPath);

		            Part p;
		            p.Name = absPath;
		            p.Type = type;
		            p.Format = format;
		            p.Description = generate_description(absPath);
		            p.Checksum = ComputeSHA256(absPath);
		            Parts.push_back(p);
		        }
		    }
		}
	}
	
	std::string MetaData::ComputeSHA256(const std::string& path)
	{
	#ifdef OP_HAS_OPENSSL

		std::ifstream file(path, std::ios::binary);
		if (!file) return "";

		EVP_MD_CTX* ctx = EVP_MD_CTX_new();
		if (!ctx) return "";

		if (EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) != 1)
		{
		    EVP_MD_CTX_free(ctx);
		    return "";
		}

		char buffer[8192];
		while (file.read(buffer, sizeof(buffer)))
		{
		    if (EVP_DigestUpdate(ctx, buffer, file.gcount()) != 1)
		    {
		        EVP_MD_CTX_free(ctx);
		        return "";
		    }
		}

		// Handle remaining bytes after last read
		if (file.gcount() > 0)
		{
		    if (EVP_DigestUpdate(ctx, buffer, file.gcount()) != 1)
		    {
		        EVP_MD_CTX_free(ctx);
		        return "";
		    }
		}

		unsigned char hash[EVP_MAX_MD_SIZE];
		unsigned int hash_len = 0;

		if (EVP_DigestFinal_ex(ctx, hash, &hash_len) != 1)
		{
		    EVP_MD_CTX_free(ctx);
		    return "";
		}

		EVP_MD_CTX_free(ctx);

		std::stringstream ss;
		for (unsigned int i = 0; i < hash_len; ++i)
		{
		    ss << std::hex << std::setw(2) << std::setfill('0')
		       << static_cast<int>(hash[i]);
		}

		return ss.str();

	#else
		(void)path;
		return "";   // Checksums disabled
	#endif
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
        
        #ifdef _WIN32
		Program.OperatingSystem = "Windows";
		#elif __linux__
		Program.OperatingSystem = "Linux";
		#elif __APPLE__
		Program.OperatingSystem = "macOS";
		#else
		Program.OperatingSystem = "Unknown";
		#endif

		#ifdef __GNUC__
		Program.Compiler = "GCC " + std::to_string(__GNUC__) + "." +
				           std::to_string(__GNUC_MINOR__);
		#elif defined(__clang__)
		Program.Compiler = "Clang";
		#elif defined(_MSC_VER)
		Program.Compiler = "MSVC";
		#else
		Program.Compiler = "Unknown";
		#endif       
		
		#ifdef MPI_PARALLEL
		Program.MpiRanks = MPI_SIZE;
		#else
		Program.MpiRanks = 1;
		#endif               
    }

    void MetaData::WriteJSON(std::string FileName, const Settings& settings, const RunTimeControl& RTC)
    {
    	OutPutDir  = settings.OutPutDir;
		VTKDir     = settings.VTKDir;
		RawDataDir = settings.RawDataDir;
		TextDir    = settings.TextDir;
		
		#ifdef __linux__
		char result[PATH_MAX];
		ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
		Program.ExecutablePath = std::string(result, (count > 0) ? count : 0);
		#else
		Program.ExecutablePath = "unknown";
		#endif
    
    	Parts.clear();
    
    	IndexOutputFiles();
    	    	    
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
		jApp["@type"] = "SoftwarejApp";
		jApp["@id"] = Program.ExecutablePath;
		jApp["name"] = Program.Title;
		jApp["softwareVersion"] = Program.VersionCore;

		jApp["identifier"] = "git:" + Program.GitHashCore;

		if (!Program.GitHashStudio.empty())
			jApp["studioVersion"] = Program.VersionStudio;

		jApp["operatingSystem"] = Program.OperatingSystem;
		jApp["compiler"] = Program.Compiler;
		jApp["mpiRanks"] = Program.MpiRanks;
		#ifdef OP_HAS_OPENSSL
		jApp["checksumSupport"] = "SHA256";
		#else
		jApp["checksumSupport"] = "disabled";
		#endif

		if (Program.RandomSeed != 0)
			jApp["randomSeed"] = Program.RandomSeed;
        j["@graph"].push_back(jApp);
        j["@graph"][0]["datePublished"] = to_iso8601(timestamp);
        j["@graph"][0]["description"] = Simulation.Description;
        j["@graph"][0]["hasPart"] = json::array();
		j["@graph"][1]["result"]  = json::array();
		for (size_t i = 0; i < Parts.size(); ++i)
		{
		    j["@graph"][0]["hasPart"][i]["@id"] = make_relative_path(Parts[i].Name,dir.string());
		    j["@graph"][1]["result"][i]["@id"] = make_relative_path(Parts[i].Name,dir.string());
		    json j1;
		    j1["@id"] = make_relative_path(Parts[i].Name,dir.string());
		    j1["@type"] = Parts[i].Type;
		    j1["description"] = Parts[i].Description;
		    j1["encodingFormat"] = Parts[i].Format;
		    j1["contentchecksum"] = Parts[i].Checksum;
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
        if (!RelatedWorks.empty())
		{
			j["@graph"][1]["mentions"] = json::array();

			for (size_t i = 0; i < RelatedWorks.size(); ++i)
			{
				const auto& work = RelatedWorks[i];

				json jWork;
				std::string workId = "./related_" + std::to_string(i) + ".json";

				jWork["@id"] = workId;
				jWork["@type"] = work.Type.empty() ? "ScholarlyArticle" : work.Type;
				jWork["name"] = work.Title;

				if (work.Year)
				    jWork["datePublished"] = std::to_string(work.Year);

				if (!work.DOI.empty())
				    jWork["identifier"] = work.DOI;

				if (!work.Journal.empty())
				    jWork["isPartOf"] = work.Journal;

				if (!work.URL.empty())
				    jWork["url"] = work.URL;

				if (!work.Description.empty())
				    jWork["description"] = work.Description;

				if (!work.Authors.empty())
				{
					jWork["author"] = json::array();

					for (size_t aIdx = 0; aIdx < work.Authors.size(); ++aIdx)
					{
						std::string personId =
							"#related_" + std::to_string(i) + "_author_" + std::to_string(aIdx);

						// Link from article
						jWork["author"].push_back({ {"@id", personId} });

						// Create separate Person entity
						json jPerson;
						jPerson["@id"] = personId;
						jPerson["@type"] = "Person";
						jPerson["name"] = work.Authors[aIdx];

						j["@graph"].push_back(jPerson);
					}
				}

				j["@graph"].push_back(jWork);

				// SAFE append
				j["@graph"][1]["mentions"].push_back({ {"@id", workId} });
			}
		}
        json simParams;

		std::string simId = "./simulation-parameters.json";

		simParams["@id"] = simId;
		simParams["@type"] = "PropertyValue";

		simParams["name"] = "Simulation Parameters";
		simParams["description"] = "Numerical and physical simulation settings";

		// Put everything as stringified JSON (safe approach)
		simParams["value"] = AddSimulationSettings(settings, RTC).dump();

		j["@graph"].push_back(simParams);

		// Link it
		j["@graph"][1]["hasPart"].push_back({ {"@id", simId} });
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
		std::string module;
		int moduleLocation = FileInterface::FindModuleLocation(inp, "AuthorInformation");
		if (moduleLocation != 0)
		{
		    int i = 0;
		    while (true)
		    {
		        module = "Author_" + std::to_string(i);
		        std::string tmpMd = FileInterface::ReadParameterFW(inp, moduleLocation, module, false, "NOAUTHOR");          
		        if (tmpMd != "NOAUTHOR")
		        {
		            AuthorInformation NewAuthor;
		            NewAuthor.Name = trim(tmpMd);

		            module = "Affiliation_" + std::to_string(i);
		            tmpMd = FileInterface::ReadParameterFW(inp, moduleLocation, module, false, "");
		            NewAuthor.Affiliation = trim(tmpMd);

		            module = "ORCID_" + std::to_string(i);
		            tmpMd = FileInterface::ReadParameterFW(inp, moduleLocation, module, false, "");
		            NewAuthor.ORCID = trim(tmpMd);

		            Authors.push_back(NewAuthor);
		            ++i;
		        }
		        else
		        {
		            break; // no more authors
		        }
		    }
		}
	}
    
    void MetaData::ReadSimulationInformation(std::stringstream& inp)
    {
        std::string module = " ";
        int moduleLocation = FileInterface::FindModuleLocation(inp, "SimulationInformation");
        if (moduleLocation != 0)
        {
            module = "Title";
            Simulation.Title = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"NOTITLE");      
            module = "Description";
            Simulation.Description = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"NODESCRIPTION");     
            module = "GUIGitHash";
            Program.GitHashStudio = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"0000");      
            module = "GUIVersion";
            Program.VersionStudio = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"0000");   
            
            int i = 0;
            do {
                module = "Keyword_" + std::to_string(i);
                std::string tmpMd  = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"NOKEYWORD");          
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
                std::string tmpMd  = FileInterface::ReadParameterFW(inp,moduleLocation,module,false,"NOARTICLE");          
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
    
    void MetaData::ReadRelatedWorks(std::stringstream& inp)
	{
		int moduleLocation = FileInterface::FindModuleLocation(inp, "RelatedWorks");
		if (moduleLocation == 0) return; // No RelatedWorks module found

		int i = 0;
		while (true)
		{
		    std::string keyBase = "Work_" + std::to_string(i);

		    std::string title = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_Title", false, "");
		    if (title.empty()) break; // No more entries

		    RelatedWork work;
		    work.Title = title;

		    // Authors as comma-separated string
		    std::string authorsStr = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_Authors", false, "");
		    if (!authorsStr.empty())
		    {
		        std::stringstream ss(authorsStr);
		        std::string author;
		        while (std::getline(ss, author, ','))
		        {
		            author.erase(0, author.find_first_not_of(" \t")); // trim leading spaces
		            author.erase(author.find_last_not_of(" \t") + 1); // trim trailing spaces
		            if (!author.empty())
		                work.Authors.push_back(author);
		        }
		    }

		    work.DOI     = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_DOI", false, "");
		    work.Journal = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_Journal", false, "");
		    work.Year    = FileInterface::ReadParameterI(inp, moduleLocation, keyBase + "_Year", false, 0); // <- integer
		    work.URL     = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_URL", false, "");
		    work.Description = FileInterface::ReadParameterFW(inp, moduleLocation, keyBase + "_Description", false, "");

		    RelatedWorks.push_back(work);
		    ++i;
		}
	}
    
    json MetaData::AddSimulationSettings(const Settings& settings, const RunTimeControl& RTC)
	{
		json simParams;

		// Phases and elements
		simParams["phases"] = settings.PhaseNames;
		simParams["elements"] = settings.ElementNames;

		// Grid parameters
		simParams["grid"]["Nx"] = settings.Grid.Nx;
		simParams["grid"]["Ny"] = settings.Grid.Ny;
		simParams["grid"]["Nz"] = settings.Grid.Nz;
		simParams["grid"]["dx"] = settings.Grid.dx;
		simParams["grid"]["gridBcells"] = settings.Grid.Bcells;
		simParams["grid"]["Eta"] = settings.Grid.Eta;
		simParams["resolution"] = (settings.Grid.Resolution == Resolutions::Single) ? "single" : "double";

		// Units
		simParams["units"]["length"] = RTC.UnitsOfLength;
		simParams["units"]["mass"]   = RTC.UnitsOfMass;
		simParams["units"]["time"]   = RTC.UnitsOfTime;
		simParams["units"]["energy"] = RTC.UnitsOfEnergy;

		// RTC
		simParams["RTC"]["dt"] = RTC.dt;
		simParams["RTC"]["simulationTime"] = RTC.SimulationTime;
		simParams["RTC"]["timeStep"] = RTC.TimeStep;
		simParams["RTC"]["maxTimeStep"] = RTC.MaxTimeStep;
		simParams["RTC"]["restartEnabled"] = RTC.RestartSwitch;

		// Parallel config
		simParams["parallel"]["openMPThreads"] = RTC.OpenMPThreads;

		// Output intervals
		simParams["outputIntervals"]["VTK"] = RTC.VTKOutputInterval;
		simParams["outputIntervals"]["console"] = RTC.ConsoleOutputInterval;
		simParams["outputIntervals"]["checkpoint"] = RTC.CheckpointInterval;

		return simParams;
	}

}// namespace openphase

