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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Philipp Engels; Johannes Goerler
 *
 */

#ifndef FILEINTERFACE_H
#define FILEINTERFACE_H

#include "Includes.h"
#include "../external/json.hpp"
#include <variant>
#include <vector>
#include <string>
#include <optional>

namespace openphase
{

using json = nlohmann::json;

// Define a path element: either a key (string) or an index (size_t)
using JsonPathElement = std::variant<std::string, size_t>;

// Traverses a json object safely based on the path

class OP_EXPORTS FileInterface                                                  ///< Input file interface class
{
 public:

    static bool PRINT_DEFAULTS;                                                 ///< Enables printing default values if no user input is provided for non-mandatory parameters
    
    static std::string getFileExtension(const std::string& filename);			///< Returns FileExtension or an empty string 

    static std::string MakeFileName(std::string Directory,
                                    std::string NameBase,
                                    int Index,
                                    std::string FileExtension);                 ///< Creates full filename string using the location directory, name base, running index and file extension

    // Methods to read input parameters from OpenPhase input files.
    // They search for $KEY in the entire file using the following syntax:
    // $KEY    commment    :   value

    static bool ParameterPresent(const std::stringstream& Inp,                  ///< Checks if Key is present in the specified module, returns true if parameter is found and false otherwise
                                 const int location, std::string Key);

    static int FindParameter(const std::stringstream& Inp,                      ///< Checks if Key is present in the specified module and returns location right before the key
                             const int location, std::string Key);

    static int FindParameterLocation(const std::stringstream& Inp,              ///< Checks if Key is present in the specified module and returns location right after the key
                             const int location, std::string Key);

    static int FindModuleLocation(const std::stringstream& Inp,                 ///< Returns location of the module
                                const std::string module);

    static double ReadParameterD(const std::stringstream& Inp,                  ///< Read double precision floating point parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const double defaultval = 0.0);
        
    static double ReadParameterD(const json& j, const std::string module, std::string Key,
    const bool mandatory, const double defaultval);
        
     static double ReadParameterD(const std::stringstream& Inp,                 ///< Read double precision floating point parameter value
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        const double defaultval = 0.0);

    static dVector3 ReadParameterV3(const std::stringstream& Inp,               ///< Reads dVector3 vector
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const dVector3 defaultval = dVector3::ZeroVector());
        
    static dVector3 ReadParameterV3(const std::stringstream& Inp,               ///< Reads dVector3 vector
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        const dVector3 defaultval = dVector3::ZeroVector());

    static iVector3 ReadParameterV3I(const std::stringstream& Inp,              ///< Reads iVector3 vector
                int currentLocation,
                std::string Key,
                const bool mandatory = true,
                const iVector3 defaultval = iVector3::ZeroVector());

    static dMatrix3x3 ReadParameterM3x3(const std::stringstream& Inp,           ///< Reads dMatrix3x3 tensor
            int currentLocation,
            const std::string Key,
            const bool mandatory = true,
            const dMatrix3x3 defaultval = dMatrix3x3::ZeroTensor());

    static dMatrix6x6 ReadParameterT6(const std::stringstream& Inp,             ///< Read dMatrix6x6 tensor
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const dMatrix6x6 defaultval = dMatrix6x6::UnitTensor());

    static dMatrix6x6 ReadParameterM6x6(const std::stringstream& Inp,           ///< Reads dMatrix6x6 tensor
                int currentLocation,
                const std::string Key,
                const bool mandatory = true,
                const dMatrix6x6 defaultval = dMatrix6x6::ZeroTensor());

    static int ReadParameterI(const std::stringstream& Inp,                     ///< Reads integer parameter value
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        int const defaultval = 0);
        
    static int ReadParameterI(const std::stringstream& Inp,                     ///< Reads integer parameter value
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        int const defaultval = 0);

    static std::string ReadParameterS(const std::stringstream& Inp,             ///< Reads string parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterS(const std::stringstream& Inp,             ///< Reads string parameter value
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterK(const std::stringstream& Inp,             ///< Reads keyword (removes spaces and converts to upper case)
            int currentLocation,
            const std::string Key,
            const bool mandatory = true,
            const std::string defaultval = "NN");

    static std::string ReadParameterK(const std::stringstream& Inp,             ///< Reads keyword (removes spaces and converts to upper case)
            int currentLocation,
            const std::vector<std::string> Key,
            const bool mandatory = true,
            const std::string defaultval = "NN");

    static std::vector<std::string> ReadParameterVS(const std::stringstream& Inp, ///< Reads list of strings as a vector separated by any punctuation marks
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const std::vector<std::string> defaultval = {""});

    static bool ReadParameterB(const std::stringstream& Inp,                    ///< Read boolean parameter value
        int currentLocation,
        std::string Key,
        const bool mandatory = true,
        const bool defaultval = false);
        
    static bool ReadParameterB(const std::stringstream& Inp,                    ///< Read boolean parameter value
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        const bool defaultval = false);

    static std::string ReadParameterF(const std::stringstream& Inp,             ///< Reads filename string (assumes no spaces in the filename)
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterF(const std::stringstream& Inp,             ///< Reads filename string (assumes no spaces in the filename)
        int currentLocation,
        const std::vector<std::string> Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");

    static std::string ReadParameterFW(const std::stringstream& Inp,            ///< Reads filename string (assumes no spaces in the filename)
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");
        
    static std::string ReadParameterFW(const std::stringstream& Inp,            ///< Reads filename string (assumes no spaces in the filename)
        int currentLocation,
       	const std::vector<std::string> Key,
        const bool mandatory = true,
        const std::string defaultval = "NN");    

    static char ReadParameterC(const std::stringstream& Inp,                    ///< Reads char parameter value
        int currentLocation,
        const std::string Key,
        const bool mandatory = true,
        const char defaultval = 'X');

    static void removeCarriageReturn(std::string& str);                         ///< Removes carriage return '\r' from string
    
    template<typename T>
    static std::optional<T> try_eval_json(const json& j, const std::vector<JsonPathElement>& path)
    {
        const json* current = &j;

        for (const auto& elem : path) 
        {
            if (std::holds_alternative<std::string>(elem)) {
                const auto& key = std::get<std::string>(elem);
                if (!current->is_object() || !current->contains(key)) return std::nullopt;
                std::cout << "/" << key;
                current = &(*current)[key];
            } else {
                const auto& key = std::to_string((int)std::get<size_t>(elem));
                if (!current->is_object() || !current->contains(key)) return std::nullopt;
                std::cout << "/" << key;
                current = &(*current)[key];
            }
        }

        try
        {
            return current->get<T>();
        }
        catch (...)
        {
            return std::nullopt;
        }
    }

    // Version 1: With default fallback
    template<typename T>
    static T ReadParameter(const json& j, const std::vector<JsonPathElement>& path, const T& defaultval) {
        auto result = try_eval_json<T>(j, path);
        std::cout << " = " << result.value_or(defaultval) << std::endl;
        return result.value_or(defaultval);
    }

    // Version 2: Without default — throws if not found
    template<typename T>
    static T ReadParameter(const json& j, const std::vector<JsonPathElement>& path) {
        auto result = try_eval_json<T>(j, path);
        if (!result.has_value()) {
            std::cerr << "No Entry available for";
            for (const auto& elem : path) 
            {
                if (std::holds_alternative<std::string>(elem)) {
                    const auto& key = std::get<std::string>(elem);
                    std::cerr << " " << key;
                } else {
                    size_t index = std::get<size_t>(elem);
                    std::cerr << " " << index;
                }
            }
            std::cerr << " and no default value provided." << std::endl;
            OP_Exit(1);
        }
        std::cout << " = " << result.value() << std::endl;
        return result.value();
    }

};

inline void FileInterface::removeCarriageReturn(std::string& str)
{
    str.erase(std::remove(str.begin(), str.end(), '\r'), str.end());
}

}// namespace openphase
#endif // FILEINTERFACE_H
