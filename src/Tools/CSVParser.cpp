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

 *   File created :   2015
 *   Main contributors :   Rajendran Mohan Kumar, Philipp Engels
 *
 */

#include "Tools/CSVParser.h"

namespace openphase
{
using namespace std;

void CSVParser::WriteHeader(const std::string fileName,
    const std::vector<std::string> headerArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str());
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteHeader()");
        OP_Exit(EXIT_FAILURE);
    }
    if (headerArray.size() == 0)
    {
        ConsoleOutput::WriteWarning("Tried to write empty header array", "CSVParser", "WriteHeader()");
    }
    else
    {
        for (auto iter = headerArray.cbegin(); iter != headerArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << headerArray.back() << endl;
    }
    file.close();
}

void CSVParser::ClearContent(const std::string fileName)
{
    std::ofstream file(fileName.c_str(), ios_base::trunc);
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "ClearContent()");
        OP_Exit(EXIT_FAILURE);
    }
    file.close();
}

void CSVParser::WriteData(const std::string fileName,
        const std::vector<int> dataArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str(), ios_base::app);
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteData()");
        OP_Exit(EXIT_FAILURE);
    }
    if (dataArray.size() == 0)
    {
        ConsoleOutput::WriteWarning("Tried to write empty data array", "CSVParser", "WriteData()");
    }
    else
    {
        for (auto iter = dataArray.cbegin(); iter != dataArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << dataArray.back() << endl;
    }
    file.close();
}

void CSVParser::WriteData(const std::string fileName,
        const std::vector<double> dataArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str(), ios_base::app);
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteData()");
        OP_Exit(EXIT_FAILURE);
    }
    if (dataArray.size() == 0)
    {
        ConsoleOutput::WriteWarning("Tried to write empty data array", "CSVParser", "WriteData()");
    }
    else
    {
        for (auto iter = dataArray.cbegin(); iter != dataArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << dataArray.back() << endl;
    }
    file.close();
}

void CSVParser::readFile(std::string& fileName,
    std::vector<std::vector<double>>& dataArray,
    std::vector<std::string>& headerArray, const std::string seperator)
{
    std::ifstream file(fileName.c_str());
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "readFile()");
        OP_Exit(EXIT_FAILURE);
    };

    const char csep = seperator[0];
    bool headerFlag = false;
    std::string line;
    while (getline(file, line))
    {
        if (!headerFlag) //populate headers
        {
            std::stringstream lineStream(line);
            std::string cell;
            headerArray.clear();
            while (std::getline(lineStream, cell, csep))
            {
                headerArray.push_back(cell);
            }
            headerFlag = true;
        }
        else //populate data
        {
            std::vector<double> row_data;
            std::stringstream lineStream(line);
            std::string cell;
            row_data.clear();
            while (std::getline(lineStream, cell, csep))
            {
                row_data.push_back(atof(cell.c_str()));
            }
            dataArray.push_back(row_data);
        }
    }
    file.close();
}

void CSVParser::readFile(std::string& fileName,
    std::vector<std::vector<std::string>>& dataArray,
    std::vector<std::string>& headerArray, std::string seperator)
{
    std::ifstream file(fileName.c_str());
    if (!file)
    {
        ConsoleOutput::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "readFile()");
        OP_Exit(EXIT_FAILURE);
    };

    const char csep = seperator[0];
    bool headerFlag = false;
    std::string line;
    while (getline(file, line))
    {
        if (!headerFlag) //populate headers
        {
            std::stringstream lineStream(line);
            std::string cell;
            headerArray.clear();
            while (std::getline(lineStream, cell, csep))
            {
                headerArray.push_back(cell);
            }
            headerFlag = true;
        }
        else //populate data
        {
            std::vector<std::string> row_data;
            std::stringstream lineStream(line);
            std::string cell;
            row_data.clear();
            while (std::getline(lineStream, cell, csep))
            {
                row_data.push_back(cell);
            }
            dataArray.push_back(row_data);
        }
    }
    file.close();
}

}// namespace openphase
