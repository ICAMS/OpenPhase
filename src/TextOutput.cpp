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
 *  File created :   2014
 *  Main contributors :   Matthias Stratmann; Johannes Goerler
 *
 */

#include "TextOutput.h"

namespace openphase
{

std::string default_separator = ", ";
std::string TextOutput::separator = default_separator;


void TextOutput::SetSeparator(std::string new_separator)
{
    separator = new_separator;
}

void TextOutput::ResetSeparator(void)
{
    separator = default_separator;
}

void TextOutput::WriteValue(std::string filename,
                            double time,
                            double value,
                            std::ios_base::fmtflags format)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            std::ofstream file(filename, std::ios::out);
            file << "Time" << separator << "Value\n";
            file.close();
        }

        std::ofstream file(filename, std::ios::app);
        file << time << separator; 
        file.setf(format, std::ios_base::floatfield);
        file << value << "\n";
        file.close();
    }
}

void TextOutput::WriteNamedValues(std::string file_name,
                                  std::vector<std::string> names,
                                  std::vector<double> values,
                                  std::ios_base::fmtflags format)
{
     assert(names.size() == values.size() and "The number of names is not equal to the number of values");
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(file_name))
        {
            std::ofstream file(file_name, std::ios::out);
            for (size_t n = 0; n < names.size() - 1; n++)
            {
                file << names[n] << separator;
            }
            file << names.back() << "\n";
            file.close();
        }

        std::ofstream file(file_name, std::ios::app);
        file.setf(format, std::ios_base::floatfield); 
        for (size_t n = 0; n < values.size() - 1; n++)
        {
            file << values[n] << separator;
        }
        file << values.back() << "\n";
        file.close();
    }
}

void TextOutput::WritedVectorNValues(std::string filename,
                                     std::vector<std::string> names,
                                     double time,
                                     dVectorN values,
                                     std::ios_base::fmtflags format)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            std::ofstream file(filename, std::ios::out);
            file << "Time";
            for (auto& N : names)
            {
                file << separator << N;
            }
            file << "\n";
            file.close();
        }

        std::ofstream file(filename, std::ios::app);
        file << time;
        file.setf(format, std::ios_base::floatfield); 
        for (size_t i = 0; i < values.size(); i++ )
        {
            file << separator << values[i];
        }
        file << "\n";
        file.close();
    }
}

}// namespace openphase
