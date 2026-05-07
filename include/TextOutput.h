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

#ifndef TEXTOUTPUT_H
#define TEXTOUTPUT_H

#include "Includes.h"

namespace openphase
{

class TextOutput                                                                ///< Output of simulation data in ASCI format
{
 public:
    static void SetSeparator(std::string new_separator);                        ///< Sets new text separator (has global effect)
    static void ResetSeparator(void);                                           ///< Resets text separator to its default value (has global effect)

    static void WriteValue(std::string filename,
                           double time,
                           double value,
                           std::ios_base::fmtflags format =
                           std::ios_base::scientific);                          ///< Writes double value in a column over time
    static void WriteNamedValues(std::string file_name,
                                 std::vector<std::string> names,
                                 std::vector<double> values,
                                 std::ios_base::fmtflags format =
                                 std::ios_base::scientific);                    ///< Writes vector of values in the named columns using specified format
    static void WritedVectorNValues(std::string filename,
                                    std::vector<std::string> names,
                                    double time,
                                    dVectorN values,
                                    std::ios_base::fmtflags format =
                                    std::ios_base::scientific);                 ///< Write dVectorN values and names of values in columns over time
 private:
    static std::string separator;                                               ///< Text separator (default ", ")
};

} // namespace openphase
#endif
