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
 *  File created :   2026
 *  Main contributors :   Muhammad Adil Ali
 *
 */

#ifndef HISTOGRAMSTATISTICS_H
#define HISTOGRAMSTATISTICS_H

#include "Includes.h"

namespace openphase
{

class OP_EXPORTS HistogramStatistics
{
 public:
    /// @brief Prints histogram and percentile statistics for the given data.
    /// numBins         = 0   -> to skip the histogram,
    /// percentilesteps = 0   -> to skip percentiles.
    static void Print(std::vector<double>& data, int numBins, int percentilesteps);
};

} // namespace openphase

#endif // HISTOGRAMSTATISTICS_H
