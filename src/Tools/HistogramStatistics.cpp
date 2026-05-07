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

#include "Tools/HistogramStatistics.h"

namespace openphase
{
using namespace std;

void HistogramStatistics::Print(vector<double>& data, int numBins, int percentilesteps)
{
    if (data.empty())
    {
        ConsoleOutput::WriteLineInsert("HistogramStatistics::Print -> No data.");
        return;
    }
    auto CreateBins = [](const vector<double>& data, int numBins)
    {
        vector<pair<double, double>> bins;
        double minValue = *min_element(data.begin(), data.end());
        double maxValue = *max_element(data.begin(), data.end());
        double binWidth = (maxValue - minValue) / numBins;

        bins.reserve(numBins);
        for (int i = 0; i < numBins; ++i)
        {
            double binStart = minValue + i * binWidth;
            double binEnd   = binStart + binWidth;
            bins.emplace_back(binStart, binEnd);
        }
        return bins;
    };

    // Print histogram bins with ASCII bar
    auto PrintBins = [](const vector<double>& data,
                        const vector<pair<double, double>>& bins)
    {
        vector<int> counts(bins.size(), 0);
        int total_count = static_cast<int>(data.size());

        for (double value : data)
        {
            for (size_t i = 0; i < bins.size(); ++i)
            {
                if (value >= bins[i].first && value < bins[i].second)
                {
                    counts[i]++;
                    break;
                }
            }
        }

        ConsoleOutput::WriteLineInsert("Histogram");
        ConsoleOutput::WriteLine();

        ostringstream oss;
        oss << left     << setw(9)  << "N-Bin "
            << setw(10) << "[Start"  <<   " -"
            << right    << setw(12)  << "End] "
            << setw(10) << " Frequency \n";

        for (size_t i = 0; i < bins.size(); ++i)
        {
            double percentage = (static_cast<double>(counts[i]) / total_count) * 100.0;
            int num_dots = static_cast<int>(round(percentage));

            oss << left << "Bin " << setw(5) << i + 1
                << scientific << setprecision(2) << "["
                << setw(10) << bins[i].first << "-  "
                << setw(8)  << bins[i].second << "]  "
                << setw(8)  << counts[i]
                << fixed << setprecision(2) << setw(5) << percentage << "[%]  | "
                << setw(num_dots) << setfill('.') << ""
                << setfill(' ') << "\n";
        }
        cout << oss.str();
    };

    // Calculate percentiles at given step intervals
    auto CalculatePercentile = [](vector<double> data, int step)
    {
        sort(data.begin(), data.end());
        dVectorN percentiles;
        percentiles.push_back(data.front());
        for (int percentile = step; percentile <= 100 - step; percentile += step)
        {
            int index = static_cast<int>((percentile / 100.0) * data.size());
            percentiles.push_back(data[index]);
        }
        percentiles.push_back(data.back());
        return percentiles;
    };

    if (numBins > 0)
    {
        auto bins = CreateBins(data, numBins);
        PrintBins(data, bins);
        ConsoleOutput::WriteLine();
    }

    if (percentilesteps > 0)
    {
        ConsoleOutput::WriteStandard("Total Count", data.size());
        ConsoleOutput::WriteLine("=");

        dVectorN PercentileVector = CalculatePercentile(data, percentilesteps);
        ConsoleOutput::WriteLineInsert("Percentile Steps : " + to_string(percentilesteps)
                                       + " > < No of counts : "
                                       + to_string(data.size() / percentilesteps));
        ConsoleOutput::WriteLine();

        for (size_t index = 0; index < PercentileVector.size(); index++)
        {
            string indexNumber = "[" + to_string(index) + "]";
            cout << right << setfill(' ') << setprecision(2)
                 << setw(4) << indexNumber
                 << setw(8) << PercentileVector[index] << " | ";
            if (!(index % 5))
                cout << endl;
        }
        ConsoleOutput::WriteLine("=");
    }
}

} // namespace openphase
