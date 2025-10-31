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

 *   File created :   2013
 *   Main contributors :
 *
 */

#include "Tools/TimeInfo.h"
#include "Settings.h"

namespace openphase
{

using namespace std;

TimeInfo::TimeInfo(const Settings& locSettings, const std::string Name, bool verbose_in)
{
    Initialize(locSettings, Name, verbose_in);
}

TimeInfo::~TimeInfo(void)
{
    ConsoleOutput::WriteStandard(thisclassname, "Exited normally");
}

void TimeInfo::Initialize(const Settings& locSettings, const std::string Name, bool verbose_in)
{
    thisclassname = "TimeInfo";
    counter = 0;
    ClockStart = 0;
    TimerName = Name;
    verbose = verbose_in;
}
void TimeInfo::Reset(void)
{
    counter = 0;
    ClockStart = 0;
}

void TimeInfo::SetStart(void)
{
    counter++;
    ClockStart = GetTime();
}

void TimeInfo::SetTimeStamp(const string Message)
{
    const double CurrentTime = GetTime();
    TimeMap[Message] += CurrentTime - ClockStart;
    ClockStart = CurrentTime;
    if(verbose)
    {
        ConsoleOutput::WriteSimple(Message);
    }
}

void TimeInfo::SkipToHere(void)
{
    const double CurrentTime = GetTime();
    TimeMap["#"] += CurrentTime - ClockStart;
    ClockStart = CurrentTime;
}

void TimeInfo::PrintWallClockSummary(void)
{
    if(counter)
    {
        ConsoleOutput::WriteLine("=");
        ConsoleOutput::WriteSimple(TimerName);
        ConsoleOutput::WriteLine("-");
        double TotalTime = 0;
        for (auto const& [message,time] : TimeMap) TotalTime += time;
        double TotalConsumedTime = TotalTime/double(counter);
        std::vector<std::pair<std::string, double>> tmp(TimeMap.begin(), TimeMap.end());
        std::sort(tmp.begin(), tmp.end(), [](const std::pair<std::string, double>& a, const std::pair<std::string, double>& b) { return a.second > b.second; });
        for (auto const& [message,time] : tmp)
        {
            if (message != "#")
            {
                std::stringstream showtime;

                double SectionConsumedTime = time/double(counter);
                    showtime << std::fixed << std::setprecision(2) << std::scientific << SectionConsumedTime << "  "
                             << std::fixed << std::setprecision(2) << ((SectionConsumedTime/TotalConsumedTime)*100.0) << " %";

                ConsoleOutput::WriteStandard(message, showtime.str());
            }
        }
        counter = 0;
        for (auto & [message,time] : TimeMap) time = 0;
        ConsoleOutput::WriteLine("-");
        ConsoleOutput::WriteStandard("Total time per time step [s]", TotalConsumedTime);
        ConsoleOutput::WriteLine("=");
    }
}

}// namespace openphase
