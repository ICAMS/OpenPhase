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
 *   Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann
 *
 */

#ifndef MOVINGFRAME_H
#define MOVINGFRAME_H

#include "Includes.h"

namespace openphase
{
class Settings;
class BoundaryConditions;

class OP_EXPORTS MovingFrame : public OPObject                                  ///< Moving frame class
{
 public:
    MovingFrame();                                                              ///< Constructor
    MovingFrame(Settings& locSettings,
                std::string InputFileName = DefaultInputFileName)               ///< Constructor. Sets internal variables and reads input.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates memory, initializes the settings);
    void ReadInput(const std::string FileName) override;                        ///< Reads input parameters
    void ReadInput(std::stringstream& FileName) override;                       ///< Reads input parameters

    size_t trigger_phase_idx;                                                   ///< Phase index which triggers the moving frame
    size_t trigger_position;                                                    ///< Position which should be reached by the trigger_phase to trigger the moving frame
    std::string moving_direction;                                               ///< Direction to move the frame into

 protected:
 private:

    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC,
              const int tStep) override                                         ///< Reads raw data from a file
    {
        // Does nothing and should not be used
        return false;
    }
    bool Write(const Settings& locSettings, const int tStep) const override     ///< Writes raw data to a file
    {
        // Does nothing and should not be used
        return false;
    }
};

} // namespace openphase
#endif
