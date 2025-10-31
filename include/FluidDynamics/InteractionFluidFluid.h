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

 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev
 *
 */

#ifndef INTERACTIONFLUIDFLUID_H
#define INTERACTIONFLUIDFLUID_H

#include "Includes.h"

namespace openphase
{

class Settings;

class OP_EXPORTS InteractionFluidFluid : public OPObject                                   ///  Handles the interaction between immiscible liquids or liquid and gas phases
{
 public:

    void Initialize(Settings& locSettings) override;                            ///  Allocates memory and initiates global settings, structure described in OPSettings.h file
    void ReadInput(std::string InputFileName) override;                         ///  Reads input parameters from a file

 protected:
 private:

};

} //namespace openphase
#endif
