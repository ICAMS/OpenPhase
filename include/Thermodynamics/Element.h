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
 *   Main contributors :   Oleg Shchyglo; Matthias Stratmann
 *
 */

#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>

namespace openphase
{

class Element                                                                   ///< Stores properties of chemical element
{
 public:
    int             Index;                                                      ///< Element's index in components list.
    double          Value;                                                      ///< Element's amount
    double          Min;                                                        ///< Minimum mole fraction of the component
    double          Max;                                                        ///< Maximum mole fraction of the component
    std::string     Name;                                                       ///< Element's name (2 characters long max).
    bool            isVacancy;                                                  ///< True if the element is a vacancy (VA)
    bool            isStoichiometric;                                           ///< True if solubility range of the element is extremely narrow
    bool            isReference;                                                ///< True if the element is used as a reference element
    bool            isFastDiffusor;                                             ///< True if the element can be considered as diffusing extremely fast
};

}
#endif//ELEMENT_H
