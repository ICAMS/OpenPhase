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

#ifndef SUBLATTICEMODEL_H
#define SUBLATTICEMODEL_H

#include "Includes.h"
#include "Thermodynamics/Element.h"

namespace openphase
{

class OP_EXPORTS SublatticeModel                                                           ///< Stores the sublattice properties
{
 public:
    size_t Ncons;
    bool hasVacancies;
    void Initialize(void);
    SublatticeModel(double n_sites, unsigned int n_cons = 0);                   ///< Constructor
    bool isElementPresent(int i);                                               ///< Returns boolean values, if constituent is present on sublattice
    std::vector<Element> Constituent;                                           ///< Constituents composing the sublattice
    double Site;                                                                ///< Sublattice contribution to the chemical formular unit
    int Reference;

 protected:
 private:
    size_t getNcons();
    bool gethasVacancies();
};

}
#endif//SUBLATTICEMODEL_H
