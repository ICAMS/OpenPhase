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

#include "Thermodynamics/SublatticeModel.h"

namespace openphase
{

using namespace std;

void SublatticeModel::Initialize(void)
{
    Ncons = getNcons();
    hasVacancies = gethasVacancies();
}

size_t SublatticeModel::getNcons()
{
    /**This will return the number of constituents of this sublattice.*/
    return Constituent.size();
}

bool SublatticeModel::gethasVacancies()
{
    /**This function will loop over all constituents and checks, if one of them
    is considered a vacancy.*/
    bool temp = false;
    for(size_t i = 0; i < Ncons; i++)
    if(Constituent[i].Index < 0)
    {
        temp = true;
    }
    return temp;
}

bool SublatticeModel::isElementPresent(int i)
{
    /**This function will loop over all constituents and checks, if element i
    is present on this sublattice. "i" has to be the index, not the component
    number!!*/
    bool temp = false;
    for(size_t n = 0; n < Ncons; n++)
    {
        if(Constituent[n].Index == i)
        temp = true;
    }
    return temp;
}

SublatticeModel::SublatticeModel(double n_sites, unsigned int n_cons)
{
    /**Constructor*/
    Site = n_sites;
}

}
