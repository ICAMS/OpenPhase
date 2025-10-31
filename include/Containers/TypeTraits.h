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

 *   File created :   2022
 *   Main contributors : Raphael Schiedung
 *
 *
 */

#ifndef TYPETRAITS_H
#define TYPETRAITS_H

#include <type_traits>

#include "NodeA.h"
#include "NodeAB.h"
#include "NodeAdv.h"
#include "NodePF.h"
#include "NodeDF.h"
#include "NodeIP.h"
#include "NodeVectorN.h"
#include "dMatrix3x3.h"
#include "dMatrix6x6.h"
#include "dMatrixNxN.h"
#include "dVector3.h"
#include "dVector6.h"
#include "dVectorN.h"
#include "iVector3.h"
#include "vStrain.h"
#include "vStress.h"

namespace openphase
{

template <typename T>
using is_vector = std::integral_constant<bool,
      std::is_same<T, dVector3>::value |
      std::is_same<T, iVector3>::value |
      std::is_same<T, dVector6>::value |
      std::is_same<T, dVectorN>::value |
      std::is_same<T, vStrain >::value |
      std::is_same<T, vStress >::value >;

template <typename T>
using is_matrix = std::integral_constant<bool,
      std::is_same<T, dMatrix3x3>::value |
      std::is_same<T, dMatrix6x6>::value |
      std::is_same<T, dMatrixNxN>::value >;

template <typename T, typename C1 = double, typename C2 = double>
using is_node = std::integral_constant<bool,
      std::is_same<T,  NodeA<C1>     >::value |
      std::is_same<T,  NodeAB<C1,C2> >::value |
      std::is_same<T,  NodeAdv       >::value |
      std::is_same<T,  NodePF        >::value |
      std::is_same<T,  NodeDF        >::value |
      std::is_same<T,  NodeIP        >::value |
      std::is_same<T,  NodeVectorN   >::value >;

template <typename T>
using has_datapointer = std::integral_constant<bool,
      //std::is_same<T,  std::vector<TT>    >::value | //TODO std::vector<> has data() method
      std::is_same<T,  dVector3     >::value |
      std::is_same<T,  iVector3     >::value >;

template <typename T>
class has_clear
{
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::clear) ) ;
    template <typename C> static two test(...);

public:
    enum
    {
        value = sizeof(test<T>(0)) == sizeof(char)
    };
};

template <typename T>
class has_pack
{
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::pack) ) ;
    template <typename C> static two test(...);

public:
    enum
    {
        value = sizeof(test<T>(0)) == sizeof(char)
    };
};

} // end of namespace ophenphase
#endif
