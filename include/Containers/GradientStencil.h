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

 *   File created :   2020
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef GRADIENTSTENCIL_H
#define GRADIENTSTENCIL_H

namespace openphase
{

enum class GradientStencils                                                     ///< List of available gradient stencils
{
    Simple,
    Isotropic,
    LB                                                                          ///< Isotropic stencils based on lattice Boltzmann stencils
};

// 1D, 2D and 3D Gradient stencils.
// Note that 1D and 2D stencils for convenience defined in 3D arrays
// but should only be used with inactive dimensions suppressed.
// Gradient stencils have no negative coefficients which should be multiplied
// by (-1) manually at the point of use.

const double GradientStencil1D[3][3][3] = {{{ 0.0, 0.0, 0.0},
                                            { 0.0, 0.5, 0.0},
                                            { 0.0, 0.0, 0.0}},

                                           {{ 0.0, 0.5, 0.0},
                                            { 0.5, 0.0, 0.5},
                                            { 0.0, 0.5, 0.0}},

                                           {{ 0.0, 0.0, 0.0},
                                            { 0.0, 0.5, 0.0},
                                            { 0.0, 0.0, 0.0}}};                 ///< 2 point simple 1D/2D/3D Gradient stencil (standard finite differences)

const double GradientStencil2D[3][3][3] = {{{      0.0, 0.25/3.0,      0.0},
                                            { 0.25/3.0,  1.0/3.0, 0.25/3.0},
                                            {      0.0, 0.25/3.0,      0.0}},

                                           {{ 0.25/3.0,  1.0/3.0, 0.25/3.0},
                                            {  1.0/3.0,      0.0,  1.0/3.0},
                                            { 0.25/3.0,  1.0/3.0, 0.25/3.0}},

                                           {{      0.0, 0.25/3.0,      0.0},
                                            { 0.25/3.0,  1.0/3.0, 0.25/3.0},
                                            {      0.0, 0.25/3.0,      0.0}}};  ///< 8 point 2D Gradient stencil from "M.Alfaraj, Y. Wang and Y. Luo, Geophysical prospecting 62 (2014) 507-517"

const double GradientStencil2D_2[3][3][3] = {{{   0.0, 0.125,   0.0},
                                              { 0.125,  0.25, 0.125},
                                              {   0.0, 0.125,   0.0}},

                                             {{ 0.125,  0.25, 0.125},
                                              {  0.25,   0.0,  0.25},
                                              { 0.125,  0.25, 0.125}},

                                             {{   0.0, 0.125,   0.0},
                                              { 0.125,  0.25, 0.125},
                                              {   0.0, 0.125,   0.0}}};         ///< 8 point 2D Gradient stencil by Sobel

const double GradientStencil3D[3][3][3] = {{{ 0.085/4.64, 0.245/4.64, 0.085/4.64},
                                            { 0.245/4.64,   1.0/4.64, 0.245/4.64},
                                            { 0.085/4.64, 0.245/4.64, 0.085/4.64}},

                                           {{ 0.245/4.64, 1.0/4.64, 0.245/4.64},
                                            {   1.0/4.64,      0.0,   1.0/4.64},
                                            { 0.245/4.64, 1.0/4.64, 0.245/4.64}},

                                           {{ 0.085/4.64, 0.245/4.64, 0.085/4.64},
                                            { 0.245/4.64,   1.0/4.64, 0.245/4.64},
                                            { 0.085/4.64, 0.245/4.64, 0.085/4.64}}};///< 27 point 3D Gradient stencil from "M.Alfaraj, Y. Wang and Y. Luo, Geophysical prospecting 62 (2014) 507-517"

/// Gradient stencils based on lattice Boltzmann stencils

const double GradientStencil2D_LB[3][3][3] = {{{     0.0, 1.0/12.0,      0.0},
                                               {1.0/12.0, 1.0/3.0,  1.0/12.0},
                                               {     0.0, 1.0/12.0,      0.0}},

                                              {{1.0/12.0,   1.0/3.0, 1.0/12.0},
                                               {1.0/3.0,        0.0, 1.0/3.0},
                                               {1.0/12.0,   1.0/3.0, 1.0/12.0}},

                                               {{     0.0, 1.0/12.0,      0.0},
                                                {1.0/12.0, 1.0/3.0,  1.0/12.0},
                                                {     0.0, 1.0/12.0,      0.0}}};///< Isotropic gradient stencil based on D2Q9 lattice Boltzmann stencil. It is the same as the 8 point 2D Gradient stencil from "M.Alfaraj, Y. Wang and Y. Luo, Geophysical prospecting 62 (2014) 507-517"

const double GradientStencil3D_LB[3][3][3] = {{{1.0/72.0, 1.0/18.0, 1.0/72.0},
                                               {1.0/18.0, 2.0/9.0,  1.0/18.0},
                                               {1.0/72.0, 1.0/18.0, 1.0/72.0}},

                                              {{1.0/18.0, 2.0/9.0, 1.0/18.0},
                                               {2.0/9.0,      0.0, 2.0/9.0},
                                               {1.0/18.0, 2.0/9.0, 1.0/18.0}},

                                              {{1.0/72.0, 1.0/18.0, 1.0/72.0},
                                               {1.0/18.0, 2.0/9.0,  1.0/18.0},
                                               {1.0/72.0, 1.0/18.0, 1.0/72.0}}};///< Isotropic gradient stencil based on D3Q27 lattice Boltzmann stencil. It is very close to the 27 point 2D Gradient stencil from "M.Alfaraj, Y. Wang and Y. Luo, Geophysical prospecting 62 (2014) 507-517"

class GradientStencil                                                           ///< Gradient stencil class (uses user specified stencil as the basis). Allows replacing the loop over array elements by the iterator which is beneficial for compact stencils
{
 public:
    struct GradientStencilEntry
    {
        int di;                                                                 ///< x coordinate of stencil element
        int dj;                                                                 ///< y coordinate of stencil element
        int dk;                                                                 ///< z coordinate of stencil element
        double weightX;                                                         ///< Weight associated with the stencil element
        double weightY;                                                         ///< Weight associated with the stencil element
        double weightZ;                                                         ///< Weight associated with the stencil element
    };
    void Set(const double UserStencil[3][3][3], double dx,
             int dNx = 1, int dNy = 1, int dNz = 1)                             ///< Sets nonzero stencil entries for gradient components using user specified stencil
    {
        if(not StencilElements.empty())
        {
            StencilElements.clear();
        }

        for(int x = -dNx; x <= dNx; ++x)
        for(int y = -dNy; y <= dNy; ++y)
        for(int z = -dNz; z <= dNz; ++z)
        if (UserStencil[x+1][y+1][z+1] != 0.0)
        {
            GradientStencilEntry locElement;
            locElement.di = x;
            locElement.dj = y;
            locElement.dk = z;

            locElement.weightX = x*UserStencil[x+1][y+1][z+1]/dx;
            locElement.weightY = y*UserStencil[x+1][y+1][z+1]/dx;
            locElement.weightZ = z*UserStencil[x+1][y+1][z+1]/dx;

            StencilElements.push_back(locElement);
        }
    }
    void Set(const double UserStencil[3][3][3], GridParameters& Dimensions)
    {
        Set(UserStencil, Dimensions.dx, Dimensions.dNx, Dimensions.dNy, Dimensions.dNz);
    }
    size_t size() const {return StencilElements.size();};                       ///< Returns the size of stencil.
    typedef std::vector<GradientStencilEntry>::iterator iterator;               ///< Iterator over stencil
    typedef std::vector<GradientStencilEntry>::const_iterator citerator;        ///< Constant iterator over stencil
    iterator  begin() {return StencilElements.begin();};                        ///< Iterator to the begin of stencil
    iterator  end()   {return StencilElements.end();};                          ///< Iterator to the end of stencil
    citerator cbegin() const {return StencilElements.cbegin();};                ///< Constant iterator to the begin of stencil
    citerator cend()   const {return StencilElements.cend();};                  ///< Constant iterator to the end of stencil
    iterator erase(iterator it) {it = StencilElements.erase(it); return it;}    ///< Erases entry
 protected:
 private:
    std::vector< GradientStencilEntry > StencilElements;                        ///< Stencil storage
};
}
#endif//GRADIENTSTENCIL_H
