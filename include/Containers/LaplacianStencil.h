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
 *   Main contributors :   Oleg Shchyglo
 *
 */

#ifndef LAPLACIANSTENCIL_H
#define LAPLACIANSTENCIL_H

#include <cstddef>
#include <vector>

namespace openphase
{

enum class LaplacianStencils                                                    ///< List of available Laplacian stencils
{
    Simple,
    Isotropic,
    LB                                                                          ///< Isotropic stencils based on lattice Boltzmann stencils
};

// 1D, 2D and 3D Laplacian stencils.
// Note that 1D and 2D stencils for convenience defined in 3D arrays
// but should only be used with inactive dimensions suppressed.

/// Standard Laplacian stencils
const double LaplacianStencil1D_3[3][3][3] = {{{0.0,   0.0,   0.0},
                                               {0.0,   1.0,   0.0},
                                               {0.0,   0.0,   0.0}},

                                              {{0.0,   1.0,   0.0},
                                               {1.0,  -2.0,   1.0},
                                               {0.0,   1.0,   0.0}},

                                              {{0.0,   0.0,   0.0},
                                               {0.0,   1.0,   0.0},
                                               {0.0,   0.0,   0.0}}};           ///< 3 point 1D Laplacian stencil

const double LaplacianStencil2D_5[3][3][3] = {{{0.0,   0.0,   0.0},
                                               {0.0,   1.0,   0.0},
                                               {0.0,   0.0,   0.0}},

                                              {{0.0,   1.0,   0.0},
                                               {1.0,  -4.0,   1.0},
                                               {0.0,   1.0,   0.0}},

                                              {{0.0,   0.0,   0.0},
                                               {0.0,   1.0,   0.0},
                                               {0.0,   0.0,   0.0}}};           ///< 5 point simple 2D Laplacian stencil

const double LaplacianStencil3D_7[3][3][3] = {{{0.0,   0.0,  0.0},
                                               {0.0,   1.0,  0.0},
                                               {0.0,   0.0,  0.0}},

                                              {{0.0,   1.0,  0.0},
                                               {1.0,  -6.0,  1.0},
                                               {0.0,   1.0,  0.0}},

                                              {{0.0,   0.0,  0.0},
                                               {0.0,   1.0,  0.0},
                                               {0.0,   0.0,  0.0}}};            ///< 7 point simple 3D Laplacian stencil

///Isotropic Laplacian stencils
const double LaplacianStencil2D_9[3][3][3] = {{{    0.0,   1.0/6.0,     0.0},
                                               {1.0/6.0,   2.0/3.0, 1.0/6.0},
                                               {    0.0,   1.0/6.0,     0.0}},

                                              {{1.0/6.0,   2.0/3.0, 1.0/6.0},
                                               {2.0/3.0, -10.0/3.0, 2.0/3.0},
                                               {1.0/6.0,   2.0/3.0, 1.0/6.0}},

                                              {{    0.0,   1.0/6.0,     0.0},
                                               {1.0/6.0,   2.0/3.0, 1.0/6.0},
                                               {    0.0,   1.0/6.0,     0.0}}}; ///< 9 point 2D Laplacian stencil by Dave Hale

const double LaplacianStencil3D_19[3][3][3] = {{{     0.0,   1.0/6.0,      0.0},
                                                { 1.0/6.0,   1.0/3.0,  1.0/6.0},
                                                {     0.0,   1.0/6.0,      0.0}},

                                               {{ 1.0/6.0,   1.0/3.0,  1.0/6.0},
                                                { 1.0/3.0,      -4.0,  1.0/3.0},
                                                { 1.0/6.0,   1.0/3.0,  1.0/6.0}},

                                               {{     0.0,   1.0/6.0,      0.0},
                                                { 1.0/6.0,   1.0/3.0,  1.0/6.0},
                                                {     0.0,   1.0/6.0,      0.0}}};///< 19 point 3D Laplacian stencil from Patra and Karttunnen paper

const double LaplacianStencil3D_27a[3][3][3] = {{{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                                 {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                                 {1.0/30.0,   1.0/10.0, 1.0/30.0}},

                                                {{1.0/10.0,   7.0/15.0, 1.0/10.0},
                                                 {7.0/15.0, -64.0/15.0, 7.0/15.0},
                                                 {1.0/10.0,   7.0/15.0, 1.0/10.0}},

                                                {{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                                 {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                                 {1.0/30.0,   1.0/10.0, 1.0/30.0}}};///< 27 point 3D Laplacian stencil by Spotz and Carey (1995)

const double LaplacianStencil3D_27b[3][3][3] = {{{1.0/48.0,   1.0/8.0, 1.0/48.0},
                                                 {1.0/8.0,   5.0/12.0, 1.0/8.0},
                                                 {1.0/48.0,   1.0/8.0, 1.0/48.0}},

                                                {{1.0/8.0,   5.0/12.0, 1.0/8.0},
                                                 {5.0/12.0, -25.0/6.0, 5.0/12.0},
                                                 {1.0/8.0,   5.0/12.0, 1.0/8.0}},

                                                {{1.0/48.0,   1.0/8.0, 1.0/48.0},
                                                 { 1.0/8.0,  5.0/12.0, 1.0/8.0},
                                                 {1.0/48.0,   1.0/8.0, 1.0/48.0}}};///< 27 point 3D Laplacian stencil by Dave Hale

///Laplacian stencils based on lattice Boltzmann stencils
const double LaplacianStencil2D_LB[3][3][3] = {{{    0.0, 1.0/6.0,      0.0},
                                                {1.0/6.0, 2.0/3.0,  1.0/6.0},
                                                {    0.0, 1.0/6.0,      0.0}},

                                               {{1.0/6.0,   2.0/3.0, 1.0/6.0},
                                                {2.0/3.0, -10.0/3.0, 2.0/3.0},
                                                {1.0/6.0,   2.0/3.0, 1.0/6.0}},

                                               {{    0.0, 1.0/6.0,     0.0},
                                                {1.0/6.0, 2.0/3.0, 1.0/6.0},
                                                {    0.0, 1.0/6.0,     0.0}}};  ///< Laplacian stencil based on D2Q9 lattice Boltzmann stencil. It is the same as 2D stencil by Dave Hale

const double LaplacianStencil3D_LB[3][3][3] = {{{1.0/36.0, 1.0/9.0, 1.0/36.0},
                                                {1.0/9.0,  4.0/9.0, 1.0/9.0 },
                                                {1.0/36.0, 1.0/9.0, 1.0/36.0}},

                                               {{1.0/9.0,   4.0/9.0, 1.0/9.0},
                                                {4.0/9.0, -38.0/9.0, 4.0/9.0},
                                                {1.0/9.0,   4.0/9.0, 1.0/9.0}},

                                               {{1.0/36.0, 1.0/9.0, 1.0/36.0},
                                                {1.0/9.0,  4.0/9.0, 1.0/9.0 },
                                                {1.0/36.0, 1.0/9.0, 1.0/36.0}}};///< Laplacian stencil based on D3Q27 lattice Boltzmann stencil

class LaplacianStencil                                                          ///< Diffusion stencil class (uses user specified Laplacian stencil as the basis). Allows replacing the loop over Laplacian elements by the iterator which is beneficial for compact stencils
{
 public:
    struct LaplacianStencilEntry
    {
        int di;                                                                 ///< x coordinate of stencil element
        int dj;                                                                 ///< y coordinate of stencil element
        int dk;                                                                 ///< z coordinate of stencil element
        double weight;                                                          ///< Weight associated with the stencil element
    };

    void SetNoCenter(const double UserStencil[3][3][3], double dx,
                     int dNx = 1, int dNy = 1, int dNz = 1)                     ///< Sets the diffusion stencil without center element using user specified Laplacian stencil
    {
        if(not StencilElements.empty())
        {
            StencilElements.clear();
        }

        for(int x = -dNx; x <= dNx; ++x)
        for(int y = -dNy; y <= dNy; ++y)
        for(int z = -dNz; z <= dNz; ++z)
        if(x != 0 or y != 0 or z!= 0)
        if (UserStencil[x+1][y+1][z+1] != 0.0)
        {
            LaplacianStencilEntry locElement;
            locElement.di = x;
            locElement.dj = y;
            locElement.dk = z;

            locElement.weight = UserStencil[x+1][y+1][z+1]/(dx*dx);
            StencilElements.push_back(locElement);
        }
        weight_zero = UserStencil[1][1][1];
    }
    void SetNoCenter(const double UserStencil[3][3][3], GridParameters& Dimensions)
    {
        SetNoCenter(UserStencil,Dimensions.dx, Dimensions.dNx, Dimensions.dNy, Dimensions.dNz);
    }
    void Set(const double UserStencil[3][3][3], double dx,
             int dNx = 1, int dNy = 1, int dNz = 1)                             ///< Sets the diffusion stencil using user specified Laplacian stencil
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
            LaplacianStencilEntry locElement;
            locElement.di = x;
            locElement.dj = y;
            locElement.dk = z;

            locElement.weight = UserStencil[x+1][y+1][z+1]/(dx*dx);
            StencilElements.push_back(locElement);
        }
        weight_zero = UserStencil[1][1][1];
    }
    void Set(const double UserStencil[3][3][3], GridParameters& Dimensions)
    {
        Set(UserStencil,Dimensions.dx, Dimensions.dNx, Dimensions.dNy, Dimensions.dNz);
    }

    size_t size() const {return StencilElements.size();};                       ///< Returns the size of stencil.
    typedef std::vector<LaplacianStencilEntry>::iterator iterator;              ///< Iterator over stencil
    typedef std::vector<LaplacianStencilEntry>::const_iterator citerator;       ///< Constant iterator over stencil
    iterator  begin() {return StencilElements.begin();};                        ///< Iterator to the begin of stencil
    iterator  end()   {return StencilElements.end();};                          ///< Iterator to the end of stencil
    citerator cbegin() const {return StencilElements.cbegin();};                ///< Constant iterator to the begin of stencil
    citerator cend()   const {return StencilElements.cend();};                  ///< Constant iterator to the end of stencil
    iterator erase(iterator it) {it = StencilElements.erase(it); return it;}    ///< Erases entry

    double weight_zero;                                                         ///< Extra copy of stencil center element
 protected:
 private:
    std::vector< LaplacianStencilEntry > StencilElements;                       ///< Stencil storage
};
}
#endif//LAPLACIANSTENCIL_H
