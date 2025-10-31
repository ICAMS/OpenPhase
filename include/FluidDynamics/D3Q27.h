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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#ifndef D3Q27_H
#define D3Q27_H

#include <vector>
#include <array>
#include <cassert>

namespace openphase
{

// 1D, 2D and 3D lattice Boltzmann stencils.
// Note that 1D and 2D stencils are defined in 3D arrays
// but should only be used with inactive dimensions suppressed.

const double LBStencil1D[3][3][3] = {{{    0.0,     0.0,     0.0},
                                      {    0.0, 1.0/6.0,     0.0},
                                      {    0.0,     0.0,     0.0}},

                                     {{    0.0, 1.0/6.0,     0.0},
                                      {1.0/6.0, 4.0/6.0, 1.0/6.0},
                                      {    0.0, 1.0/6.0,     0.0}},

                                     {{    0.0,     0.0,     0.0},
                                      {    0.0, 1.0/6.0,     0.0},
                                      {    0.0,     0.0,     0.0}}};            ///< D1Q3 lattice Boltzmann stencil

const double LBStencil2D[3][3][3] = {{{     0.0, 1.0/36.0,      0.0},
                                      {1.0/36.0, 1.0/9.0,  1.0/36.0},
                                      {     0.0, 1.0/36.0,      0.0}},

                                     {{1.0/36.0, 1.0/9.0, 1.0/36.0},
                                      {1.0/9.0,  4.0/9.0, 1.0/9.0},
                                      {1.0/36.0, 1.0/9.0, 1.0/36.0}},

                                     {{     0.0, 1.0/36.0,     0.0},
                                      {1.0/36.0, 1.0/9.0, 1.0/36.0},
                                      {     0.0, 1.0/36.0,     0.0}}};          ///< D2Q9 lattice Boltzmann stencil

const double LBStencil3D[3][3][3] = {{{1.0/216.0, 1.0/54.0, 1.0/216.0},
                                      {1.0/54.0,  2.0/27.0, 1.0/54.0},
                                      {1.0/216.0, 1.0/54.0, 1.0/216.0}},

                                     {{1.0/54.0, 2.0/27.0, 1.0/54.0},
                                      {2.0/27.0, 8.0/27.0, 2.0/27.0},
                                      {1.0/54.0, 2.0/27.0, 1.0/54.0}},

                                     {{1.0/216.0, 1.0/54.0, 1.0/216.0},
                                      {1.0/54.0,  2.0/27.0, 1.0/54.0},
                                      {1.0/216.0, 1.0/54.0, 1.0/216.0}}};       ///< D3Q27 lattice Boltzmann stencil

class OP_EXPORTS D3Q27                                                                     ///< Lattice Boltzmann populations storage/manipulator
{
 public:
    D3Q27()
    {
        set_to_zero();
    };

    D3Q27(const D3Q27& rhs)
    {
        storage = rhs.storage;
    };

    double& operator()(int x, int y, int z)
    {
        assert(std::abs(x) < 2);
        assert(std::abs(y) < 2);
        assert(std::abs(z) < 2);

        return storage[idx(x,y,z)];
    };

    double operator()(int x, int y, int z) const
    {
        assert(std::abs(x) < 2);
        assert(std::abs(y) < 2);
        assert(std::abs(z) < 2);

        return storage[idx(x,y,z)];
    };

    D3Q27& operator=(const D3Q27& rhs)
    {
        storage = rhs.storage;
        return *this;
    };

    void set_to_zero()
    {
        storage.fill(0.0);
    };

    D3Q27 operator* (const double value) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,y,z)] * value;
        }
        return locPopulations;
    }

    D3Q27& operator*=(const double value)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[idx(x,y,z)] *= value;
        }
        return *this;
    }

    D3Q27 operator/ (const double value) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,y,z)] / value;
        }
        return locPopulations;
    }

    D3Q27& operator/=(const double value)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[idx(x,y,z)] /= value;
        }
        return *this;
    }

    D3Q27 operator+ (const D3Q27& rhs) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,y,z)] + rhs(x,y,z);
        }
        return locPopulations;
    };

    D3Q27& operator+=(const D3Q27& rhs)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[idx(x,y,z)] += rhs(x,y,z);
        }
        return *this;
    };

    D3Q27 operator- (const D3Q27& rhs) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,y,z)] - rhs(x,y,z);
        }
        return locPopulations;
    };

    D3Q27& operator-=(const D3Q27& rhs)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[idx(x,y,z)] -= rhs(x,y,z);
        }
        return *this;
    };

    D3Q27 inverted(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(-x,-y,-z)];
        }
        return locPopulations;
    };

    D3Q27 inverted(int x, int y, int z)
    {
        D3Q27 locPopulations = *this;
        locPopulations(x,y,z) = storage[idx(-x,-y,-z)];

        return locPopulations;
    };

    D3Q27 Xreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(-x,y,z)];
        }
        return locPopulations;
    };

    D3Q27 Yreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,-y,z)];
        }
        return locPopulations;
    };

    D3Q27 Zreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[idx(x,y,-z)];
        }
        return locPopulations;
    };

    double* data()
    {
        return storage.data();
    };

    const double* const_data() const
    {
        return storage.data();
    };

    /// Appends data to a data buffer which can called by a storage class e.g.
    /// Storage3D for MPI-node communication
    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 27; ++i)
        {
            buffer.push_back(storage[i]);
        }
    };

    /// Extracts data from a data buffer which can called by a storage class e.g.
    /// Storage3D for MPI-node communication
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 27; ++i)
        {
            storage[i] = buffer[it];
            ++it;
        }
    };

 protected:
    std::array<double, 27> storage;

 private:
    size_t idx(int x, int y, int z) const
    {
        size_t index = 13 + 9*x + 3*y + z;
        assert(index < 27);
        return index;
    }
};

} //namespace openphase
#endif
