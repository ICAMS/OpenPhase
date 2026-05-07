/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 */

#ifndef STORAGE3D_H
#define STORAGE3D_H

#include <array>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef MPI_PARALLEL
#include "mpi_wrapper.h"
#endif

#include "Macros.h"
#include "GridParameters.h"
#include "Tensor.h"
#include "TypeTraits.h"
#include "Storage3DTensor.h"

namespace openphase
{
// Auxiliary method which helps detecting memory leaks
/*static double getMemoryUsageMB()
{
    std::ifstream file("/proc/self/status");
    std::string line;
    while (getline(file, line))
    {
        if (line.substr(0, 6) == "VmRSS:")
        {
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stol(value) / 1024.0;  // Convert kB to MB
        }
    }
    return -1.0;
}*/

template <typename A, typename T, typename = void>
class ClearCaller3D;

template <typename A, typename T>
class ClearCaller3D<A, T, typename std::enable_if<has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        for (size_t i = 0; i < self.size(); ++i)
        {
            self[i].clear();
        }
        return T();
    }
};

template <typename A, typename T>
class ClearCaller3D<A, T, typename std::enable_if<!has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        std::cerr << "ERROR: Storage3D.h: ClearCaller::Clear() called for non-POD with no clear() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <typename A, typename T, typename = void>
class PackCaller3D;

template <typename A, typename T>
class PackCaller3D<A, T, typename std::enable_if<!std::is_class<T>::value && std::is_standard_layout<T>::value && std::is_trivial<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            buffer.push_back(self(i,j,k));
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            self(i,j,k) = buffer[it]; ++it;
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller3D<A, T, typename std::enable_if<std::is_class<T>::value && has_pack<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            self(i,j,k).pack(buffer);
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            self(i,j,k).unpack(buffer, it);
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller3D<A, T, typename std::enable_if<std::is_class<T>::value && !has_pack<T>::value>::type>
{
 public:
    static T call(A& self, std::vector<long int> window)
    {
        std::stringstream message;
        std::cerr << "ERROR: Storage3D.h: PackCaller::pack() called for non-POD with no pack() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <class T>
class Storage3D<T,0>                                                            ///< 3D storage template class specification for Rank == 0. Can handle any types of data.
{
 public:
    friend class ClearCaller3D< Storage3D<T,0>, T>;
    friend class PackCaller3D< Storage3D<T,0>, T>;

    ~Storage3D()                                                                ///< Destructor.
    {

    }

    Storage3D()                                                                 ///< Constructor. Creates empty storage object.
    {
        DX = 0;
        DY = 0;
        DZ = 0;

        b_cells = 0;

        Size_X = 0;
        Size_Y = 0;
        Size_Z = 0;

        Size_X_BC = 0;
        Size_Y_BC = 0;
        Size_Z_BC = 0;

        Total_X = 0;
        Total_Y = 0;
        Total_Z = 0;

        Offset_X = 0;
        Offset_Y = 0;
        Offset_Z = 0;
    }

    Storage3D(const Storage3D<T,0>& rhs)                                        ///< Copy constructor. Creates the storage object and initializes it with the content of rhs.
    {
        if (rhs.IsAllocated())
        {
            DX = rhs.DX;
            DY = rhs.DY;
            DZ = rhs.DZ;

            Size_X = rhs.Size_X;
            Size_Y = rhs.Size_Y;
            Size_Z = rhs.Size_Z;

            Total_X = rhs.Total_X;
            Total_Y = rhs.Total_Y;
            Total_Z = rhs.Total_Z;

            Offset_X = rhs.Offset_X;
            Offset_Y = rhs.Offset_Y;
            Offset_Z = rhs.Offset_Z;

            b_cells = rhs.b_cells;

            Size_X_BC = Size_X + 2*b_cells*DX;
            Size_Y_BC = Size_Y + 2*b_cells*DY;
            Size_Z_BC = Size_Z + 2*b_cells*DZ;
            const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

            locData.resize(Size);

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,b_cells,)
            {
                locData[Index(i,j,k)] = rhs(i,j,k);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        }
        else
        {
            locData.resize(0);
        }
    }

    Storage3D(const long int nx, const long int ny, const long int nz,
              const long int dnx, const long int dny, const long int dnz,
              const long int bc)                                                ///< Constructor. Creates and allocates memory for the storage of the given size.
    {
        DX = dnx;
        DY = dny;
        DZ = dnz;

        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Total_X = nx;
        Total_Y = ny;
        Total_Z = nz;

        Offset_X = 0;
        Offset_Y = 0;
        Offset_Z = 0;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        locData.resize(Size);
    }

    Storage3D(const GridParameters Dimensions, const long int bc)               ///< Constructor. Creates and allocates memory for the storage of the given size.
    {
        DX = Dimensions.dNx;
        DY = Dimensions.dNy;
        DZ = Dimensions.dNz;

        Size_X = Dimensions.Nx*DX + 1 - DX;
        Size_Y = Dimensions.Ny*DY + 1 - DY;
        Size_Z = Dimensions.Nz*DZ + 1 - DZ;

        Total_X = Dimensions.TotalNx;
        Total_Y = Dimensions.TotalNy;
        Total_Z = Dimensions.TotalNz;

        Offset_X = Dimensions.OffsetX;
        Offset_Y = Dimensions.OffsetY;
        Offset_Z = Dimensions.OffsetZ;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        locData.resize(Size);
    }

    void Allocate(const Storage3D<T,0>& rhs)                                    ///< Allocates previously created empty object to the dimensions of the rhs.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of allocating of a nonempty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                DX = rhs.DX;
                DY = rhs.DY;
                DZ = rhs.DZ;

                Size_X = rhs.Size_X;
                Size_Y = rhs.Size_Y;
                Size_Z = rhs.Size_Z;

                Total_X = rhs.Total_X;
                Total_Y = rhs.Total_Y;
                Total_Z = rhs.Total_Z;

                Offset_X = rhs.Offset_X;
                Offset_Y = rhs.Offset_Y;
                Offset_Z = rhs.Offset_Z;

                b_cells = rhs.b_cells;

                Size_X_BC = Size_X + 2*b_cells*DX;
                Size_Y_BC = Size_Y + 2*b_cells*DY;
                Size_Z_BC = Size_Z + 2*b_cells*DZ;
                const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

                locData.resize(Size);
            }
        }
    }

    void AllocateCopy(const Storage3D<T,0>& rhs)                                ///< Allocates previously created empty object to the dimensions of the rhs and initializes it with rhs's content.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                DX = rhs.DX;
                DY = rhs.DY;
                DZ = rhs.DZ;

                Size_X = rhs.Size_X;
                Size_Y = rhs.Size_Y;
                Size_Z = rhs.Size_Z;

                Total_X = rhs.Total_X;
                Total_Y = rhs.Total_Y;
                Total_Z = rhs.Total_Z;

                Offset_X = rhs.Offset_X;
                Offset_Y = rhs.Offset_Y;
                Offset_Z = rhs.Offset_Z;

                b_cells = rhs.b_cells;

                Size_X_BC = Size_X + 2*b_cells*DX;
                Size_Y_BC = Size_Y + 2*b_cells*DY;
                Size_Z_BC = Size_Z + 2*b_cells*DZ;
                const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

                locData.resize(Size);

                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,b_cells,)
                {
                    locData[Index(i,j,k)] = rhs(i,j,k);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                locData.resize(0);
            }
        }
    }

    T& operator()(const long int x, const long int y, const long int z)         ///< Random access operator. Returns the reference to the data stored in a given grid location.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locData[Index(x,y,z)];
    }

    T const& operator()(const long int x, const long int y, const long int z) const///< Random access operator. Returns const reference to the data stored in a given grid location.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locData[Index(x,y,z)];
    }

    T at(const double x, const double y, const double z) const                  ///< Returns the interpolated data in the given location between the grid points.
    {
        long int x0 = floor(x)*DX;
        long int y0 = floor(y)*DY;
        long int z0 = floor(z)*DZ;
        double dx = fabs(x - x0)*DX;
        double dy = fabs(y - y0)*DY;
        double dz = fabs(z - z0)*DZ;

        T tempValue = locData[Index(x0, y0, z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz));

        if(DX) tempValue += locData[Index(x0+1, y0  , z0  )]*(dx*(1.0 - dy)*(1.0 - dz));
        if(DY) tempValue += locData[Index(x0  , y0+1, z0  )]*((1.0 - dx)*dy*(1.0 - dz));
        if(DZ) tempValue += locData[Index(x0  , y0  , z0+1)]*((1.0 - dx)*(1.0 - dy)*dz);

        if(DX and DY) tempValue += locData[Index(x0+1, y0+1, z0  )]*(dx*dy*(1.0 - dz));
        if(DX and DZ) tempValue += locData[Index(x0+1, y0  , z0+1)]*(dx*(1.0 - dy)*dz);
        if(DY and DZ) tempValue += locData[Index(x0  , y0+1, z0+1)]*((1.0 - dx)*dy*dz);

        if(DX and DY and DZ) tempValue += locData[Index(x0+1, y0+1, z0+1)]*(dx*dy*dz);

        return tempValue;
    }

    T& operator[](size_t idx)                                                   ///< Random access operator. Returns the reference to the data stored in a given flattened storage location.
    {
        assert(idx < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locData[idx];
    }

    T const& operator[](const size_t idx) const                                 ///< Random access operator. Returns const reference to the data stored in a given flattened storage location.
    {
        assert(idx < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locData[idx];
    }

    size_t Allocate(const long int nx, const long int ny, const long int nz,
                    const long int dnx, const long int dny, const long int dnz,
                    const long int bc)                                          ///< Allocates previously created empty object to the given dimensions.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of reallocation of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        DX = dnx;
        DY = dny;
        DZ = dnz;

        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Total_X = nx;
        Total_Y = ny;
        Total_Z = nz;

        Offset_X = 0;
        Offset_Y = 0;
        Offset_Z = 0;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;
        const size_t AllocatedMemory = sizeof(T)*Size;

        locData.resize(Size);
        return AllocatedMemory;
    }

    size_t Allocate(const GridParameters Dimensions, const long int bc)         ///< Allocates previously created empty object to the given dimensions.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of reallocation of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        DX = Dimensions.dNx;
        DY = Dimensions.dNy;
        DZ = Dimensions.dNz;

        Size_X = Dimensions.Nx*DX + 1 - DX;
        Size_Y = Dimensions.Ny*DY + 1 - DY;
        Size_Z = Dimensions.Nz*DZ + 1 - DZ;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;
        const size_t AllocatedMemory = sizeof(T)*Size;

        locData.resize(Size);
        return AllocatedMemory;
    }

    void Reallocate(const long int nx, const long int ny, const long int nz)    ///< Reallocates memory of already allocated container to new dimensions. Keeps the number of boundary cells the same. Scrambles the stored data.
    {
        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        locData.resize(Size);
    }

    bool IsNotAllocated() const                                                 ///< Returns true if current container is not allocated.
    {
        return locData.size() == 0;
    }

    bool IsAllocated() const                                                    ///< Returns true if current container is allocated.
    {
        return locData.size() != 0;
    }

    bool IsSize(const long int Nx, const long int Ny, const long int Nz)        ///< Returns true if current container has given dimensions.
    {
        return (Size_X == Nx and Size_Y == Ny and Size_Z == Nz);
    }

    long int ActiveDimensions() const                                           ///< Returns the number of active dimensions: 3 for all three dimensions active, 2 for any two dimensions active, and 1 for a single dimension.
    {
        return DX + DY + DZ;
    }

    void Remesh(const long int nX, const long int nY, const long int nZ)        ///< Changes the dimensions of the container while keeping the data. Uses linear interpolation to recover the data at the new grid locations.
    {
        long int nx = nX*DX + 1 - DX;
        long int ny = nY*DY + 1 - DY;
        long int nz = nZ*DZ + 1 - DZ;

        std::vector<T> tempArray((nx + 2*b_cells*DX)*(ny + 2*b_cells*DY)*(nz + 2*b_cells*DZ));

        double Xscale = double(Size_X)/double(nx);
        double Yscale = double(Size_Y)/double(ny);
        double Zscale = double(Size_Z)/double(nz);

        #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE, OMP_CHUNKSIZE)
        for(long int x = 0; x < nx; x++)
        for(long int y = 0; y < ny; y++)
        for(long int z = 0; z < nz; z++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5)*DX;
            long int y0 = floor((y - ny*0.5)*Yscale + Size_Y * 0.5)*DY;
            long int z0 = floor((z - nz*0.5)*Zscale + Size_Z * 0.5)*DZ;

            double dx = (x*Xscale - x0)*DX;
            double dy = (y*Yscale - y0)*DY;
            double dz = (z*Zscale - z0)*DZ;

            T tmpData = locData[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz));

            if(DX) tmpData += locData[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz));
            if(DY) tmpData += locData[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz));
            if(DZ) tmpData += locData[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz);

            if(DX and DY) tmpData += locData[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz));
            if(DX and DZ) tmpData += locData[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz);
            if(DY and DZ) tmpData += locData[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz);

            if(DX and DY and DZ) tmpData += locData[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

            tempArray[(((x + b_cells*DX)*(ny + 2*b_cells*DY) + y + b_cells*DY)*(nz + 2*b_cells*DZ) + z + b_cells*DZ)] = tmpData;
        }

        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;

        locData = tempArray;
    }

    bool rotate(const long int newdimx, const long int newdimy, const long int newdimz)///< Switching storage dimensions, e.g. changes storage diemsnions arrangement.
    {
        if ((newdimx == newdimy) or (newdimx == newdimz) or (newdimy == newdimz))
        {
            std::cerr << "flip_storage: Double rotation axis.";
            return false;
        }
        if (newdimx < 0 or newdimx > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimx.";
            return false;
        }
        if (newdimy < 0 or newdimy > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimy.";
            return false;
        }
        if (newdimz < 0 or newdimz > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimz.";
            return false;
        }

        long int newdim[3] = {newdimx, newdimy, newdimz};
        long int oldSize[3] = {Size_X, Size_Y, Size_Z};

        // newdimx == 0; newdimy == 1; newdimz == 2
        if (newdimx == 0 and newdimy == 1 and newdimz == 2) return false;

        std::vector<T> tempArray((Size_X + 2*b_cells*DX)*(Size_Y + 2*b_cells*DY)*(Size_Z + 2*b_cells*DZ));

        // dimx == 0; dimy == 2; dimz == 1 (rotates around positive x)
        if (newdimx == 0 and newdimy == 2 and newdimz == 1)
        {
            for (long int i = 0; i < Size_X; i++)
            for (long int j = 0; j < Size_Y; j++)
            for (long int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_Z_BC*(i + b_cells*DX) + Size_Z - 1 - k + b_cells*DZ)*Size_Y_BC + j + b_cells*DY]
                          = locData[Index(i,j,k)];
            }
        }
        // dimx == 2; dimy == 1; dimz == 0 (rotates around positive y)
        if (newdimx == 2 and newdimy == 1 and newdimz == 0)
        {
            for (long int i = 0; i < Size_X; i++)
            for (long int j = 0; j < Size_Y; j++)
            for (long int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_Y_BC*(k + b_cells*DZ) + j + b_cells*DY)*Size_X_BC + (Size_X - i - 1) + b_cells*DX]
                          = locData[Index(i,j,k)];
            }
        }
        // dimx == 1; dimy == 0; dimz == 2 (rotates around positive z)
        if (newdimx == 1 and newdimy == 0 and newdimz == 2)
        {
            for (long int i = 0; i < Size_X; i++)
            for (long int j = 0; j < Size_Y; j++)
            for (long int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_X_BC*(j + b_cells*DY) + i + b_cells*DX)*Size_Z_BC + k + b_cells*DZ]
                          = locData[Index(i,j,k)];
            }
        }

        Size_X = oldSize[newdim[0]];
        Size_Y = oldSize[newdim[1]];
        Size_Z = oldSize[newdim[2]];

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;

        locData = tempArray;

        return true;
    }

    Storage3D<T,0>& operator=(const Storage3D<T,0>& rhs)                        ///< Assignment operator. Assigns the copy of the rhs to the current container.
    {
        if(locData.size() != 0 and rhs.IsAllocated())
        {
            if (rhs.Size_X_BC == Size_X_BC and
                rhs.Size_Y_BC == Size_Y_BC and
                rhs.Size_Z_BC == Size_Z_BC)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,rhs.Bcells(),)
                {
                    locData[Index(i,j,k)] = rhs(i,j,k);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                if (rhs.IsAllocated())
                {
                    DX = rhs.DX;
                    DY = rhs.DY;
                    DZ = rhs.DZ;

                    Size_X = rhs.Size_X;
                    Size_Y = rhs.Size_Y;
                    Size_Z = rhs.Size_Z;

                    Total_X = rhs.Total_X;
                    Total_Y = rhs.Total_Y;
                    Total_Z = rhs.Total_Z;

                    Offset_X = rhs.Offset_X;
                    Offset_Y = rhs.Offset_Y;
                    Offset_Z = rhs.Offset_Z;

                    b_cells = rhs.b_cells;

                    Size_X_BC = Size_X + 2*b_cells*DX;
                    Size_Y_BC = Size_Y + 2*b_cells*DY;
                    Size_Z_BC = Size_Z + 2*b_cells*DZ;
                    const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

                    locData.resize(Size);

                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,b_cells,)
                    {
                        locData[Index(i,j,k)] = rhs(i,j,k);
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                }
            }
        }
        else if (locData.size() == 0 and rhs.IsAllocated())
        {
            AllocateCopy(rhs);
        }
        return *this;
    }

    void Clear(void)                                                            ///< Clears the data while keeping the allocated storage dimension unchanged.
    {
        if(locData.size() != 0)
        {
            if constexpr (std::is_arithmetic<T>::value)
            {
                std::fill(locData.begin(), locData.end(), 0);
            }
            else
            {
                ClearCaller3D< Storage3D<T,0> ,T> K;
                K.call(*this);
            }
        }
    }

    std::vector<double> pack(std::vector<long int> window)                      ///< Packs (writes) selected container content to the MPI communication buffer.
    {
        std::vector<double> buffer;
        PackCaller3D< Storage3D<T,0>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer, std::vector<long int> window)      ///< Unpacks (reads) selected container content from the MPI communication buffer.
    {
        PackCaller3D< Storage3D<T,0>, T> K;
        K.unpack(*this, buffer, window);
    }
    
    void WriteToFile(std::string FileName)                                      ///< Writes container content to the file in MPI-parallel mode.
    {
        #ifdef MPI_PARALLEL
        std::vector<std::string> sbuffer(Size_Y*Size_Z);
        for (long int j = 0; j < Size_Y; j++)
        for (long int k = 0; k < Size_Z; k++)
        {
            std::stringstream ssbuffer;
            std::vector<double> buffer;
            std::vector<long int> window(6);
            window[0] = 0;
            window[1] = Size_X;
            window[2] = j;
            window[3] = j+1;
            window[4] = k;
            window[5] = k+1;
            PackCaller3D< Storage3D<T,0>, T> K;
            K.pack(*this, buffer, window);
            for (size_t i = 0; i < buffer.size(); ++i)
            {
                ssbuffer << std::any_cast<double>(buffer[i]);
            }
            sbuffer[k+Size_Z*j] = ssbuffer.str();
        }
        op_mpi_write_data(FileName, sbuffer, Size_X, Size_Y, Size_Z);
        #endif
    }

    std::vector<double> pack()                                                  ///< Packs (writes) the entire container content to the MPI communication buffer.
    {
        std::vector<double> buffer;
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        window[2] = 0;
        window[3] = Size_Y;
        window[4] = 0;
        window[5] = Size_Z;
        PackCaller3D< Storage3D<T,0>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer)                                    ///< Unpacks (reads) the entire container content from the MPI communication buffer.
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        window[2] = 0;
        window[3] = Size_Y;
        window[4] = 0;
        window[5] = Size_Z;
        PackCaller3D< Storage3D<T,0>, T> K;
        K.unpack(*this, buffer, window);
    }

    void set_to_value(T val)                                                    ///< Sets all entries in the container to the given value.
    {
        if(locData.size() != 0)
        {
            const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;
            for(size_t idx = 0; idx < Size; idx++)
            {
                locData[idx] = val;
            }
        }
        else
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name()
                      << ", 0>::set_to_value() operation on the nonallocated storage!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
    }

    T* data(void)
    {
        return locData.data();
    }

    T* data(void) const
    {
        return locData.data();
    }

    long int sizeX() const                                                      ///< Returns X dimension of the storage.
    {
        return Size_X;
    }

    long int sizeY() const                                                      ///< Returns Y dimension of the storage.
    {
        return Size_Y;
    }

    long int sizeZ() const                                                      ///< Returns Z dimension of the storage.
    {
        return Size_Z;
    }

    long int Bcells() const                                                     ///< Returns the number of boundary cells
    {
        return b_cells;
    }

    long int BcellsX() const                                                    ///< Returns the number of boundary cells in X directions
    {
        return b_cells*DX;
    }

    long int BcellsY() const                                                    ///< Returns the number of boundary cells in Y directions
    {
        return b_cells*DY;
    }

    long int BcellsZ() const                                                    ///< Returns the number of boundary cells in Z directions
    {
        return b_cells*DZ;
    }

    long int dNx() const                                                        ///< Returns 1 is X dimension is active, 0 otherwise.
    {
        return DX;
    }

    long int dNy() const                                                        ///< Returns 1 is Y dimension is active, 0 otherwise.
    {
        return DY;
    }

    long int dNz() const                                                        ///< Returns 1 is Z dimension is active, 0 otherwise.
    {
        return DZ;
    }

    size_t size() const                                                         ///< Returns the size of the internal storage
    {
        return Size_X_BC*Size_Y_BC*Size_Z_BC;
    }

    void SetNewBcells(const int new_b_cells)                                    ///< Sets new number of boundary cells. Destroys the stored data.
    {
        b_cells = new_b_cells;

        Reallocate(Size_X, Size_Y, Size_Z);
    }

    ////NOTE: Method is currently not used in OpenPhase and the << operator is
    ////not defined for all OpenPhase file types (Matrix3x3)!
    /*std::string print(void) const
    {
        std::stringstream out;
        for (long int k = 0; k < Size_Z; ++k)
        for (long int j = 0; j < Size_Y; ++j)
        for (long int i = 0; i < Size_X; ++i)
        {
            out << " " << locData[Index(i,j,k)] << "\n";
        }
        return out.str();
    };*/

    ////NOTE: Method is currently not used in OpenPhase and the << operator is
    ////not defined for all OpenPhase file types (Matrix3x3)!
    /*std::string print_slice(const std::string dim, const long int scoord) const
    {
        std::stringstream out;

        if (dim == "z" or dim == "Z")
        if (scoord <= Size_Z)
        {
            for (long int j = 0; j < Size_Y; ++j)
            {
                for (long int i = 0; i < Size_X; ++i)
                {
                    out << " " << locData[Index(i,j,scoord)];
                }
                out << "\n";
            }
        }
        if (dim == "y" or dim == "Y")
        if (scoord <= Size_Y)
        {
            for (long int k = 0; k < Size_Z; ++k)
            {
                for (long int i = 0; i < Size_X; ++i)
                {
                    out << " " << locData[Index(i,scoord,k)];
                }
                out << "\n";
            }
        }
        if (dim == "x" or dim == "X")
        if (scoord <= Size_X)
        {
            for (long int k = 0; k < Size_Z; ++k)
            {
                for (long int j = 0; j < Size_Y; ++j)
                {
                    out << " " << locData[Index(scoord,j,k)];
                }
                out << "\n";
            }
        }
        return out.str();
    };*/

    ////Raphael: Method is currently not used in OpenPhase and the << operator is
    ////not defined for all OpenPhase file types (Matrix3x3)!
    /*std::string print_subset(const long int xmin, const long int xmax,
                             const long int ymin, const long int ymax,
                             const long int zmin, const long int zmax) const
    {
        std::stringstream out;
        if(xmin < -b_cells*DX or xmin > xmax or xmin > sizeX() + b_cells*DX) {out << "Storage3D.print() - bad index"; return out.str();}
        if(ymin < -b_cells*DY or ymin > ymax or ymin > sizeY() + b_cells*DY) {out << "Storage3D.print() - bad index"; return out.str();}
        if(zmin < -b_cells*DZ or zmin > xmax or zmin > sizeZ() + b_cells*DZ) {out << "Storage3D.print() - bad index"; return out.str();}
        if(xmax < -b_cells*DX or xmax < xmin or xmax > sizeX() + b_cells*DX) {out << "Storage3D.print() - bad index"; return out.str();}
        if(ymax < -b_cells*DY or ymax < ymin or ymax > sizeY() + b_cells*DY) {out << "Storage3D.print() - bad index"; return out.str();}
        if(zmax < -b_cells*DZ or zmax < zmin or zmax > sizeZ() + b_cells*DZ) {out << "Storage3D.print() - bad index"; return out.str();}

        for (long int k = zmin; k < zmax; ++k)
        {
            for (long int j = ymin; j < ymax; ++j)
            {
                for (long int i = xmin; i < xmax; ++i)
                {
                    out << " " << locData[Index(i,j,k)];
                }
                out << "\n";
            }
            out << "\n";
        }
        return out.str();
    };*/

    bool InLimits(const long int x, const long int y, const long int z)         ///< Returns true if the given location is within the dimensions limits of the current container.
    {
        if (x < -b_cells*DX) return false;
        if (y < -b_cells*DY) return false;
        if (z < -b_cells*DZ) return false;
        if (x >= Size_X + b_cells*DX) return false;
        if (y >= Size_Y + b_cells*DY) return false;
        if (z >= Size_Z + b_cells*DZ) return false;
        return true;
    }

 protected:
 private:
    long int Size_X;
    long int Size_Y;
    long int Size_Z;

    long int Size_X_BC;
    long int Size_Y_BC;
    long int Size_Z_BC;
    
    long int Total_X;
    long int Total_Y;
    long int Total_Z;
    
    long int Offset_X;
    long int Offset_Y;
    long int Offset_Z;

    long int b_cells;

    long int DX;
    long int DY;
    long int DZ;

    std::vector<T> locData;

    size_t Index(const long int x, const long int y, const long int z) const    ///< Converts 3D coordinates into flattened storage index.
    {
        assert(x < Size_X + b_cells*DX && "Access beyond storage range");
        assert(y < Size_Y + b_cells*DY && "Access beyond storage range");
        assert(z < Size_Z + b_cells*DZ && "Access beyond storage range");
        assert(x >= - b_cells*DX && "Access beyond storage range");
        assert(y >= - b_cells*DY && "Access beyond storage range");
        assert(z >= - b_cells*DZ && "Access beyond storage range");

        return ((Size_Y_BC*(x + b_cells*DX) + y + b_cells*DY)*Size_Z_BC + z + b_cells*DZ);
    }
};

}// namespace openphase
#endif
