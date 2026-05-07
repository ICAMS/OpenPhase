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

#ifndef STORAGE3DTENSOR_H
#define STORAGE3DTENSOR_H

#include <array>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef MPI_PARALLEL
#include "mpi_wrapper.h"
#endif

#include "Macros.h"
#include "Tensor.h"
#include "GridParameters.h"
#include "TypeTraits.h"

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
class ClearCallerTensor3D;

template <typename A, typename T>
class ClearCallerTensor3D<A, T, typename std::enable_if<has_clear<T>::value>::type>
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
class ClearCallerTensor3D<A, T, typename std::enable_if<!has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        std::cerr << "ERROR: Storage3DTensor.h: ClearCaller::Clear() called for non-POD with no clear() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <typename A, typename T, typename = void>
class PackCallerTensor3D;

template <typename A, typename T>
class PackCallerTensor3D<A, T, typename std::enable_if<!std::is_class<T>::value && std::is_standard_layout<T>::value && std::is_trivial<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            for (size_t n = 0; n < self(i,j,k).size(); ++n)
            {
                buffer.push_back(self(i,j,k)[n]);
            }
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
            for (size_t n = 0; n < self(i,j,k).size(); ++n)
            {
                self(i,j,k)[n] = buffer[it]; ++it;
            }
        }
        return T();
    }
};

template <typename A, typename T>
class PackCallerTensor3D<A, T, typename std::enable_if<std::is_class<T>::value && has_pack<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            for (size_t n = 0; n < self(i,j,k).size(); ++n)
            {
                self(i,j,k)[n].pack(buffer);
            }
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t  it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        for (long int j = window[2]; j < window[3]; ++j)
        for (long int k = window[4]; k < window[5]; ++k)
        {
            for (size_t n = 0; n < self(i,j,k).size(); ++n)
            {
                self(i,j,k)[n].unpack(buffer, it);
            }
        }
        return T();
    }
};

template <typename A, typename T>
class PackCallerTensor3D<A, T, typename std::enable_if<std::is_class<T>::value && !has_pack<T>::value>::type>
{
 public:
    static T call(A& self, std::vector<long int> window)
    {
        std::cerr << "ERROR: Storage3DTensor.h: PackCallerTensor::pack() called for non-POD with no pack() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <class T, size_t Rank>
class Storage3D                                                                 ///< 3D storage template class. Can handle any types of data.
{
 public:
    friend class ClearCallerTensor3D< Storage3D<T, Rank> , T>;
    friend class PackCallerTensor3D< Storage3D<T, Rank> , T>;

    ~Storage3D()                                                                ///< Destructor.
    {

    }

    Storage3D()                                                                 ///< Constructor. Creates empty storage object.
    {
        DX = 0;
        DY = 0;
        DZ = 0;

        Size_X = 0;
        Size_Y = 0;
        Size_Z = 0;

        Total_X = 0;
        Total_Y = 0;
        Total_Z = 0;

        Offset_X = 0;
        Offset_Y = 0;
        Offset_Z = 0;

        b_cells = 0;

        Size_X_BC = 0;
        Size_Y_BC = 0;
        Size_Z_BC = 0;
        Size_D    = 0;

        locTensors.resize(0);
        locData.resize(0);
    }

    Storage3D(const long int nx,  const long int ny,  const long int nz,
              const long int dnx, const long int dny, const long int dnz,
              const std::array<size_t,Rank> nn,
              const long int bc)                                                ///< Constructor. Creates and allocates memory for the storage of the given size.
    {
        TensorDimensions = nn;

        DX = dnx;
        DY = dny;
        DZ = dnz;

        Size_X = nx*DX + 1 - DX;
        Size_Y = ny*DY + 1 - DY;
        Size_Z = nz*DZ + 1 - DZ;

        Total_X = Size_X;
        Total_Y = Size_Y;
        Total_Z = Size_Z;

        Offset_X = 0;
        Offset_Y = 0;
        Offset_Z = 0;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData.resize(Size*Size_D);

        locTensors.resize(Size);

        for(size_t i = 0; i < Size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    Storage3D(const GridParameters Dimensions,
              const std::array<size_t,Rank> nn,
              const long int bc)                                                ///< Constructor. Creates and allocates memory for the storage of the given size.
    {
        TensorDimensions = nn;

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

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData.resize(Size*Size_D);
        locTensors.resize(Size);

        for(size_t i = 0; i < Size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    Storage3D(const Storage3D<T,Rank>& rhs)                                     ///< Copy constructor. Creates the storage object and initializes it with the content of rhs.
    {
        if (rhs.IsAllocated())
        {
            TensorDimensions = rhs.TensorDimensions;

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

            Size_D = 1;
            for(size_t n = 0; n < TensorDimensions.size(); n++)
            {
                Size_D *= TensorDimensions[n];
            }

            locData.resize(Size*Size_D);
            locTensors.resize(Size);

            for(size_t i = 0; i < Size; i++)
            {
                locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
            }

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,0,)
                {
                    locTensors[Index(i,j,k)] = rhs(i,j,k);
                }
            OMP_PARALLEL_STORAGE_LOOP_END
        }
        else
        {
            locData.resize(0);
            locTensors.resize(0);
        }
    }

    Tensor<T, Rank>& operator()(const long int x, const long int y, const long int z) ///< Random access operator. Returns the reference to the tensor object stored in a given grid location.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locTensors[Index(x,y,z)];
    }

    Tensor<T, Rank> const& operator()(const long int x, const long int y,
            const long int z) const                                             ///< Random access operator. Returns const reference to the tensor object stored in a given grid location.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locTensors[Index(x,y,z)];
    }

    T& operator()(const long int x, const long int y, const long int z, const std::array<size_t, Rank> Position)///< Random access operator. Returns the reference to the tensor element stored in a given grid location and position in the tensor.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        assert(IndexT(Position) < size_t(Size_D) && "Access beyond storage range");

        return locData[Index(x,y,z)*Size_D + IndexT(Position)];
    }

    const T& operator()(const long int x, const long int y, const long int z, const std::array<size_t, Rank> Position) const///< Random access operator. Returns const reference to the tensor element stored in a given grid location and position in the tensor.
    {
        assert(Index(x,y,z) < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        assert(IndexT(Position) < size_t(Size_D) && "Access beyond storage range");

        return locData[Index(x,y,z)*Size_D + IndexT(Position)];
    }

    Tensor<T, Rank> at(const double x, const double y, const double z) const    ///< Returns the interpolated tensor in the given location between the grid points.
    {
        long int x0 = floor(x)*DX;
        long int y0 = floor(y)*DY;
        long int z0 = floor(z)*DZ;
        double dx = fabs(x - x0)*DX;
        double dy = fabs(y - y0)*DY;
        double dz = fabs(z - z0)*DZ;

        Tensor<T, Rank> tempTensor = locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz));

        if(dx) tempTensor += locTensors[Index(x0+DX,y0   ,z0   )]*(dx*(1.0 - dy)*(1.0 - dz));
        if(dy) tempTensor += locTensors[Index(x0   ,y0+DY,z0   )]*((1.0 - dx)*dy*(1.0 - dz));
        if(dz) tempTensor += locTensors[Index(x0   ,y0   ,z0+DZ)]*((1.0 - dx)*(1.0 - dy)*dz);

        if(dx and dy) tempTensor += locTensors[Index(x0+DX,y0+DY,z0   )]*(dx*dy*(1.0 - dz));
        if(dx and dz) tempTensor += locTensors[Index(x0+DX,y0   ,z0+DZ)]*(dx*(1.0 - dy)*dz);
        if(dy and dz) tempTensor += locTensors[Index(x0   ,y0+DY,z0+DZ)]*((1.0 - dx)*dy*dz);

        if(dx and dy and dz) tempTensor += locTensors[Index(x0+DX,y0+DY,z0+DZ)]*(dx*dy*dz);

        return tempTensor;
    }

    Tensor<T, Rank>& operator[](const size_t idx)                               ///< Random access operator. Returns the reference to the tensor object stored in a given flattened storage location.
    {
        assert(idx < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locTensors[idx];
    }

    Tensor<T, Rank>const& operator[](const size_t idx) const                    ///< Random access operator. Returns const reference to the tensor object stored in a given flattened storage location.
    {
        assert(idx < size_t(Size_X_BC*Size_Y_BC*Size_Z_BC) && "Access beyond storage range");
        return locTensors[idx];
    }

    void Allocate(const long int nx, const long int ny, const long int nz,
                  const long int dnx, const long int dny, const long int dnz,
                  const std::array<size_t,Rank> nn, const long int bc)          ///< Allocates previously created empty object to the given dimensions.
    {
        if(locTensors.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank
                      << ">::Allocate()\n"
                      << "Attempt of allocating of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        std::copy(nn.begin(), nn.end(), TensorDimensions.begin());

        DX = dnx;
        DY = dny;
        DZ = dnz;

        Size_X = nx*DX + 1 - DX;
        Size_Y = ny*DY + 1 - DY;
        Size_Z = nz*DZ + 1 - DZ;

        b_cells = bc;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;
        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData.resize(Size*Size_D);
        locTensors.resize(Size);

        for(size_t i = 0; i < Size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    void Allocate(const GridParameters Dimensions,
                  const std::array<size_t,Rank> nn, const long int bc)          ///< Allocates previously created empty object to the given dimensions.
    {
        if(locTensors.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank
                      << ">::Allocate()\n"
                      << "Attempt of allocating of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        std::copy(nn.begin(), nn.end(), TensorDimensions.begin());

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
        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData.resize(Size*Size_D);
        locTensors.resize(Size);

        for(size_t i = 0; i < Size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    void Allocate(const Storage3D<T,Rank>& rhs)                                 ///< Allocates previously created empty object to the dimensions of the rhs.
    {
        if(locTensors.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", "
                      <<  Rank << ">::Allocate()\n"
                      << "Attempt of allocating of a nonempty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                TensorDimensions = rhs.TensorDimensions;

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

                Size_D = 1;
                for(size_t n = 0; n < TensorDimensions.size(); n++)
                {
                    Size_D *= TensorDimensions[n];
                }

                locData.resize(Size*Size_D);
                locTensors.resize(Size);

                for(size_t i = 0; i < Size; i++)
                {
                    locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
                }
            }
        }
    }

    void AllocateCopy(const Storage3D<T,Rank>& rhs)                             ///< Allocates previously created empty object to the dimensions of the rhs and initializes it with rhs's content.
    {
        if(locTensors.size() != 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank
                      << ">::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                TensorDimensions = rhs.TensorDimensions;

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

                Size_D = 1;
                for(size_t n = 0; n < TensorDimensions.size(); n++)
                {
                    Size_D *= TensorDimensions[n];
                }

                locData.resize(Size*Size_D);
                locTensors.resize(Size);

                for(size_t i = 0; i < Size; i++)
                {
                    locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
                }

                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,0,)
                    {
                        locTensors[Index(i,j,k)] = rhs(i,j,k);
                    }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
        }
    }

    Storage3D<T,Rank>& operator=(const Storage3D<T,Rank>& rhs)                  ///< Assignment operator. Assigns the copy of the rhs to the current container.
    {
        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                TensorDimensions = rhs.TensorDimensions;

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

                Size_D = 1;
                for(size_t n = 0; n < TensorDimensions.size(); n++)
                {
                    Size_D *= TensorDimensions[n];
                }

                locData.resize(Size*Size_D);
                locTensors.resize(Size);

                for(size_t i = 0; i < Size; i++)
                {
                    locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
                }

                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,rhs,b_cells,)
                    {
                        locTensors[Index(i,j,k)] = rhs(i,j,k);
                    }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                AllocateCopy(rhs);
            }
        }
        return *this;
    }

    void Reallocate(const long int nx, const long int ny, const long int nz)    ///< Reallocates memory of already allocated container to new dimensions. Keeps tensor dimensions and the number of boundary cells the same. Destroys the stored data.
    {
        Size_X = nx*DX + 1 - DX;
        Size_Y = ny*DY + 1 - DY;
        Size_Z = nz*DZ + 1 - DZ;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;

        const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;

        locTensors.clear();
        locTensors.resize(Size);
        locData.clear();
        locData.resize(Size*Size_D);

        for(size_t i = 0; i < Size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    bool IsNotAllocated() const                                                 ///< Returns true if current container is not allocated.
    {
        return locTensors.size() == 0;
    }

    bool IsAllocated() const                                                    ///< Returns true if current container is allocated.
    {
        return locTensors.size() != 0;
    }

    bool IsSize(const long int Nx, const long int Ny, const long int Nz)        ///< Returns true if current container has given dimensions.
    {
        return (Size_X == Nx and Size_Y == Ny and Size_Z == Nz);
    }

    void Remesh(const long int nX, const long int nY, const long int nZ)        ///< Changes the dimensions of the container while keeping the data. Uses linear interpolation to recover the data at the new grid locations.
    {
        long int nx = nX*DX + 1 - DX;
        long int ny = nY*DY + 1 - DY;
        long int nz = nZ*DZ + 1 - DZ;

        size_t new_size = (nx + 2*b_cells*DX)*(ny + 2*b_cells*DY)*(nz + 2*b_cells*DZ);

        std::vector<T> tempData(new_size*Size_D);
        std::vector<Tensor<T, Rank>> tempTensors(new_size);

        for(size_t i = 0; i < new_size; i++)
        {
            tempTensors[i].Assign(&tempData[i*Size_D],TensorDimensions);
        }

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

            Tensor<T, Rank> tempTensor = locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz));

            if(DX) tempTensor += locTensors[Index(x0+1,y0,z0)]*(dx*(1.0 - dy)*(1.0 - dz));
            if(DY) tempTensor += locTensors[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz));
            if(DZ) tempTensor += locTensors[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz);

            if(DX and DY) tempTensor += locTensors[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz));
            if(DX and DZ) tempTensor += locTensors[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz);
            if(DY and DZ) tempTensor += locTensors[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz);

            if(DX and DY and DZ) tempTensor += locTensors[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

            tempTensors[(((ny + 2*b_cells*DY)*(x + b_cells*DX) + y + b_cells*DY)*(nz + 2*b_cells*DZ) + z + b_cells*DZ)] = tempTensor;
        }

        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Size_X_BC = Size_X + 2*b_cells*DX;
        Size_Y_BC = Size_Y + 2*b_cells*DY;
        Size_Z_BC = Size_Z + 2*b_cells*DZ;

        locTensors.clear();
        locTensors.resize(new_size);

        locData = tempData;

        for(size_t i = 0; i < new_size; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    auto* data(void)                                                            ///< Returns the pointer to the raw data storage.
    {
        if constexpr (has_datapointer<T>::value)
        {
            return locData->data();
        }
        else return locData;
    }

    void Clear(void)                                                            ///< Clears the data while keeping the allocated storage dimension unchanged.
    {
        ClearCallerTensor3D< Storage3D<T, Rank>, T> K;
        K.call(*this);
    }

    std::vector<double> pack(std::vector<long int> window)                      ///< Packs (writes) selected container content to the MPI communication buffer.
    {
        std::vector<double> buffer;
        PackCallerTensor3D< Storage3D<T,Rank>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer, std::vector<long int> window)      ///< Unpacks (reads) selected container content from the MPI communication buffer.
    {
        PackCallerTensor3D< Storage3D<T,Rank>, T> K;
        K.unpack(*this, buffer, window);
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
        PackCallerTensor3D< Storage3D<T,Rank>, T> K;
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
        PackCallerTensor3D< Storage3D<T,Rank>, T> K;
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
            window[3] = j;
            window[4] = k;
            window[5] = k;
            PackCallerTensor3D< Storage3D<T,Rank>, T> K;
            K.pack(*this, buffer, window);
            for (long int i = 0; i < buffer.size(); ++i)
            {
                ssbuffer << std::any_cast<double>(buffer[i]);
            }
            sbuffer[k+Size_Z*j] = ssbuffer.str();
        }
        op_mpi_write_data(FileName, sbuffer, Size_X, Size_Y, Size_Z);
        #endif
    }

    void ReadFromFile(std::string FileName)                                     ///< Reads container content from the file
    {
    
    }

    void set_to_value(Tensor<T, Rank>& val)                                     ///< Sets all entries in the container to the given value.
    {
        if(locTensors.size() != 0)
        {
            if(TensorDimensions == val.TensorDimensions)
            {
                const size_t Size = Size_X_BC*Size_Y_BC*Size_Z_BC;
                for(size_t idx = 0; idx < Size; idx++)
                {
                    locTensors[idx] = val;
                }
            }
            else
            {
                std::cerr << "ERROR: Storage3DTensor<" << typeid(T).name() << ", "
                          <<  Rank << ">::set_to_value(): Wrong tensor size!\n"
                          << "Terminating!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
        }
        else
        {
            std::cerr << "ERROR: Storage3DTensor<" << typeid(T).name() << ", " <<  Rank
                      << ">::set_to_value() operation on a not allocated storage!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
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

    constexpr size_t tensor_rank(void) const                                    ///< Returns the rank of the stored tensors.
    {
        return Rank;
    }

    size_t tensor_size(size_t n) const                                          ///< Returns the size of the tensor dimension n.
    {
        if(n < Rank)
        {
            return TensorDimensions[n];
        }
        else
        {
            return 0;
        }
    }

    size_t tensor_size() const                                                  ///< Returns the size of the storage occupied by the single tensor
    {
        return Size_D;
    }

    size_t size() const                                                         ///< Returns the size of the storage of tensors.
    {
        return Size_X_BC*Size_Y_BC*Size_Z_BC;
    }

    int ActiveDimensions() const                                                ///< Returns the number of active dimensions: 3 for all three dimensions active, 2 for any two dimensions active, and 1 for a single dimension.
    {
        return DX + DY + DZ;
    }

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

    void SetNewBcells(const int new_b_cells)                                    ///< Sets new number of boundary cells. Destroys the stored data.
    {
        b_cells = new_b_cells;

        Reallocate(Size_X, Size_Y, Size_Z);
    }

 protected:
    long int Size_X;
    long int Size_Y;
    long int Size_Z;

    long int Size_X_BC;
    long int Size_Y_BC;
    long int Size_Z_BC;
    long int Size_D;
    
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

    std::array<size_t, Rank> TensorDimensions;
    std::vector<Tensor<T, Rank>> locTensors;
    std::vector<T> locData;

 private:
    size_t Index(const long int x, const long int y, const long int z) const    ///< Converts 3D coordinates into flattened storage index.
    {
        assert(x < Size_X + b_cells*DX && "Access beyond storage range");
        assert(y < Size_Y + b_cells*DY && "Access beyond storage range");
        assert(z < Size_Z + b_cells*DZ && "Access beyond storage range");
        assert(x >= - b_cells*DX && "Access beyond storage range");
        assert(y >= - b_cells*DY && "Access beyond storage range");
        assert(z >= - b_cells*DZ && "Access beyond storage range");

        return (((x + b_cells*DX)*Size_Y_BC + y + b_cells*DY)*Size_Z_BC + z + b_cells*DZ);
    }

    size_t IndexT(const std::array<size_t,Rank> position) const                 ///< Converts tensor indices into flattened tensor storage index.
    {
        switch (Rank)
        {
            case 1:
                return position[0];
            case 2:
                return position[0]*TensorDimensions[1] + position[1];
            default:
                size_t locIndex = position[0];
                for(size_t n = 1; n < Rank; n++)
                {
                    locIndex *= TensorDimensions[n];
                    locIndex += position[n];
                }
                return locIndex;
        }
    }
};

}// namespace openphase
#endif
