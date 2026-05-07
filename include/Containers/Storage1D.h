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

#ifndef STORAGE1D_H
#define STORAGE1D_H

#include <array>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "Tensor.h"
#include "TypeTraits.h"

namespace openphase
{

template <typename A, typename T, typename = void>
class ClearCaller1D;

template <typename A, typename T>
class ClearCaller1D<A, T, typename std::enable_if<has_clear<T>::value>::type>
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
class ClearCaller1D<A, T, typename std::enable_if<!has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        std::cerr << "ERROR: Storage1D.h: ClearCaller1D::Clear() called for non-POD with no clear() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <typename A, typename T, typename = void>
class PackCaller1D;

template <typename A, typename T>
class PackCaller1D<A, T, typename std::enable_if<!std::is_class<T>::value && std::is_standard_layout<T>::value && std::is_trivial<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        {
            buffer.push_back(self(i));
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        {
            self(i) = buffer[it]; ++it;
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller1D<A, T, typename std::enable_if<std::is_class<T>::value && has_pack<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        {
            self(i).pack(buffer);
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t  it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        {
            self(i).unpack(buffer, it);
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller1D<A, T, typename std::enable_if<std::is_class<T>::value && !has_pack<T>::value>::type>
{
 public:
    static T call(A& self, std::vector<long int> window)
    {
        std::stringstream message;
        std::cerr << "ERROR: Storage1D.h: PackCaller1D::pack() called for non-POD with no pack() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <class T>
class Storage1D                                                                 ///< 1D storage template class. Can handle any type of numerical values
{
 public:
    friend class ClearCaller1D< Storage1D<T>, T>;
    friend class PackCaller1D< Storage1D<T>, T>;

    Storage1D():                                                                ///< Default constructor. Creates empty container.
        Size_X(0), b_cells(0), locData(0)
    {
    }

    Storage1D(const Storage1D<T>& rhs):                                         ///< Copy constructor. Creates and initializes the container with the copy of the rhs.
                Size_X(rhs.Size_X),
                b_cells(rhs.b_cells),
                locData(rhs.locData)
    {
    }

    Storage1D(const long int nx, const long int bc)                             ///< Constructor. Creates a container of the specified size.
    {
        Size_X = nx;
        b_cells = bc;

        locData.resize(Size_X + 2*b_cells);
    }

    void Allocate(const long int nx, const long int bc)                         ///< Allocates previously created empty storage container.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ">::Allocate()\n"
                      << "Attempt of reallocation of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        Size_X = nx;
        b_cells = bc;

        locData.resize(Size_X + 2*b_cells);
    }

    void Reallocate(const long int nx)                                          ///< Reallocates non-empty container to new size while keeping the b_cells.
    {
        if(locData.size() == 0)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ">::Rellocate()\n"
                      << "Attempt of reallocation of an empty storage!\n"
                      << "If it is intended, use Allocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        Size_X = nx;

        locData.resize(Size_X + 2*b_cells);
    }

    void AllocateCopy(const Storage1D<T>& rhs)                                  ///< Allocates current container and initializes it with the data from the rhs.
    {
        if(locData.size() != 0)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ">::AllocateCopy()\n"
                      << "Attempt of allocating and copying into a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &rhs)
        {
            if (rhs.IsAllocated())
            {
                Size_X = rhs.Size_X;
                b_cells = rhs.b_cells;

                locData = rhs.locData;
            }
        }
    }

    T& operator()(const long int x)                                             ///< Random access operator. Returns the reference to the element pointed to by the index x.
    {
        return locData[Index(x)];
    }

    T const& operator()(const long int x) const                                 ///< Random access operator. Returns const reference to the element pointed to by the index x.
    {
        return locData[Index(x)];
    }

    T at(const double x) const                                                  ///< Arbitrary position access operator. Returns interpolated value at arbitrary location inside of the container dimensions.
    {
        long int x0 = floor(x);
        double dx = fabs(x - x0);

        T tempValue = locData[Index(x0)]*(1.0 - dx)
                    + locData[Index(x0+1)]*dx;

        return tempValue;
    }

    T& operator[](size_t idx)                                                   ///< Random access operator to the raw storage. Returns the reference to the element pointed to by the index x.
    {
        assert(idx < size_t(Size_X + 2*b_cells) && "Access beyond storage range");
        return locData[idx];
    }

    T const& operator[](const size_t idx) const                                 ///< Random access operator to the raw storage. Returns const reference to the element pointed to by the index x.
    {
        assert(idx < size_t(Size_X + 2*b_cells) && "Access beyond storage range");
        return locData[idx];
    }

    bool IsNotAllocated() const                                                 ///< Returns true if the storage is not allocated.
    {
        return locData.size() == 0;
    }

    bool IsAllocated() const                                                    ///< Returns true if the storage is allocated.
    {
        return locData.size() != 0;
    }

    bool IsSize(const long int Nx)                                              ///< Returns true if the current storage is of the specified size Nx.
    {
        return (Size_X == Nx);
    }

    void Remesh(const long int nx)                                              ///< Remeshes the storage to new dimensions while keeping the data. Uses linear interpolation.
    {
        std::vector<T> tempArray(nx + 2*b_cells);

        double Xscale = double(Size_X)/double(nx);

        for(long int x = 0; x < nx; x++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            double dx = (x*Xscale - x0);

            tempArray[(x + b_cells)] = locData[Index(x0)]*(1.0 - dx) +
                                      locData[Index(x0+1)]*dx;
        }

        Size_X = nx;
        locData = tempArray;
    }

    Storage1D<T>& operator=(const Storage1D<T>& rhs)                            ///< Assignment operator. Assigns the copy of the rhs to the current container.
    {
        Size_X = rhs.Size_X;
        b_cells = rhs.b_cells;
        locData = rhs.locData;

        return *this;
    }

    void Clear(void)                                                            ///< Clears the container while keeping its storage size.
    {
        if constexpr (std::is_arithmetic<T>::value)
        {
            for(size_t i = 0; i < locData.size(); i++)
            {
                locData[i] = (T)0;
            }
        }
        else
        {
            ClearCaller1D< Storage1D<T> ,T> K;
            K.call(*this);
        }
    }

//    std::vector<double> pack(std::vector<long int> window)
//    {
//        std::vector<double> buffer;
//        PackCaller1D< Storage1D<T>, T> K;
//        K.pack(*this, buffer, window);
//        return buffer;
//    }
//
//    void unpack(std::vector<double>& buffer, std::vector<long int> window)
//    {
//        PackCaller1D< Storage1D<T>, T> K;
//        K.unpack(*this, buffer, window);
//    }

    std::vector<double> pack()                                                  ///< Packs (writes) the content of the current container into the MPI communication buffer.
    {
        std::vector<double> buffer;
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1D< Storage1D<T>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer)                                    ///< Unpacks (reads) the content of the current container from the MPI communication buffer.
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1D< Storage1D<T>, T> K;
        K.unpack(*this, buffer, window);
    }

    void set_to_value(T value)                                                  ///< Sets all entries of the container to the specified value.
    {
        for(size_t idx = 0; idx < locData.size(); idx++)
        {
            locData[idx] = value;
        }
    }

    T* data(void)                                                               ///< Returns pointer to the raw data storage.
    {
        return locData.data();
    }

    const T* data(void) const                                                   ///< Returns const pointer to the raw data storage.
    {
        return locData.data();
    }

    size_t sizeX() const                                                        ///< Returns the size of the current container without boundary cells.
    {
        return Size_X;
    }

    long int Bcells() const                                                     ///< Returns the number of boundary cells.
    {
        return b_cells;
    }

    size_t size() const                                                         ///< Returns total size of the container including boundary cells.
    {
        return Size_X + 2*b_cells;
    }

    void SetNewBcells(const int new_b_cells)                                    ///< Sets new number of boundary cells. Scrambles stored data.
    {
        b_cells = new_b_cells;

        Reallocate(Size_X);
    }

    bool InLimits(const long int idx)                                           ///< Returns true if the specified idx is inside of the container dimensions.
    {
        if (idx < -b_cells) return false;
        if (idx >= Size_X + b_cells) return false;
        return true;
    }

 protected:
 private:
    long int Size_X;                                                            ///< Storage size without boundary cells.
    long int b_cells;                                                            ///< Number of boundary cells.

    std::vector<T> locData;                                                     ///< Internal raw data storage.

    size_t Index(const long int x) const                                        ///< Converts specified storage position x to the raw storage index.
    {
        assert(x < Size_X + b_cells && "Access beyond storage range");
        assert(x >= - b_cells && "Access beyond storage range");

        return (x + b_cells);
    }
};

}// namespace openphase
#endif
