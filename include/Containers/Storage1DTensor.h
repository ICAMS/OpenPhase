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

#ifndef STORAGE1DTENSOR_H
#define STORAGE1DTENSOR_H

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
class ClearCaller1DTensor;

template <typename A, typename T>
class ClearCaller1DTensor<A, T, typename std::enable_if<has_clear<T>::value>::type>
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
class ClearCaller1DTensor<A, T, typename std::enable_if<!has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        std::cerr << "ERROR: Storage1DTensor.h: ClearCaller1DTensor::Clear() called for non-POD with no clear() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <typename A, typename T, typename = void>
class PackCaller1DTensor;

template <typename A, typename T>
class PackCaller1DTensor<A, T, typename std::enable_if<!std::is_class<T>::value && std::is_standard_layout<T>::value && std::is_trivial<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        {
            for (size_t n = 0; n < self(i).size(); ++n)
            {
                buffer.push_back(self(i)[n]);
            }
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        {
            for (size_t n = 0; n < self(i).size(); ++n)
            {
                self(i)[n] = buffer[it]; ++it;
            }
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller1DTensor<A, T, typename std::enable_if<std::is_class<T>::value && has_pack<T>::value>::type>
{
 public:
    static T pack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        buffer.clear();
        for (long int i = window[0]; i < window[1]; ++i)
        {
            for (size_t n = 0; n < self(i).size(); ++n)
            {
                self(i)[n].pack(buffer);
            }
        }
        return T();
    }
    static T unpack(A& self, std::vector<double>& buffer, std::vector<long int> window)
    {
        size_t  it = 0;
        for (long int i = window[0]; i < window[1]; ++i)
        {
            for (size_t n = 0; n < self(i).size(); ++n)
            {
                self(i)[n].unpack(buffer, it);
            }
        }
        return T();
    }
};

template <typename A, typename T>
class PackCaller1DTensor<A, T, typename std::enable_if<std::is_class<T>::value && !has_pack<T>::value>::type>
{
 public:
    static T call(A& self, std::vector<long int> window)
    {
        std::cerr << "ERROR: Storage1DTensor.h: PackCallerTensor1DTensor::pack() called for non-POD with no pack() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

template <class T, size_t Rank>
class Storage1DTensor                                                           ///< 1D storage template class of tensors. Can handle any type of values
{
 public:
    friend class ClearCaller1DTensor< Storage1DTensor<T, Rank> , T>;
    friend class PackCaller1DTensor< Storage1DTensor<T, Rank> , T>;

    ~Storage1DTensor()                                                          ///< Destructor.
    {
        delete[] locTensors;
        delete[] locData;
    }

    Storage1DTensor()                                                           ///< Default constructor. Creates empty container.
    {
        Size_X = 0;
        Size_D = 0;
        b_cells = 0;

        locTensors = nullptr;
        locData    = nullptr;
    }

    Storage1DTensor(const long int nx,
               const std::array<size_t,Rank> nn,
               const long int bc)                                               ///< Constructor. Creates a container of the specified size.
    {
        TensorDimensions = nn;

        Size_X = nx;
        b_cells = bc;

        size_t Size_X_BC = Size_X + 2*b_cells;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData = new T[Size_X_BC*Size_D]();

        locTensors = new Tensor<T, Rank>[Size_X_BC] ();

        for(size_t i = 0; i < Size_X_BC; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    Storage1DTensor(const Storage1DTensor<T,Rank>& rhs)                         ///< Copy constructor. Creates and initializes the container with the copy of the rhs.
    {
        if (rhs.IsAllocated())
        {
            TensorDimensions = rhs.TensorDimensions;

            Size_X = rhs.Size_X;
            b_cells = rhs.b_cells;

            size_t Size_X_BC = Size_X + 2*b_cells;

            Size_D = 1;
            for(size_t n = 0; n < TensorDimensions.size(); n++)
            {
                Size_D *= TensorDimensions[n];
            }

            locData = new T[Size_X_BC*Size_D]();

            locTensors = new Tensor<T, Rank>[Size_X_BC] ();

            for(size_t i = 0; i < Size_X_BC; i++)
            {
                locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
            }

            for(size_t i = 0; i < Size_X_BC; i++)
            {
                locTensors[i] = rhs[i];
            }
        }
        else
        {
            locTensors = nullptr;
        }
    }

    void Allocate(const long int nx,
                  const std::array<size_t,Rank> nn,
                  const long int bc)                                            ///< Allocates previously created empty storage container.
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1DTensor<" << typeid(T).name() << ", " <<  Rank
                      << ">::Allocate()\n"
                      << "Attempt of allocating of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        std::copy(nn.begin(), nn.end(), TensorDimensions.begin());

        Size_X = nx;
        b_cells = bc;

        size_t Size_X_BC = Size_X + 2*b_cells;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData = new T[Size_X_BC*Size_D]();

        locTensors = new Tensor<T, Rank>[Size_X_BC] ();

        for(size_t i = 0; i < Size_X_BC; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    void Allocate(const Storage1DTensor<T,Rank>& rhs)                           ///< Allocates current container to the dimensions of the rhs.
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1DTensor<" << typeid(T).name() << ", "
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

                Size_X = rhs.Size_X;
                b_cells = rhs.b_cells;

                size_t Size_X_BC = Size_X + 2*b_cells;

                Size_D = 1;
                for(size_t n = 0; n < TensorDimensions.size(); n++)
                {
                    Size_D *= TensorDimensions[n];
                }

                locData = new T[Size_X_BC*Size_D]();

                locTensors = new Tensor<T, Rank>[Size_X_BC] ();

                for(size_t i = 0; i < Size_X_BC; i++)
                {
                    locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
                }
            }
        }
    }

    void AllocateCopy(const Storage1DTensor<T,Rank>& rhs)                       ///< Allocates current container and initializes it with the data from the rhs.
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1DTensor<"
                      << typeid(T).name() << ", " <<  Rank
                      << ">::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        if (this != &rhs)
        {
            if(rhs.IsAllocated())
            {
                TensorDimensions = rhs.TensorDimensions;

                Size_X = rhs.Size_X;
                b_cells = rhs.b_cells;

                size_t Size_X_BC = Size_X + 2*b_cells;

                Size_D = 1;
                for(size_t n = 0; n < TensorDimensions.size(); n++)
                {
                    Size_D *= TensorDimensions[n];
                }

                locData = new T[Size_X_BC*Size_D]();

                locTensors = new Tensor<T, Rank>[Size_X_BC] ();

                for(size_t i = 0; i < Size_X_BC; i++)
                {
                    locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
                }

                for(size_t i = 0; i < Size_X_BC; i++)
                {
                    locTensors[i] = rhs[i];
                }
            }
        }
    }

    void Reallocate(const long int nx)                                          ///< Reallocates non-empty container to new size while keeping the b_cells.
    {
        delete[] locTensors;
        delete[] locData;

        Size_X = nx;

        size_t Size_X_BC = Size_X + 2*b_cells;

        Size_D = 1;
        for(size_t n = 0; n < TensorDimensions.size(); n++)
        {
            Size_D *= TensorDimensions[n];
        }

        locData = new T[Size_X_BC*Size_D]();

        locTensors = new Tensor<T, Rank>[Size_X_BC] ();

        for(size_t i = 0; i < Size_X_BC; i++)
        {
            locTensors[i].Assign(&locData[i*Size_D],TensorDimensions);
        }
    }

    Tensor<T, Rank>& operator()(const long int x)                               ///< Random access operator. Returns the reference to the tensor pointed to by the index x.
    {
        return locTensors[Index(x)];
    }

    Tensor<T, Rank> const& operator()(const long int x) const                   ///< Random access operator. Returns const reference to the tensor pointed to by the index x.
    {
        return locTensors[Index(x)];
    }

    T& operator()(const long int x, const std::array<size_t, Rank> Position)    ///< Random access operator. Returns the reference to the storage element pointed to by the index x and position inside of tensor.
    {
        return locData[Index(x)*Size_D + IndexT(Position)];
    }

    const T& operator()(const long int x, const std::array<size_t, Rank> Position) const///< Random access operator. Returns const reference to the storage element pointed to by the index x and position inside of tensor.
    {
        return locData[Index(x)*Size_D + IndexT(Position)];
    }

    Tensor<T, Rank> at(const double x) const                                    ///< Arbitrary position access operator. Returns interpolated tensor at arbitrary location inside of the container dimensions.
    {
        long int x0 = floor(x);
        double dx = fabs(x - x0);

        Tensor<T, Rank> tempTensor = locTensors[Index(x0)]*(1.0 - dx)
                                   + locTensors[Index(x0+1)]*dx;

        return tempTensor;
    }

    Tensor<T, Rank>& operator[](const size_t idx)                               ///< Random access operator to the raw storage. Returns the reference to the element pointed to by the index x.
    {
        assert(idx < size() && "Access beyond storage range");
        return locTensors[idx];
    }

    Tensor<T, Rank>const& operator[](const size_t idx) const                    ///< Random access operator to the raw storage. Returns const reference to the element pointed to by the index x.
    {
        assert(idx < size() && "Access beyond storage range");
        return locTensors[idx];
    }

    bool IsNotAllocated() const                                                 ///< Returns true if the storage is not allocated.
    {
        return (locTensors == nullptr);
    }

    bool IsAllocated() const                                                    ///< Returns true if the storage is allocated.
    {
        return !(locTensors == nullptr);
    }

    bool IsSize(const long int Nx)                                              ///< Returns true if the current storage is of the specified size Nx.
    {
        return (Size_X == Nx);
    }

    void Remesh(const long int nx)                                              ///< Remeshes the storage to new dimensions while keeping the data. Uses linear interpolation.
    {
        T* tempData = new T [(nx + 2*b_cells)*Size_D] ();
        Tensor<T, Rank>* tempTensors = new Tensor<T, Rank>[(nx + 2*b_cells)] ();

        for(long int i = 0; i < (nx + 2*b_cells); i++)
        {
            tempTensors[i].Assign(&tempData[i],TensorDimensions);
        }

        double Xscale = double(Size_X)/double(nx);
        for(long int x = 0; x < nx; x++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            double dx = (x*Xscale - x0);

            tempTensors[(x + b_cells)] =
                  locTensors[Index(x0)]*(1.0 - dx) +
                  locTensors[Index(x0+1)]*dx;
        }
        delete[] locTensors;
        delete[] locData;

        Size_X = nx;

        locTensors = tempTensors;
        locData    = tempData;
    }

    Storage1DTensor<T,Rank>& operator=(const Storage1DTensor<T,Rank>& rhs)      ///< Assignment operator. Assigns the copy of the rhs to the current container.
    {
        if (this != &rhs)
        {
            if (IsAllocated())
            {
                for(long int i = 0; i < (Size_X + 2*b_cells); i++)
                {
                    locTensors[i] = rhs[i];
                }
            }
            else
            {
                AllocateCopy(rhs);
            }
        }
        return *this;
    }

    void Clear(void)                                                            ///< Clears the container while keeping its storage size.
    {
        ClearCaller1DTensor< Storage1DTensor<T, Rank>, T> K;
        K.call(*this);
    }

//    std::vector<double> pack(std::vector<long int> window)                      ///< Packs (writes) the content of the current container into the MPI communication buffer.
//    {
//        std::vector<double> buffer;
//        PackCaller1DT< Storage1DT<T,Rank>, T> K;
//        K.pack(*this, buffer, window);
//        return buffer;
//    }
//
//    void unpack(std::vector<double>& buffer, std::vector<long int> window)      ///< Unpacks (reads) the content of the current container from the MPI communication buffer.
//    {
//        PackCaller1DT< Storage1DT<T,Rank>, T> K;
//        K.unpack(*this, buffer, window);
//    }

    std::vector<double> pack()                                                  ///< Packs (writes) the content of the current container into the MPI communication buffer.
    {
        std::vector<double> buffer;
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1DTensor< Storage1DTensor<T,Rank>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer)                                    ///< Unpacks (reads) the content of the current container from the MPI communication buffer.
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1DTensor< Storage1DTensor<T,Rank>, T> K;
        K.unpack(*this, buffer, window);
    }

    void set_to_value(Tensor<T, Rank>& val)                                     ///< Sets all entries of the container to the specified value.
    {
        if(locTensors != nullptr)
        {
            if(TensorDimensions == val.TensorDimensions)
            {
                for(long int i = 0; i < (Size_X + 2*b_cells); i++)
                {
                    locTensors[i] = val;
                }
            }
            else
            {
                std::cerr << "ERROR: Storage1DTensor<" << typeid(T).name() << ", "
                          <<  Rank << ">::set_to_value(): Wrong tensor size!\n"
                          << "Terminating!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
        }
        else
        {
            std::cerr << "ERROR: Storage1DTensor<" << typeid(T).name() << ", " <<  Rank
                      << ">::set_to_value() operation on a not allocated storage!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
    }

    auto* data(void)                                                            ///< Returns pointer to the raw data storage.
    {
        if constexpr (has_datapointer<T>::value)
        {
            return locData->data();
        }
        else return locData;
    }

    size_t sizeX() const                                                        ///< Returns the size of the current container without boundary cells.
    {
        return Size_X;
    }

    long int Bcells() const                                                     ///< Returns the number of boundary cells.
    {
        return b_cells;
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
    long int Size_X;
    size_t   Size_D;
    long int b_cells;

    std::array<size_t, Rank> TensorDimensions;
    Tensor<T, Rank>* locTensors;
    T* locData;

 private:
    size_t Index(const long int x) const
    {
        assert(x < Size_X + b_cells && "Access beyond storage range");
        assert(x >= -b_cells && "Access beyond storage range");

        return (x + b_cells);
    }

    size_t IndexT(const std::array<size_t,Rank> position) const
    {
        size_t locIndexT = 0;
        switch (Rank)
        {
            case 1:
            {
                locIndexT = position[0];
                break;
            }
            case 2:
            {
                locIndexT = position[0]*TensorDimensions[1] + position[1];
                break;
            }
            default:
            {
                locIndexT = position[0];
                for(size_t n = 1; n < Rank; n++)
                {
                    locIndexT *= TensorDimensions[n];
                    locIndexT += position[n];
                }
                break;
            }
        }

        assert(locIndexT < size_t(Size_D) && "Access beyond storage range");

        return locIndexT;
    }
};

}// namespace openphase
#endif
