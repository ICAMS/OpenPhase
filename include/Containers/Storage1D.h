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

 *   File created :   2024
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Marvin Tegeler;
 *                         Dmitry Medvedev
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
        for (long int i = -self.b_cells(); i < self.size() + self.b_cells(); ++i)
        {
            self(i).clear();
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
class PackCallerTensor1D;

template <typename A, typename T, typename = void>
class PackCaller1D;

template <typename A, typename T>
class PackCallerTensor1D<A, T, typename std::enable_if<!std::is_class<T>::value && std::is_standard_layout<T>::value && std::is_trivial<T>::value>::type>
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
class PackCallerTensor1D<A, T, typename std::enable_if<std::is_class<T>::value && has_pack<T>::value>::type>
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
class PackCallerTensor1D<A, T, typename std::enable_if<std::is_class<T>::value && !has_pack<T>::value>::type>
{
 public:
    static T call(A& self, std::vector<long int> window)
    {
        std::cerr << "ERROR: Storage1D.h: PackCallerTensor1D::pack() called for non-POD with no pack() method.\n";
        OP_Exit(EXIT_FAILURE);
        return T();
    }
};

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

template <class T, size_t Rank>
class Storage1D                                                                 /// 3D storage template class of vector values. Can handle any type of values
{
 public:
    friend class ClearCaller1D< Storage1D<T, Rank> , T>;
    friend class PackCallerTensor1D< Storage1D<T, Rank> , T>;
    Storage1D()
    {
        Size_X    = 0;
        Size_X_BC = 0;
        Size_D    = 0;
        Bcells    = 0;

        locTensors = nullptr;
        locData    = nullptr;
    }

    Storage1D(const long int nx,
              const std::array<size_t,Rank> nn,
              const long int bc)
    {
        TensorDimensions = nn;

        Size_X = nx;
        Bcells = bc;

        Size_X_BC = Size_X + 2*Bcells;

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

    Storage1D(const Storage1D<T,Rank>& Field)
    {
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                TensorDimensions = Field.TensorDimensions;

                Size_X = Field.Size_X;
                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;

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

                const size_t totsize = Field(0).size();
                for (long int i = -Bcells; i < Size_X + Bcells; i++)
                for (size_t n = 0; n < totsize; ++n)
                {
                    locTensors[Index(i)][n] = Field(i)[n];
                }
            }
            else
            {
                locTensors = nullptr;
            }
        }
    }

    Tensor<T, Rank>& operator()(const long int x)
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        return locTensors[Index(x)];
    }

    Tensor<T, Rank> const& operator()(const long int x) const
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        return locTensors[Index(x)];
    }

    T& operator()(const long int x, const std::array<size_t, Rank> Position)
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        assert(IndexT(Position) < size_t(Size_D) && "Access beyond storage range");

        return locData[Index(x)*Size_D + IndexT(Position)];
    }

    const T& operator()(const long int x, const std::array<size_t, Rank> Position) const
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        assert(IndexT(Position) < size_t(Size_D) && "Access beyond storage range");

        return locData[Index(x)*Size_D + IndexT(Position)];
    }

    Tensor<T, Rank> at(const double x) const
    {
        long int x0 = floor(x);
        double dx = fabs(x - x0);

        Tensor<T, Rank> tempTensor = locTensors[Index(x0)]*(1.0 - dx)
                                   + locTensors[Index(x0+1)]*dx;

        return tempTensor;
    }

    Tensor<T, Rank>& operator[](const size_t idx)
    {
        assert(idx < size_t(Size_X_BC) && "Access beyond storage range");
        return locTensors[idx];
    }

    Tensor<T, Rank>const& operator[](const size_t idx) const
    {
        assert(idx < size_t(Size_X_BC) && "Access beyond storage range");
        return locTensors[idx];
    }

    void Allocate(const long int nx,
                  const std::array<size_t,Rank> nn,
                  const long int bc)
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", " <<  Rank
                      << ">::Allocate()\n"
                      << "Attempt of allocating of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        std::copy(nn.begin(), nn.end(), TensorDimensions.begin());

        Size_X = nx;
        Bcells = bc;

        Size_X_BC = Size_X + 2*Bcells;

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

    void Allocate(const Storage1D<T,Rank>& Field)
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", "
                      <<  Rank << ">::Allocate()\n"
                      << "Attempt of allocating of a nonempty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                TensorDimensions = Field.TensorDimensions;

                Size_X = Field.Size_X;

                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;
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

    void AllocateCopy(const Storage1D<T,Rank>& Field)
    {
        if(IsAllocated())
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", " <<  Rank
                      << ">::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        if (this != &Field)
        {
            if(Field.IsAllocated())
            {
                TensorDimensions = Field.TensorDimensions;

                Size_X = Field.Size_X;

                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;

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

                const size_t totsize = Field(0).size();
                for(long int i = -Bcells; i < Size_X + Bcells; i++)
                for (size_t n = 0; n < totsize; ++n)
                {
                    locTensors[Index(i)][n] = Field(i)[n];
                }
            }
        }
    }

    Storage1D<T,Rank>& operator=(const Storage1D<T,Rank>& Field)
    {
        if (this != &Field)
        {
            if (IsAllocated())
            {
                const size_t totsize = Field(0).size();
                for(long int i = -Bcells; i < Size_X + Bcells; i++)
                for (size_t n = 0; n < totsize; ++n)
                {
                    locTensors[Index(i)][n] = Field(i)[n];
                }
            }
            else
            {
                AllocateCopy(Field);
            }
        }
        return *this;
    }

    void Reallocate(const long int nx)
    {
        delete[] locTensors;
        delete[] locData;

        Size_X = nx;

        Size_X_BC = Size_X + 2*Bcells;

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

    bool IsNotAllocated() const
    {
        return (locTensors == nullptr);
    }

    bool IsAllocated() const
    {
        return !(locTensors == nullptr);
    }

    bool IsSize(const long int Nx)
    {
        return (Size_X == Nx);
    }

    void Remesh(const long int nx)
    {
        T* tempData = new T [(nx + 2*Bcells)*Size_D] ();
        Tensor<T, Rank>* tempTensors = new Tensor<T, Rank>[(nx + 2*Bcells)] ();

        for(long int i = 0; i < (nx + 2*Bcells); i++)
        {
            tempTensors[i].Assign(&tempData[i],TensorDimensions);
        }

        double Xscale = double(Size_X)/double(nx);
        for(long int x = 0; x < nx; x++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            double dx = (x*Xscale - x0);

            tempTensors[(x + Bcells)] =
                  locTensors[Index(x0)]*(1.0 - dx) +
                  locTensors[Index(x0+1)]*dx;
        }
        delete[] locTensors;
        delete[] locData;

        Size_X = nx;

        Size_X_BC = Size_X + 2*Bcells;

        locTensors = tempTensors;
        locData    = tempData;
    }

    ~Storage1D()
    {
        delete[] locTensors;
        delete[] locData;
    }

    auto* data(void)
    {
        if constexpr (has_datapointer<T>::value)
        {
            return locData->data();
        }
        else return locData;
    }

    void Clear(void)
    {
        ClearCaller1D< Storage1D<T, Rank>, T> K;
        K.call(*this);
    }

    std::vector<double> pack(std::vector<long int> window)
    {
        std::vector<double> buffer;
        PackCallerTensor1D< Storage1D<T,Rank>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer, std::vector<long int> window)
    {
        PackCallerTensor1D< Storage1D<T,Rank>, T> K;
        K.unpack(*this, buffer, window);
    }

    std::vector<double> pack()
    {
        std::vector<double> buffer;
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCallerTensor1D< Storage1D<T,Rank>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }
    void unpack(std::vector<double>& buffer)
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCallerTensor1D< Storage1D<T,Rank>, T> K;
        K.unpack(*this, buffer, window);
    }
    void set_to_value(Tensor<T, Rank>& val)
    {
        if(locTensors != nullptr)
        {
            if(TensorDimensions == val.TensorDimensions)
            {
                for(size_t idx = 0; idx < Size_X_BC; idx++)
                {
                    locTensors[idx] = val;
                }
            }
            else
            {
                std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", "
                          <<  Rank << ">::set_to_value(): Wrong tensor size!\n"
                          << "Terminating!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
        }
        else
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank
                      << ">::set_to_value() operation on a not allocated storage!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
    }

    size_t size() const
    {
        return Size_X;
    }

    size_t total_size() const
    {
        return Size_X_BC;
    }

    long int b_cells() const
    {
        return Bcells;
    }

    bool InLimits(const long int x)
    {
        if (x < -Bcells) return false;
        if (x >= Size_X + Bcells) return false;
        return true;
    }

    void SetNewBcells(const int new_b_cells)
    {
        Bcells = new_b_cells;

        Reallocate(Size_X);
    }
 protected:
    long int Size_X;
    size_t   Size_X_BC;
    size_t   Size_D;
    long int Bcells;

    std::array<size_t, Rank> TensorDimensions;
    Tensor<T, Rank>* locTensors;
    T* locData;

 private:
    size_t Index(const long int x) const
    {
        assert(x < Size_X + Bcells && "Access beyond storage range");
        assert(x >= - Bcells && "Access beyond storage range");

        return (x + Bcells);
    }

    size_t IndexT(const std::array<size_t,Rank> position) const
    {
        switch (Rank)
        {
            case 1:
            {
                return position[0];
            }
            case 2:
            {
                return position[0]*TensorDimensions[1] + position[1];
            }
            default:
            {
                size_t locIndex = position[0];
                for(size_t n = 1; n < Rank; n++)
                {
                    locIndex *= TensorDimensions[n];
                    locIndex += position[n];
                }
                return locIndex;
            }
        }
    }
};

template <class T>
class Storage1D<T,0>                                                            /// 3D storage template class specification for Rank == 0. Can handle any type of numerical values
{
 public:
    friend class ClearCaller1D< Storage1D<T,0>, T>;
    friend class PackCaller1D< Storage1D<T,0>, T>;

    Storage1D()
    {
        Size_X = 0;
        Bcells = 0;
        Size_X_BC = 0;
        locData = nullptr;
    }

    Storage1D(const Storage1D<T,0>& Field)
    {
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Size_X = Field.Size_X;

                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;

                locData = new T[Size_X_BC] ();

                for(long int i = -Bcells; i < Size_X + Bcells; i++)
                {
                    locData[Index(i)] = Field(i);
                }
            }
            else
            {
                locData = nullptr;
            }
        }
    }

    Storage1D(const long int nx, const long int bc)
    {
        Size_X = nx;
        Bcells = bc;

        Size_X_BC = Size_X + 2*Bcells;

        locData = new T[Size_X_BC] ();
    }

    void Allocate(const Storage1D<T,0>& Field)
    {
        if(locData != nullptr)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of allocating of a nonempty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Size_X = Field.Size_X;
                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;

                locData = new T [Size_X_BC] ();
            }
        }
    }

    void AllocateCopy(const Storage1D<T,0>& Field)
    {
        if(locData != nullptr)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", 0>::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignment operator instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Size_X = Field.Size_X;
                Bcells = Field.Bcells;

                Size_X_BC = Size_X + 2*Bcells;

                locData = new T[Size_X_BC] ();

                for(long int i = -Bcells; i < Size_X + Bcells; i++)
                {
                    locData[Index(i)] = Field(i);
                }
            }
            else
            {
                locData = nullptr;
            }
        }
    }

    T& operator()(const long int x)
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        return locData[Index(x)];
    }

    T const& operator()(const long int x) const
    {
        assert(Index(x) < size_t(Size_X_BC) && "Access beyond storage range");
        return locData[Index(x)];
    }

    T at(const double x) const
    {
        long int x0 = floor(x);
        double dx = fabs(x - x0);

        T tempValue = locData[Index(x0)]*(1.0 - dx)
                    + locData[Index(x0+1)]*dx;

        return tempValue;
    }

    T& operator[](size_t idx)
    {
        assert(idx < size_t(Size_X_BC) && "Access beyond storage range");
        return locData[idx];
    }

    T const& operator[](const size_t idx) const
    {
        assert(idx < size_t(Size_X_BC) && "Access beyond storage range");
        return locData[idx];
    }

    size_t Allocate(const long int nx, const long int bc)
    {
        if(locData != nullptr)
        {
            std::cerr << "ERROR: Storage1D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of reallocation of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        Size_X = nx;
        Bcells = bc;

        Size_X_BC = Size_X + 2*Bcells;

        const size_t AllocatedMemory = sizeof(T)*Size_X_BC;

        locData = new T[Size_X_BC] ();
        return AllocatedMemory;
    }

    void Reallocate(const long int nx)
    {
        delete[] locData;

        Size_X = nx;

        Size_X_BC = Size_X + 2*Bcells;
        locData = new T[Size_X_BC] ();
    }

    bool IsNotAllocated() const
    {
        return (locData == nullptr);
    }

    bool IsAllocated() const
    {
        return !(locData == nullptr);
    }

    bool IsSize(const long int Nx)
    {
        return (Size_X == Nx);
    }

    void Remesh(const long int nx)
    {
        T* tempArray = new T[(nx + 2*Bcells)] ();

        double Xscale = double(Size_X)/double(nx);

        for(long int x = 0; x < nx; x++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            double dx = (x*Xscale - x0);

            tempArray[(x + Bcells)] = locData[Index(x0)]*(1.0 - dx) +
                                      locData[Index(x0+1)]*dx;
        }
        delete[] locData;
        Size_X = nx;

        Size_X_BC = Size_X + 2*Bcells;

        locData = tempArray;
    }

    ~Storage1D()
    {
        delete[] locData;
    }

    Storage1D<T,0>& operator=(const Storage1D<T,0>& locStorage1D)
    {
        if(locData != nullptr and locStorage1D.IsAllocated())
        {
            if (locStorage1D.Size_X_BC == Size_X_BC)
            {
                for(long int i = -Bcells; i < Size_X + Bcells; i++)
                {
                    locData[Index(i)] = locStorage1D(i);
                }
            }
            else
            {
                std::cerr << "ERROR: Wrong storage size in Storage1D<"
                          << typeid(T).name() << ", 0>::operator=() !\n"
                          << "Size_X_BC: " << locStorage1D.Size_X_BC << " != " << Size_X_BC << "\n"
                          << "Terminating!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
        }
        else if (locData == nullptr and locStorage1D.IsAllocated())
        {
            AllocateCopy(locStorage1D);
        }
        return *this;
    }

    void Clear(void)
    {
        if(locData != nullptr)
        {
            if constexpr (std::is_arithmetic<T>::value)
            {
                memset(locData, 0, sizeof(T)*Size_X_BC);
            }
            else
            {
                ClearCaller1D< Storage1D<T,0> ,T> K;
                K.call(*this);
            }
        }
    }

    std::vector<double> pack(std::vector<long int> window)
    {
        std::vector<double> buffer;
        PackCaller1D< Storage1D<T,0>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer, std::vector<long int> window)
    {
        PackCaller1D< Storage1D<T,0>, T> K;
        K.unpack(*this, buffer, window);
    }

    std::vector<double> pack()
    {
        std::vector<double> buffer;
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1D< Storage1D<T,0>, T> K;
        K.pack(*this, buffer, window);
        return buffer;
    }

    void unpack(std::vector<double>& buffer)
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = Size_X;
        PackCaller1D< Storage1D<T,0>, T> K;
        K.unpack(*this, buffer, window);
    }

    void set_to_value(T val)
    {
        if(locData != nullptr)
        {
            for(size_t idx = 0; idx < Size_X_BC; idx++)
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
        return locData;
    }

    T* data(void) const
    {
        return locData;
    }

    size_t size() const
    {
        return Size_X;
    }

    long int b_cells() const
    {
        return Bcells;
    }

    size_t total_size() const
    {
        return Size_X_BC;
    }

    void SetNewBcells(const int new_b_cells)
    {
        Bcells = new_b_cells;

        Reallocate(Size_X);
    }

    bool InLimits(const long int x)
    {
        if (x < -Bcells) return false;
        if (x >= Size_X + Bcells) return false;
        return true;
    }

 protected:
 private:
    long int Size_X;
    size_t   Size_X_BC;
    long int Bcells;

    T* locData;

    size_t Index(const long int x) const
    {
        assert(x < Size_X + Bcells && "Access beyond storage range");
        assert(x >= - Bcells && "Access beyond storage range");

        return (x + Bcells);
    }
};

}// namespace openphase
#endif
