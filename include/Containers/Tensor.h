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

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Raphael Schiedung
 *
 */

#ifndef TENSOR_H
#define TENSOR_H

#include <array>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Globals.h"

namespace openphase
{
template <class T, size_t Rank>
class Tensor                                                                    /// Tensor template class. It can handle any type of values but POD is preferred
{
 public:
    ~Tensor()
    {
        if(allocated)
        {
            delete[] locData;
        }
    }
    Tensor() :
        Dimensions(),
        totSize(0),
        allocated(false),
        assigned(false),
        locData(nullptr)
    {

    }
    Tensor(std::array<size_t, Rank> locDimensions) :
        Dimensions(locDimensions),
        totSize(1),
        allocated(true),
        assigned(false)
    {
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }
        if(totSize)
        {
            locData = new T[totSize] ();
        }
        else
        {
            totSize = 0;
            locData = nullptr;
            allocated = false;
        }
    }
    Tensor(T* data_ptr, std::array<size_t, Rank> locDimensions) :
        Dimensions(locDimensions),
        totSize(1),
        allocated(false),
        assigned(true),
        locData(data_ptr)
    {
        for(size_t n = 0; n < Rank; n++)
        {
            totSize *= Dimensions[n];
        }
    }
    Tensor(const Tensor<T, Rank>& locTensor) :
        Dimensions(locTensor.Dimensions),
        totSize(locTensor.totSize),
        assigned(false)
    {
        if(totSize)
        {
            locData = new T[totSize] ();
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] = locTensor[n];
            }
            allocated = true;
        }
        else
        {
            locData = nullptr;
            allocated = false;
        }
    }
    T& operator()(const std::array<size_t,Rank> locPosition)
    {
        for(size_t n = 0; n < Rank; n++)
            assert(locPosition[n] < Dimensions[n] && "Access beyond storage range");

        return locData[Index(locPosition)];
    }
    T const& operator()(const std::array<size_t,Rank> locPosition) const
    {
        for(size_t n = 0; n < Rank; n++)
            assert(locPosition[n] < Dimensions[n] && "Access beyond storage range");

        return locData[Index(locPosition)];
    }
    template<size_t n>
    Tensor<T, n> get(const std::array<size_t,Rank-n> locPosition) const
    {
        for(size_t m = 0; m < Rank-n; m++)
            assert(locPosition[m] < Dimensions[m] && "Access beyond storage range");

        std::array<size_t,Rank> tmpPosition;
        tmpPosition.fill(0);
        for(size_t m = 0; m < Rank-n; m++)
        {
            tmpPosition[m] = locPosition[m];
        }

        std::array<size_t,n> locDimensions;
        for(size_t m = Rank - n, p = 0; m < Rank; m++,p++)
        {
            locDimensions[p] = Dimensions[m];
        }

        Tensor<T, n> myReturn(locDimensions);
        for(size_t p = Index(tmpPosition), m = 0; p < totSize; p++, m++)
        {
            myReturn[m] = locData[p];
        }
        return myReturn;
    }
    T& operator[](const size_t position)
    {
        assert(position < totSize && "Access beyond storage range");

        return locData[position];
    }
    T const& operator[](const size_t position) const
    {
        assert(position < totSize && "Access beyond storage range");

        return locData[position];
    }
    void Assign(T* data_ptr, const std::array<size_t, Rank> locDimensions)
    {
        if(not allocated)
        {
            locData = data_ptr;
            allocated = false;
            assigned  = true;
            Dimensions = locDimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Assignment attempt of already allocated Tensor!\n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }

    }
    size_t Allocate(const std::array<size_t, Rank> locDimensions)
    {
        if(not (allocated or assigned))
        {
            Dimensions = locDimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
                return sizeof(T)*totSize;
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Allocation attempt of already allocated or assigned Tensor!\n"
                    << "Reallocate() should be used instead!\n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return 0;
    }

    void Reallocate(const std::array<size_t, Rank> locDimensions)
    {
        if(allocated)
        {
            delete[] locData;
        }
        if(not assigned)
        {
            Dimensions = locDimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Reallocation attempt of already assigned Tensor!\n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
    }

    void set_to_value(T value)
    {
        if(allocated or assigned)
        for (size_t i = 0; i < totSize; i++)
        {
            locData[i] = value;
        }
    }
    void set_to_zero(void)
    {
        if(allocated or assigned)
        for (size_t i = 0; i < totSize; i++)
        {
            locData[i] = T();
        }
    }
    size_t size(size_t n) const
    {
        if(n < Rank)
        {
            return Dimensions[n];
        }
        else
        {
            return 0;
        }
    }
    size_t size(void) const
    {
        return totSize;
    }
    size_t rank(void) const
    {
        return Rank;
    }
    const T* data(void) const
    {
        return locData;
    }
    T* data(void)
    {
        return locData;
    }
    Tensor<T, Rank>& operator=(const Tensor<T, Rank>& locTensor)
    {
        if (this == &locTensor)
        {
            return *this;
        }

        if ((not allocated) and (not assigned))
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] = locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in assignment operator!\n";
            throw std::logic_error(message.str());
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator+=(const Tensor<T2, Rank>& locTensor)
    {
        if ((not allocated) and (not assigned))
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] += (T)locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                << ">: Different tensor size in += operator!\n";

            message << "Dimensions lside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << Dimensions[i] << " ";
            message << "\n";
            message << "Dimensions rside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << locTensor.Dimensions[i] << " \n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator-=(const Tensor<T2, Rank>& locTensor)
    {
        if ((not allocated) and (not assigned))
        {
            Dimensions = locTensor.Dimensions;
            totSize = 1;
            for(size_t n = 0; n < Rank; n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
        }
        if (locTensor.Dimensions == Dimensions)
        {
            for(size_t n = 0; n < totSize; n++)
            {
                locData[n] += (T)locTensor[n];
            }
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                << ">: Different tensor size in -= operator!"
                << std::endl;

            message << "Dimensions lside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << Dimensions[i] << " ";
            message << std::endl;
            message << "Dimensions rside: ";
                for (size_t i = 0; i < Rank; i++)
                    message << locTensor.Dimensions[i] << " \n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator/=(const T2 number)
    {
        if (allocated or assigned)
        for(size_t n = 0; n < totSize; n++)
        {
            locData[n] /= (T)number;
        }
        return *this;
    }
    template<typename T2>
    Tensor<T, Rank>& operator*=(const T2 number)
    {
        if (allocated or assigned)
        for(size_t n = 0; n < totSize; n++)
        {
            locData[n] *= (T)number;
        }
        return *this;
    }
    Tensor<T, Rank> operator+(const Tensor<T, Rank>& locTensor) const
    {
        if (locTensor.Dimensions == Dimensions)
        {
            Tensor<T,Rank> myReturn(locTensor);

            for(size_t n = 0; n < totSize; n++)
            {
                myReturn[n] += locData[n];
            }
            return myReturn;
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in summation operator!\n";
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return *this;
    }
    Tensor<T, Rank> operator-(const Tensor<T, Rank>& locTensor) const
    {
        if (locTensor.Dimensions == Dimensions)
        {
            Tensor<T,Rank> myReturn(locTensor);

            for(size_t n = 0; n < totSize; n++)
            {
                myReturn[n] -= locData[n];
            }
            return myReturn;
        }
        else
        {
            std::stringstream message;
            message << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in subtraction operator!\n";
            throw std::logic_error(message.str());
            std::cerr << message.str();
            OP_Exit(EXIT_FAILURE);
        }
        return *this;
    }
    Tensor<T, Rank> operator*(double val) const
    {
        Tensor<T, Rank> myReturn(*this);
        for(size_t n = 0; n < totSize; n++)
        {
            myReturn[n] *= val;
        }
        return myReturn;
    }
    Tensor<T, Rank> operator/(double val) const
    {
        Tensor<T, Rank> myReturn(*this);
        for(size_t n = 0; n < totSize; n++)
        {
            myReturn[n] /= val;
        }
        return myReturn;
    }
    void normalize()
    {
        double sum = 0.0;
        for (size_t n = 0; n < totSize; n++)
        {
            sum += locData[n];
        }
        for (size_t n = 0; n < totSize; n++)
        {
            locData[n] /= sum;
        }
    }
    bool IsAllocated() const
    {
        return allocated;
    }
    bool IsAssigned() const
    {
        return assigned;
    }
    void DeAllocate()
    {
        if(allocated)
        {
            delete[] locData;
        }
        allocated = false;
        assigned = false;
        totSize = 0;
        for(size_t n = 0; n < Rank; n++)
        {
            Dimensions[n] = 0;
        }
    }
 protected:
    std::array<size_t, Rank> Dimensions;
    size_t totSize;
    bool allocated;
    bool assigned;
    T*  locData;

 private:
    size_t Index(const std::array<size_t,Rank> position) const
    {
        switch (Rank)
        {
            case 1:
                return position[0];
            case 2:
                return position[0]*Dimensions[1] + position[1];
            default:
                size_t locIndex = position[0];
                for(size_t n = 1; n < Rank; n++)
                {
                    locIndex *= Dimensions[n];
                    locIndex += position[n];
                }
                return locIndex;
        }
    }
};

}// namespace openphase
#endif
