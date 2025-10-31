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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef STORAGE_H
#define STORAGE_H

#include <cassert>

namespace openphase
{

template <class T>
class Storage  /// 1D storage template class. Can handle any type of values
{
 public:
    Storage()
    {
        Array = nullptr;
        Size_X = 0;
    }
    Storage(const Storage<T>& rhs)
    {
        Size_X = rhs.size();
        Array = new T[Size_X] ();
        for(size_t n = 0; n != Size_X; n++)
        {
            Array[n] = rhs[n];
        }
    }
    Storage<T>& operator=(const Storage<T>& rhs)
    {
        if(Size_X == rhs.size())
        {
            Reallocate(rhs.size());
        }
        for(size_t n = 0; n != Size_X; n++)
        {
            Array[n] = rhs[n];
        }
        return *this;
    }

    T& operator[](const size_t x)
    {
        assert(x < Size_X && "Access beyond storage range");
        return Array[x];
    }
    T const& operator[](const size_t x) const
    {
        assert(x < Size_X && "Access beyond storage range");
        return Array[x];
    }
    void Allocate(const size_t nx)
    {
        Size_X = nx;
        Array = new T[Size_X] ();
    }
    void Reallocate(const size_t nx)
    {
        delete[] Array;
        Size_X = nx;
        Array = new T[Size_X] ();
    }
    size_t size() const
    {
        return Size_X;
    }
    bool IsNotAllocated() const
    {
        return (Array == nullptr);
    }
    bool IsAllocated() const
    {
        return !(Array == nullptr);
    }
    ~Storage()
    {
        delete[] Array;
    }

 protected:
 private:
    T* Array;
    size_t Size_X;
};

}// namespace openphase
#endif
