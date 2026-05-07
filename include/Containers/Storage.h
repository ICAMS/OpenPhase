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

#ifndef STORAGE_H
#define STORAGE_H

#include <cassert>

namespace openphase
{

template <class T>
class Storage                                                                   ///< Storage template class. Can handle any type of values. Wraps around std::vector for consistency of interface with other storages.
{
 public:
    Storage() : Array()                                                         ///< Default constructor. Creates empty storage.
    {
    }

    Storage(const Storage<T>& rhs) : Array(rhs.Array)                           ///< Copy constructor. Initializes current storage with the copy of the rhs.
    {
    }

    Storage<T>& operator=(const Storage<T>& rhs)                                ///< Assignment operator. Assigns the content of the rhs to the current Storage.
    {
        Array = rhs.Array;
        return *this;
    }

    T& operator[](const size_t idx)                                             ///< Random access operator. Returns the reference to the value pointed to by the idx.
    {
        assert(idx < Array.size() && "Access beyond storage range");
        return Array[idx];
    }

    T const& operator[](const size_t idx) const                                 ///< Random access operator. Returns const reference to the value pointed to by the idx.
    {
        assert(idx < Array.size() && "Access beyond storage range");
        return Array[idx];
    }

    void Allocate(const size_t size_in)                                         ///< Allocates storage to the given size_in.
    {
        Array.resize(size_in);
    }

    void Reallocate(const size_t new_size)                                      ///< Reallocates storage to the new_size. Erases old data.
    {
        Array.clear();
        Array.resize(new_size);
    }

    size_t size() const                                                         ///< Returns storage size.
    {
        return Array.size();
    }

    bool IsAllocated() const                                                    ///< Returns true if storage is allocated.
    {
        return Array.size() != 0;
    }

    bool IsNotAllocated() const                                                 ///< Returns true if storage is not allocated.
    {
        return Array.size() == 0;
    }
 protected:
 private:
    std::vector<T> Array;
};

template <>
class Storage<bool>                                                             ///< Bool type specialization of the Storage template class. Wraps around std::vector for consistency of interface with other storages. Uses std::vector<int> as internal sotrage due to std::vector peculiarity.
{
 public:
    Storage() : Array()                                                         ///< Default constructor. Creates empty storage.
    {
    }

    Storage(const Storage<bool>& rhs) : Array(rhs.Array)                        ///< Copy constructor. Initializes current storage with the copy of the rhs.
    {
    }

    Storage<bool>& operator=(const Storage<bool>& rhs)                          ///< Assignment operator. Assigns the content of rhs to the current Storage.
    {
        Array = rhs.Array;
        return *this;
    }

    int& operator[](const size_t idx)                                           ///< Random access operator. Returns the reference to the value pointed to by the idx.
    {
        assert(idx < Array.size() && "Access beyond storage range");
        return Array[idx];
    }

    int const& operator[](const size_t idx) const                               ///< Random access operator. Returns const reference to the value pointed to by the idx.
    {
        assert(idx < Array.size() && "Access beyond storage range");
        return Array[idx];
    }

    void Allocate(const size_t size_in)                                         ///< Allocates storage to the given size_in.
    {
        Array.resize(size_in);
    }

    void Reallocate(const size_t new_size)                                      ///< Reallocates storage to the new_size. Erases old data.
    {
        Array.clear();
        Array.resize(new_size);
    }

    size_t size() const                                                         ///< Returns storage size.
    {
        return Array.size();
    }

    bool IsAllocated() const                                                    ///< Returns true if storage is allocated.
    {
        return Array.size() != 0;
    }

    bool IsNotAllocated() const                                                 ///< Returns true if storage is not allocated.
    {
        return Array.size() == 0;
    }
 protected:
 private:
    std::vector<int> Array;
};

}// namespace openphase
#endif
