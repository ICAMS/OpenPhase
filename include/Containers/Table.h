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


#ifndef TABLE_H
#define TABLE_H

#include <cassert>
#include <vector>

namespace openphase
{

template <class T>
class OP_EXPORTS Table                                                          ///< Table template class. Can handle any type of values. Bool type specialization is below
{
 public:
    Table()                                                                     ///< Default constructor. Constructs empty table.
    {
        Size_N = 0;
        Size_M = 0;
        allocated = false;
    }

    Table(const Table<T>& rhs)                                                  ///< Copy constructor. Initializes current container with a copy of rhs
    {
        Size_N = rhs.Size_N;
        Size_M = rhs.Size_M;

        storage = rhs.storage;

        allocated = true;
    }

    Table(const size_t N, const size_t M)                                       ///< Constructor. Allocates internal storage to the specified dimensions (N,M).
    {
        Allocate(N, M);
    }

    T const& operator()(const size_t n, const size_t m) const                   ///< Random access operator. Returns constant reference to the table element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    T& operator()(const size_t n, const size_t m)                               ///< Random access operator. Returns reference to the table element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    T const& operator[](const size_t idx) const                                 ///< Random access operator. Returns constant reference to the internal storage element with the given idx.
    {
        assert(idx < storage.size() && "Access beyond storage range");

        return storage[idx];
    }

    T& operator[](const size_t idx)                                             ///< Random access operator. Returns the reference to the internal storage element with the given idx.
    {
        assert(idx < storage.size() && "Access beyond storage range");

        return storage[idx];
    }

    void Allocate(const size_t N, const size_t M)                               ///< Allocates internal storage to the specified table dimensions (N,M)
    {
        assert(storage.size() == 0 && "Attempt of allocating of already allocated Table");
        Size_N = N;
        Size_M = M;

        storage.resize(Size_N*Size_M);
        allocated = true;
    }

    void Reallocate(const size_t N, const size_t M)                             ///< Reallocates internal storage to new dimensions. Previously stored data is erased.
    {
        Size_N = N;
        Size_M = M;
        storage.clear();
        storage.resize(Size_N*Size_M);
        allocated = true;
    }

    void set(const size_t n, const size_t m, const T value)                     ///< Sets table element (n,m) to the specified value
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] = value;
    }

    void set_to_value(const T value)                                            ///< Sets all table elements to the specified value
    {
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            storage[Index(n,m)] = value;
        }
    }

    T get(const size_t n, const size_t m) const                                 ///< Returns the value of Table element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    size_t sizeN() const                                                        ///< Returns first dimension of the stored table
    {
        return Size_N;
    }

    size_t sizeM() const                                                        ///< Returns the second dimension of the stored table
    {
        return Size_M;
    }

    size_t size() const                                                         ///< Returns the size of the internal flattened storage
    {
        return Size_N*Size_M;
    }

    bool IsNotAllocated() const                                                 ///< Returns true if Table is not allocated (internal storage size is zero).
    {
        return (storage.size() == 0);
    }

    bool IsAllocated() const                                                    ///< Returns true if Table is allocated (internal storage size is not zero).
    {
        return (storage.size() != 0);
    }

    Table<T>& operator=(const Table<T>& rhs)                                    ///< Assignment operator, assumes equal dimensions tables for lhs and rhs, or empty lhs. Does not check dimensions.
    {
        if(allocated)
        {
            assert(Size_N == rhs.sizeN() && "Tables have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Tables have incompatible dimensions!");

            storage = rhs.storage;
        }
        else
        {
            Size_N = rhs.Size_N;
            Size_M = rhs.Size_M;
            storage = rhs.storage;
            allocated = true;
        }
        return *this;
    }

 protected:
 private:
    std::vector<T> storage;
    size_t Size_N;
    size_t Size_M;
    bool allocated;
    size_t Index(const size_t n, const size_t m) const
    {
        return n*Size_M + m;
    }
};

template <>
class OP_EXPORTS Table<bool>                                                    ///< Table template class specialization for bool type (uses int instead of bool due to std::vector<bool> peculiarity).
{
 public:
    Table()                                                                     ///< Default constructor. Constructs empty table.
    {
        Size_N = 0;
        Size_M = 0;
        allocated = false;
    }

    Table(const Table<bool>& rhs)                                               ///< Copy constructor. Initializes current container with a copy of rhs
    {
        Size_N = rhs.Size_N;
        Size_M = rhs.Size_M;

        storage = rhs.storage;

        allocated = true;
    }

    Table(const size_t N, const size_t M)                                       ///< Constructor. Allocates internal storage to the specified dimensions (N,M).
    {
        Allocate(N, M);
    }

    int const& operator()(const size_t n, const size_t m) const                 ///< Random access operator. Returns constant reference to the table element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    int& operator()(const size_t n, const size_t m)                             ///< Random access operator. Returns reference to the table element (n,m)

    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    int const& operator[](const size_t idx) const                               ///< Random access operator. Returns constant reference to the internal storage element with the given idx.
    {
        assert(idx < storage.size() && "Access beyond storage range");

        return storage[idx];
    }

    int& operator[](const size_t idx)                                           ///< Random access operator. Returns the reference to the internal storage element with the given idx.
    {
        assert(idx < storage.size() && "Access beyond storage range");

        return storage[idx];
    }

    void Allocate(const size_t N, const size_t M)                               ///< Allocates internal storage to the specified table dimensions (N,M)
    {
        assert(storage.size() == 0 && "Attempt of allocating of already allocated Table");
        Size_N = N;
        Size_M = M;

        storage.resize(Size_N*Size_M);
        allocated = true;
    }

    void Reallocate(const size_t N, const size_t M)                             ///< Reallocates internal storage to new dimensions. Previously stored data is erased.
    {
        Size_N = N;
        Size_M = M;
        storage.clear();
        storage.resize(Size_N*Size_M);
        allocated = true;
    }

    void set(const size_t n, const size_t m, const bool value)                  ///< Sets table element (n,m) to the specified value
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] = value;
    }

    void set_to_value(const bool value)                                         ///< Sets all table elements to the specified value
    {
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            storage[Index(n,m)] = value;
        }
    }

    bool get(const size_t n, const size_t m) const                              ///< Returns the value of Table element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    size_t sizeN() const                                                        ///< Returns first dimension of the stored table
    {
        return Size_N;
    }

    size_t sizeM() const                                                        ///< Returns the second dimension of the stored table
    {
        return Size_M;
    }

    size_t size() const                                                         ///< Returns the size of the internal flattened storage
    {
        return storage.size();
    }

    bool IsNotAllocated() const                                                 ///< Returns true if Table is not allocated (internal storage size is zero).
    {
        return (storage.size() == 0);
    }

    bool IsAllocated() const                                                    ///< Returns true if Table is allocated (internal storage size is not zero).
    {
        return (storage.size() != 0);
    }

    Table<bool>& operator=(const Table<bool>& rhs)                              ///< Assignment operator, assumes equal dimensions tables for lhs and rhs, or empty lhs. Does not check dimensions.
    {
        if(allocated)
        {
            assert(Size_N == rhs.sizeN() && "Tables have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Tables have incompatible dimensions!");

            storage = rhs.storage;
        }
        else
        {
            Size_N = rhs.Size_N;
            Size_M = rhs.Size_M;
            storage = rhs.storage;
            allocated = true;
        }
        return *this;
    }

 protected:
 private:
    std::vector<int> storage;
    size_t Size_N;
    size_t Size_M;
    bool allocated;
    size_t Index(const size_t n, const size_t m) const
    {
        return n*Size_M + m;
    }
};

}// namespace openphase
#endif
