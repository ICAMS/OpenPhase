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


#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

namespace openphase
{

template <class T>
class OP_EXPORTS Matrix                                                         ///< Matrix template class. Recommended for numerical values. Can handle any type of values except "bool".
{
 public:
    Matrix()                                                                    ///< Default constructor. Constructs empty matrix.
    {
        Size_N = 0;
        Size_M = 0;
        allocated = false;
    }

    Matrix(const Matrix<T>& rhs)                                                ///< Copy constructor. Initializes current container with the copy of rhs
    {
        Size_N = rhs.sizeN();
        Size_M = rhs.sizeM();

        storage = rhs.storage;

        allocated = true;
    }

    Matrix(const size_t N, const size_t M)                                      ///< Constructor. Allocates internal storage to the specified dimensions (N,M).
    {
        Allocate(N, M);
    }

    T const& operator()(const size_t n, const size_t m) const                   ///< Random access operator. Returns constant reference to the matrix element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    T& operator()(const size_t n, const size_t m)                               ///< Random access operator. Returns reference to the matrix element (n,m)
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
        assert(storage.size() == 0 && "Attempt of allocating of already allocated Matrix");
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

    void set(const size_t n, const size_t m, const T value)                     ///< Sets matrix element (n,m) to the specified value
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] = value;
    }

    void set_to_value(const T value)                                            ///< Sets all matrix elements to the specified value
    {
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            storage[Index(n,m)] = value;
        }
    }

    void add(const size_t n, const size_t m, const T value)                     ///< Adds value to the matrix element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] += value;
    }

    T get(const size_t n, const size_t m) const                                 ///< Returns the value of the matrix element (n,m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }

    T get_min() const                                                           ///< Returns the value of the minimum element of the matrix
    {
        T min = std::numeric_limits<T>::max();
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        if (storage[Index(n,m)] < min)
        {
             min = storage[Index(n,m)];
        }
        return min;
    }

    T get_max() const                                                           ///< Returns the value of the maximum element of the matrix
    {
        T max = std::numeric_limits<T>::min();
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        if (storage[Index(n,m)] > max)
        {
             max = storage[Index(n,m)];
        }
        return max;
    }

    static Matrix<T> max(const Matrix<T>& M1, const Matrix<T>& M2)              ///< Returns matrix containing maximum elements between the two matrices
    {
        assert(M1.sizeN() == M2.sizeN() && "Matrices have incompatible dimensions!");
        assert(M1.sizeM() == M2.sizeM() && "Matrices have incompatible dimensions!");

        Matrix<T> Out(M1.sizeN(), M1.sizeM());
        for(size_t n = 0; n < M1.sizeN(); n++)
        for(size_t m = 0; m < M1.sizeM(); m++)
        {
            Out(n,m) = std::max(M1(n,m), M2(n,m));
        }
        return Out;
    }

    static Matrix<T> min(const Matrix<T>& M1, const  Matrix<T>& M2)             ///< Returns matrix containing minimum elements between the two matrices
    {
        assert(M1.sizeN() == M2.sizeN() && "Matrices have incompatible dimensions!");
        assert(M1.sizeM() == M2.sizeM() && "Matrices have incompatible dimensions!");

        Matrix<T> Out(M1.sizeN(), M1.sizeM());
        for(size_t n = 0; n < M1.sizeN(); n++)
        for(size_t m = 0; m < M1.sizeM(); m++)
        {
            Out(n,m) = std::min(M1(n,m), M2(n,m));
        }
        return Out;
    }

    T* data()                                                                   ///< Returns pointer to the internally stored data array.
    {
        return storage.data();
    }

    size_t size()                                                               ///< Total size of the internal flattened storage.
    {
        return storage.size();
    }

    size_t sizeN() const                                                        ///< Returns first dimension of the stored matrix.
    {
        return Size_N;
    }

    size_t sizeM() const                                                        ///< Returns second dimension of the stored matrix.
    {
        return Size_M;
    }

    bool IsNotAllocated() const                                                 ///< Returns true if Table is not allocated (internal storage size is zero).
    {
        return (storage.size() == 0);
    }

    bool IsAllocated() const                                                    ///< Returns true if Table is allocated (internal storage size is not zero).
    {
        return (storage.size() != 0);
    }

    Matrix<T> operator*(const Matrix<T>& rhs) const                             ///< Returns matrix dot product between the current matrix and rhs.
    {
        assert(Size_N == rhs.sizeM() && "Matrices have incompatible dimensions!");
        assert(Size_M == rhs.sizeN() && "Matrices have incompatible dimensions!");

        size_t sizeN = rhs.sizeM();
        size_t sizeM = Size_N;

        Matrix<T> Out(sizeN, sizeM);

        for (size_t n = 0; n < sizeN; n++)
        for (size_t m = 0; m < sizeM; m++)
        for (size_t p = 0; p < sizeM; p++)
        {
            Out(n,m) += get(n,p)*rhs.get(p,m);
        }

        return Out;
    }

    Matrix<T> operator+(const Matrix<T>& rhs)                                   ///< Returns a matrix containing the sum of two matrices.
    {
        assert(Size_N == rhs.sizeN() && "Matrices have incompatible dimensions!");
        assert(Size_M == rhs.sizeM() && "Matrices have incompatible dimensions!");

        Matrix<T> Out(*this);

        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            Out(n,m) += rhs(n,m);
        }
        return Out;
    }

    Matrix<T>& operator+=(const Matrix<T>& rhs)                                 ///< Adds rhs to the current matrix. For empty matrices assumes zeros.
    {
        if(allocated and rhs.IsAllocated())
        {
            assert(Size_N == rhs.sizeN() && "Matrices have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Matrices have incompatible dimensions!");

            for (size_t n = 0; n < Size_N; n++)
            for (size_t m = 0; m < Size_M; m++)
            {
                storage[Index(n,m)] += rhs(n,m);
            }
        }
        else if(rhs.IsAllocated())
        {
            Size_N = rhs.sizeN();
            Size_M = rhs.sizeM();

            storage = rhs.storage;

            allocated = true;
        }
        return *this;
    }

    Matrix<T>& operator=(const Matrix<T>& rhs)                                  ///< Assignment operator. Assigned the content of rhs to the current matrix.
    {
        if(allocated and rhs.IsAllocated())
        {
            assert(Size_N == rhs.sizeN() && "Matrices have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Matrices have incompatible dimensions!");

            storage = rhs.storage;
        }
        else if (rhs.IsAllocated())
        {
            Size_N = rhs.Size_N;
            Size_M = rhs.Size_M;

            storage = rhs.storage;
            allocated = true;
        }
        return *this;
    }

    std::string print(void) const                                               ///< Returns formatted string with matrix content.
    {
        std::stringstream out;
        out << std::setprecision(6)
            << std::right
            << std::setw(8);
        for(size_t n = 0; n < Size_N; n++)
        {
            out << "||";
            for(size_t m = 0; m < Size_M; m++)
            {
                out << storage[Index(n,m)] << " ";
            }
            out << "||\n";
        }
        return out.str();
    }

 protected:
 private:
    std::vector<T> storage;                                                     ///< Internal storage.
    size_t Size_N;                                                              ///< First dimension of the matrix.
    size_t Size_M;                                                              ///< Second dimension of the matrix.
    bool allocated;                                                             ///< True if internal storage is allocated.
    size_t Index(const size_t n, const size_t m) const                          ///< Matrix to internal storage index conversion
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return n*Size_M + m;
    }
};

}// namespace openphase
#endif
