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
 *                         Philipp Engels; Raphael Schiedung
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
class OP_EXPORTS Matrix                                                         ///< Matrix template class. Can handle any type of values except "bool"
{
 public:
    Matrix()
    {
        Size_N = 0;
        Size_M = 0;
        allocated = false;
    }
    Matrix(const Matrix<T>& rhs)
    {
        Size_N = rhs.sizeN();
        Size_M = rhs.sizeM();

        storage = rhs.storage;

        allocated = true;
    }
    Matrix(const size_t N, const size_t M)
    {
        Allocate(N, M);
    }
    T const& operator()(const size_t n, const size_t m) const
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }
    T& operator()(const size_t n, const size_t m)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }
    void Allocate(const size_t N, const size_t M)
    {
        Size_N = N;
        Size_M = M;
        storage.resize(Size_N*Size_M);
        allocated = true;
    }
    void Reallocate(const size_t N, const size_t M)
    {
        Size_N = N;
        Size_M = M;
        storage.resize(Size_N*Size_M);
        allocated = true;
    }
    void set(const size_t n, const size_t m, const T value)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] = value;
    }
    void set_to_value(const T value)
    {
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            storage[Index(n,m)] = value;
        }
    }
    void add(const size_t n, const size_t m, const T value)
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        storage[Index(n,m)] += value;
    }
    T get(const size_t n, const size_t m) const
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return storage[Index(n,m)];
    }
    size_t sizeN() const
    {
        return Size_N;
    }
    size_t sizeM() const
    {
        return Size_M;
    }
    T get_min() const
    {
        T min = std::numeric_limits<T>::max();
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            if (storage[Index(n,m)] < min) min = storage[Index(n,m)];
        }
        return min;
    }
    T get_max() const
    {
        T max = std::numeric_limits<T>::min();
        for (size_t n = 0; n < Size_N; n++)
        for (size_t m = 0; m < Size_M; m++)
        {
            if (storage[Index(n,m)] > max) max = storage[Index(n,m)];
        }
        return max;
    }
    bool IsNotAllocated() const
    {
        return (storage.size() == 0);
    }
    bool IsAllocated() const
    {
        return (storage.size() != 0);
    }
    std::string print(void) const
    {
        std::stringstream out;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_M; m++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(8) << storage[Index(n,m)];
            if (m == Size_M-1)
            {
                out << "||\n";
            }
            else
            {
                out << " ";
            }
        }
        return out.str();
    }
    Matrix<T> operator*(const Matrix<T>& rhs) const
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
    Matrix<T> operator+(const Matrix<T>& rhs)
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
    Matrix<T>& operator+=(const Matrix<T>& rhs)
    {
        if(allocated)
        {
            assert(Size_N == rhs.sizeN() && "Matrices have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Matrices have incompatible dimensions!");

            for (size_t n = 0; n < Size_N; n++)
            for (size_t m = 0; m < Size_M; m++)
            {
                storage[Index(n,m)] += rhs(n,m);
            }
        }
        else
        {
            Size_N = rhs.sizeN();
            Size_M = rhs.sizeM();

            storage = rhs.storage;

            allocated = true;
        }
        return *this;
    }
    Matrix<T>& operator=(const Matrix<T>& rhs)
    {
        if(allocated)
        {
            assert(Size_N == rhs.sizeN() && "Matrices have incompatible dimensions!");
            assert(Size_M == rhs.sizeM() && "Matrices have incompatible dimensions!");

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
    T* data()
    {
        return storage.data();
    }
    static Matrix<T> max(const Matrix<T>& M1, const Matrix<T>& M2)
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
    static Matrix<T> min(const Matrix<T>& M1, const  Matrix<T>& M2)
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

 protected:
 private:
    std::vector<T> storage;
    size_t Size_N;
    size_t Size_M;
    bool allocated;
    size_t Index(const size_t n, const size_t m) const
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_M && "Access beyond storage range");

        return n*Size_M + m;
    };
};

}// namespace openphase
#endif
