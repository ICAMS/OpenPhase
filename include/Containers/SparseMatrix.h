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

 *   File created :   2019
 *   Main contributors :   Marvin Tegeler, Raphael Schiedung
 *
 */

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <map>
#include <fenv.h>
#include <vector>
#include <math.h>

#include "Globals.h"

namespace openphase
{

struct OP_EXPORTS SparseMatrix
{
    SparseMatrix(){};
    SparseMatrix(const size_t n){Allocate(n);};                                 ///< Allocates n rows
    SparseMatrix(const size_t n, const size_t m){Allocate(n,m);};               ///< Allocates n rows and reserves memory for m entries

    void Allocate(const size_t n);                                              ///< Allocates n rows
    void Allocate(const size_t n, const size_t m);                              ///< Allocates n rows and reserves memory for m entries

    SparseMatrix  operator*(const double& scalar) const;                        ///< Multiplies matrix by scalar value
    SparseMatrix& operator*(const double& scalar);                              ///< Multiplies matrix by scalar value

    void Sort();                                                                ///< Sorts storage by index
    void clear();
protected:
    struct entry
    {
        size_t index;
        double value;
    };
    std::vector<std::vector<entry>> storage;
    size_t size;
    size_t reserve;
};

struct SparseMatrixCSC;
struct OP_EXPORTS SparseMatrixCSR: public SparseMatrix                                     ///< Compressed sparse row
{
    friend SparseMatrixCSC;

    SparseMatrixCSR(){};
    SparseMatrixCSR(const size_t n): SparseMatrix(n){};                         ///< Allocates n rows
    SparseMatrixCSR(const size_t n, const size_t m): SparseMatrix(n,m){};       ///< Allocates n rows and reserves memory of m entries per row

    using SparseMatrix::Allocate;

    double  operator()(const size_t i, const size_t j) const;                   ///< returns matrix element ij
    double& operator()(const size_t i, const size_t j);                         ///< returns matrix element ij

    using SparseMatrix::operator*;
    SparseMatrixCSR operator*(const SparseMatrixCSR& rMatrix) const;            ///< Multiplies matrix by scalar value
    template<class Vector> Vector operator*(const Vector& vector) const         ///< Multiplies matrix by vector
    {
        Vector tmp(size,0);
        #pragma omp parallel for
        for (size_t i = 0; i < size; ++i)
        {
            for (auto j = storage[i].cbegin(); j != storage[i].cend(); ++j)
            {
                tmp[i] += j->value*vector[j->index];
            }
        }
        return tmp;
    }
    SparseMatrixCSR& operator=(const SparseMatrixCSC& rhs);
    friend SparseMatrixCSR operator*(const SparseMatrixCSR& lhs, const SparseMatrixCSR& rhs); ///< Multiplies matrix by scalar value
    friend bool        operator==(const SparseMatrixCSR& lhs, const SparseMatrixCSC& rhs);
    friend bool inline operator==(const SparseMatrixCSC& lhs, const SparseMatrixCSR& rhs){return   rhs == lhs;};
    friend bool inline operator!=(const SparseMatrixCSR& lhs, const SparseMatrixCSC& rhs){return !(lhs == rhs);};
    friend bool inline operator!=(const SparseMatrixCSC& lhs, const SparseMatrixCSR& rhs){return !(lhs == rhs);};

    SparseMatrixCSR transposed() const;                                         ///< Returns Transposed Matrix
    SparseMatrixCSR& transpose();                                               ///< Transposes Matrix
    size_t MaxRow(const size_t i, const size_t j);                              ///< Returns row index of max value in column j starting for row i
    size_t NonZeroRow(const size_t i, const size_t j);                          ///< Returns row index of first non zero value in column j starting for row i
    void DivideRowBy(size_t i, double scalar);                                  ///< Divides row i by scalar
    void SwapRows(const size_t i, const size_t j);                              ///< Swaps rows i and j
};

struct OP_EXPORTS SparseMatrixCSC: public SparseMatrix                                     ///< Compressed sparse column
{
    friend SparseMatrixCSR;

    SparseMatrixCSC(){};
    SparseMatrixCSC(const size_t n): SparseMatrix(n){};                         ///< Allocates n columns
    SparseMatrixCSC(const size_t n, const size_t m): SparseMatrix(n,m){};       ///< Allocates n columns and reserves memory of m entries per columns

    using SparseMatrix::Allocate;

    double  operator()(const size_t i, const size_t j) const;                   ///< returns matrix element ij
    double& operator()(const size_t i, const size_t j);                         ///< returns matrix element ij

    using SparseMatrix::operator*;
    template<class Vector> Vector operator*(const Vector& vector) const         ///< Multiplies matrix by vector
    {
        Vector tmp(size,0);
        #pragma omp parallel for
        for (size_t j = 0; j < size; ++j)
        {
            for (auto i = storage[j].cbegin(); i != storage[j].cend(); ++i)
            {
                tmp[i->index] += i->value*vector[j];
            }
        }
        return tmp;
    }
    SparseMatrixCSC& operator=(const SparseMatrixCSR& rhs);
    friend SparseMatrixCSC  operator*(const SparseMatrixCSC& lhs, const SparseMatrixCSC& rhs); ///< Multiplies matrix by scalar value
    friend bool        operator==(const SparseMatrixCSR& lhs, const SparseMatrixCSC& rhs);
    friend bool inline operator==(const SparseMatrixCSC& lhs, const SparseMatrixCSR& rhs);
    friend bool inline operator!=(const SparseMatrixCSR& lhs, const SparseMatrixCSC& rhs);
    friend bool inline operator!=(const SparseMatrixCSC& lhs, const SparseMatrixCSR& rhs);

    SparseMatrixCSC transposed() const;                                         ///< Returns Transposed Matrix
    SparseMatrixCSC& transpose();                                               ///< Transposes Matrix
    size_t MaxRow(const size_t i, const size_t j);                              ///< Returns row index of max value in column j starting for row i
    size_t NonZeroRow(const size_t i, const size_t j);                          ///< Returns row index of first non zero value in column j starting for row i
    void DivideRowBy(size_t i, double scalar);                                  ///< Divides row i by scalar
    void SwapRows(const size_t i, const size_t j);                              ///< Swaps rows i and j
};

// NOTE commented class because it does not work
//class BiCGStab
//{
//public:
//    typedef SparseMatrix Matrix;
//    typedef std::vector<double> Vector;
//
//    static void BiCGstab_Jacobi(const Matrix& A, const Vector& b, Vector& x, const size_t l, const double TOL);
//};

}// namespace openphase
#endif
