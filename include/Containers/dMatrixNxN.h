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

#ifndef DMATRIXNXN_H
#define DMATRIXNXN_H

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dVectorN.h"

namespace openphase
{

class dMatrixNxN
{
 public:
    dMatrixNxN() : Size_N(0)                                                    ///< Default constructor. Creates empty matrix.
    {
    }

    dMatrixNxN(const size_t N): storage(N*N,0.0), Size_N(N)                     ///< Constructor. Allocates matrix container of given size NxN.
    {
    }

    dMatrixNxN(const size_t N, const double value): storage(N*N,value), Size_N(N)///< Constructor. Allocates matrix container of given size NxN and initializes all of its elements with the given value.
    {
    }

    dMatrixNxN(const dMatrixNxN& rhs): storage(rhs.storage), Size_N(rhs.Size_N) ///< Copy constructor. Creates the matrix and initializes it with the copy of the rhs.
    {
    }

    void Allocate(const size_t N)                                               ///< Allocates storage of empty container to the given size NxN.
    {
        Size_N = N;
        storage.resize(N*N,0.0);
    }

    void Allocate(const size_t N, const double value)                           ///< Allocates storage of empty container to given size NxN and initializes all of its elements with the given value.
    {
        Size_N = N;
        storage.resize(N*N,value);
    }

    double& operator()(const size_t n, const size_t m)                          ///< Random access operator. Returns the reference to the element pointed to by the indices n,m.
    {
        return storage[Index(n,m)];
    }

    double const& operator()(const size_t n, const size_t m) const              ///< Random access operator. Returns const reference to the element pointed to by the indices n,m.
    {
        return storage[Index(n,m)];
    }

    double& operator[](const size_t idx)                                        ///< Random access operator. Returns the reference to the internal flattened storage element pointed to by the idx.
    {
        return storage[idx];
    }

    double const& operator[](const size_t idx) const                            ///< Random access operator. Returns const reference to the internal flattened storage element pointed to by the idx.
    {
        return storage[idx];
    }

    double* data(void)                                                          ///< Returns the pointer to the stored data.
    {
        return storage.data();
    }

    const double* const_data(void) const                                        ///< Returns const pointer to the stored data.
    {
        return storage.data();
    }

    double norm(void) const                                                     ///< Returns Frobenius norm of the matrix.
    {
        double tmp = 0.0;
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp += storage[p]*storage[p];
        }
        return sqrt(tmp);
    }

    double determinant(void) const                                              ///< Returns determinant of the matrix.
    {
        double det = 0.0;

        for(size_t i = 0; i < Size_N; i++)
        {
            double line_product = 1.0;
            for(size_t j = 0; j < Size_N; j++)
            {
                line_product *= storage[Index((i+j)%Size_N,j)];
            }
            det += line_product;
        }

        for(size_t i = 0; i < Size_N; i++)
        {
            double line_product = 1.0;
            for(size_t j = 0; j < Size_N; j++)
            {
                line_product *= storage[Index((i-j+Size_N)%Size_N,j)];
            }
            det -= line_product;
        }
        return det;
    }

    bool is_singular(void) const                                                ///< Returns true if the current matrix is singular.
    {
        return fabs(determinant()) < DBL_EPSILON;
    }

    dMatrixNxN operator*(const double rhs) const                                ///< Returns the copy of the current matrix with all its elements multiplied by rhs.
    {
        dMatrixNxN tmp(*this);
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp.storage[p] *= rhs;
        }
        return tmp;
    }

    dMatrixNxN& operator*=(const double rhs)                                    ///< Multiplies all entries of the current matrix by the rhs.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            storage[p] *= rhs;
        }
        return *this;
    }

    dMatrixNxN operator*(const dMatrixNxN& rhs) const                           ///< Returns the product of the current matrix and the rhs matrix.
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp;
        tmp.set_to_zero();
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        for(size_t p = 0; p < Size_N; p++)
        {
            tmp(n,m) += storage[Index(n,p)]*rhs(p,m);
        }
        return tmp;
    }

    dVectorN operator*(const dVectorN& rhs) const                               ///< Returns the product of the current matrix and the rhs vector.
    {
        assert(Size_N == rhs.size());

        dVectorN tmp(Size_N, 0.0);

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp[n] += storage[Index(n,m)]*rhs[m];
        }
        return tmp;
    }

    dMatrixNxN operator+(const dMatrixNxN& rhs) const                           ///< Returns the sum of the current matrix and the rhs matrix.
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp(*this);
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp.storage[p] += rhs.storage[p];
        }
        return tmp;
    }

    dMatrixNxN operator-(const dMatrixNxN& rhs) const                           ///< Returns the difference of the current matrix and the rhs matrix.
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp(*this);
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp.storage[p] -= rhs.storage[p];
        }
        return tmp;
    }

    dMatrixNxN& operator+=(const dMatrixNxN& rhs)                               ///< Adds rhs matrix to the current matrix.
    {
        assert(Size_N == rhs.size());

        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            storage[p] += rhs.storage[p];
        }
        return *this;
    }

    dMatrixNxN& operator-=(const dMatrixNxN& rhs)                               ///< Subtracts rhs matrix from the current matrix.
    {
        assert(Size_N == rhs.size());

        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            storage[p] -= rhs.storage[p];
        }
        return *this;
    }

    dMatrixNxN operator/(const double rhs) const                                ///< Returns the copy of the current matrix with all its elements divides by the scalar rhs.
    {
        dMatrixNxN tmp(*this);
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp.storage[p] /= rhs;
        }
        return tmp;
    }

    dMatrixNxN& operator/=(const double rhs)                                    ///< Divides all elements of the current matrix by the scalar rhs.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            storage[p] /= rhs;
        }
        return *this;
    }

    dMatrixNxN& operator=(const dMatrixNxN& rhs)                                ///< Assignment operator. Replaces the content of the current matrix with the copy of the rhs.
    {
        assert(Size_N == rhs.size());

        storage = rhs.storage;
        return *this;
    }

    dMatrixNxN Hadamard_product(const dMatrixNxN& rhs) const                    ///< Returns Hadamard (component-wise) product of the current matrix and the rhs.
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp(*this);
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            tmp.storage[p] *= rhs.storage[p];
        }
        return tmp;
    }

    dMatrixNxN inverted(void) const                                             ///< Returns inverted copy of the current matrix.
    {
        dMatrixNxN Out(*this);

        std::vector<size_t> indxc(Size_N);
        std::vector<size_t> indxr(Size_N);
        std::vector<size_t> ipiv(Size_N, 0);
        size_t icol = 0;
        size_t irow = 0;
        double pivinv;
        double dum;
        //dMatrixNxN Uni(Size_N);
        //Uni.set_to_unity();

        for(size_t i = 0; i < Size_N; i++)
        {
            double big = 0.0;
            for(size_t j = 0; j < Size_N; j++)
            if(ipiv[j] != 1)
            for(size_t k = 0; k < Size_N; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out(j,k)) >= big)
                    {
                        big = fabs(Out(j,k));
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cerr << "dMatrixNxN: Can Not Compute Inverse Matrix."
                              << this->print() << "Matrix is Singular 1!!!\n";
                    OP_Exit(EXIT_FAILURE);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (size_t l = 0; l < Size_N; l++)
                {
                    double temp = Out(irow,l);
                    Out(irow,l) = Out(icol,l);
                    Out(icol,l) = temp;
                };
//                for (size_t l = 0; l < Size_N; l++)
//                {
//                    double temp = Uni(irow,l);
//                    Uni(irow,l) = Uni(icol,l);
//                    Uni(icol,l) = temp;
//                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out(icol,icol)) <= DBL_EPSILON)
            {
                std::cerr << "dMatrixNxN: Can Not Compute Inverse Matrix.\n"
                          << this->print() << "Matrix is Singular 2!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
            pivinv = 1.0/Out(icol,icol);
            Out(icol,icol) = 1.0;
            for(size_t l = 0; l < Size_N; l++) Out(icol,l) *= pivinv;
            //for(size_t l = 0; l < Size_N; l++) Uni(icol,l) *= pivinv;
            for(size_t ll = 0; ll < Size_N; ll++)
            if(ll != icol)
            {
                dum = Out(ll,icol);
                Out(ll,icol) = 0.0;
                for(size_t l = 0; l < Size_N; l++) Out(ll,l) -= Out(icol,l)*dum;
                //for(size_t l = 0; l < Size_N; l++) Uni(ll,l) -= Uni(icol,l)*dum;
            }
        }
        for(int l = Size_N - 1; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(size_t k = 0; k < Size_N; k++)
            {
                double temp = Out(k,indxr[l]);
                Out(k,indxr[l]) = Out(k,indxc[l]);
                Out(k,indxc[l]) = temp;
            };
        }
        return Out;
    }

    dMatrixNxN& invert(void)                                                    ///< Inverts the current matrix.
    {
        dMatrixNxN tmp = inverted();
        storage = tmp.storage;
        return *this;
    }

    dMatrixNxN& set_to_value(double value)                                      ///< Sets all elements of the current matrix to the given value
    {
        std::fill(storage.begin(),storage.end(),value);
        return *this;
    }

    dMatrixNxN& set_to_zero(void)                                               ///< Sets all elements of the current matrix to zero.
    {
        std::fill(storage.begin(),storage.end(),0.0);
        return *this;
    }

    dMatrixNxN& set_to_unity(void)                                              ///< Sets all elements of the current matrix to one.
    {
        std::fill(storage.begin(),storage.end(),0.0);
        for(size_t p = 0; p < Size_N; p++)
        {
            storage[Index(p,p)] = 1.0;
        }
        return *this;
    }

    dMatrixNxN transposed(void) const                                           ///< Returns transposed copy of the current matrix.
    {
        dMatrixNxN tmp;

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(m,n)];
        }
        return tmp;
    }

    dMatrixNxN& transpose(void)                                                 ///< Transposes the current matrix.
    {
        dMatrixNxN tmp = transposed();
        storage = tmp.storage;
        return *this;
    }

    size_t size(void) const                                                     ///< Returns the size of the internal flattened storage.
    {
        return storage.size();
    }

    size_t sizeX() const                                                        ///< Returns the size of the first matrix dimension.
    {
        return Size_N;
    }

    size_t sizeY() const                                                        ///< Returns the size of the second matrix dimension.
    {
        return Size_N;
    }

    void pack(std::vector<double>& buffer) const                                ///< Packs (writes) matrix content into the MPI communication buffer.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            buffer.push_back(storage[p]);
        }
    }

    void unpack(std::vector<double>& buffer, size_t& it)                        ///< Unpacks (reads) matrix content from the MPI communication buffer.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            storage[p] = buffer[it]; it++;
        }
    }

    std::string print(void) const                                               ///< Returns formatted string containing matrix elements for printing.
    {
        std::stringstream out;
        out << std::setprecision(6);
        for(size_t n = 0; n < Size_N; n++)
        {
            out << "||";
            for(size_t m = 0; m < Size_N; m++)
            {
                out << std::setw(8) << storage[Index(n,m)] << " ";
            }
            out << "||\n";
        }
        return out.str();
    }

    void read_binary(std::istream& inp)                                         ///< Reads matrix content from the binary input stream inp.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            inp >> storage[p];
        }
    }

    void read_ASCII(std::istream& inp)                                          ///< Reads matrix content from the ASCII input stream inp.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            inp >> storage[p];
        }
    }

    void write_binary(std::ostream& outp) const                                 ///< Writes matrix content to the binary output stream outp.
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            outp << storage[p];
        }
    }

    void write_ASCII(std::ostream& outp, const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const///< Writes matrix content to the ASCII output stream outp.
    {
        outp << std::setprecision(precision) << std::defaultfloat;
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            outp << storage[p] << sep;
        }
        outp << std::endl;
    }

    std::string get_output_string(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const ///< Returns formatted string with the matrix content
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
           out << storage[p] << sep;
        }
        return out.str();
    }

    std::vector<double> get_vector() const                                      ///< Returns vector with flattened matrix content.
    {
        return storage;
    }

    std::vector<float> get_vector_float() const                                 ///< Returns vector with float accuracy with flattened matrix content.
    {
        std::vector<float> out;
        for(size_t n = 0; n < Size_N*Size_N; n++)
        {
            out.push_back((float)storage[n]);
        }
        return out;
    }

    //============================= Deprecated methods =========================
    [[deprecated]]
    void read(std::fstream& inp)
    {
        for(size_t p = 0; p < Size_N*Size_N; p++)
        {
            inp >> storage[p];
        }
    }

    [[deprecated]]
    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::scientific;
        for(size_t n = 0; n < Size_N; n++)
        {
           for(size_t m = 0; m < Size_N; m++)
           {
               out << storage[Index(n,m)] << sep;
           }
           out << "\n";
        }
        return out.str();
    }

 protected:
 private:

    std::vector<double> storage;
    size_t Size_N;
    size_t Index(const size_t n, const size_t m) const                          ///< Converts matrix indices into flattened linear storage index.
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_N && "Access beyond storage range");

        return n*Size_N + m;
    };
};

}// namespace openphase
#endif
