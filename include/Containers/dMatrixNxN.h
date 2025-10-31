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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev; Philipp Engels
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
    dMatrixNxN()
    {
        Size_N = 0;
    };
    dMatrixNxN(const size_t N): storage(N*N,0.0)
    {
        Size_N = N;
    };
    dMatrixNxN(const size_t N, const double value): storage(N*N,value)
    {
        Size_N = N;
    };
    dMatrixNxN(const dMatrixNxN& other): storage(other.storage)
    {
        Size_N = other.Size_N;
    };
    void Allocate(const size_t N)
    {
        Size_N = N;
        storage.resize(N*N,0.0);
    };
    void Allocate(const size_t N, const double value)
    {
        Size_N = N;
        storage.resize(N*N,value);
    };
    double& operator()(const size_t n, const size_t m)
    {
        return storage[Index(n,m)];
    };
    double const& operator()(const size_t n, const size_t m) const
    {
        return storage[Index(n,m)];
    };
    void pack(std::vector<double>& buffer)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            buffer.push_back(storage[Index(n,m)]);
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            storage[Index(n,m)] = buffer[it];
            ++it;
        }
    }
    double norm(void) const
    {
        double tmp = 0.0;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp += storage[Index(n,m)]*storage[Index(n,m)];
        }
        return sqrt(tmp);
    };
    double det(void) const
    {
        double determinant = 0.0;

        for(size_t i = 0; i < Size_N; i++)
        {
            double line_product = 1.0;
            for(size_t j = 0; j < Size_N; j++)
            {
                line_product *= storage[Index((i+j)%Size_N,j)];
            }
            determinant += line_product;
        }

        for(size_t i = 0; i < Size_N; i++)
        {
            double line_product = 1.0;
            for(size_t j = 0; j < Size_N; j++)
            {
                line_product *= storage[Index((i-j+Size_N)%Size_N,j)];
            }
            determinant -= line_product;
        }
        std::cout << determinant << "\n";
        return determinant;
    };

    bool is_singular(void) const
    {
        if(fabs(det()) > DBL_EPSILON)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    dMatrixNxN operator*(const double rhs) const
    {
        dMatrixNxN tmp;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(n,m)]*rhs;
        }
        return tmp;
    };
    dMatrixNxN& operator*=(const double rhs)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            storage[Index(n,m)] *= rhs;
        }
        return *this;
    };
    dMatrixNxN operator*(const dMatrixNxN& rhs) const
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
    };
    dVectorN operator*(const dVectorN& rhs) const
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
    dMatrixNxN operator+(const dMatrixNxN& rhs) const
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(n,m)] + rhs(n,m);
        }
        return tmp;
    };
    dMatrixNxN operator-(const dMatrixNxN& rhs) const
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(n,m)] - rhs(n,m);
        }
        return tmp;
    };
    dMatrixNxN& operator+=(const dMatrixNxN& rhs)
    {
        assert(Size_N == rhs.size());

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            storage[Index(n,m)] += rhs(n,m);
        }
        return *this;
    };
    dMatrixNxN operator/(const double rhs) const
    {
        dMatrixNxN tmp;
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(n,m)]/rhs;
        }
        return tmp;
    };
    dMatrixNxN& operator/=(const double rhs)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            storage[Index(n,m)] /= rhs;
        }
        return *this;
    };
    dMatrixNxN& operator-=(const dMatrixNxN& rhs)
    {
        assert(Size_N == rhs.size());

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            storage[Index(n,m)] -= rhs(n,m);
        }
        return *this;
    };
    dMatrixNxN& operator=(const dMatrixNxN& rhs)
    {
        assert(Size_N == rhs.size());

        memmove(data(), rhs.const_data(), Size_N*Size_N*sizeof(double));
        return *this;
    };
    dMatrixNxN H_product(const dMatrixNxN& rhs) const
    {
        assert(Size_N == rhs.size());

        dMatrixNxN tmp(Size_N);

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(n,m)]*rhs(n,m);
        }
        return tmp;
    };
    dMatrixNxN& invert(void)
    {
        dMatrixNxN Out(Size_N);

        for(size_t i = 0; i < Size_N; i++)
        for(size_t j = 0; j < Size_N; j++)
        {
            Out(i,j) = storage[Index(i,j)];
        }

        std::vector<size_t> indxc(Size_N);
        std::vector<size_t> indxr(Size_N);
        std::vector<size_t> ipiv(Size_N, 0);
        size_t icol = 0;
        size_t irow = 0;
        double pivinv;
        double dum;
        dMatrixNxN Uni(Size_N);
        Uni.set_to_unity();

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
                for (size_t l = 0; l < Size_N; l++)
                {
                    double temp = Uni(irow,l);
                    Uni(irow,l) = Uni(icol,l);
                    Uni(icol,l) = temp;
                };
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
            for(size_t l = 0; l < Size_N; l++) Uni(icol,l) *= pivinv;
            for(size_t ll = 0; ll < Size_N; ll++)
            if(ll != icol)
            {
                dum = Out(ll,icol);
                Out(ll,icol) = 0.0;
                for(size_t l = 0; l < Size_N; l++) Out(ll,l) -= Out(icol,l)*dum;
                for(size_t l = 0; l < Size_N; l++) Uni(ll,l) -= Uni(icol,l)*dum;
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

        for(size_t i = 0; i < Size_N; i++)
        for(size_t j = 0; j < Size_N; j++)
        {
            storage[Index(i,j)] = Out(i,j);
        }
        return *this;
    };

    dMatrixNxN inverted(void) const
    {
        dMatrixNxN Out(Size_N);

        for(size_t i = 0; i < Size_N; i++)
        for(size_t j = 0; j < Size_N; j++)
        {
            Out(i,j) = storage[Index(i,j)];
        }

        std::vector<size_t> indxc(Size_N);
        std::vector<size_t> indxr(Size_N);
        std::vector<size_t> ipiv(Size_N, 0);
        size_t icol = 0;
        size_t irow = 0;
        double pivinv;
        double dum;
        dMatrixNxN Uni(Size_N);
        Uni.set_to_unity();

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
                for (size_t l = 0; l < Size_N; l++)
                {
                    double temp = Uni(irow,l);
                    Uni(irow,l) = Uni(icol,l);
                    Uni(icol,l) = temp;
                };
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
            for(size_t l = 0; l < Size_N; l++) Uni(icol,l) *= pivinv;
            for(size_t ll = 0; ll < Size_N; ll++)
            if(ll != icol)
            {
                dum = Out(ll,icol);
                Out(ll,icol) = 0.0;
                for(size_t l = 0; l < Size_N; l++) Out(ll,l) -= Out(icol,l)*dum;
                for(size_t l = 0; l < Size_N; l++) Uni(ll,l) -= Uni(icol,l)*dum;
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
    };

    dMatrixNxN& set_to_value(double value)
    {
        for(size_t x = 0; x < Size_N*Size_N; x++)
        {
            storage[x] = value;
        }
        return *this;
    };

    dMatrixNxN& set_to_zero(void)
    {
        for(size_t x = 0; x < Size_N*Size_N; x++)
        {
            storage[x] = 0.0;
        }
        return *this;
    };

    dMatrixNxN& set_to_unity(void)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            if(n != m)
            {
                storage[Index(n,m)] = 0.0;
            }
            else
            {
                storage[Index(n,m)] = 1.0;
            }
        }
        return *this;
    };

    dMatrixNxN& transpose(void)
    {
        for(size_t n = 0; n < Size_N-1; n++)
        for(size_t m = n+1; m < Size_N; m++)
        {
            double tmp = storage[Index(n,m)];
            storage[Index(n,m)] = storage[Index(m,n)];
            storage[Index(m,n)] = tmp;
        }
        return *this;
    };

    dMatrixNxN transposed(void) const
    {
        dMatrixNxN tmp;

        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            tmp(n,m) = storage[Index(m,n)];
        }
        return tmp;
    };

    std::string print(void) const
    {
        std::stringstream out;
        for(size_t n = 0; n < Size_N; n++)
        {
            out << "|| " << std::setprecision(6);
            for(size_t m = 0; m < Size_N; m++)
            {
                out << std::setw(8) << storage[Index(n,m)] << " ";
            }
            out << "||\n";
        }
        return out.str();
    };
    void read(std::fstream& inp)
    {
        for(size_t n = 0; n < Size_N; n++)
        for(size_t m = 0; m < Size_N; m++)
        {
            inp >> storage[Index(n,m)];
        }
    };
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
    };
    double* data(void)
    {
        return storage.data();
    };
    const double* const_data(void) const
    {
        return storage.data();
    };
    size_t size(void) const
    {
        return Size_N;
    }
    size_t sizeX() const
    {
        return Size_N;
    }
    size_t sizeY() const
    {
        return Size_N;
    }
 protected:
 private:

    std::vector<double> storage;
    size_t Size_N;
    size_t Index(const size_t n, const size_t m) const
    {
        assert(n < Size_N && "Access beyond storage range");
        assert(m < Size_N && "Access beyond storage range");

        return n*Size_N + m;
    };
};

}// namespace openphase
#endif
