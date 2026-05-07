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

#ifndef DMATRIX6X6_H
#define DMATRIX6X6_H

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

#include "dMatrix3x3.h"
#include "dVector3.h"

namespace openphase
{

class dMatrix6x6
{
 public:
    dMatrix6x6() : storage()
    {
    }

    dMatrix6x6(const dMatrix6x6& rhs) : storage(rhs.storage)
    {
    }

    double& operator()(const size_t i, const size_t j)
    {
        assert(i < 6 && "Access beyond storage range");
        assert(j < 6 && "Access beyond storage range");

        return storage[Index(i,j)];
    }

    const double& operator()(const size_t i, const size_t j) const
    {
        assert(i < 6 && "Access beyond storage range");
        assert(j < 6 && "Access beyond storage range");

        return storage[Index(i,j)];
    }

    double& operator[](const size_t idx)
    {
        assert(idx < 36 && "Access beyond storage range");

        return storage[idx];
    }

    const double& operator[](const size_t idx) const
    {
        assert(idx < 36 && "Access beyond storage range");

        return storage[idx];
    }

    dMatrix6x6& set_to_zero(void)
    {
        storage.fill(0.0);
        return *this;
    }

    dMatrix6x6& set_to_unity(void)
    {
        set_to_zero();
        for(int i = 0; i < 6; i++)
        {
            storage[Index(i,i)] = 1.0;
        }
        return *this;
    }

    dMatrix6x6& set_to_value(double value)
    {
        for (int n = 0; n < 36; n++)
        {
            storage[n] = value;
        }
        return *this;
    }

    double norm(void) const
    {
        double tmp = 0.0;
        for (int n = 0; n < 36; n++)
        {
            tmp += std::pow(storage[n], 2.0);
        }
        return sqrt(tmp);
    }

    double determinant(void) const
    {
        double determinant = 0.0;

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[Index((i+j)%6,j)];
            }
            determinant += line_product;
        }

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[Index((i-j+6)%6,j)];
            }
            determinant -= line_product;
        }
        return determinant;
    }

//    bool is_singular(void)
//    {
//        bool result = true;
//        if(fabs(determinant()) > DBL_EPSILON)
//        {
//            result = false;
//        }
//        return result;
//    }

    dMatrix6x6 operator*(const double m) const
    {
        dMatrix6x6 tmp(*this);
        for(int n = 0; n < 36; n++)
        {
            tmp.storage[n] *= m;
        }
        return tmp;
    }

    dMatrix6x6& operator*=(const double m)
    {
        for(int n = 0; n < 36; n++)
        {
            storage[n] *= m;
        }
        return *this;
    }

    dMatrix6x6 operator*(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        for(int k = 0; k < 6; k++)
        {
            tmp(i,j) += storage[Index(i,k)]*rhs(k,j);
        }
        return tmp;
    }

    dMatrix6x6 operator+(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp(*this);
        for(int n = 0; n < 36; n++)
        {
            tmp.storage[n] += rhs.storage[n];
        }
        return tmp;
    }

    dMatrix6x6& operator+=(const dMatrix6x6& rhs)
    {
        for(int n = 0; n < 36; n++)
        {
            storage[n] += rhs.storage[n];
        }
        return *this;
    }

    dMatrix6x6 operator-(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp(*this);
        for(int n = 0; n < 36; n++)
        {
            tmp.storage[n] -= rhs.storage[n];
        }
        return tmp;
    }

    dMatrix6x6& operator-=(const dMatrix6x6& rhs)
    {
        for(int n = 0; n < 36; n++)
        {
            storage[n] -= rhs.storage[n];
        }
        return *this;
    }

    dMatrix6x6 operator/(const double m) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[Index(i,j)]/m;
        }
        return tmp;
    }

    dMatrix6x6& operator/=(const double m)
    {
        for(int n = 0; n < 36; n++)
        {
            storage[n] /= m;
        }
        return *this;
    }

    dMatrix6x6& operator=(const dMatrix6x6& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    bool operator==(dMatrix6x6 rhs) const
    {
        bool result = true;
        for(int n = 0; n < 36; n++)
        {
            if(fabs(storage[n] - rhs.storage[n]) > DBL_EPSILON)
            {
                result = false;
                break;
            }
        }
        return result;
    }

    bool operator!=(dMatrix6x6 rhs) const
    {
        return !(*this == rhs);
    }

    dMatrix6x6 Hadamard_product(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp(*this);

        for(int n = 0; n < 36; n++)
        {
            tmp.storage[n] *= rhs.storage[n];
        }
        return tmp;
    }

    dMatrix6x6 inverted(void) const
    {
        dMatrix6x6 Out(*this);

        std::array<int,6> indxc{};
        std::array<int,6> indxr{};
        std::array<int,6> ipiv {};

        int icol = 0;
        int irow = 0;
        double pivinv;
        double dum;

        for(int i = 0; i < 6; i++)
        {
            double big = 0.0;
            for(int j = 0; j < 6; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < 6; k++)
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
                    std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix."
                              << this->print() << "Matrix is Singular 1!!!\n";
                    OP_Exit(EXIT_FAILURE);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    double temp = Out(irow,l);
                    Out(irow,l) = Out(icol,l);
                    Out(icol,l) = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out(icol,icol)) <= DBL_EPSILON)
            {
                std::cerr << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                          << this->print() << "Matrix is Singular 2!!!\n";
                OP_Exit(EXIT_FAILURE);
            }
            pivinv = 1.0/Out(icol,icol);
            Out(icol,icol) = 1.0;
            for(int l = 0; l < 6; l++) Out(icol,l) *= pivinv;
            for(int ll = 0; ll < 6; ll++)
            if(ll != icol)
            {
                dum = Out(ll,icol);
                Out(ll,icol) = 0.0;
                for(int l = 0; l < 6; l++) Out(ll,l) -= Out(icol,l)*dum;
            }
        }

        for(int l = 5; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < 6; k++)
            {
                double temp = Out(k,indxr[l]);
                Out(k,indxr[l]) = Out(k,indxc[l]);
                Out(k,indxc[l]) = temp;
            };
        }
        return Out;
    }

    dMatrix6x6& invert(void)
    {
        dMatrix6x6 tmp = inverted();
        storage = tmp.storage;
        return *this;
    }

    dMatrix6x6 transposed(void) const
    {
        dMatrix6x6 Out;

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            Out(i,j) = storage[Index(j,i)];
        }
        return Out;
    }
    dMatrix6x6& transpose(void)
    {
        dMatrix6x6 tmp = transposed();
        storage = tmp.storage;
        return *this;
    }

    dMatrix6x6 rotated(const dMatrix3x3& RotationMatrix)
    {
        dMatrix6x6 Out;
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};

        for(int a = 0; a < 6; a++)
        for(int b = 0; b < 6; b++)
        {
            int m = VoigtIndex[a][0];
            int n = VoigtIndex[a][1];
            int p = VoigtIndex[b][0];
            int q = VoigtIndex[b][1];

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
                // active rotation
                Out(a,b) += RotationMatrix(m,i)*
                            RotationMatrix(n,j)*
                            tensor(i,j,k,l)*
                            RotationMatrix(p,k)*
                            RotationMatrix(q,l);
            }
        }

        return Out;
    }

    /*
    dMatrix6x6 rotated(const dMatrix3x3& RotationMatrix)
    {
        double Out[3][3][3][3];

        for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++)
        for(int p = 0; p < 3; p++)
        for(int q = 0; q < 3; q++)
        {
            Out[m][n][p][q] = 0.0;

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
                // active rotation
                Out[m][n][p][q] += RotationMatrix(m,i)*
                                   RotationMatrix(n,j)*
                                   tensor(i,j,k,l)*
                                   RotationMatrix(p,k)*
                                   RotationMatrix(q,l);
            }
        }
        dMatrix6x6 OUT;
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};

        for(int m = 0; m < 6; m++)
        for(int n = 0; n < 6; n++)
        {
            int i = VoigtIndex[m][0];
            int j = VoigtIndex[m][1];
            int k = VoigtIndex[n][0];
            int l = VoigtIndex[n][1];

            OUT(m,n) = Out[i][j][k][l];
        }
        return OUT;
    }*/

    dMatrix6x6& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix6x6 tmp = rotated(RotationMatrix);
        storage = tmp.storage;
        return *this;
    }

    dMatrix3x3 project(const dVector3& n)
    {
        dMatrix3x3 Out;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
        {
            Out(j,k) += n[i]*tensor(i,j,k,l)*n[l];
        }
        return Out;
    }

    double& tensor(const int i, const int j, const int k, const int l)
    {
        return storage[Index((i==j)?(i):(6-(i+j)),
                             (k==l)?(k):(6-(k+l)))];
    }

    const double& tensor(const int i, const int j, const int k, const int l) const
    {
        return storage[Index((i==j)?(i):(6-(i+j)),
                             (k==l)?(k):(6-(k+l)))];
    }

    static dMatrix6x6 UnitTensor(void)
    {
        dMatrix6x6 myUnitTensor;
        return myUnitTensor.set_to_unity();
    }

    static dMatrix6x6 ZeroTensor(void)
    {
        dMatrix6x6 myZeroTensor;
        return myZeroTensor.set_to_zero();
    }

    constexpr size_t sizeX() const
    {
        return 6u;
    }

    constexpr size_t sizeY() const
    {
        return 6u;
    }

    constexpr size_t size() const                                               ///< Returns the size of the internal flattened storage.
    {
        return 36u;
    }

    double* data(void)
    {
        return storage.data();
    }

    const double* data(void) const
    {
        return storage.data();
    }

    void pack(std::vector<double>& buffer)
    {
        for(int n = 0; n < 36; n++)
        {
            buffer.push_back(storage[n]);
        }
    }

    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for(int n = 0; n < 36; n++)
        {
            storage[n] = buffer[it]; ++it;
        }
    }

    std::string print(void) const
    {
        std::stringstream out;
        out << std::setprecision(6);
        for(int i = 0; i < 6; i++)
        {
            out << "||";
            for(int j = 0; j < 6; j++)
            {
                out << std::setw(10) << storage[Index(i,j)] << " ";
            }
            out << "||\n";
        }
        return out.str();
    }

    void read_binary(std::istream& inp)
    {
        for(int n = 0; n < 36; n++)
        {
            inp >> storage[n];
        }
    }

    void read_ASCII(std::istream& inp)
    {
        for(int n = 0; n < 36; n++)
        {
            inp >> storage[n];
        }
    }

    void write_binary(std::ostream& outp) const
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            outp << storage[Index(i,j)];
        }
    }

    void write_ASCII(std::ostream& outp, const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        outp << std::setprecision(precision) << std::defaultfloat;
        for(int n = 0; n < 36; n++)
        {
            outp << storage[n] << sep;
        }
        outp << std::endl;
    }

    std::string get_output_string(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int n = 0; n < 36; n++)
        {
           out << storage[n] << sep;
        }
        return out.str();
    }

    std::vector<double> get_vector() const
    {
        std::vector<double> out;
        for(int n = 0; n < 36; n++)
        {
            out.push_back(storage[n]);
        }
        return out;
    }

    std::vector<float> get_vector_float() const
    {
        std::vector<float> out;
        for(int n = 0; n < 36; n++)
        {
            out.push_back((float)storage[n]);
        }
        return out;
    }

    //============================ Deprecated methods ==========================
    [[deprecated]]
    void read(std::istream& inp)
    {
        for(int n = 0; n < 36; n++)
        {
            inp >> storage[n];
        }
    }

    [[deprecated]]
    std::string write(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 6; i++)
        {
           for(int j = 0; j < 6; j++)
           {
               out << storage[Index(i,j)] << sep;
           }
           out << "\n";
        }
        return out.str();
    }

    [[deprecated]]
    void read(std::fstream& inp)
    {
        for(int n = 0; n < 36; n++)
        {
            inp >> storage[n];
        }
    }

    [[deprecated]]
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        for(int n = 0; n < 36; n++)
        {
            out.push_back(storage[n]);
        }
        return out;
    }

    [[deprecated]]
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        for(int n = 0; n < 36; n++)
        {
            out.push_back((float)storage[n]);
        }
        return out;
    }

    [[deprecated]]
    void read(std::stringstream& inp)
    {
        for(int n = 0; n < 36; n++)
        {
            inp >> storage[n];
        }
    };
    //==========================================================================

 protected:
 private:
    std::array<double,36> storage;
    size_t Index(const size_t x, const size_t y) const
    {
        return 6*x + y;
    }
};

}// namespace openphase
#endif
