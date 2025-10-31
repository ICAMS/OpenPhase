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
 *                         Philipp Engels
 *
 */

#ifndef DMATRIX3X3_H
#define DMATRIX3X3_H

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

#include "Globals.h"

namespace openphase
{

class OP_EXPORTS dMatrix3x3
{
 public:
    dMatrix3x3():
        storage({0,0,0,0,0,0,0,0,0})
    {

    }
    dMatrix3x3(const dMatrix3x3& rhs):
        storage(rhs.storage)
    {

    }
    dMatrix3x3(const std::initializer_list<double> matinit)
    {
        assert(matinit.size() == 9 && "Initialization list size is not equal to storage range.");

        int ii = 0;
        for (auto it = matinit.begin(); it != matinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }
    double& operator()(const size_t i, const size_t j)
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[idx(i,j)];
    }
    double const& operator()(const size_t i, const size_t j) const
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[idx(i,j)];
    }
    dMatrix3x3& set_to_zero(void)
    {
        storage.fill(0);
        return *this;
    }
    dMatrix3x3& set_to_unity(void)
    {
        storage.fill(0);
        storage[idx(0,0)] = 1.0;
        storage[idx(1,1)] = 1.0;
        storage[idx(2,2)] = 1.0;
        return *this;
    }
    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            buffer.push_back(storage[idx(i,j)]);
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            storage[idx(i,j)] = buffer[it]; ++it;
        }
    }
    void set(double r00, double r01, double r02,
             double r10, double r11, double r12,
             double r20, double r21, double r22)
    {
        storage[idx(0,0)] = r00;
        storage[idx(0,1)] = r01;
        storage[idx(0,2)] = r02;
        storage[idx(1,0)] = r10;
        storage[idx(1,1)] = r11;
        storage[idx(1,2)] = r12;
        storage[idx(2,0)] = r20;
        storage[idx(2,1)] = r21;
        storage[idx(2,2)] = r22;
    }
    dMatrix3x3& operator=(const dMatrix3x3& rhs)
    {
        storage = rhs.storage;
        return *this;
    }
    bool operator==(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(fabs(storage[idx(i,j)]-rhs(i,j)) > DBL_EPSILON)
            {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(fabs(storage[idx(i,j)]-rhs(i,j)) > DBL_EPSILON)
            {
                return true;
            }
        }
        return false;
    }
    dMatrix3x3 operator*(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[idx(i,j)]*m;
        }
        return tmp;
    }
    dMatrix3x3 operator/(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[idx(i,j)]/m;
        }
        return tmp;
    }
    dMatrix3x3 operator*(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        {
            tmp(i,j) += storage[idx(i,k)]*rhs(k,j);
        }
        return tmp;
    }
    dMatrix3x3 operator+(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[idx(i,j)] + rhs(i,j);
        }
        return tmp;
    }
    dMatrix3x3 operator-(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[idx(i,j)] - rhs(i,j);
        }
        return tmp;
    }
    dMatrix3x3& operator*=(const dMatrix3x3& rhs)
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        {
            tmp(i,j) += storage[idx(i,k)]*rhs(k,j);
        }

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] = tmp(i,j);
        }
        return *this;
    }
    dMatrix3x3& operator+=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] += rhs(i,j);
        }
        return *this;
    }
    dMatrix3x3& operator-=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] -= rhs(i,j);
        }
        return *this;
    }
    dMatrix3x3& operator*=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] *= m;
        }
        return *this;
    }
    dMatrix3x3& operator/=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] /= m;
        }
        return *this;
    }
    static constexpr size_t size()
    {
        return 9u;
    }
    double& operator[](const size_t i)
    {
        assert(i < size() && "Access beyond storage range");

        int column = i % 3;
        int row = int((i - column)/3);
        return storage[idx(row,column)];
    }
    const double& operator[](const size_t i) const
    {
        assert(i < size() && "Access beyond storage range");

        int column = i % 3;
        int row = int((i - column)/3);
        return storage[idx(row,column)];
    }
    dMatrix3x3 Hadamard_product(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[idx(i,j)]*rhs(i,j);
        }
        return tmp;
    }
    inline double determinant(void) const
    {
        return (storage[idx(0,0)]*storage[idx(1,1)]*storage[idx(2,2)] +
                storage[idx(0,1)]*storage[idx(1,2)]*storage[idx(2,0)] +
                storage[idx(0,2)]*storage[idx(1,0)]*storage[idx(2,1)] -
                storage[idx(0,2)]*storage[idx(1,1)]*storage[idx(2,0)] -
                storage[idx(0,1)]*storage[idx(1,0)]*storage[idx(2,2)] -
                storage[idx(0,0)]*storage[idx(1,2)]*storage[idx(2,1)]);
    }
    dMatrix3x3& invert(void)
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0)
        {
            detInv = 1.0/detInv;
        }
        else
        {
            std::cerr << "dMatrix3x3: Can Not Compute Inverse Matrix.\n"
                      << this->print() << "Matrix is Singular !!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        tmp(0,0) = (storage[idx(1,1)]*storage[idx(2,2)] -
                    storage[idx(1,2)]*storage[idx(2,1)])*detInv;
        tmp(1,0) =-(storage[idx(1,0)]*storage[idx(2,2)] -
                    storage[idx(1,2)]*storage[idx(2,0)])*detInv;
        tmp(2,0) = (storage[idx(1,0)]*storage[idx(2,1)] -
                    storage[idx(1,1)]*storage[idx(2,0)])*detInv;
        tmp(0,1) =-(storage[idx(0,1)]*storage[idx(2,2)] -
                    storage[idx(0,2)]*storage[idx(2,1)])*detInv;
        tmp(1,1) = (storage[idx(0,0)]*storage[idx(2,2)] -
                    storage[idx(0,2)]*storage[idx(2,0)])*detInv;
        tmp(2,1) =-(storage[idx(0,0)]*storage[idx(2,1)] -
                    storage[idx(0,1)]*storage[idx(2,0)])*detInv;
        tmp(0,2) = (storage[idx(0,1)]*storage[idx(1,2)] -
                    storage[idx(1,1)]*storage[idx(0,2)])*detInv;
        tmp(1,2) =-(storage[idx(0,0)]*storage[idx(1,2)] -
                    storage[idx(0,2)]*storage[idx(1,0)])*detInv;
        tmp(2,2) = (storage[idx(0,0)]*storage[idx(1,1)] -
                    storage[idx(0,1)]*storage[idx(1,0)])*detInv;

        memmove(storage.data(), tmp.data(), 9*sizeof(double));
        return *this;
    }
    dMatrix3x3 inverted(void) const
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0)
        {
            detInv = 1.0/detInv;
        }
        else
        {
            std::cerr << "dMatrix3x3: Can Not Compute Inverse Matrix.\n"
                      << this->print() << "Matrix is Singular !!!\n";
            OP_Exit(EXIT_FAILURE);
        }

        tmp(0,0) = (storage[idx(1,1)]*storage[idx(2,2)] -
                    storage[idx(1,2)]*storage[idx(2,1)])*detInv;
        tmp(1,0) =-(storage[idx(1,0)]*storage[idx(2,2)] -
                    storage[idx(1,2)]*storage[idx(2,0)])*detInv;
        tmp(2,0) = (storage[idx(1,0)]*storage[idx(2,1)] -
                    storage[idx(1,1)]*storage[idx(2,0)])*detInv;
        tmp(0,1) =-(storage[idx(0,1)]*storage[idx(2,2)] -
                    storage[idx(0,2)]*storage[idx(2,1)])*detInv;
        tmp(1,1) = (storage[idx(0,0)]*storage[idx(2,2)] -
                    storage[idx(0,2)]*storage[idx(2,0)])*detInv;
        tmp(2,1) =-(storage[idx(0,0)]*storage[idx(2,1)] -
                    storage[idx(0,1)]*storage[idx(2,0)])*detInv;
        tmp(0,2) = (storage[idx(0,1)]*storage[idx(1,2)] -
                    storage[idx(1,1)]*storage[idx(0,2)])*detInv;
        tmp(1,2) =-(storage[idx(0,0)]*storage[idx(1,2)] -
                    storage[idx(0,2)]*storage[idx(1,0)])*detInv;
        tmp(2,2) = (storage[idx(0,0)]*storage[idx(1,1)] -
                    storage[idx(0,1)]*storage[idx(1,0)])*detInv;

        return tmp;
    }
    dMatrix3x3& transpose(void)
    {
        for(int i =   0; i < 2; i++)
        for(int j = i+1; j < 3; j++)
        {
            double tmp = storage[idx(i,j)];
            storage[idx(i,j)] = storage[idx(j,i)];
            storage[idx(j,i)] = tmp;
        }
        return *this;
    }
    dMatrix3x3 transposed(void) const
    {
        dMatrix3x3 TempMat;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[idx(j,i)];
        }
        return TempMat;
    }
    dMatrix3x3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat;

        memmove(TempMat.data(), storage.data(), 9*sizeof(double));

        TempMat = RotationMatrix * (TempMat * RotationMatrix.transposed());

        memmove(storage.data(), TempMat.data(), 9*sizeof(double));
        return *this;
    }
    dMatrix3x3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 Out;
        memmove(Out.data(), storage.data(), 9*sizeof(double));

        Out = RotationMatrix * (Out * RotationMatrix.transposed());

        return Out;
    }
    dMatrix3x3& rotateU(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat;

        memmove(TempMat.data(), storage.data(), 9*sizeof(double));

        TempMat = RotationMatrix * TempMat;

        memmove(storage.data(), TempMat.data(), 9*sizeof(double));
        return *this;
    }
    dMatrix3x3 rotatedU(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 Out;
        memmove(Out.data(), storage.data(), 9*sizeof(double));

        Out = RotationMatrix * Out;

        return Out;
    }
    dMatrix3x3& power_Elements(double &n)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[idx(i,j)] = pow(storage[idx(i,j)],n);
        }
        return *this;
    }
    dMatrix3x3 powered_Elements(double &n) const
    {
        dMatrix3x3 TempMat;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            TempMat(i,j) = pow(storage[idx(i,j)],n);
        }
        return TempMat;
    }
    double double_contract(const dMatrix3x3& rhs) const
    {
        double tmp = 0.0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp += storage[idx(i,j)]*rhs(i,j);
        }
        return tmp;
    }
    double trace(void) const
    {
        return storage[idx(0,0)] + storage[idx(1,1)] + storage[idx(2,2)];
    }
    double max(void) const
    {
        double value = -DBL_MAX;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            value = std::max(value,storage[idx(i,j)]);
        }
        return value;
    }
    double min(void) const
    {
        double value = DBL_MAX;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            value = std::min(value,storage[idx(i,j)]);
        }
        return value;
    }
    dMatrix3x3& symmetrize(void)
    {
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        if(i != j)
        {
            double tmp = (storage[idx(i,j)] + storage[idx(j,i)]) * 0.5;
            storage[idx(i,j)] = tmp;
            storage[idx(j,i)] = tmp;
        }
        return *this;
    }
    dMatrix3x3 symmetrized(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[idx(i,j)] + storage[idx(j,i)]) * 0.5;
        }

        return TempMat;
    }
    dMatrix3x3 get_skew(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[idx(i,j)] - storage[idx(j,i)]) * 0.5;
        }

        return TempMat;
    }
    double norm(void) const                                                     ///< Frobenius norm
    {
        double tmp = 0.0;
        for(int i = 0; i < 3; i++)
        {
            tmp += storage[idx(i,0)] * storage[idx(i,0)]
                 + storage[idx(i,1)] * storage[idx(i,1)]
                 + storage[idx(i,2)] * storage[idx(i,2)];
        }

        return sqrt(tmp);
    }
    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < 3; i++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(12) << storage[idx(i,0)] << " "
                        << std::setw(12) << storage[idx(i,1)] << " "
                        << std::setw(12) << storage[idx(i,2)] << "||\n";
        }
        return out.str();
    }
    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 3; i++)
        {
           for(int j = 0; j < 3; j++)
           {
               out << storage[idx(i,j)] << sep;
           }
           out << "\n";
        }
        return out.str();
    }
    void read(std::fstream& inp)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            inp >> storage[idx(i,j)];
        }
    }
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            out.push_back((double)storage[idx(i,j)]);
        }
        return out;
    }
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            out.push_back((float)storage[idx(i,j)]);
        }
        return out;
    }
    void read(std::stringstream& inp)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            inp >> storage[idx(i,j)];
        }
    }
    double* data(void)
    {
        return storage.data();
    }
    const double* data(void) const
    {
        return storage.data();
    }
    static dMatrix3x3 UnitTensor(void)
    {
        dMatrix3x3 unity;
        return unity.set_to_unity();
    }
    static dMatrix3x3 ZeroTensor()
    {
        dMatrix3x3 myZeroTensor;
        return myZeroTensor.set_to_zero();
    }
    constexpr size_t sizeX() const
    {
        return 3u;
    }
    constexpr size_t sizeY() const
    {
        return 3u;
    }
    double p_norm(double p)
    {
        //NOTE this is the simpler "entry-wise" norm and not supremum version
        double norm = 0.0;
        for (size_t i = 0; i < sizeX(); i++)
        for (size_t j = 0; j < sizeY(); j++)
        {
            norm += std::pow(std::abs(storage[idx(i,j)]),p);
        }
        return std::pow(norm, 1.0/p);
    }
 protected:
 private:
    std::array<double,9> storage;
    size_t idx(const size_t x, const size_t y) const
    {
        return 3*x + y;
    }
};

}// namespace openphase
#endif //DMATRIX3X3_H
