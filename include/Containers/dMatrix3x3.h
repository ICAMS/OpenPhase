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

#ifndef DMATRIX3X3_H
#define DMATRIX3X3_H

#include "Globals.h"

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
#include <array>

namespace openphase
{

class OP_EXPORTS dMatrix3x3
{
 public:
    dMatrix3x3() : storage()
    {
    }

    dMatrix3x3(const dMatrix3x3& rhs) : storage(rhs.storage)
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

    static dMatrix3x3 ZeroTensor()
    {
        dMatrix3x3 myZeroTensor;
        return myZeroTensor.set_to_zero();
    }

    static dMatrix3x3 UnitTensor(void)
    {
        dMatrix3x3 myUnitTensor;
        return myUnitTensor.set_to_unity();
    }

    double& operator()(const size_t i, const size_t j)
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[Index(i,j)];
    }

    double const& operator()(const size_t i, const size_t j) const
    {
        assert(i < 3 && "Access beyond storage range");
        assert(j < 3 && "Access beyond storage range");

        return storage[Index(i,j)];
    }

    double& operator[](const size_t idx)
    {
        assert(idx < 9 && "Access beyond storage range");

        return storage[idx];
    }

    const double& operator[](const size_t idx) const
    {
        assert(idx < 9 && "Access beyond storage range");

        return storage[idx];
    }

    dMatrix3x3& set_to_zero(void)
    {
        storage.fill(0.0);
        return *this;
    }

    dMatrix3x3& set_to_unity(void)
    {
        storage.fill(0.0);
        for(int i = 0; i < 3; i++)
        {
            storage[Index(i,i)] = 1.0;
        }
        return *this;
    }

    void set(double r00, double r01, double r02,
             double r10, double r11, double r12,
             double r20, double r21, double r22)
    {
        storage[Index(0,0)] = r00;
        storage[Index(0,1)] = r01;
        storage[Index(0,2)] = r02;
        storage[Index(1,0)] = r10;
        storage[Index(1,1)] = r11;
        storage[Index(1,2)] = r12;
        storage[Index(2,0)] = r20;
        storage[Index(2,1)] = r21;
        storage[Index(2,2)] = r22;
    }

    dMatrix3x3& operator=(const dMatrix3x3& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    bool operator==(const dMatrix3x3& rhs)
    {
        bool result = true;
        for(int n = 0; n < 9; n++)
        if(fabs(storage[n] - rhs.storage[n]) > DBL_EPSILON)
        {
            result = false;
            break;
        }
        return result;
    }

    bool operator!=(const dMatrix3x3& rhs)
    {
        return !(*this == rhs);
    }

    dMatrix3x3 operator*(const double m) const
    {
        dMatrix3x3 tmp(*this);
        for(int n = 0; n < 9; n++)
        {
            tmp.storage[n] *= m;
        }
        return tmp;
    }

    dMatrix3x3 operator/(const double m) const
    {
        dMatrix3x3 tmp(*this);
        for(int n = 0; n < 9; n++)
        {
            tmp.storage[n] /= m;
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
            tmp(i,j) += storage[Index(i,k)]*rhs(k,j);
        }
        return tmp;
    }

    dMatrix3x3 operator+(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp(*this);
        for(int n = 0; n < 9; n++)
        {
            tmp.storage[n] += rhs.storage[n];
        }
        return tmp;
    }

    dMatrix3x3 operator-(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp(*this);
        for(int n = 0; n < 9; n++)
        {
            tmp.storage[n]-= rhs.storage[n];
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
            tmp(i,j) += storage[Index(i,k)]*rhs(k,j);
        }
        storage = tmp.storage;
        return *this;
    }

    dMatrix3x3& operator+=(const dMatrix3x3& rhs)
    {
        for(int n = 0; n < 9; n++)
        {
            storage[n] += rhs.storage[n];
        }
        return *this;
    }

    dMatrix3x3& operator-=(const dMatrix3x3& rhs)
    {
        for(int n = 0; n < 9; n++)
        {
            storage[n] -= rhs.storage[n];
        }
        return *this;
    }

    dMatrix3x3& operator*=(const double m)
    {
        for(int n = 0; n < 9; n++)
        {
            storage[n] *= m;
        }
        return *this;
    }

    dMatrix3x3& operator/=(const double m)
    {
        for(int n = 0; n < 9; n++)
        {
            storage[n] /= m;
        }
        return *this;
    }

    dMatrix3x3 Hadamard_product(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp(*this);

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) *= rhs(i,j);
        }
        return tmp;
    }

    double determinant(void) const
    {
        return (storage[Index(0,0)]*storage[Index(1,1)]*storage[Index(2,2)] +
                storage[Index(0,1)]*storage[Index(1,2)]*storage[Index(2,0)] +
                storage[Index(0,2)]*storage[Index(1,0)]*storage[Index(2,1)] -
                storage[Index(0,2)]*storage[Index(1,1)]*storage[Index(2,0)] -
                storage[Index(0,1)]*storage[Index(1,0)]*storage[Index(2,2)] -
                storage[Index(0,0)]*storage[Index(1,2)]*storage[Index(2,1)]);
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

        tmp(0,0) = (storage[Index(1,1)]*storage[Index(2,2)] -
                    storage[Index(1,2)]*storage[Index(2,1)])*detInv;
        tmp(1,0) =-(storage[Index(1,0)]*storage[Index(2,2)] -
                    storage[Index(1,2)]*storage[Index(2,0)])*detInv;
        tmp(2,0) = (storage[Index(1,0)]*storage[Index(2,1)] -
                    storage[Index(1,1)]*storage[Index(2,0)])*detInv;
        tmp(0,1) =-(storage[Index(0,1)]*storage[Index(2,2)] -
                    storage[Index(0,2)]*storage[Index(2,1)])*detInv;
        tmp(1,1) = (storage[Index(0,0)]*storage[Index(2,2)] -
                    storage[Index(0,2)]*storage[Index(2,0)])*detInv;
        tmp(2,1) =-(storage[Index(0,0)]*storage[Index(2,1)] -
                    storage[Index(0,1)]*storage[Index(2,0)])*detInv;
        tmp(0,2) = (storage[Index(0,1)]*storage[Index(1,2)] -
                    storage[Index(1,1)]*storage[Index(0,2)])*detInv;
        tmp(1,2) =-(storage[Index(0,0)]*storage[Index(1,2)] -
                    storage[Index(0,2)]*storage[Index(1,0)])*detInv;
        tmp(2,2) = (storage[Index(0,0)]*storage[Index(1,1)] -
                    storage[Index(0,1)]*storage[Index(1,0)])*detInv;

        return tmp;
    }

    dMatrix3x3& invert(void)
    {
        dMatrix3x3 TempMat = inverted();
        storage = TempMat.storage;
        return *this;
    }

    dMatrix3x3 transposed(void) const
    {
        dMatrix3x3 TempMat;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[Index(j,i)];
        }
        return TempMat;
    }

    dMatrix3x3& transpose(void)
    {
        dMatrix3x3 TempMat = transposed();
        storage = TempMat.storage;
        return *this;
    }

    dMatrix3x3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 TempMat(*this);
        TempMat = RotationMatrix * (TempMat * RotationMatrix.transposed());
        return TempMat;
    }

    dMatrix3x3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat = rotated(RotationMatrix);
        storage = TempMat.storage;
        return *this;
    }

    dMatrix3x3 rotatedU(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 TempMat(*this);
        TempMat = RotationMatrix * TempMat;
        return TempMat;
    }

    dMatrix3x3& rotateU(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat = rotatedU(RotationMatrix);
        storage = TempMat.storage;
        return *this;
    }

    dMatrix3x3& power_elements(double &p)
    {
        for(int n = 0; n < 9; n++)
        {
            storage[n] = pow(storage[n],p);
        }
        return *this;
    }

    dMatrix3x3 powered_elements(double &p) const
    {
        dMatrix3x3 TempMat;
        for(int n = 0; n < 9; n++)
        {
            TempMat.storage[n] = pow(storage[n],p);
        }
        return TempMat;
    }

    double double_contract(const dMatrix3x3& rhs) const
    {
        double tmp = 0.0;
        for (int n = 0; n < 9; n++)
        {
            tmp += storage[n]*rhs.storage[n];
        }
        return tmp;
    }

    double trace(void) const
    {
        return storage[Index(0,0)] + storage[Index(1,1)] + storage[Index(2,2)];
    }

    double max(void) const
    {
        double value = -DBL_MAX;
        for(int n = 0; n < 9; n++)
        {
            value = std::max(value,storage[n]);
        }
        return value;
    }

    double min(void) const
    {
        double value = DBL_MAX;
        for(int n = 0; n < 9; n++)
        {
            value = std::min(value,storage[n]);
        }
        return value;
    }

    dMatrix3x3 symmetrized(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[Index(i,j)] + storage[Index(j,i)]) * 0.5;
        }

        return TempMat;
    }

    dMatrix3x3& symmetrize(void)
    {
        dMatrix3x3 TempMat = symmetrized();
        storage = TempMat.storage;
        return *this;
    }

    dMatrix3x3 get_skew(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = (storage[Index(i,j)] - storage[Index(j,i)]) * 0.5;
        }

        return TempMat;
    }

    double p_norm(double p) const
    {
        //NOTE this is the simpler "entry-wise" norm and not supremum version
        double norm = 0.0;
        for (size_t n = 0; n < 9; n++)
        {
            norm += std::pow(std::fabs(storage[n]),p);
        }
        return std::pow(norm, 1.0/p);
    }

    double norm(void) const                                                     ///< Frobenius norm
    {
        return p_norm(2.0);
    }

    constexpr size_t sizeX() const
    {
        return 3u;
    }

    constexpr size_t sizeY() const
    {
        return 3u;
    }

    constexpr size_t size() const                                               ///< Returns the size of the internal flattened storage.
    {
        return 9u;
    }

    double* data(void)
    {
        return storage.data();
    }

    const double* data(void) const
    {
        return storage.data();
    }

    void pack(std::vector<double>& buffer) const
    {
        for (int n = 0; n < 9; ++n)
        {
            buffer.push_back(storage[n]);
        }
    }

    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int n = 0; n < 9; ++n)
        {
            storage[n] = buffer[it]; ++it;
        }
    }

    std::string print(void) const
    {
        std::stringstream out;
        out << std::setprecision(6);
        for(int i = 0; i < 3; i++)
        {
            out << "||";
            for(int j = 0; j < 3; j++)
            {
                out << std::setw(10) << storage[Index(i,j)] << " ";
            }
            out << "||\n";
        }
        return out.str();
    }

    void read_binary(std::istream& inp)
    {
        for(size_t n = 0; n < 9; n++)
        {
            inp >> storage[n];
        }
    }

    void read_ASCII(std::istream& inp)
    {
        for(size_t n = 0; n < 9; n++)
        {
            inp >> storage[n];
        }
    }

    void write_binary(std::ostream& outp) const
    {
        for(size_t n = 0; n < 9; n++)
        {
            outp << storage[n];
        }
    }

    void write_ASCII(std::ostream& outp, const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        outp << std::setprecision(precision) << std::defaultfloat;
        for(size_t n = 0; n < 9; n++)
        {
            outp << storage[n] << sep;
        }
        outp << std::endl;
    }

    std::string get_output_string(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int n = 0; n < 9; n++)
        {
            out << storage[n] << sep;
        }
        return out.str();
    }

    std::vector<double> get_vector() const
    {
        std::vector<double> out(9);
        for(int n = 0; n < 9; n++)
        {
            out[n] = storage[n];
        }
        return out;
    }

    std::vector<float> get_vector_float() const
    {
        std::vector<float> out(9);
        for(int n = 0; n < 9; n++)
        {
            out[n] = (float)storage[n];
        }
        return out;
    }

    //============================= Deprecated methods =========================
    [[deprecated]]
    std::string write(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int n = 0; n < 9; n++)
        {
            out << storage[n] << sep;
        }
        return out.str();
    }

    [[deprecated]]
    void read(std::fstream& inp)
    {
        for(int n = 0; n < 9; n++)
        {
            inp >> storage[n];
        }
    }

    [[deprecated]]
    std::vector<double> writeBinary() const
    {
        std::vector<double> out(9);
        for(int n = 0; n < 9; n++)
        {
            out[n] = storage[n];
        }
        return out;
    }

    [[deprecated]]
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out(9);
        for(int n = 0; n < 9; n++)
        {
            out[n] = (float)storage[n];
        }
        return out;
    }

    [[deprecated]]
    void read(std::stringstream& inp)
    {
        for(int n = 0; n < 9; n++)
        {
            inp >> storage[n];
        }
    }
    //==========================================================================

 protected:
 private:
    std::array<double,9> storage;
    size_t Index(const size_t x, const size_t y) const
    {
        return 3*x + y;
    }
};

}// namespace openphase
#endif //DMATRIX3X3_H
