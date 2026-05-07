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

#ifndef DVECTOR3_H
#define DVECTOR3_H

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

#include "dMatrix3x3.h"

namespace openphase
{

class dVector3
{
 public:
    dVector3()
    {
        storage.fill(0.0);
    }

    dVector3(const dVector3& vecinit) :
        storage(vecinit.storage)
    {

    }

    dVector3(std::array<double,3> vecinit) :
        storage(vecinit)
    {

    }

    dVector3(std::array<int,3> vecinit)
    {
        storage[0] = vecinit[0];
        storage[1] = vecinit[1];
        storage[2] = vecinit[2];
    }

    dVector3(std::initializer_list<double> vecinit)
    {
        assert(vecinit.size() == 3 && "Initialization list size is not equal to storage range");
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    static dVector3 ZeroVector(void)
    {
        return dVector3();
    }

    double& operator[](const size_t i)
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    }

    double const& operator[](const size_t i) const
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    }

//    bool operator<(dVector3 rhs)
//    {
//        return length() < rhs.length();
//    }
//
//    bool operator>(dVector3 rhs)
//    {
//        return length() > rhs.length();
//    }

    double getX(void) const
    {
        return storage[0];
    }

    double getY(void) const
    {
        return storage[1];
    }

    double getZ(void) const
    {
        return storage[2];
    }

    void setX(const double newX)
    {
        storage[0] = newX;
    }

    void setY(const double newY)
    {
        storage[1] = newY;
    }

    void setZ(const double newX)
    {
        storage[2] = newX;
    }

    void set_to_zero(void)
    {
        storage.fill(0.0);
    }

//    void set_to_unit_X(void)/
//    {
//        storage[0] = 1.0;
//        storage[1] = 0.0;
//        storage[2] = 0.0;
//    }
//
//    void set_to_unit_Y(void)
//    {
//        storage[0] = 0.0;
//        storage[1] = 1.0;
//        storage[2] = 0.0;
//    }
//
//    void set_to_unit_Z(void)
//    {
//        storage[0] = 0.0;
//        storage[1] = 0.0;
//        storage[2] = 1.0;
//    }

    dVector3 operator*(const double m) const
    {
        dVector3 tmp(*this);
        tmp[0] *= m;
        tmp[1] *= m;
        tmp[2] *= m;
        return tmp;
    }

    dVector3 operator/(const double m) const
    {
        assert(std::fabs(m) > std::numeric_limits<double>::epsilon() && "Division by zero in dVector3::operator/.");
        dVector3 tmp(*this);
        tmp[0] /= m;
        tmp[1] /= m;
        tmp[2] /= m;
        return tmp;
    }

    double operator*(const dVector3& rhs) const
    {
        return storage[0]*rhs[0] + storage[1]*rhs[1] + storage[2]*rhs[2];
    }

    double abs() const
    {
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    }

    double length(void) const
    {
        return abs();
    }

    double p_norm(double p)
    {
        assert(p >= 1.0 && "p should be >= 1.0 in dVector3::p_norm to be a valid norm.");
        double norm = 0.0;
        for (size_t i = 0; i < size(); i++)
        {
            norm += std::pow(std::fabs(storage[i]),p);
        }
        return std::pow(norm, 1.0/p);
    }

    dVector3 cross(const dVector3& rhs) const
    {
        dVector3 tmp;
        tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
        tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
        tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
        return tmp;
    }

    dMatrix3x3 dyadic(const dVector3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int n = 0; n < 3; n++)
        for(int m = 0; m < 3; m++)
        {
            tmp(n,m) = storage[n]*rhs[m];
        }
        return tmp;
    }

    dVector3 operator+(const dVector3& rhs) const
    {
        dVector3 tmp(*this);
        tmp[0] += rhs[0];
        tmp[1] += rhs[1];
        tmp[2] += rhs[2];
        return tmp;
    }

    dVector3 operator-(const dVector3& rhs) const
    {
        dVector3 tmp(*this);
        tmp[0] -= rhs[0];
        tmp[1] -= rhs[1];
        tmp[2] -= rhs[2];
        return tmp;
    }

    dVector3& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        return *this;
    }

    dVector3& operator/=(const double m)
    {
        assert(std::fabs(m) > std::numeric_limits<double>::epsilon() && "Division by zero in dVector3::operator/=.");
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        return *this;
    }

    dVector3& operator-=(const dVector3& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        return *this;
    }

    dVector3& operator+=(const dVector3& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        return *this;
    }

    dVector3& operator=(const dVector3& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    dVector3& operator=(const double rhs[3])
    {
        storage[0] = rhs[0];
        storage[1] = rhs[1];
        storage[2] = rhs[2];
        return *this;
    }

    dVector3& normalize(void)
    {
        double norm = length();
        if(norm > DBL_EPSILON)
        {
            double norm_inv = 1.0/norm;
            storage[0] *= norm_inv;
            storage[1] *= norm_inv;
            storage[2] *= norm_inv;
        }
        else
        {
            set_to_zero();
        }
        return *this;
    }

    dVector3 normalized(void) const
    {
        dVector3 tmp(*this);
        tmp.normalize();
        return tmp;
    }

    dVector3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dVector3 tmp;
        for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            tmp[i] += RotationMatrix(i,j) * storage[j];
        }
        storage = tmp.storage;
        return *this;
    }

    dVector3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dVector3 Out(*this);
        Out.rotate(RotationMatrix);
        return Out;
    }

    dVector3 Xreflected(void) const
    {
        dVector3 Out(*this);
        Out[0] *= -1.0;
        return Out;
    }

    dVector3 Yreflected(void) const
    {
        dVector3 Out(*this);
        Out[1] *= -1.0;
        return Out;
    }

    dVector3 Zreflected(void) const
    {
        dVector3 Out(*this);
        Out[2] *= -1.0;
        return Out;
    }

    constexpr size_t size(void) const
    {
        return 3;
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
        for (int i = 0; i < 3; ++i)
        {
            buffer.push_back(storage[i]);
        }
    }

    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 3; ++i)
        {
            storage[i] = buffer[it]; ++it;
        }
    }

    std::string print(void) const
    {
        std::stringstream out;

        out << "(" << storage[0] << ", "
                   << storage[1] << ", "
                   << storage[2] << ")";
        return out.str();
    }

    void read_binary(std::istream& inp)
    {
        for(size_t i = 0; i < 3; i++)
        {
            inp >> storage[i];
        }
    }

    void read_ASCII(std::istream& inp)
    {
        for(size_t i = 0; i < 3; i++)
        {
            inp >> storage[i];
        }
    }

    void write_binary(std::ostream& outp) const
    {
        for(size_t i = 0; i < 3; i++)
        {
            outp << storage[i];
        }
    }

    void write_ASCII(std::ostream& outp, const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        outp << std::setprecision(precision) << std::defaultfloat;
        for(size_t i = 0; i < 3; i++)
        {
            outp << storage[i] << sep;
        }
        outp << std::endl;
    }

    std::string get_output_string(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 3; i++)
        {
            out << storage[i] << sep;
        }
        return out.str();
    }

    std::vector<double> get_vector() const
    {
        std::vector<double> out(3);
        for(int i = 0; i < 3; i++)
        {
            out[i] = storage[i];
        }
        return out;
    }

    std::vector<float> get_vector_float() const
    {
        std::vector<float> out(3);
        for(int i = 0; i < 3; i++)
        {
            out[i] = (float)storage[i];
        }
        return out;
    }

    //============================= Deprecated methods==========================
    [[deprecated]]
    void read(std::istream& inp)
    {
        for(size_t i = 0; i < 3; i++)
        {
            inp >> storage[i];
        }
    }

    [[deprecated]]
    std::string write(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 3; i++)
        {
            out << storage[i] << sep;
        }
        return out.str();
    }

    [[deprecated]]
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        for(int i = 0; i < 3; i++)
        {
            out.push_back((float)storage[i]);
        }
        return out;
    }

    [[deprecated]]
    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        for(int i = 0; i < 3; i++)
        {
            out.push_back((double)storage[i]);
        }
        return out;
    }
    //==========================================================================

 protected:
 private:
    std::array<double, 3> storage;
};

}// namespace openphase
#endif
