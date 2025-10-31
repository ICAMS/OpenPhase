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

#ifndef IVECTOR3_H
#define IVECTOR3_H

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

#include "dVector3.h"

namespace openphase
{

class iVector3
{
 public:

    iVector3()
    {
        storage.fill(0);
    }

    iVector3(const iVector3& vecinit) :
        storage(vecinit.storage)
    {

    }

    iVector3(std::array<long int,3> vecinit) :
        storage(vecinit)
    {

    }

    iVector3(std::initializer_list<long int> vecinit)
    {
        assert(vecinit.size() == 3 && "Initialization list size is not equal to storage range");
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    long int& operator[](const size_t i)
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    }

    const long int & operator[](const size_t i) const
    {
        assert(i < 3 && "Access beyond storage range");
        return storage[i];
    }

    bool operator<(iVector3 rhs) const
    {
        return length() < rhs.length();
    }

    bool operator>(iVector3 rhs) const
    {
        return length() > rhs.length();
    }

    bool operator==(iVector3 rhs) const
    {
        return length_sqr() == rhs.length_sqr();
    }

    long int get_x(void) const
    {
        return storage[0];
    }

    void set_x(const long int new_x)
    {
        storage[0] = new_x;
    }

    long int get_y(void) const
    {
        return storage[1];
    }

    void set_y(const long int new_y)
    {
        storage[1] = new_y;
    }

    long int get_z(void) const
    {
        return storage[2];
    }

    void set_z(const long int new_z)
    {
        storage[2] = new_z;
    }

    void pack(std::vector<long int>& buffer)
    {
        for (int i = 0; i < 3; ++i)
        {
            buffer.push_back(storage[i]);
        }
    }

    void unpack(std::vector<long int>& buffer, size_t& it)
    {
        for (int i = 0; i < 3; ++i)
        {
            storage[i] = buffer[it]; ++it;
        }
    }

    void set_to_zero(void)
    {
        storage.fill(0);
    }

    iVector3 operator*(const long int m) const
    {
        iVector3 tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        return tmp;
    }

    iVector3 operator/(const long int m) const
    {
        iVector3 tmp;
        tmp[0] = storage[0]/m;
        tmp[1] = storage[1]/m;
        tmp[2] = storage[2]/m;
        return tmp;
    }

    long int operator*(const iVector3& rhs) const
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
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    }

    long int length_sqr(void) const
    {
        return (storage[0]*storage[0] +
                storage[1]*storage[1] +
                storage[2]*storage[2]);
    }

    iVector3 cross(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
        tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
        tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
        return tmp;
    }

    iVector3 operator+(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        return tmp;
    }

    iVector3 operator-(const iVector3& rhs) const
    {
        iVector3 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        return tmp;
    }

    iVector3& operator*=(const long int m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        return *this;
    }

    iVector3& operator/=(const long int m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        return *this;
    }

    iVector3& operator-=(const iVector3& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        return *this;
    }

    iVector3& operator+=(const iVector3& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        return *this;
    }

    iVector3& operator=(const iVector3& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    iVector3& operator=(const long int rhs[3])
    {
        memmove(storage.data(), rhs, 3*sizeof(long int));
        return *this;
    }

    dVector3 normalized(void) const
    {
        double norm = length();
        dVector3 tmp;
        if(norm != 0.0)
        {
            double norm_inv = 1.0/norm;
            tmp[0] = storage[0]*norm_inv;
            tmp[1] = storage[1]*norm_inv;
            tmp[2] = storage[2]*norm_inv;
        }
        else
        {
            tmp.set_to_zero();
        }
        return tmp;
    }

    iVector3 Xreflected(void) const
    {
        iVector3 Out(*this);
        Out[0] *= -1.0;
        return Out;
    }

    iVector3 Yreflected(void) const
    {
        iVector3 Out(*this);
        Out[1] *= -1.0;
        return Out;
    }

    iVector3 Zreflected(void) const
    {
        iVector3 Out(*this);
        Out[2] *= -1.0;
        return Out;
    }

    std::string print(void) const
    {
        std::stringstream out;

        out << "(" << storage[0] << ", "
                   << storage[1] << ", "
                   << storage[2] << ")";
        return out.str();
    }

    void write(std::fstream& out) const
    {
        for(int i = 0; i < 3; i++)
        {
            out << storage[i] << " ";
        }
        out << std::endl;
    }

    void write(std::stringstream& out) const
    {
        for(int i = 0; i < 3; i++)
        {
            out << storage[i] << " ";
        }
        out << std::endl;
    }

    void read(std::fstream& inp)
    {
        for(int i = 0; i < 3; i++)
        {
            inp >> storage[i];
        }
    }

    void read(std::stringstream& inp)
    {
        for(int i = 0; i < 3; i++)
        {
            inp >> storage[i];
        }
    }

    static iVector3 ZeroVector(void)
    {
        return iVector3({0,0,0});
    }

    constexpr size_t size(void) const
    {
        return 3;
    }

    long int* data(void)
    {
        return storage.data();
    }

    const long int* data(void) const
    {
        return storage.data();
    }
 protected:
 private:
    std::array<long int, 3> storage;
};

}// namespace openphase
#endif
