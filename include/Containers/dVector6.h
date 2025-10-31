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

#ifndef DVECTOR6_H
#define DVECTOR6_H

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

namespace openphase
{

class dVector6
{
 public:
    dVector6()
    {
        storage.fill(0.0);
    }

    dVector6(const dVector6& rhs) :
        storage(rhs.storage)
    {

    }

    dVector6(std::array<double,6> vecinit) :
        storage(vecinit)
    {

    }

    dVector6(std::initializer_list<double> vecinit)
    {
        assert(vecinit.size() == 6 && "Initialization list size is not equal to storage range.");

        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    dVector6& operator=(const dVector6& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    double& operator[](const size_t i)
    {
        assert(i < 6 && "Access beyond storage range");
        return storage[i];
    }

    double const& operator[](const size_t i) const
    {
        assert(i < 6 && "Access beyond storage range");
        return storage[i];
    }

    dVector6& set_to_zero(void)
    {
        storage.fill(0.0);
        return *this;
    }

    void pack(std::vector<double>& buffer)
    {
        for (int i = 0; i < 6; ++i)
        {
            buffer.push_back(storage[i]);
        }
    }

    void unpack(std::vector<double>& buffer, size_t& it)
    {
        for (int i = 0; i < 6; ++i)
        {
            storage[i] = buffer[it]; ++it;
        }
    }

    dVector6& set_to_unity(void)
    {
        storage[0] = 1.0;
        storage[1] = 1.0;
        storage[2] = 1.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
        return *this;
    }

    /*bool operator==(const vStrain& rhs)
    {
        for (int i = 0; i < 6; i++)
        {
            if (storage[i] != rhs[i]){return false;};
        }
        return true;
    };*/

    dVector6 operator*(const double m) const
    {
        dVector6 tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        tmp[3] = storage[3]*m;
        tmp[4] = storage[4]*m;
        tmp[5] = storage[5]*m;
        return tmp;
    }

    dVector6& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        storage[3] *= m;
        storage[4] *= m;
        storage[5] *= m;
        return *this;
    }

    dVector6 operator/(const double m) const
    {
        dVector6 tmp;
        double num = 1.0/m;
        tmp[0] = storage[0]*num;
        tmp[1] = storage[1]*num;
        tmp[2] = storage[2]*num;
        tmp[3] = storage[3]*num;
        tmp[4] = storage[4]*num;
        tmp[5] = storage[5]*num;
        return tmp;
    }

    dVector6& operator/=(const double m)
    {
        double num = 1.0/m;
        storage[0] *= num;
        storage[1] *= num;
        storage[2] *= num;
        storage[3] *= num;
        storage[4] *= num;
        storage[5] *= num;
        return *this;
    }

    dVector6 operator+(const double m) const
    {
        dVector6 tmp;
        tmp[0] = storage[0]+m;
        tmp[1] = storage[1]+m;
        tmp[2] = storage[2]+m;
        tmp[3] = storage[3]+m;
        tmp[4] = storage[4]+m;
        tmp[5] = storage[5]+m;
        return tmp;
    }

    dVector6 operator-(const double m) const
    {
        dVector6 tmp;
        tmp[0] = storage[0]-m;
        tmp[1] = storage[1]-m;
        tmp[2] = storage[2]-m;
        tmp[3] = storage[3]-m;
        tmp[4] = storage[4]-m;
        tmp[5] = storage[5]-m;
        return tmp;
    }

    dVector6 operator-(const dVector6& rhs) const
    {
        dVector6 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    }

    dVector6& operator-=(const dVector6& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        storage[3] -= rhs[3];
        storage[4] -= rhs[4];
        storage[5] -= rhs[5];
        return *this;
    }

    dVector6 operator+(const dVector6& rhs) const
    {
        dVector6 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        tmp[3] = storage[3] + rhs[3];
        tmp[4] = storage[4] + rhs[4];
        tmp[5] = storage[5] + rhs[5];
        return tmp;
    }

    dVector6& operator+=(const dVector6& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        storage[3] += rhs[3];
        storage[4] += rhs[4];
        storage[5] += rhs[5];
        return *this;
    }

    double norm(void) const  /// Frobenius norm
    {
        /*
         * multiplication by 0.5 of the off diagonal elements before taking
         * square and doubling the result due to double appearance of the off
         * diagonal elements => 0.5^2 * 2 = 0.5
         */
        double tmp = storage[0] * storage[0]
                   + storage[1] * storage[1]
                   + storage[2] * storage[2]
                   + storage[3] * storage[3] * 0.5
                   + storage[4] * storage[4] * 0.5
                   + storage[5] * storage[5] * 0.5;
        return sqrt(tmp);
    }

    double max_abs(void) const
    {
        // Returns maximum absolute value
        double tempmax = 0.0;
        for (int i = 0; i < 6; i++)
        {
            if (std::abs(storage[i]) > tempmax)
            {
                tempmax = std::abs(storage[i]);
            }
        }
        return tempmax;
    }

    double trace() const
    {
        return storage[0] + storage[1] + storage[2];
    }

    std::string print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for(int i = 0; i < 6; i++)
        {
            out << storage[i]<< " " << " | ";
        }
        out << " >";
        return out.str();
    }

    std::string write(const int precision = 16, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        out << storage[0] << sep;
        out << storage[1] << sep;
        out << storage[2] << sep;
        out << storage[5] << sep;
        out << storage[3] << sep;
        out << storage[4] << sep;
        return out.str();
    }

    std::vector<float> writeCompressed() const
    {
        std::vector<float> out;
        out.push_back((float)storage[0]);
        out.push_back((float)storage[1]);
        out.push_back((float)storage[2]);
        out.push_back((float)storage[5]);
        out.push_back((float)storage[3]);
        out.push_back((float)storage[4]);
        return out;
    }

    std::vector<double> writeBinary() const
    {
        std::vector<double> out;
        out.push_back((double)storage[0]);
        out.push_back((double)storage[1]);
        out.push_back((double)storage[2]);
        out.push_back((double)storage[5]);
        out.push_back((double)storage[3]);
        out.push_back((double)storage[4]);
        return out;
    };
    double* data(void)
    {
        return storage.data();
    };
    const double* data(void) const
    {
        return storage.data();
    };

    constexpr size_t size(void) const
    {
        return 6u;
    }
 protected:
 private:
    std::array<double, 6> storage;
};

}// namespace openphase
#endif
