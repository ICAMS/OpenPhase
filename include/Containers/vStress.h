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

#ifndef VSTRESS_H
#define VSTRESS_H

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
#include "dVector3.h"

namespace openphase
{

class vStress                                                                   ///< Stress Voigt vector
{
    /*
     *  This is a special version of the six component Voigt vector
     *  to store stresses. It implements a number of stress specific
     *  methods.
     */
 public:
    vStress()
    {
        storage.fill(0.0);
    }

    vStress(const vStress& rhs) :
        storage(rhs.storage)
    {

    }

    double& operator[](const int i)
    {
        assert(i < 6 && "Access beyond storage range");
        return storage[i];
    }

    double const& operator[](const int i) const
    {
        assert(i < 6 && "Access beyond storage range");
        return storage[i];
    }

    vStress& set_to_zero()
    {
        storage.fill(0.0);
        return *this;
    }

    vStress& operator=(const vStress& rhs)
    {
        storage = rhs.storage;
        return *this;
    }

    vStress operator*(const double m) const
    {
        vStress tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        tmp[3] = storage[3]*m;
        tmp[4] = storage[4]*m;
        tmp[5] = storage[5]*m;
        return tmp;
    }

    vStress operator/(const double m) const
    {
        vStress tmp;
        double num = 1.0/m;
        tmp[0] = storage[0]*num;
        tmp[1] = storage[1]*num;
        tmp[2] = storage[2]*num;
        tmp[3] = storage[3]*num;
        tmp[4] = storage[4]*num;
        tmp[5] = storage[5]*num;
        return tmp;
    }

    vStress& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        storage[3] *= m;
        storage[4] *= m;
        storage[5] *= m;
        return *this;
    }

    vStress& operator/=(const double m)
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

    vStress operator+(const vStress& rhs) const
    {
        vStress tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        tmp[3] = storage[3] + rhs[3];
        tmp[4] = storage[4] + rhs[4];
        tmp[5] = storage[5] + rhs[5];
        return tmp;
    }

    vStress& operator+=(const vStress& rhs)
    {
        storage[0] += rhs[0];
        storage[1] += rhs[1];
        storage[2] += rhs[2];
        storage[3] += rhs[3];
        storage[4] += rhs[4];
        storage[5] += rhs[5];
        return *this;
    }

    vStress operator-(const vStress& rhs) const
    {
        vStress tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    }

    vStress& operator-=(const vStress& rhs)
    {
        storage[0] -= rhs[0];
        storage[1] -= rhs[1];
        storage[2] -= rhs[2];
        storage[3] -= rhs[3];
        storage[4] -= rhs[4];
        storage[5] -= rhs[5];
        return *this;
    }

    double norm(void) const                                                     ///< Frobenius norm
    {
        double tmp = storage[0] * storage[0]
                   + storage[1] * storage[1]
                   + storage[2] * storage[2]
                   + storage[3] * storage[3] * 2.0
                   + storage[4] * storage[4] * 2.0
                   + storage[5] * storage[5] * 2.0;
        return sqrt(tmp);
    }

    double double_contract(const vStress& Bstress) const  // "double-dot product"
    {
        double tmp =
        storage[0] * Bstress[0] +
        storage[1] * Bstress[1] +
        storage[2] * Bstress[2] +
        2.0*storage[3] * Bstress[3] +
        2.0*storage[4] * Bstress[4] +
        2.0*storage[5] * Bstress[5];
        return tmp;
    }

    double Pressure() const
    {
        return -(storage[0] + storage[1] + storage[2])/3.0;
    }

    double trace() const
    {
        return storage[0] + storage[1] + storage[2];
    }

    double determinant() const
    {
        return storage[0]*storage[1]*storage[2] -
               storage[0]*storage[3]*storage[3] -
               storage[1]*storage[4]*storage[4] -
               storage[2]*storage[5]*storage[5] +
               storage[3]*storage[4]*storage[5]*2.0;
    }

    dVector3 invariants(void) const
    {
        const double I1 = trace();
        const double I2 = 0.5*(trace()*trace() -
                (storage[0]*storage[0] +
                 storage[1]*storage[1] +
                 storage[2]*storage[2] +
                 storage[3]*storage[3]*2 +
                 storage[4]*storage[4]*2 +
                 storage[5]*storage[5]*2));
        const double I3 = determinant();
        return dVector3({I1,I2,I3});
    }

    double Mises(void) const
    {
        double vMises = 0.0;
        vMises = (storage[0]-storage[1])*(storage[0]-storage[1]) +
                 (storage[1]-storage[2])*(storage[1]-storage[2]) +
                 (storage[2]-storage[0])*(storage[2]-storage[0]) +
                 6.0*(storage[3]*storage[3]+storage[4]*storage[4]+storage[5]*storage[5]);

        return sqrt(0.5*vMises);
    }

    vStress rotated(const dMatrix3x3& RotationMatrix) const
    {
        const dMatrix3x3& locStressTensor = tensor().rotated(RotationMatrix);

        vStress locStress;

        locStress[0] = locStressTensor(0,0);
        locStress[5] = locStressTensor(0,1);
        locStress[4] = locStressTensor(0,2);
        locStress[1] = locStressTensor(1,1);
        locStress[3] = locStressTensor(1,2);
        locStress[2] = locStressTensor(2,2);

        return locStress;
    }

    vStress& rotate(const dMatrix3x3& RotationMatrix)
    {
        const dMatrix3x3& locStressTensor = tensor().rotated(RotationMatrix);

        storage[0] = locStressTensor(0,0);
        storage[5] = locStressTensor(0,1);
        storage[4] = locStressTensor(0,2);
        storage[1] = locStressTensor(1,1);
        storage[3] = locStressTensor(1,2);
        storage[2] = locStressTensor(2,2);

        return *this;
    }

    double get_tensor(const int i, const int j) const
    {
        return storage[(i==j)?(i):(6-(i+j))];
    }

    dMatrix3x3 tensor(void) const
    {
        dMatrix3x3 tmp;
        tmp(0,0) = storage[0];
        tmp(0,1) = storage[5];
        tmp(0,2) = storage[4];
        tmp(1,0) = storage[5];
        tmp(1,1) = storage[1];
        tmp(1,2) = storage[3];
        tmp(2,0) = storage[4];
        tmp(2,1) = storage[3];
        tmp(2,2) = storage[2];
        return tmp;
    }

    constexpr size_t size(void) const
    {
        return 6u;
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

    std::string print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for(int i = 0; i < 6; i++)
        {
            out << storage[i] << " " << " | ";
        }
        out << " >";
        return out.str();
    }

    void read_binary(std::istream& inp)
    {
        for(size_t i = 0; i < 6; i++)
        {
            inp >> storage[i];
        }
    }

    void read_ASCII(std::istream& inp)
    {
        for(size_t i = 0; i < 6; i++)
        {
            inp >> storage[i];
        }
    }

    void write_binary(std::ostream& outp) const
    {
        for(size_t i = 0; i < 6; i++)
        {
            outp << storage[i];
        }
    }

    void write_ASCII(std::ostream& outp, const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        outp << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 6; i++)
        {
            outp << storage[i] << sep;
        }
        outp << std::endl;
    }

    std::string get_output_string(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 6; i++)
        {
            out << storage[i] << sep;
        }
        return out.str();
    }

    std::vector<double> get_vector() const
    {
        std::vector<double> out(6);
        for(int i = 0; i < 6; i++)
        {
            out[i] = storage[i];
        }
        return out;
    }

    std::vector<float> get_vector_float() const
    {
        std::vector<float> out(6);
        for(int i = 0; i < 6; i++)
        {
            out[i] = (float)storage[i];
        }
        return out;
    }

    //============================= Deprecated methods==========================
    [[deprecated]]
    std::string print_line(void) const
    {
        std::stringstream out;
        out << " ";
        for(int i = 0; i < 6; i++)
        {
            out << storage[i] << ", ";
        }
        out << " ";

        return out.str();
    }

    [[deprecated]]
    std::string write(const int precision = std::numeric_limits<double>::digits10 + 1, const char sep = ' ') const
    {
        std::stringstream out;
        out << std::setprecision(precision) << std::defaultfloat;
        for(int i = 0; i < 6; i++)
        {
            out << storage[i] << sep;
        }
        return out.str();
    }

    [[deprecated]]
    std::vector<float> writeCompressed() const
    {
        std::vector<float> out(6);
        for(int i = 0; i < 6; i++)
        {
            out[i] = (float)storage[i];
        }
        return out;
    }

    [[deprecated]]
    std::vector<double> writeBinary() const
    {
        std::vector<double> out(6);
        for(int i = 0; i < 6; i++)
        {
            out[i] = storage[i];
        }
        return out;
    }
    //==========================================================================

protected:
private:
    std::array<double,6> storage;
};

}// namespace openphase
#endif
