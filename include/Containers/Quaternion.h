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

#ifndef QUATERNION_H
#define QUATERNION_H

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

#include "dVector3.h"
#include "dMatrix3x3.h"

namespace openphase
{

class OP_EXPORTS Quaternion                                                     ///< Stores and manages Quaternions
{
 public:

    Quaternion();                                                               ///< Default constructor
    Quaternion(const Quaternion& rhS);                                          ///< Copy constructor
    Quaternion(std::array<double,4> vecinit);                                   ///< Constructs quaternion using array of four values

    void set(const double sIn, const double x, const double y, const double z); ///< Sets quaternion using individual entries
    void set(const double sIn, const dVector3 vIn);                             ///< Sets quaternion using first entry and imaginary part vector
    void set(const dMatrix3x3& RotMatrix);                                      ///< Sets quaternion from rotation matrix
    void set_entry(const int idx, const double val);                            ///< Sets individual entry
    void set(dVector3 Axis, const double Angle);                                ///< Sets quaternion from rotation axis and angle
    void set_to_zero_rotation();                                                ///< Sets zero rotation quaternion
    void set_to_zero();                                                         ///< Sets all quaternion components to zero

    double& operator[](const size_t index);                                     ///< Bi-directional access operator for accessing individual entries
    const double& operator[](const size_t index) const;                         ///< Bi-directional access operator for accessing individual entries

    Quaternion  operator+(const Quaternion rhS) const;
    Quaternion  operator-(const Quaternion rhS) const;
    Quaternion  operator*(const Quaternion rhS) const;
    Quaternion  operator*(const double scalar) const;
    Quaternion  operator/(const double divisor) const;

    Quaternion& operator=(const Quaternion rhS);
    Quaternion& operator+=(const Quaternion rhS);
    Quaternion& operator-=(const Quaternion rhS);
    Quaternion& operator*=(const Quaternion rhS);
    Quaternion& operator*=(const double scalar);
    Quaternion& operator/=(const double divisor);

    double length() const;                                                      ///< Returns length of the quaternion
    Quaternion& normalize();                                                    ///< Normalizes quaternion
    Quaternion  normalized() const;                                             ///< Returns normalized quaternion
    Quaternion& conjugate();                                                    ///< Conjugates quaternion (reverse rotation)
    Quaternion  conjugated() const;                                             ///< Returns conjugated quaternion (reverse rotation)
    Quaternion& invert();                                                       ///< Conjugates quaternion (reverse rotation if quaternion is normalized)
    Quaternion  inverted() const;                                               ///< Returns inverted quaternion (reverse rotation if quaternion is normalized)

    void set_rotation_matrix(void);                                             ///< Sets rotation matrix for the quaternion
    dMatrix3x3 get_rotation_matrix(const bool Active = true);                   ///< Returns rotation matrix (default: active rotation)

    static Quaternion lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2,
                                                            const double t);    ///< Linear interpolation between two quaternions, where 0.0 < t < 1.0
    static Quaternion slerp(const Quaternion& a, const Quaternion& b,
                                                            const double t);    ///< Linear interpolation between two quaternions, where 0.0 < t < 1.0

    constexpr size_t size() const {return 4;};                                  ///< Returns the size of the container

    void pack(std::vector<double>& buffer);                                     ///< Packs the quaternion into communication buffer
    void unpack(std::vector<double>& buffer, size_t& it);                       ///< Unpacks the quaternion from communication buffer

    std::string print(void) const;                                              ///< Returns the quaternion as formated string

    void read_binary(std::istream& inp);                                        ///< Reads quaternion from a file stream in binary format
    void read_ASCII(std::istream& inp);                                         ///< Reads quaternion from a file stream in ASCII format

    void write_binary(std::ostream& out) const;                                 ///< Writes quaternion to a file stream in binary format
    void write_ASCII(std::ostream& out,
                     const int precision = 16, const char sep = ' ') const;     ///< Writes quaternion to a file stream in ASCII

    std::string get_output_string(
            const int precision = std::numeric_limits<double>::digits10 + 1,
            const char sep = ' ') const;                                        ///< Returns a string of quaternion values with given precision separated by "sep"
    std::vector<double> get_vector() const;                                     ///< Returns quaternion values as a 4-component vector
    std::vector<float>  get_vector_float() const;                               ///< Returns quaternion values as a 4-component vector of floats
//=========================== Deprecated methods ===============================
//    [[deprecated]]
//    void write(std::fstream& out) const;                                        ///< Writes quaternion to a file stream
//    [[deprecated]]
//    void write(std::stringstream& out) const;                                   ///< Writes quaternion to a string stream
//
//    [[deprecated]]
//    void read(std::fstream& inp);                                               ///< Reads quaternion from a file stream
//    [[deprecated]]
//    void read(std::stringstream& inp);                                          ///< Reads quaternion from a string stream
//==============================================================================
    dMatrix3x3 RotationMatrix;                                                  ///< Rotation matrix corresponding to the quaternion
 protected:
 private:
    double s;                                                                   ///< Real component
    dVector3 v;                                                                 ///< Imaginary components
};
}
#endif
