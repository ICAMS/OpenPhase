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

 *   File created :   2014
 *   Main contributors :   Philipp Engels; Efim Borukhovich; Hesham Salama
 *
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#include "Includes.h"

namespace openphase
{
class EulerAngles;

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
    void set_to_zero();                                                         ///< Sets zero rotation quaternion

    double& operator[](const size_t index);                                     ///< Bi-directional access operator for accessing individual entries
    double const& operator[](const size_t index) const;                         ///< Bi-directional access operator for accessing individual entries
    static constexpr size_t size() {return 4u;};

    Quaternion& operator=(const Quaternion rhS);
    Quaternion  operator+(const Quaternion rhS) const;
    Quaternion& operator+=(const Quaternion rhS);
    Quaternion  operator-(const Quaternion rhS) const;
    Quaternion& operator-=(const Quaternion rhS);
    Quaternion  operator*(const double scalar) const;
    Quaternion& operator*=(const double scalar);
    Quaternion  operator*(const Quaternion rhS) const;
    Quaternion& operator*=(const Quaternion rhS);
    Quaternion  operator/(const double divisor);
    Quaternion& operator/=(const double divisor);

    double length() const;                                                      ///< Returns length of the quaternion
    Quaternion& normalize();                                                    ///< Normalizes quaternion
    Quaternion  normalized() const;                                             ///< Returns normalized quaternion
    Quaternion& conjugate();                                                    ///< Conjugates quaternion (reverse rotation)
    Quaternion  conjugated() const;                                             ///< Returns conjugated quaternion (reverse rotation)
    Quaternion& invert();                                                       ///< Conjugates quaternion (reverse rotation if quaternion is normalized)
    Quaternion  inverted() const;                                               ///< Returns inverted quaternion (reverse rotation if quaternion is normalized)

    void pack(std::vector<double>& buffer);                                     ///< Packs the quaternion into communication buffer
    void unpack(std::vector<double>& buffer, size_t& it);                       ///< Unpacks the quaternion from communication buffer

    void setRotationMatrix(void);                                               ///< Sets rotation matrix for the quaternion
    dMatrix3x3 getRotationMatrix(const bool Active = true);                     ///< Returns rotation matrix (default: active rotation)

    static Quaternion lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2,
                                                            const double t);    ///< Linear interpolation between two quaternions, where 0.0 < t < 1.0
    static Quaternion slerp(const Quaternion& a, const Quaternion& b,
                                                            const double t);    ///< Linear interpolation between two quaternions, where 0.0 < t < 1.0
    std::string print(void) const;                                              ///< Returns the quaternion as formated string

    void write(std::fstream& out) const;                                        ///< Writes quaternion to a file stream
    void write(std::stringstream& out) const;                                   ///< Writes quaternion to a string stream

    void read(std::fstream& inp);                                               ///< Reads quaternion from a file stream
    void read(std::stringstream& inp);                                          ///< Reads quaternion from a string stream

    dMatrix3x3 RotationMatrix;                                                  ///< Rotation matrix corresponding to the quaternion
 protected:
 private:
    double s;                                                                   ///< Real component
    dVector3 v;                                                                 ///< Imaginary components
};

inline Quaternion::Quaternion()
{
    s = 1.0;
    v.set_to_zero();
    RotationMatrix.set_to_unity();
}

inline Quaternion::Quaternion(const Quaternion& rhS)
{
    s = rhS.s;
    v = rhS.v;
    setRotationMatrix();
}

inline Quaternion::Quaternion(std::array<double,4> vecinit)
{
    //assert(vecinit.size() == 4 && "Initialization list size is not equal to storage range");

    int ii = 0;
    s = *vecinit.begin();
    for (auto it = vecinit.begin()+1; it != vecinit.end(); it++)
    {
        v[ii] = *it;
        ii += 1;
    }
    setRotationMatrix();
}

inline double& Quaternion::operator[](const size_t index)
{
    assert(index < 4 && "Access beyond storage range");

    if (index == 0)
    {
        return s;
    };
    return v[index-1];
}

inline double const& Quaternion::operator[](const size_t index) const
{
    assert(index < 4 && "Access beyond storage range");

    if (index == 0){return s;};
    return v[index-1];
}

inline void Quaternion::set(const double sIn, const double x, const double y, const double z)
{
    s = sIn;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    setRotationMatrix();
}

inline void Quaternion::set(const double sIn, const dVector3 vIn)
{
    s = sIn;
    v = vIn;
    setRotationMatrix();
}

inline void Quaternion::set_entry(const int idx, const double val)
{
    if (idx == 0)
    {
        s = val;
    }
    else
    {
        v[idx] = val;
    }
    setRotationMatrix();
}

inline void Quaternion::set(dVector3 Axis, const double Angle)
{
    Axis.normalize();
    // https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
    s = cos(Angle/2);
    v[0] = Axis[0] * sin(Angle/2);
    v[1] = Axis[1] * sin(Angle/2);
    v[2] = Axis[2] * sin(Angle/2);
    normalize();
    setRotationMatrix();
}

inline void Quaternion::set_to_zero()
{
    s    = 1.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    RotationMatrix.set_to_unity();
}

inline Quaternion& Quaternion::operator=(const Quaternion rhS)
{
    s = rhS.s;
    v = rhS.v;
    setRotationMatrix();
    return *this;
}

inline Quaternion Quaternion::operator+(const Quaternion rhS) const
{
    Quaternion result;

    result.s = s + rhS.s;
    result.v = v + rhS.v;
    //result.setRotationMatrix();                                               // RotationMatrix will be set in assignment operator
    return result;
}

inline Quaternion& Quaternion::operator+=(const Quaternion rhS)
{
    s = s + rhS.s;
    v = v + rhS.v;
    setRotationMatrix();
    return *this;
}

inline Quaternion Quaternion::operator-(const Quaternion rhS) const
{
    Quaternion result;
    result.s = s - rhS.s;
    result.v = v - rhS.v;
    //result.setRotationMatrix();                                               // RotationMatrix will be set in assignment operator
    return result;
}

inline Quaternion& Quaternion::operator-=(const Quaternion rhS)
{
    s = s - rhS.s;
    v = v - rhS.v;
    setRotationMatrix();
    return *this;
}

inline Quaternion Quaternion::operator*(const double scalar) const
{
    Quaternion result;
    result.s = s * scalar;
    result.v = v * scalar;
    //result.setRotationMatrix();                                               // RotationMatrix will be set in assignment operator
    return result;
}

inline Quaternion& Quaternion::operator*=(const double scalar)
{
    s *= scalar;
    v *= scalar;
    setRotationMatrix();
    return *this;
}

inline Quaternion Quaternion::operator*(const Quaternion rhS) const
{
    Quaternion result;
    result.v[0] =  v[0]*rhS.s    + v[1]*rhS.v[2] - v[2]*rhS.v[1] + s*rhS.v[0];
    result.v[1] = -v[0]*rhS.v[2] + v[1]*rhS.s    + v[2]*rhS.v[0] + s*rhS.v[1];
    result.v[2] =  v[0]*rhS.v[1] - v[1]*rhS.v[0] + v[2]*rhS.s    + s*rhS.v[2];
    result.s    = -v[0]*rhS.v[0] - v[1]*rhS.v[1] - v[2]*rhS.v[2] + s*rhS.s;
    //result.setRotationMatrix();                                               // RotationMatrix will be set in assignment operator
    return result;
}

inline Quaternion& Quaternion::operator*=(const Quaternion rhS)
{
    s = s*rhS.s - v*rhS.v;
    v = v.cross(rhS.v) + rhS.v*s + v*rhS.s;
    setRotationMatrix();
    return *this;
}

inline Quaternion Quaternion::operator/(const double divisor)
{
    Quaternion result;
    result.s = s/divisor;
    result.v = v/divisor;
    //result.setRotationMatrix();                                               // RotationMatrix will be set in assignment operator
    return result;
}

inline Quaternion& Quaternion::operator/=(const double divisor)
{
    s = s/divisor;
    v = v/divisor;
    setRotationMatrix();
    return *this;
}

inline double Quaternion::length() const
{
    return sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline Quaternion& Quaternion::normalize()
{
    *this /= sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return *this;
}

inline Quaternion Quaternion::normalized() const
{
    Quaternion result;
    double norm = length();
    result.s = s/norm;
    result.v = v/norm;
    return result;
}

inline Quaternion& Quaternion::conjugate()
{
    v = v*(-1.0);
    setRotationMatrix();
    return *this;
}
inline Quaternion Quaternion::conjugated() const
{
    Quaternion result;
    result.s = s;
    result.v = v*(-1.0);
    return result;
}

inline Quaternion& Quaternion::invert()
{
    // if |q|=1 then q inverse = q conjugate
    v = v*(-1.0);
    *this /= (s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return *this;
}

inline Quaternion Quaternion::inverted() const
{
    Quaternion result = *this;
    return result = result.conjugate()/(s*s+v*v);
}

inline void Quaternion::pack(std::vector<double>& buffer)
{
    buffer.push_back(s);
    v.pack(buffer);
}

inline void Quaternion::unpack(std::vector<double>& buffer, size_t& it)
{
    s = buffer[it]; ++it;
    v.unpack(buffer, it);
}

inline std::string Quaternion::print(void) const
{
    std::stringstream out;
    out << "{ "   << s
        << ", < " << v[0]
        << ", "   << v[1]
        << ", "   << v[2] << " > }";
    return out.str();
}

inline void Quaternion::write(std::fstream& out) const
{
    out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    out << s << " ";
    for(int i = 0; i < 3; i++)
    {
        out << v[i] << " ";
    }
    out << std::endl;
}

inline void Quaternion::write(std::stringstream& out) const
{
    out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    out << s << " ";
    for(int i = 0; i < 3; i++)
    {
        out << v[i] << " ";
    }
    out << std::endl;
}

inline void Quaternion::read(std::fstream& inp)
{
    inp >> s;
    for(int i = 0; i < 3; i++)
    {
        inp >> v[i];
    }
}

inline void Quaternion::read(std::stringstream& inp)
{
    inp >> s;
    for(int i = 0; i < 3; i++)
    {
        inp >> v[i];
    }
}

}
#endif
