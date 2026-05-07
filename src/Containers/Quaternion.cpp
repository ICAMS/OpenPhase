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
 *  File created :   2014
 *  Main contributors :   Philipp Engels; Efim Borukhovich; Hesham Salama
 *
 */

#include "Containers/Quaternion.h"

namespace openphase
{

Quaternion::Quaternion()
{
    s = 1.0;
    v.set_to_zero();
    RotationMatrix.set_to_unity();
}

Quaternion::Quaternion(const Quaternion& rhS)
{
    s = rhS.s;
    v = rhS.v;
    set_rotation_matrix();
}

Quaternion::Quaternion(std::array<double,4> vecinit)
{
    int ii = 0;
    s = *vecinit.begin();
    for (auto it = vecinit.begin()+1; it != vecinit.end(); it++)
    {
        v[ii] = *it;
        ii += 1;
    }
    set_rotation_matrix();
}

double& Quaternion::operator[](const size_t index)
{
    assert(index < 4 && "Access beyond storage range");

    if (index == 0)
    {
        return s;
    };
    return v[index-1];
}

double const& Quaternion::operator[](const size_t index) const
{
    assert(index < 4 && "Access beyond storage range");

    if (index == 0){return s;};
    return v[index-1];
}

void Quaternion::set(const double sIn, const double x, const double y, const double z)
{
    s = sIn;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    set_rotation_matrix();
}

void Quaternion::set(const double sIn, const dVector3 vIn)
{
    s = sIn;
    v = vIn;
    set_rotation_matrix();
}

void Quaternion::set_entry(const int idx, const double val)
{
    if (idx == 0)
    {
        s = val;
    }
    else
    {
        v[idx] = val;
    }
    set_rotation_matrix();
}

void Quaternion::set(dVector3 Axis, const double Angle)
{
    Axis.normalize();
    // https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
    s = cos(Angle/2);
    v[0] = Axis[0] * sin(Angle/2);
    v[1] = Axis[1] * sin(Angle/2);
    v[2] = Axis[2] * sin(Angle/2);
    normalize();
    set_rotation_matrix();
}

void Quaternion::set_to_zero_rotation()
{
    s    = 1.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    RotationMatrix.set_to_unity();
}

void Quaternion::set_to_zero()
{
    s    = 0.0;
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    RotationMatrix.set_to_unity();
}

Quaternion& Quaternion::operator=(const Quaternion rhS)
{
    s = rhS.s;
    v = rhS.v;
    normalize();
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::operator+(const Quaternion rhS) const
{
    Quaternion result;

    result.s = s + rhS.s;
    result.v = v + rhS.v;
    //result.setRotationMatrix();                                               // RotationMatrix will be set upon assignment
    return result;
}

Quaternion& Quaternion::operator+=(const Quaternion rhS)
{
    s = s + rhS.s;
    v = v + rhS.v;
    normalize();
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::operator-(const Quaternion rhS) const
{
    Quaternion result;
    result.s = s - rhS.s;
    result.v = v - rhS.v;
    //result.setRotationMatrix();                                               // RotationMatrix will be set upon assignment
    return result;
}

Quaternion& Quaternion::operator-=(const Quaternion rhS)
{
    s = s - rhS.s;
    v = v - rhS.v;
    normalize();
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::operator*(const double scalar) const
{
    Quaternion result;
    result.s = s * scalar;
    result.v = v * scalar;
    //result.setRotationMatrix();                                               // RotationMatrix will be set upon assignment
    return result;
}

Quaternion& Quaternion::operator*=(const double scalar)
{
    s *= scalar;
    v *= scalar;
    normalize();
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::operator*(const Quaternion rhS) const
{
    Quaternion result;
    result.v[0] =  v[0]*rhS.s    + v[1]*rhS.v[2] - v[2]*rhS.v[1] + s*rhS.v[0];
    result.v[1] = -v[0]*rhS.v[2] + v[1]*rhS.s    + v[2]*rhS.v[0] + s*rhS.v[1];
    result.v[2] =  v[0]*rhS.v[1] - v[1]*rhS.v[0] + v[2]*rhS.s    + s*rhS.v[2];
    result.s    = -v[0]*rhS.v[0] - v[1]*rhS.v[1] - v[2]*rhS.v[2] + s*rhS.s;
    //result.setRotationMatrix();                                               // RotationMatrix will be set upon assignment
    return result;
}

Quaternion& Quaternion::operator*=(const Quaternion rhS)
{
    s = s*rhS.s - v*rhS.v;
    v = v.cross(rhS.v) + rhS.v*s + v*rhS.s;
    normalize();
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::operator/(const double divisor) const
{
    Quaternion result;
    result.s = s/divisor;
    result.v = v/divisor;
    //result.setRotationMatrix();                                               // RotationMatrix will be set upon assignment
    return result;
}

Quaternion& Quaternion::operator/=(const double divisor)
{
    s = s/divisor;
    v = v/divisor;
    normalize();
    set_rotation_matrix();
    return *this;
}

double Quaternion::length() const
{
    return sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

Quaternion& Quaternion::normalize()
{
    double norm = length();
    s /= norm;
    v /= norm;
    return *this;
}

Quaternion Quaternion::normalized() const
{
    Quaternion result;
    double norm = length();
    result.s = s/norm;
    result.v = v/norm;
    return result;
}

Quaternion& Quaternion::conjugate()
{
    v *= -1.0;
    set_rotation_matrix();
    return *this;
}

Quaternion Quaternion::conjugated() const
{
    Quaternion result;
    result.s = s;
    result.v = v*(-1.0);
    return result;
}

Quaternion& Quaternion::invert()
{
    // if |q|=1 then q inverse = q conjugate
    v *= -1.0;
    double norm = length();
    s /= norm;
    v /= norm;
    return *this;
}

Quaternion Quaternion::inverted() const
{
    Quaternion result;
    result.s = s;
    result.v = v*(-1.0);
    double norm = result.length();
    result.s /= norm;
    result.v /= norm;
    return result;
}

void Quaternion::pack(std::vector<double>& buffer)
{
    buffer.push_back(s);
    v.pack(buffer);
}

void Quaternion::unpack(std::vector<double>& buffer, size_t& it)
{
    s = buffer[it]; ++it;
    v.unpack(buffer, it);
}

std::string Quaternion::print(void) const
{
    std::stringstream out;
    out << "{ "   << s
        << ", < " << v[0]
        << ", "   << v[1]
        << ", "   << v[2] << " > }";
    return out.str();
}

void Quaternion::read_binary(std::istream& inp)
{
    inp >> s;
    for(int i = 0; i < 3; i++)
    {
        inp >> v[i];
    }
}

void Quaternion::read_ASCII(std::istream& inp)
{
    inp >> s;
    for(int i = 0; i < 3; i++)
    {
        inp >> v[i];
    }
}

void Quaternion::write_binary(std::ostream& out) const
{
    out << s;
    for(int i = 0; i < 3; i++)
    {
        out << v[i];
    }
}

void Quaternion::write_ASCII(std::ostream& out,
                                    const int precision, const char sep) const
{
    out << std::setprecision(precision) << std::defaultfloat;
    out << s << sep;
    for(int i = 0; i < 3; i++)
    {
        out << v[i] << sep;
    }
    out << std::endl;
}

//void Quaternion::write(std::fstream& out) const
//{
//    out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
//
//    out << s << " ";
//    for(int i = 0; i < 3; i++)
//    {
//        out << v[i] << " ";
//    }
//    out << std::endl;
//}
//
//void Quaternion::write(std::stringstream& out) const
//{
//    out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
//
//    out << s << " ";
//    for(int i = 0; i < 3; i++)
//    {
//        out << v[i] << " ";
//    }
//    out << std::endl;
//}
//
//void Quaternion::read(std::fstream& inp)
//{
//    inp >> s;
//    for(int i = 0; i < 3; i++)
//    {
//        inp >> v[i];
//    }
//}
//
//void Quaternion::read(std::stringstream& inp)
//{
//    inp >> s;
//    for(int i = 0; i < 3; i++)
//    {
//        inp >> v[i];
//    }
//}

std::string Quaternion::get_output_string(const int precision, const char sep) const
{
    std::stringstream out;
    out << std::setprecision(precision) << std::defaultfloat;
    out << s << sep;
    for(int i = 0; i < 3; i++)
    {
        out << v[i] << sep;
    }
    return out.str();
}

std::vector<double> Quaternion::get_vector() const
{
    std::vector<double> out(4);
    out[0] = s;
    for(int i = 0; i < 3; i++)
    {
        out[i+1] = v[i];
    }
    return out;
}

std::vector<float> Quaternion::get_vector_float() const
{
    std::vector<float> out(4);
    out[0] = (float)s;
    for(int i = 0; i < 3; i++)
    {
        out[i+1] = (float)v[i];
    }
    return out;
}


void Quaternion::set(const dMatrix3x3& RotMatrix)
{
    assert(std::fabs((RotMatrix.transposed()*RotMatrix - dMatrix3x3::UnitTensor()).determinant()) <= std::numeric_limits<float>::epsilon() && "matrix is not orthogonal");

    // Found at http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    RotationMatrix = RotMatrix;

    double trace = RotMatrix.trace();

    if(1.0 + trace > 0)
    {
        double ss = sqrt(trace + 1.0)*2.0;
        set(0.25*ss,
           (RotMatrix(2,1) - RotMatrix(1,2))/ss,
           (RotMatrix(0,2) - RotMatrix(2,0))/ss,
           (RotMatrix(1,0) - RotMatrix(0,1))/ss);
    }
    else
    {
        if (RotMatrix(0,0) > RotMatrix(1,1) && RotMatrix(0,0) > RotMatrix(2,2))
        {
            double ss = 2.0 * sqrt( 1.0 + RotMatrix(0,0) - RotMatrix(1,1) - RotMatrix(2,2));
            set((RotMatrix(2,1) - RotMatrix(1,2))/ss,
                 0.25 * ss,
                (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                (RotMatrix(0,2) + RotMatrix(2,0))/ss);
        }
        else if (RotMatrix(1,1) > RotMatrix(2,2))
        {
            double ss = 2.0 * sqrt(1.0 + RotMatrix(1,1) - RotMatrix(0,0) - RotMatrix(2,2));
            set((RotMatrix(0,2) - RotMatrix(2,0))/ss,
                (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                 0.25 * ss,
                (RotMatrix(1,2) + RotMatrix(2,1))/ss);
        }
        else
        {
            double ss = 2.0 * sqrt(1.0 + RotMatrix(2,2) - RotMatrix(0,0) - RotMatrix(1,1));
            set((RotMatrix(1,0) - RotMatrix(0,1))/ss,
                (RotMatrix(0,2) + RotMatrix(2,0))/ss,
                (RotMatrix(1,2) + RotMatrix(2,1))/ss,
                 0.25 * ss);
        }
    }
}

void Quaternion::set_rotation_matrix()
{
    if(length() < DBL_EPSILON)
    {
        RotationMatrix.set_to_unity();
    }
    else
    {
        Quaternion Qcopy = normalized();

        RotationMatrix(0,0) = 1.0 - 2.0*(Qcopy.v[1]*Qcopy.v[1] + Qcopy.v[2]*Qcopy.v[2]);
        RotationMatrix(0,1) =       2.0*(Qcopy.v[0]*Qcopy.v[1] -    Qcopy.s*Qcopy.v[2]);
        RotationMatrix(0,2) =       2.0*(Qcopy.v[0]*Qcopy.v[2] +    Qcopy.s*Qcopy.v[1]);
        RotationMatrix(1,0) =       2.0*(Qcopy.v[0]*Qcopy.v[1] +    Qcopy.s*Qcopy.v[2]);
        RotationMatrix(1,1) = 1.0 - 2.0*(Qcopy.v[0]*Qcopy.v[0] + Qcopy.v[2]*Qcopy.v[2]);
        RotationMatrix(1,2) =       2.0*(Qcopy.v[1]*Qcopy.v[2] -    Qcopy.s*Qcopy.v[0]);
        RotationMatrix(2,0) =       2.0*(Qcopy.v[0]*Qcopy.v[2] -    Qcopy.s*Qcopy.v[1]);
        RotationMatrix(2,1) =       2.0*(Qcopy.v[1]*Qcopy.v[2] +    Qcopy.s*Qcopy.v[0]);
        RotationMatrix(2,2) = 1.0 - 2.0*(Qcopy.v[0]*Qcopy.v[0] + Qcopy.v[1]*Qcopy.v[1]);
    }
}

dMatrix3x3 Quaternion::get_rotation_matrix(const bool Active)
{
    dMatrix3x3 RotationMatrix;

    Quaternion Qcopy = normalized();

    double q = Qcopy.s * Qcopy.s - (Qcopy.v[0] * Qcopy.v[0] + Qcopy.v[1] * Qcopy.v[1] + Qcopy.v[2] * Qcopy.v[2]);

    if(Active)
    {
        RotationMatrix(0,0) = q + 2.0 * Qcopy.v[0] * Qcopy.v[0];
        RotationMatrix(0,1) = 2.0 * (Qcopy.v[0] * Qcopy.v[1] - Qcopy.s * Qcopy.v[2]);
        RotationMatrix(0,2) = 2.0 * (Qcopy.v[0] * Qcopy.v[2] + Qcopy.s * Qcopy.v[1]);
        RotationMatrix(1,0) = 2.0 * (Qcopy.v[1] * Qcopy.v[0] + Qcopy.s * Qcopy.v[2]);
        RotationMatrix(1,1) = q + 2.0 *Qcopy.v[1] * Qcopy.v[1];
        RotationMatrix(1,2) = 2.0 * (Qcopy.v[1] * Qcopy.v[2] - Qcopy.s * Qcopy.v[0]);
        RotationMatrix(2,0) = 2.0 * (Qcopy.v[2] * Qcopy.v[0] - Qcopy.s * Qcopy.v[1]);
        RotationMatrix(2,1) = 2.0 * (Qcopy.v[2] * Qcopy.v[1] + Qcopy.s * Qcopy.v[0]);
        RotationMatrix(2,2) = q + 2.0 *Qcopy.v[2] * Qcopy.v[2];
    }
    else
    {
        RotationMatrix(0,0) = q + 2.0 * Qcopy.v[0] * Qcopy.v[0];
        RotationMatrix(0,1) = 2.0 * (Qcopy.v[0] * Qcopy.v[1] + Qcopy.s * Qcopy.v[2]);
        RotationMatrix(0,2) = 2.0 * (Qcopy.v[0] * Qcopy.v[2] - Qcopy.s * Qcopy.v[1]);
        RotationMatrix(1,0) = 2.0 * (Qcopy.v[1] * Qcopy.v[0] - Qcopy.s * Qcopy.v[2]);
        RotationMatrix(1,1) = q + 2.0 *Qcopy.v[1] * Qcopy.v[1];
        RotationMatrix(1,2) = 2.0 * (Qcopy.v[1] * Qcopy.v[2] + Qcopy.s * Qcopy.v[0]);
        RotationMatrix(2,0) = 2.0 * (Qcopy.v[2] * Qcopy.v[0] + Qcopy.s * Qcopy.v[1]);
        RotationMatrix(2,1) = 2.0 * (Qcopy.v[2] * Qcopy.v[1] - Qcopy.s * Qcopy.v[0]);
        RotationMatrix(2,2) = q + 2.0 *Qcopy.v[2] * Qcopy.v[2];
    }
    return RotationMatrix;
}

Quaternion Quaternion::lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2, const double t)
{
    // linear interpolation between two quaternions, where 0.0 < t < 1.0
    Quaternion result = (rhSQ1*(1.0 - t) + rhSQ2*t).normalized();
    return result;
}

Quaternion Quaternion::slerp(const Quaternion& Qa, const Quaternion& Qb, const double t)
{
    // linear interpolation between two quaternions, where 0.0 < t < 1.0
    // taken from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
    Quaternion result;

    double cosHalfTheta = Qa.s    * Qb.s
                        + Qa.v[0] * Qb.v[0]
                        + Qa.v[1] * Qb.v[1]
                        + Qa.v[2] * Qb.v[2];
    double sign = 1.0;
    if (cosHalfTheta < 0)
    {
        sign = -1.0;
    }
    if(fabs(cosHalfTheta) >= 1.0)                                               // both rotations equal
    {
        result = Qa;
        return result;
    }
    else
    {
        double halfTheta = acos(sign*cosHalfTheta);
        double sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
        if(fabs(sinHalfTheta) < 2.0*DBL_EPSILON)                                // the rotations are opposite of each other
        {
            result.s    = (Qa.s   *0.5 + sign*Qb.s   *0.5);
            result.v[0] = (Qa.v[0]*0.5 + sign*Qb.v[0]*0.5);
            result.v[1] = (Qa.v[1]*0.5 + sign*Qb.v[1]*0.5);
            result.v[2] = (Qa.v[2]*0.5 + sign*Qb.v[2]*0.5);
        }
        else                                                                    // general case
        {
            double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
            double ratioB = sin(t * halfTheta) / sinHalfTheta;

            result.s    = (Qa.s    * ratioA + sign*Qb.s    * ratioB);
            result.v[0] = (Qa.v[0] * ratioA + sign*Qb.v[0] * ratioB);
            result.v[1] = (Qa.v[1] * ratioA + sign*Qb.v[1] * ratioB);
            result.v[2] = (Qa.v[2] * ratioA + sign*Qb.v[2] * ratioB);
        }
        result.set_rotation_matrix();
    }
    return result;
}

}//namespace openphase
