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
 *  Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                        Hesham Salama
 *
 */
#include "BasicIncludes.h"
#include "Containers/EulerAngles.h"
#include "Containers/dVector3.h"
#include "Containers/dMatrix3x3.h"
#include "Containers/Quaternion.h"

namespace openphase
{

AngleModes AngleRepresentation(std::string locAngleRepresentationString)
{
    std::transform(locAngleRepresentationString.begin(),
                   locAngleRepresentationString.end(),
                   locAngleRepresentationString.begin(),
                   [](unsigned char c){ return std::toupper(c); });

    if(auto it  = AngleRepresentationTable.find(locAngleRepresentationString);
            it != AngleRepresentationTable.end())
    {
        return it->second;
    }

    return AngleModes::Radians; // default (not found)
}

EulerConventions EulerConvention(std::string locConventionString)
{
    std::transform(locConventionString.begin(),
                   locConventionString.end(),
                   locConventionString.begin(),
                   [](unsigned char c){ return std::toupper(c); });

    if (auto it  = EulerConventionTable.find(locConventionString);
             it != EulerConventionTable.end())
    {
        return it->second;
    }

    return EulerConventions::NNN; // default (not found)
}

EulerAngles::EulerAngles()
{
    is_set = false;

    Convention = EulerConventions::NNN;

    Q.fill(0.0);
    CosQ.fill(1.0);
    SinQ.fill(0.0);
}

EulerAngles::EulerAngles(std::initializer_list<double> Angles, EulerConventions locConvention)
{
    assert(Angles.size() == 3 && "Initialization list size is not equal to storage range");

    unsigned int ii = 0;
    for (auto it = Angles.begin(); it != Angles.end(); it++)
    {
        Q[ii] = *it;
        ii += 1;
    }

    Convention = locConvention;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

EulerAngles::EulerAngles(std::initializer_list<double> Angles, std::string locConventionString)
{
    assert(Angles.size() == 3 && "Initialization list size is not equal to storage range");

    unsigned int ii = 0;
    for (auto it = Angles.begin(); it != Angles.end(); it++)
    {
        Q[ii] = *it;
        ii += 1;
    }

    Convention = EulerConvention(locConventionString);

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

EulerAngles::EulerAngles(const EulerAngles& rhs)
{
    Convention = rhs.Convention;

    Q = rhs.Q;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = rhs.is_set;
}

inline bool EulerAngles::operator==(const EulerAngles& rhs)
{
    bool result = true;
    if(Convention == rhs.Convention)
    {
        for(int i = 0; i < 3; i++)
        if(fabs(Q[i] - rhs.Q[i]) > DBL_EPSILON)
        {
            result = false;
            break;
        }
    }
    else
    {
        result = false;
    }
    return result;
}

void EulerAngles::set(const double q1, const double q2, const double q3, std::string locConventionString)
{
    set(q1,q2,q3,EulerConvention(locConventionString));
}

void EulerAngles::set(const double q1, const double q2, const double q3, EulerConventions locConvention)
{
    Convention = locConvention;

    Q[0] = q1;
    Q[1] = q2;
    Q[2] = q3;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

void EulerAngles::set_trigonometric_functions()
{
    if(is_set == false)
    {
        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        is_set = true;
    }
}

void EulerAngles::set_to_zero(void)
{
    if(Convention != EulerConventions::NNN)
    {
        Q.fill(0.0);
        CosQ.fill(1.0);
        SinQ.fill(0.0);

        is_set = false;
    }
    else
    {
        std::cerr << " EulerAngles::set_to_zero(): Trying to set values of"
                  << " EulerAngles object that has no valid convention!"
                  << " Use EulerAngles::set_to_zero(q1, q2, q3, Convention)"
                  << " instead! Terminating!\n";
        OP_Exit(EXIT_FAILURE);
    }
}

void EulerAngles::set_convention(const std::string locConventionString)
{
    set_convention(EulerConvention(locConventionString));
}

void EulerAngles::set_convention(const EulerConventions locConvention)
{
    if(Convention == EulerConventions::NNN)
    {
        Convention = locConvention;
    }
    else
    {
        std::cerr << " EulerAngles::set_convention(): Trying to set convention of"
                  << " EulerAngles object that has a valid convention already!"
                  << " Use EulerAngles::set(q1, q2, q3, Convention)"
                  << " instead! Terminating!\n";
        OP_Exit(EXIT_FAILURE);
    }
}

//Conversion implemented from three.js Math //
void EulerAngles::set(const dMatrix3x3& RotMatrix, std::string locConventionString)
{
    set(RotMatrix, EulerConvention(locConventionString));
}

void EulerAngles::set(const dMatrix3x3& RotMatrix, EulerConventions locConvention)
{
    switch(locConvention)
    {
        case EulerConventions::XYZ:
        {
            Q[1] = asin(std::max(-1.0,std::min(1.0,RotMatrix(0,2))));
            if (fabs( RotMatrix(0,2) ) < 0.9999999 )
            {
                Q[0] = atan2( - RotMatrix(1,2), RotMatrix(2,2) );
                Q[2] = atan2( - RotMatrix(0,1), RotMatrix(0,0) );
            }
            else
            {
                Q[0] = atan2( RotMatrix(2,1), RotMatrix(1,1) );
                Q[2] = 0;
            }
            break;
        }
        case EulerConventions::YXZ:
        {
            Q[0] = asin(-std::max(-1.0,std::min(RotMatrix(1,2),1.0)));

            if (fabs( RotMatrix(1,2) ) < 0.9999999 )
            {
                Q[1] = atan2( RotMatrix(0,2), RotMatrix(2,2) );
                Q[2] = atan2( RotMatrix(1,0), RotMatrix(1,1) );
            }
            else
            {
                Q[1] = atan2( RotMatrix(2,0), RotMatrix(0,0) );
                Q[2] = 0;
            }
            break;
        }
        case EulerConventions::ZXY:
        {
            Q[0] = asin(std::max(-1.0,std::min(RotMatrix(2,1),1.0)));

            if (fabs( RotMatrix(2,0) ) < 0.9999999 )
            {
                Q[1] = atan2( RotMatrix(2,0), RotMatrix(2,2) );
                Q[2] = atan2( RotMatrix(0,1), RotMatrix(1,1) );
            }
            else
            {
                Q[1] = 0;
                Q[2] = atan2( RotMatrix(1,0), RotMatrix(0,0) );
            }
            break;
        }
        case EulerConventions::ZYX:
        {
            Q[1] = asin(-std::max(-1.0,std::min(RotMatrix(2,0),1.0)));

            if (fabs( RotMatrix(2,0) ) < 0.9999999 )
            {
                Q[0] = atan2( RotMatrix(2,1), RotMatrix(2,2) );
                Q[2] = atan2( RotMatrix(1,0), RotMatrix(0,0) );
            }
            else
            {
                Q[0] = 0;
                Q[2] = atan2( RotMatrix(0,1), RotMatrix(1,1) );
            }
            break;
        }
        case EulerConventions::YZX:
        {
            Q[2] = asin(std::max(-1.0,std::min(RotMatrix(1,0),1.0)));

            if (fabs( RotMatrix(1,0) ) < 0.9999999 )
            {
                Q[0] = atan2( RotMatrix(1,2), RotMatrix(1,1) );
                Q[1] = atan2( RotMatrix(2,0), RotMatrix(0,0) );
            }
            else
            {
                Q[0] = 0;
                Q[1] = atan2( RotMatrix(0,2), RotMatrix(2,2) );
            }
            break;
        }
        case EulerConventions::XZY:
        {
            Q[2] = asin(-std::max(-1.0,std::min(RotMatrix(0,1),1.0)));

            if (fabs( RotMatrix(0,1) ) < 0.9999999 )
            {
                Q[0] = atan2( RotMatrix(2,1), RotMatrix(1,1) );
                Q[1] = atan2( RotMatrix(0,2), RotMatrix(0,0) );
            }
            else
            {
                Q[0] = atan2( RotMatrix(1,2), RotMatrix(2,2) );
                Q[1] = 0;
            }
            break;
        }
        default:
        {
            std::cerr << " Wrong Euler Convention selected"
                      << " Check->EulerAngles::setFromRotationMatrix(dMatrix3x3)\n";
            OP_Exit(EXIT_FAILURE);
        }
    }

    Convention = locConvention;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

//void EulerAngles::set(dVector3 Axis, const double Angle, std::string locConventionString)
//{
//    set(Axis, Angle, EulerConvention(locConventionString));
//}
//
//void EulerAngles::set(dVector3 Axis, const double Angle, EulerConventions locConvention)
//{
//    double s = sin(Angle);
//    double c = cos(Angle);
//    double t = 1.0 - c;
//    double norm = Axis.length();
//    if(norm == 0)
//    {
//        std::cerr << "Axis length must be greater than zero"
//                  << "EulerAngles::set(axis,angle)"
//                  << "Terminating!\n";
//        OP_Exit(EXIT_FAILURE);
//    }
//    else
//    {
//        Axis.normalize();
//    }
//    // north pole singularity detected
//    if ((Axis[0]*Axis[1]*t + Axis[2]*s) > 0.998)
//    {
//        Q[0] = 2.0*atan2(Axis[0]*sin(Angle/2),cos(Angle/2));
//        Q[1] = Pi/2.0;
//        Q[2] = 0;
//    }
//    // south pole singularity detected
//    if ((Axis[0]*Axis[1]*t + Axis[2]*s) < -0.998)
//    {
//        Q[0] = -2.0*atan2(Axis[0]*sin(Angle/2),cos(Angle/2));
//        Q[1] = -Pi/2.0;
//        Q[2] = 0;
//    }
//    else
//    {
//        Q[0] = atan2(Axis[1] * s - Axis[0] * Axis[2] * t , 1.0 - (Axis[1]*Axis[1] + Axis[2]*Axis[2]) * t);
//        Q[1] = asin (Axis[0] * Axis[1] * t + Axis[2] * s);
//        Q[2] = atan2(Axis[0] * s - Axis[1] * Axis[2] * t , 1.0 - (Axis[0]*Axis[0] + Axis[2]*Axis[2]) * t);
//    }
//
//    Convention = locConvention;
//
//    SinQ[0] = sin(Q[0]);
//    SinQ[1] = sin(Q[1]);
//    SinQ[2] = sin(Q[2]);
//
//    CosQ[0] = cos(Q[0]);
//    CosQ[1] = cos(Q[1]);
//    CosQ[2] = cos(Q[2]);
//
//    is_set = true;
//}

void EulerAngles::add(const double q1, const double q2, const double q3)
{
    Q[0] += q1;
    Q[1] += q2;
    Q[2] += q3;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

EulerAngles& EulerAngles::operator=(const EulerAngles& rhs)
{
    assert(((Convention == rhs.Convention) or (Convention == EulerConventions::NNN)) and "Euler angles conventions do not coincide!");

    Convention = rhs.Convention;

    Q = rhs.Q;

    SinQ[0] = sin(rhs.Q[0]);
    SinQ[1] = sin(rhs.Q[1]);
    SinQ[2] = sin(rhs.Q[2]);

    CosQ[0] = cos(rhs.Q[0]);
    CosQ[1] = cos(rhs.Q[1]);
    CosQ[2] = cos(rhs.Q[2]);

    is_set = true;

    return *this;
}

EulerAngles EulerAngles::operator+(const EulerAngles& rhs) const
{
    assert(Convention == rhs.Convention && "Euler angle conventions do not coincide!");

    EulerAngles returnAng;

    returnAng.Convention = rhs.Convention;

    returnAng.Q[0] = Q[0] + rhs.Q[0];
    returnAng.Q[1] = Q[1] + rhs.Q[1];
    returnAng.Q[2] = Q[2] + rhs.Q[2];

    returnAng.SinQ[0] = sin(returnAng.Q[0]);
    returnAng.SinQ[1] = sin(returnAng.Q[1]);
    returnAng.SinQ[2] = sin(returnAng.Q[2]);

    returnAng.CosQ[0] = cos(returnAng.Q[0]);
    returnAng.CosQ[1] = cos(returnAng.Q[1]);
    returnAng.CosQ[2] = cos(returnAng.Q[2]);

    returnAng.is_set = true;

    return returnAng;
}

EulerAngles EulerAngles::operator-(const EulerAngles& rhs) const
{
    assert(Convention == rhs.Convention && "Euler angle conventions do not coincide!");

    EulerAngles returnAng;

    returnAng.Convention = rhs.Convention;

    returnAng.Q[0] = Q[0] - rhs.Q[0];
    returnAng.Q[1] = Q[1] - rhs.Q[1];
    returnAng.Q[2] = Q[2] - rhs.Q[2];

    returnAng.SinQ[0] = sin(returnAng.Q[0]);
    returnAng.SinQ[1] = sin(returnAng.Q[1]);
    returnAng.SinQ[2] = sin(returnAng.Q[2]);

    returnAng.CosQ[0] = cos(returnAng.Q[0]);
    returnAng.CosQ[1] = cos(returnAng.Q[1]);
    returnAng.CosQ[2] = cos(returnAng.Q[2]);

    returnAng.is_set = true;

    return returnAng;
}

EulerAngles& EulerAngles::operator+=(const EulerAngles& rhs)
{
    assert(Convention == rhs.Convention && "Euler angle conventions do not coincide!");

    Q[0] += rhs.Q[0];
    Q[1] += rhs.Q[1];
    Q[2] += rhs.Q[2];

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;

    return *this;
}

EulerAngles& EulerAngles::operator-=(const EulerAngles& rhs)
{
    assert(Convention == rhs.Convention && "Euler angle conventions do not coincide!");

    Q[0] += rhs.Q[0];
    Q[1] += rhs.Q[1];
    Q[2] += rhs.Q[2];

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;

    return *this;
}

std::string EulerAngles::get_convention(void) const
{
    std::string returnCon = EulerConventionStrings[(int)Convention];
    return returnCon;
}

EulerAngles EulerAngles::operator*(const double rhs) const
{
    EulerAngles returnAng;
    returnAng.set_convention(Convention);
    returnAng.set_to_zero();

    returnAng.Q[0] = Q[0]*rhs;
    returnAng.Q[1] = Q[1]*rhs;
    returnAng.Q[2] = Q[2]*rhs;

    returnAng.SinQ[0] = sin(returnAng.Q[0]);
    returnAng.SinQ[1] = sin(returnAng.Q[1]);
    returnAng.SinQ[2] = sin(returnAng.Q[2]);

    returnAng.CosQ[0] = cos(returnAng.Q[0]);
    returnAng.CosQ[1] = cos(returnAng.Q[1]);
    returnAng.CosQ[2] = cos(returnAng.Q[2]);

    return returnAng;
}

EulerAngles& EulerAngles::operator*=(const double rhs)
{
    Q[0] = Q[0]*rhs;
    Q[1] = Q[1]*rhs;
    Q[2] = Q[2]*rhs;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    return *this;
}

void EulerAngles::read_binary(std::istream& inp)
{
    for(size_t i = 0; i < 3; i++)
    {
        inp >> Q[i];
    }
    std::string ConventionS;
    for(size_t i = 0; i < 3; i++)
    {
        inp >> ConventionS[i];
    }

    Convention = EulerConvention(ConventionS);

    set_trigonometric_functions();
    is_set = true;
}

void EulerAngles::read_ASCII(std::istream& inp)
{
    for(size_t i = 0; i < 3; i++)
    {
        inp >> Q[i];
    }
    std::string ConventionS;
    inp >> ConventionS;

    Convention = EulerConvention(ConventionS);

    set_trigonometric_functions();
    is_set = true;
}

void EulerAngles::read_degrees(std::istream& inp)
{
    for(size_t i = 0; i < 3; i++)
    {
        inp >> Q[i];
        Q[i] *= Pi/180.0;
    }
    std::string ConventionS;
    inp >> ConventionS;

    Convention = EulerConvention(ConventionS);

    set_trigonometric_functions();
    is_set = true;
}

void EulerAngles::write_binary(std::ostream& out) const
{
    for(int i = 0; i < 3; i++)
    {
        out << Q[i];
    }
    for(int i = 0; i < 3; i++)
    {
        out << EulerConventionStrings[(int)Convention][i];
    }
}

void EulerAngles::write_ASCII(std::ostream& out,
                              const int precision, const char sep) const
{
    out << std::setprecision(precision) << std::defaultfloat;
    for(int i = 0; i < 3; i++)
    {
        out << Q[i] << sep;
    }
    out << EulerConventionStrings[(int)Convention];
    out << std::endl;
}

void EulerAngles::write_degrees(std::ostream& out,
                                const int precision, const char sep) const
{
    out << std::setprecision(precision) << std::defaultfloat;
    for(int i = 0; i < 3; i++)
    {
        out << Q[i]*180.0/Pi << sep;
    }
    out << EulerConventionStrings[(int)Convention];
    out << std::endl;
}


std::string EulerAngles::get_output_string(std::string ModeString,
                                           const int precision,
                                           const char sep) const
{
    return get_output_string(AngleRepresentation(ModeString));
}

std::string EulerAngles::get_output_string(AngleModes Mode,
                                           const int precision,
                                           const char sep) const
{
    std::stringstream out;
    out << std::setprecision(precision) << std::defaultfloat;
    for(int i = 0; i < 3; i++)
    {
        if(Mode == AngleModes::Radians)
        {
            out << Q[i] << sep;
        }
        else
        {
            out << Q[i]*180.0/Pi << sep;
        }
    }
    out << EulerConventionStrings[(int)Convention];
    return out.str();
}

std::string EulerAngles::print(void) const
{
    std::stringstream out;

    out << "["
        << Q[0]*180.0/Pi << ", "
        << Q[1]*180.0/Pi << ", "
        << Q[2]*180.0/Pi << "]["
        << EulerConventionStrings[(int)Convention]
        << "]";
    return out.str();
}

std::string EulerAngles::print_entire(void) const
{
    std::stringstream out;

    out << "Angle      [" << Q[0]*180.0/Pi << ", "
                          << Q[1]*180.0/Pi << ", "
                          << Q[2]*180.0/Pi << "]\n"
        << "Convention [" << EulerConventionStrings[(int)Convention] << "]\n"
        << "Sin        [" << SinQ[0] << ", "
                          << SinQ[1] << ", "
                          << SinQ[2] << "]\n"
        << "Cos        [" << CosQ[0] << ", "
                          << CosQ[1] << ", "
                          << CosQ[2] << "]\n"
        << "Is set:     " << is_set << "\n";
    return out.str();
}

dMatrix3x3 EulerAngles::get_rotation_matrix() const
{
    dMatrix3x3 result;
    double c1 = CosQ[0];
    double c2 = CosQ[1];
    double c3 = CosQ[2];
    double s1 = SinQ[0];
    double s2 = SinQ[1];
    double s3 = SinQ[2];

    //The rotation matrices follow notations from http://en.wikipedia.org/wiki/Euler_angles
    switch (Convention)
    {
        // Proper Euler angles
        case EulerConventions::XZX:
        {
            double XZX[3][3] =
            {{                 c2,            - c3*s2,              s2*s3},
             {              c1*s2,   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3},
             {              s1*s2,   c1*s3 + c2*c3*s1,   c1*c3 - c2*s1*s3}};

            memmove(result.data(), XZX, 9*sizeof(double));
            break;
        }
        case EulerConventions::XYX:
        {
            double XYX[3][3] =
            {{                 c2,             s2*s3,              c3*s2},
             {              s1*s2,   c1*c3 - c2*s1*s3, - c1*s3 - c2*c3*s1},
             {            - c1*s2,   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3}};

            memmove(result.data(), XYX, 9*sizeof(double));
            break;
        }
        case EulerConventions::YXY:
        {
            double YXY[3][3] =
            {{   c1*c3 - c2*s1*s3,              s1*s2,   c2*c3*s1 + c1*s3},
             {              s2*s3,                 c2,            - c3*s2},
             { - c3*s1 - c1*c2*s3,              c1*s2,   c1*c2*c3 - s1*s3}};
            memmove(result.data(), YXY, 9*sizeof(double));
            break;
        }
        case EulerConventions::YZY:
        {
            double YZY[3][3] =
            {{   c1*c2*c3 - s1*s3,            - c1*s2,   c3*s1 + c1*c2*s3},
             {              c3*s2,                 c2,              s2*s3},
             {-(c2*c3*s1) - c1*s3,              s1*s2,   c1*c3 - c2*s1*s3}};
            memmove(result.data(), YZY, 9*sizeof(double));
            break;
        }
        case EulerConventions::ZYZ:
        {
            double ZYZ[3][3] =
            {{   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3,             c1*s2},
             {   c2*c3*s1 + c1*s3,   c1*c3 - c2*s1*s3,             s1*s2},
             {            - c3*s2,              s2*s3,                c2}};
            memmove(result.data(), ZYZ, 9*sizeof(double));
            break;
        }
        case EulerConventions::ZXZ:
        {
            double ZXZ[3][3] =
            {{   c1*c3 - c2*s1*s3, - c2*c3*s1 - c1*s3,             s1*s2},
             {   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3,           - c1*s2},
             {              s2*s3,              c3*s2,                c2}};
            memmove(result.data(), ZXZ, 9*sizeof(double));
            break;
        }

        // Tait-Bryan angles
        case EulerConventions::XYZ:
        {
            double XYZ[3][3] =
            {{              c2*c3,            - c2*s3,                s2},
             {   c1*s3 + c3*s1*s2,   c1*c3 - s1*s2*s3,           - c2*s1},
             {   s1*s3 - c1*c3*s2,   c1*s2*s3 + c3*s1,             c1*c2}};
            memmove(result.data(), XYZ, 9*sizeof(double));
            break;
        }
        case EulerConventions::ZYX:
        {
            double ZYX[3][3] =
            {{              c1*c2,  c1*s2*s3 - c3*s1,   s1*s3 + c1*c3*s2},
             {              c2*s1,  c1*c3 - s1*s2*s3,   c3*s1*s2 - c1*s3},
             {               - s2,             c2*s3,              c2*c3}};
            memmove(result.data(), ZYX, 9*sizeof(double));
            break;
        }
        default:
        {
            std::cerr << "Wrong Euler convention is used!"
                      << " Check EulerAngles::getRotationMatrix()"<< std::endl;
            OP_Exit(EXIT_FAILURE);
        }
    }
    return result;
}

/*
 * Source:  https://de.mathworks.com/help/fusion/ref/quaternion.euler.html
 *          https://scholar.google.de/scholar?cluster=3204262265835591787
 *
 *          Switch between Active(point rotation) and Passive(frame rotation) with a flag
 *          --> By default: Active rotation
 */
void EulerAngles::set(Quaternion Quat, std::string locConventionString, const bool active)
{
    set(Quat, EulerConvention(locConventionString), active);
}

void EulerAngles::set(Quaternion Quat, EulerConventions locConvention, const bool active)
{
    Quat.normalize();

    double q1 = Quat[0];
    double q2 = Quat[1];
    double q3 = Quat[2];
    double q4 = Quat[3];

    switch(locConvention)
    {
        case EulerConventions::XYX:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q3-q1*q4), -2.0*(q2*q4+q1*q3));
                Q[1] = -acos(q2*q2 + q1*q1 - q4*q4 - q3*q3);
                Q[2] = atan2(-2.0*(q2*q3+q1*q4), 2.0*(q2*q4-q1*q3));
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q3+q1*q4),-2.0*(q2*q4-q1*q3));
                Q[1] = acos(q2*q2 + q1*q1 - q4*q4 - q3*q3);
                Q[2] = atan2(2.0*(q2*q3-q1*q4),2.0*(q2*q4+q1*q3));
            }
            break;
        }
        case EulerConventions::XYZ:
        {
            if(active)
            {
                Q[0] = atan2(2.0*(q3*q4+q1*q2), q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[1] = asin(-2.0*(q2*q4-q1*q3));
                Q[2] = atan2(2.0*(q2*q3+q1*q4), q1*q1 + q2*q2 - q3*q3 - q4*q4);
            }
            else
            {
                Q[0] = atan2(-2.0*(q3*q4-q1*q2), q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[1] = asin(2.0*(q2*q4+q1*q3));
                Q[2] = atan2(-2.0*(q2*q3-q1*q4), q1*q1 + q2*q2 - q3*q3 - q4*q4);
            }
            break;
        }
        case EulerConventions::XZX:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q4+q1*q3), 2.0*(q2*q3-q1*q4));
                Q[1] = -acos(q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[2] = atan2(-2.0*(q2*q4-q1*q3), -2.0*(q2*q3+q1*q4));
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q4-q1*q3), 2.0*(q2*q3+q1*q4));
                Q[1] = acos(q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[2] = atan2(2.0*(q2*q4+q1*q3),-2.0*(q2*q3-q1*q4));

            }
            break;
        }
        case EulerConventions::XZY:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q3*q4-q1*q2), q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[1] = asin(2.0*(q2*q3+q1*q4));
                Q[2] = atan2(-2.0*(q2*q4-q1*q3), q1*q1 + q2*q2 - q3*q3 - q4*q4);
            }
            else
            {
                Q[0] = atan2(2.0*(q3*q4+q1*q2), q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[1] = asin(-2.0*(q2*q3-q1*q4));
                Q[2] = atan2(2.0*(q2*q4+q1*q3), q1*q1 + q2*q2 - q3*q3 - q4*q4);
            }
            break;
        }
        case EulerConventions::YXY:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q3+q1*q4), 2.0*(q3*q4-q1*q2));
                Q[1] = -acos(q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[2] = atan2(-2.0*(q2*q3-q1*q4), -2.0*(q3*q4+q1*q2));
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q3-q1*q4), 2.0*(q3*q4+q1*q2));
                Q[1] = acos(q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[2] = atan2(2.0*(q2*q3+q1*q4),-2.0*(q3*q4-q1*q2));

            }
            break;
        }
        case EulerConventions::YXZ:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q4-q1*q3), q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[1] = asin(2.0*(q3*q4+q1*q2));
                Q[2] = atan2(-2.0*(q2*q3-q1*q4), q1*q1 - q2*q2 + q3*q3 - q4*q4);
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q4+q1*q3), q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[1] = asin(-2.0*(q3*q4-q1*q2));
                Q[2] = atan2(2.0*(q2*q3+q1*q4), q1*q1 - q2*q2 + q3*q3 - q4*q4);
            }
            break;
        }
        case EulerConventions::YZX:
        {
            if(active)
            {
                Q[0] = atan2(2.0*(q2*q4+q1*q3), q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[1] = asin(-2.0*(q2*q3-q1*q4));
                Q[2] = atan2(2.0*(q3*q4+q1*q2), q1*q1 - q2*q2 + q3*q3 - q4*q4);
            }
            else
            {
                Q[0] = atan2(-2.0*(q2*q4-q1*q3), q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[1] = asin(2.0*(q2*q3+q1*q4));
                Q[2] = atan2(-2.0*(q3*q4-q1*q2), q1*q1 - q2*q2 + q3*q3 - q4*q4);
            }
            break;
        }
        case EulerConventions::YZY:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q3*q4-q1*q2), -2.0*(q2*q3+q1*q4));
                Q[1] = -acos(q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[2] = atan2(-2.0*(q3*q4+q1*q2), 2.0*(q2*q3-q1*q4));
            }
            else
            {
                Q[0] = atan2(2.0*(q3*q4+q1*q2),-2.0*(q2*q3-q1*q4));
                Q[1] = acos(q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[2] = atan2(2.0*(q3*q4-q1*q2), 2.0*(q2*q3+q1*q4));

            }
            break;
        }
        case EulerConventions::ZXY:
        {
            if(active)
            {
                Q[2] = atan2(2.0*(q2*q4+q1*q3), q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[1] = asin(-2.0*(q3*q4-q1*q2));
                Q[0] = atan2(2.0*(q2*q3+q1*q4), q1*q1 - q2*q2 + q3*q3 - q4*q4);
            }
            else
            {
                Q[0] = atan2(-2.0*(q2*q3-q1*q4), q1*q1 - q2*q2 + q3*q3 - q4*q4);
                Q[1] = asin(2.0*(q3*q4+q1*q2));
                Q[2] = atan2(-2.0*(q2*q4-q1*q3), q1*q1 - q2*q2 - q3*q3 + q4*q4);
            }
            break;
        }
        case EulerConventions::ZXZ:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q4-q1*q3), -2.0*(q3*q4+q1*q2));
                Q[1] =-acos(q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[2] = atan2(-2.0*(q2*q4+q1*q3), 2.0*(q3*q4-q1*q2));
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q4+q1*q3), -2.0*(q3*q4-q1*q2));
                Q[1] = acos(q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[2] = atan2(2.0*(q2*q4-q1*q3), 2.0*(q3*q4+q1*q2));
            }
            break;
        }
        case EulerConventions::ZYX:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q2*q3-q1*q4), q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[1] = asin(2.0*(q2*q4+q1*q3));
                Q[2] = atan2(-2.0*(q3*q4-q1*q2), q1*q1 - q2*q2 - q3*q3 + q4*q4);
            }
            else
            {
                Q[0] = atan2(2.0*(q2*q3+q1*q4), q1*q1 + q2*q2 - q3*q3 - q4*q4);
                Q[1] = asin(-2.0*(q2*q4-q1*q3));
                Q[2] = atan2(2.0*(q3*q4+q1*q2), q1*q1 - q2*q2 - q3*q3 + q4*q4);
            }
            break;
        }
        case EulerConventions::ZYZ:
        {
            if(active)
            {
                Q[0] = atan2(-2.0*(q3*q4+q1*q2), 2.0*(q2*q4-q1*q3));
                Q[1] = -acos(q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[2] = atan2(-2.0*(q3*q4-q1*q2), -2.0*(q2*q4+q1*q3));
            }
            else
            {
                Q[0] = atan2(2.0*(q3*q4-q1*q2), 2.0*(q2*q4+q1*q3));
                Q[1] = acos(q1*q1 - q2*q2 - q3*q3 + q4*q4);
                Q[2] = atan2(2.0*(q3*q4+q1*q2), -2.0*(q2*q4-q1*q3));
            }
            break;
        }
        default:
        {
            std::cerr << " Wrong Euler Convention selected"
                      << " Check EulerAngles::set(Quaternion)" << std::endl;
            OP_Exit(EXIT_FAILURE);
        }
    }
    Convention = locConvention;

    SinQ[0] = sin(Q[0]);
    SinQ[1] = sin(Q[1]);
    SinQ[2] = sin(Q[2]);

    CosQ[0] = cos(Q[0]);
    CosQ[1] = cos(Q[1]);
    CosQ[2] = cos(Q[2]);

    is_set = true;
}

/*Source: http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
 *        https://www.astro.rug.nl/software/kapteyn/_downloads/fa29752e4cd69adcfa2fc03b1c020f4e/attitude.pdf
 *        https://de.mathworks.com/help/fusion/ref/quaternion.html#mw_54b60363-8447-46b1-b461-b50aef77772d
 *
 *        Switch between Active (point rotation) and Passive (frame rotation) with a flag
 *        Default setting: Active rotation
 */
Quaternion EulerAngles::get_quaternion(const bool Active) const
{
    double CosAng1h = cos(0.5*Q[0]);
    double SinAng1h = sin(0.5*Q[0]);

    double CosAng2h = cos(0.5*Q[1]);
    double SinAng2h = sin(0.5*Q[1]);

    double CosAng3h = cos(0.5*Q[2]);
    double SinAng3h = sin(0.5*Q[2]);

    Quaternion result;

    switch (Convention)
    {
        case EulerConventions::XYX:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0),
                            sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0));
            }
            break;
        }
        case EulerConventions::XYZ:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h);
            }
            else
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h);
            }
            break;
        }
        case EulerConventions::XZX:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0),
                           -sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0));
            }
            break;
        }
        case EulerConventions::XZY:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h);
            }
            else
            {
                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h);
            }
            break;
        }
        case EulerConventions::YXY:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0),
                           -sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0));
            }
            break;
        }
        case EulerConventions::YXZ:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h);
            }
            else
            {
                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h);
            }
            break;
        }
        case EulerConventions::YZX:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h);
            }
            else
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h);
            }
            break;
        }
        case EulerConventions::YZY:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                            sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0));
            }
            break;
        }
        case EulerConventions::ZXY:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            else
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            break;
        }
        case EulerConventions::ZXZ:
        {
            if(Active)
            {

                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0),
                            sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0));
            }
            break;
        }
        case EulerConventions::ZYX:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                            CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            else
            {

                result.set( CosAng1h*CosAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h - SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*CosAng2h*SinAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            break;
        }
        case EulerConventions::ZYZ:
        {
            if(Active)
            {
                result.set( CosAng1h*CosAng2h*CosAng3h - SinAng1h*CosAng2h*SinAng3h,
                           -CosAng1h*SinAng2h*SinAng3h + SinAng1h*SinAng2h*CosAng3h,
                            CosAng1h*SinAng2h*CosAng3h + SinAng1h*SinAng2h*SinAng3h,
                            CosAng1h*CosAng2h*SinAng3h + SinAng1h*CosAng2h*CosAng3h);
            }
            else
            {
                result.set( cos(0.5*Q[1])*cos((Q[0]+Q[2])/2.0),
                           -sin(0.5*Q[1])*sin((Q[0]-Q[2])/2.0),
                            sin(0.5*Q[1])*cos((Q[0]-Q[2])/2.0),
                            cos(0.5*Q[1])*sin((Q[0]+Q[2])/2.0));
            }
            break;
        }
        default:
        {
            std::cerr << "Wrong or not supported convention."
                      << " Check EulerAngles::getQuaternion()" << std::endl;
            OP_Exit(EXIT_FAILURE);
        }
    }
    return result;
}

void EulerAngles::get_axis_angle(dVector3& Axis, double& Angle) const
{
    //Source: https://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToAngle/index.html

    //Needs care!!! Implementation assumes only one convention for Euler angles from avionics (heading, attitude, bank).

    double c1 = CosQ[0];
    double c2 = CosQ[1];
    double c3 = CosQ[2];
    double s1 = SinQ[0];
    double s2 = SinQ[1];
    double s3 = SinQ[2];

    double w = c1*c2*c3 - s1*s2*s3;

    // When all Euler angles are zero (Angle = 0) we can set axis to arbitrary direction to avoid division by zero
    if (fabs(w) < FLT_EPSILON)
    {
          Angle   = 0.0;
          Axis[0] = 1.0;
          Axis[1] = 0.0;
          Axis[2] = 0.0;
    }
    else
    {
        Angle = 2.0 * acos(w);

        Axis[0] = c1*c2*s3 + s1*s2*c3;
        Axis[1] = s1*c2*c3 + c1*s2*s3;
        Axis[2] = c1*s2*c3 - s1*c2*s3;

        Axis.normalize();
    }
}

}//namespace openphase
