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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Hesham Salama
 *
 */

#ifndef EULERANGLES_H
#define EULERANGLES_H

#include "Includes.h"

namespace openphase
{

enum EulerConvention{ XZX, XYX, YXY, YZY, ZXZ, ZYZ, /* proper Euler angles*/
                      XYZ, YZX, ZXY, XZY, ZYX, YXZ, /* Tait–Bryan angles*/
                      NNN /*default convention -> convention not set*/};

const std::vector<std::string>
     EulerConventionS{"XZX", "XYX", "YXY", "YZY", "ZXZ", "ZYZ", /* proper Euler angles*/
                      "XYZ", "YZX", "ZXY", "XZY", "ZYX", "YXZ", /* Tait–Bryan angles*/
                      "NNN" /*default convention -> convention not set*/};

class Quaternion;

class OP_EXPORTS EulerAngles                                                    ///< Orientation angles and their Cos() and Sin().
{
//  Storages:
 public:
    union
    {
        double Q[3];
    };

    union
    {
        double CosQ[3];
    };

    union
    {
        double SinQ[3];
    };

    EulerConvention Convention;                                                 ///< Stores the convention for the Euler angles

    bool IsSet;                                                                 ///< Indicates whether sin/cos are set.

//  Methods:

    EulerAngles(): IsSet(false)                                                 ///< Constructor
    {
        Convention = NNN;

        Q[0] = 0.0;
        Q[1] = 0.0;
        Q[2] = 0.0;

        SinQ[0] = 0.0;
        SinQ[1] = 0.0;
        SinQ[2] = 0.0;

        CosQ[0] = 1.0;
        CosQ[1] = 1.0;
        CosQ[2] = 1.0;
    };

    EulerAngles(std::initializer_list<double> Angles, EulerConvention locConvention)
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

        IsSet = true;
    }

    EulerAngles(const EulerAngles& rhs)
    {
        Convention = rhs.Convention;

        Q[0] = rhs.Q[0];
        Q[1] = rhs.Q[1];
        Q[2] = rhs.Q[2];

        SinQ[0] = rhs.SinQ[0];
        SinQ[1] = rhs.SinQ[1];
        SinQ[2] = rhs.SinQ[2];

        CosQ[0] = rhs.CosQ[0];
        CosQ[1] = rhs.CosQ[1];
        CosQ[2] = rhs.CosQ[2];

        IsSet = rhs.IsSet;
    };

    bool operator==(const EulerAngles& rhs)
    {
        double epsilon = fabs(Q[0] - rhs.Q[0]) + fabs(Q[1] - rhs.Q[1]) + fabs(Q[2] - rhs.Q[2]);
        if(Convention == rhs.Convention and epsilon < 3.0*DBL_EPSILON)
        {
            return true;
        }
        return false;
    }

    void set(const double q1, const double q2, const double q3, EulerConvention locConvention)
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

        IsSet = true;
    };

    void setTrigonometricFunctions()
    {
        if(IsSet == false)
        {
            SinQ[0] = sin(Q[0]);
            SinQ[1] = sin(Q[1]);
            SinQ[2] = sin(Q[2]);

            CosQ[0] = cos(Q[0]);
            CosQ[1] = cos(Q[1]);
            CosQ[2] = cos(Q[2]);

            IsSet = true;
        }
    };

    void set_to_zero(void)
    {
        if(Convention != NNN)
        {
            Q[0] = 0.0;
            Q[1] = 0.0;
            Q[2] = 0.0;

            SinQ[0] = 0.0;
            SinQ[1] = 0.0;
            SinQ[2] = 0.0;

            CosQ[0] = 1.0;
            CosQ[1] = 1.0;
            CosQ[2] = 1.0;

            IsSet = false;
        }
        else
        {
            std::cerr << " EulerAngles::set_to_zero(): Trying to set values of"
                      << " EulerAngles object that has no valid convention!"
                      << " Use EulerAngles::set_to_zero(q1, q2, q3, Convention)"
                      << " instead! Terminating!\n";
            OP_Exit(EXIT_FAILURE);
        }
    };

    void set_convention(const EulerConvention locConvention)
    {
        if(Convention == NNN)
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
    };

    //Conversion implemented from three.js Math //
    void set(dMatrix3x3 RotMatrix, EulerConvention locConvention)
    {
        switch(locConvention)
        {
            case XYZ:
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
            case YXZ:
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
            case ZXY:
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
            case ZYX:
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
            case YZX:
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
            case XZY:
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

        IsSet = true;
    }

    void set(Quaternion Quat, EulerConvention locConvention, const bool active = true);

    void set(dVector3 Axis, const double Angle)
    {
        double s = sin(Angle);
        double c = cos(Angle);
        double t = 1.0 - c;
        double norm = Axis.length();
        if(norm == 0)
        {
            std::cerr << "Axis length must be greater than zero"
                      << "EulerAngles::set(axis,angle)"
                      << "Terminating!\n";
            OP_Exit(EXIT_FAILURE);
        }
        else
        {
            Axis.normalize();
        }
        // north pole singularity detected
        if ((Axis[0]*Axis[1]*t + Axis[2]*s) > 0.998)
        {
            Q[0] = 2.0*atan2(Axis[0]*sin(Angle/2),cos(Angle/2));
            Q[1] = Pi/2.0;
            Q[2] = 0;
        }
        // south pole singularity detected
        if ((Axis[0]*Axis[1]*t + Axis[2]*s) < -0.998)
        {
            Q[0] = -2.0*atan2(Axis[0]*sin(Angle/2),cos(Angle/2));
            Q[1] = -Pi/2.0;
            Q[2] = 0;
        }
        else
        {
            Q[0] = atan2(Axis[1] * s - Axis[0] * Axis[2] * t , 1.0 - (Axis[1]*Axis[1] + Axis[2]*Axis[2]) * t);
            Q[1] = asin (Axis[0] * Axis[1] * t + Axis[2] * s);
            Q[2] = atan2(Axis[0] * s - Axis[1] * Axis[2] * t , 1.0 - (Axis[0]*Axis[0] + Axis[2]*Axis[2]) * t);
        }

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;
    }

    void add(const double q1, const double q2, const double q3)
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

        IsSet = true;
    };

    EulerAngles& operator=(const EulerAngles& rhs)
    {

        assert(((Convention == rhs.Convention) or (Convention == NNN)) and "Euler angle conventions do not coincide!");

        Convention = rhs.Convention;

        Q[0] = rhs.Q[0];
        Q[1] = rhs.Q[1];
        Q[2] = rhs.Q[2];

        SinQ[0] = sin(rhs.Q[0]);
        SinQ[1] = sin(rhs.Q[1]);
        SinQ[2] = sin(rhs.Q[2]);

        CosQ[0] = cos(rhs.Q[0]);
        CosQ[1] = cos(rhs.Q[1]);
        CosQ[2] = cos(rhs.Q[2]);

        IsSet = true;

        return *this;
    };

    EulerAngles operator+(const EulerAngles& rhs) const
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

        returnAng.IsSet = true;

        return returnAng;
    };

    EulerAngles operator-(const EulerAngles& rhs) const
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

        returnAng.IsSet = true;

        return returnAng;
    };

    EulerAngles& operator+=(const EulerAngles& rhs)
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

        IsSet = true;

        return *this;
    };

    EulerAngles& operator-=(const EulerAngles& rhs)
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

        IsSet = true;

        return *this;
    };

    EulerAngles get_degree(void) const
    {
        EulerAngles returnAng;

        returnAng.Convention = Convention;

        returnAng.SinQ[0] = sin(Q[0]);
        returnAng.SinQ[1] = sin(Q[1]);
        returnAng.SinQ[2] = sin(Q[2]);

        returnAng.CosQ[0] = cos(Q[0]);
        returnAng.CosQ[1] = cos(Q[1]);
        returnAng.CosQ[2] = cos(Q[2]);

        returnAng.Q[0] = Q[0]*(180.0/Pi);
        returnAng.Q[1] = Q[1]*(180.0/Pi);
        returnAng.Q[2] = Q[2]*(180.0/Pi);

        returnAng.IsSet = true;

        return returnAng;
    };

    std::string get_convention(void) const
    {
        std::string returnCon = EulerConventionS[Convention];
        return returnCon;
    };

    EulerAngles operator*(const double rhs) const
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
    };

    EulerAngles& operator*=(const double rhs)
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
    };

    dMatrix3x3 getRotationMatrix() const;
    Quaternion getQuaternion(const bool Active = true) const;
    void getAxisAngle(dVector3& Axis, double& Angle);

    std::string print(void) const
    {
        std::stringstream out;

        out << "[" << Q[0] << ", "
                   << Q[1] << ", "
                   << Q[2] << "]["
                   << EulerConventionS[Convention]
                   << "]";
       return out.str();
    };
    std::string print_degree(void) const
    {
        std::stringstream out;

        out << "[" << Q[0]* (180.0/Pi) << ", "
                   << Q[1]* (180.0/Pi) << ", "
                   << Q[2]* (180.0/Pi) << "]["
                   << EulerConventionS[Convention]
                   << "]";
       return out.str();
    };
    std::string print_entire(void) const
    {
        std::stringstream out;

        out << "Angle      [" << Q[0] << ", "
                              << Q[1] << ", "
                              << Q[2] << "]\n"
            << "Convention [" << EulerConventionS[Convention]
                              << "]\n"
            << "Sin        [" << SinQ[0] << ", "
                              << SinQ[1] << ", "
                              << SinQ[2] << "]\n"
            << "Cos        [" << CosQ[0] << ", "
                              << CosQ[1] << ", "
                              << CosQ[2] << "]\n"
            << "Is set:     " << IsSet << "\n";
       return out.str();
    };
 protected:
 private:
};

}// namespace openphase
#endif
