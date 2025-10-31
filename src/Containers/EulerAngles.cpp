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

#include "Includes.h"

namespace openphase
{

dMatrix3x3 EulerAngles::getRotationMatrix() const
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
        case XZX:
        {
            double XZX[3][3] =
            {{                 c2,            - c3*s2,              s2*s3},
             {              c1*s2,   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3},
             {              s1*s2,   c1*s3 + c2*c3*s1,   c1*c3 - c2*s1*s3}};

            memmove(result.data(), XZX, 9*sizeof(double));
            break;
        }
        case XYX:
        {
            double XYX[3][3] =
            {{                 c2,             s2*s3,              c3*s2},
             {              s1*s2,   c1*c3 - c2*s1*s3, - c1*s3 - c2*c3*s1},
             {            - c1*s2,   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3}};

            memmove(result.data(), XYX, 9*sizeof(double));
            break;
        }
        case YXY:
        {
            double YXY[3][3] =
            {{   c1*c3 - c2*s1*s3,              s1*s2,   c2*c3*s1 + c1*s3},
             {              s2*s3,                 c2,            - c3*s2},
             { - c3*s1 - c1*c2*s3,              c1*s2,   c1*c2*c3 - s1*s3}};
            memmove(result.data(), YXY, 9*sizeof(double));
            break;
        }
        case YZY:
        {
            double YZY[3][3] =
            {{   c1*c2*c3 - s1*s3,            - c1*s2,   c3*s1 + c1*c2*s3},
             {              c3*s2,                 c2,              s2*s3},
             {-(c2*c3*s1) - c1*s3,              s1*s2,   c1*c3 - c2*s1*s3}};
            memmove(result.data(), YZY, 9*sizeof(double));
            break;
        }
        case ZYZ:
        {
            double ZYZ[3][3] =
            {{   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3,             c1*s2},
             {   c2*c3*s1 + c1*s3,   c1*c3 - c2*s1*s3,             s1*s2},
             {            - c3*s2,              s2*s3,                c2}};
            memmove(result.data(), ZYZ, 9*sizeof(double));
            break;
        }
        case ZXZ:
        {
            double ZXZ[3][3] =
            {{   c1*c3 - c2*s1*s3, - c2*c3*s1 - c1*s3,             s1*s2},
             {   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3,           - c1*s2},
             {              s2*s3,              c3*s2,                c2}};
            memmove(result.data(), ZXZ, 9*sizeof(double));
            break;
        }

        // Tait-Bryan angles
        case XYZ:
        {
            double XYZ[3][3] =
            {{              c2*c3,            - c2*s3,                s2},
             {   c1*s3 + c3*s1*s2,   c1*c3 - s1*s2*s3,           - c2*s1},
             {   s1*s3 - c1*c3*s2,   c1*s2*s3 + c3*s1,             c1*c2}};
            memmove(result.data(), XYZ, 9*sizeof(double));
            break;
        }
        case ZYX:
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
void EulerAngles::set(Quaternion Quat, EulerConvention locConvention, const bool active)
{
    Quat.normalize();

    double q1 = Quat[0];
    double q2 = Quat[1];
    double q3 = Quat[2];
    double q4 = Quat[3];

    switch(locConvention)
    {
        case XYX:
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
        case XYZ:
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
        case XZX:
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
        case XZY:
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
        case YXY:
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
        case YXZ:
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
        case YZX:
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
        case YZY:
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
        case ZXY:
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
        case ZXZ:
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
        case ZYX:
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
        case ZYZ:
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

    IsSet = true;
}

/*Source: http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
 *        https://www.astro.rug.nl/software/kapteyn/_downloads/fa29752e4cd69adcfa2fc03b1c020f4e/attitude.pdf
 *        https://de.mathworks.com/help/fusion/ref/quaternion.html#mw_54b60363-8447-46b1-b461-b50aef77772d
 *
 *        Switch between Active (point rotation) and Passive (frame rotation) with a flag
 *        Default setting: Active rotation
 */
Quaternion EulerAngles::getQuaternion(const bool Active) const
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
        case XYX:
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
        case XYZ:
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
        case XZX:
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
        case XZY:
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
        case YXY:
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
        case YXZ:
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
        case YZX:
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
        case YZY:
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
        case ZXY:
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
        case ZXZ:
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
        case ZYX:
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
        case ZYZ:
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

void EulerAngles::getAxisAngle(dVector3& Axis, double& Angle)
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
