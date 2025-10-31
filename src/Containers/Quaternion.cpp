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

#include "Containers/Quaternion.h"

namespace openphase
{

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

void Quaternion::setRotationMatrix()
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

dMatrix3x3 Quaternion::getRotationMatrix(const bool Active)
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
        result.setRotationMatrix();
    }
    return result;
}

}//namespace openphase
