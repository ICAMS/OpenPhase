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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Muhammad Adil Ali; Hesham Salama
 *
 */

#include "Tools.h"
#include "Crystallography.h"
#include "PhaseField.h"

namespace openphase
{

using namespace std;

std::vector<iVector3> Tools::findIndexPermutations(iVector3& locIndices)
{
    std::vector<iVector3> myReturn;
    // dropping the signs
    int x_val = std::abs(locIndices.get_x());
    int y_val = std::abs(locIndices.get_y());
    int z_val = std::abs(locIndices.get_z());

    // sorting indices
    if (x_val > y_val)
    {
        std::swap(x_val, y_val);
    }
    if (y_val > z_val)
    {
        std::swap(y_val, z_val);
    }
    if (x_val > y_val)
    {
        std::swap(x_val, y_val);
    }

    // generating permutations
    if(x_val == 0)
    {
        if(y_val == 0) // x == 0; y == 0;
        {
            if(z_val == 0) // x == 0; y == 0; z == 0; (1)
            {
                myReturn.push_back({ 0, 0, 0});
            }
            else // x == 0; y == 0; z != 0; (6)
            {
                myReturn.push_back({     0,     0, z_val});
                myReturn.push_back({     0, z_val,     0});
                myReturn.push_back({ z_val,     0,     0});

                myReturn.push_back({     0,     0,-z_val});
                myReturn.push_back({     0,-z_val,     0});
                myReturn.push_back({-z_val,     0,     0});
            }
        }
        else if(y_val == z_val) // x == 0; y == z; (12)
        {
            myReturn.push_back({     0, y_val, y_val});
            myReturn.push_back({ y_val,     0, y_val});
            myReturn.push_back({ y_val, y_val,     0});

            myReturn.push_back({     0, y_val,-y_val});
            myReturn.push_back({ y_val,     0,-y_val});
            myReturn.push_back({ y_val,-y_val,     0});

            myReturn.push_back({     0,-y_val, y_val});
            myReturn.push_back({-y_val,     0, y_val});
            myReturn.push_back({-y_val, y_val,     0});

            myReturn.push_back({     0,-y_val,-y_val});
            myReturn.push_back({-y_val,     0,-y_val});
            myReturn.push_back({-y_val,-y_val,     0});
        }
        else // x == 0; y != z; (24)
        {
            myReturn.push_back({     0, y_val, z_val});
            myReturn.push_back({ y_val,     0, z_val});
            myReturn.push_back({ y_val, z_val,     0});

            myReturn.push_back({     0, y_val,-z_val});
            myReturn.push_back({ y_val,     0,-z_val});
            myReturn.push_back({ y_val,-z_val,     0});

            myReturn.push_back({     0,-y_val, z_val});
            myReturn.push_back({-y_val,     0, z_val});
            myReturn.push_back({-y_val, z_val,     0});

            myReturn.push_back({     0,-y_val,-z_val});
            myReturn.push_back({-y_val,     0,-z_val});
            myReturn.push_back({-y_val,-z_val,     0});

            myReturn.push_back({     0, z_val, y_val});
            myReturn.push_back({ z_val,     0, y_val});
            myReturn.push_back({ z_val, y_val,     0});

            myReturn.push_back({     0, z_val,-y_val});
            myReturn.push_back({ z_val,     0,-y_val});
            myReturn.push_back({ z_val,-y_val,     0});

            myReturn.push_back({     0,-z_val, y_val});
            myReturn.push_back({-z_val,     0, y_val});
            myReturn.push_back({-z_val, y_val,     0});

            myReturn.push_back({     0,-z_val,-y_val});
            myReturn.push_back({-z_val,     0,-y_val});
            myReturn.push_back({-z_val,-y_val,     0});
        }
    }
    else if (x_val == y_val and y_val == z_val) // x == y == z; (8)
    {
        myReturn.push_back({ x_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, x_val});

        myReturn.push_back({ x_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-x_val});
        myReturn.push_back({-x_val,-x_val, x_val});

        myReturn.push_back({-x_val,-x_val,-x_val});
    }
    else if (x_val == y_val) // x == y != z; (24)
    {
        myReturn.push_back({ x_val, x_val, z_val});
        myReturn.push_back({ x_val, z_val, x_val});
        myReturn.push_back({ z_val, x_val, x_val});

        myReturn.push_back({ x_val, x_val,-z_val});
        myReturn.push_back({ x_val,-z_val, x_val});
        myReturn.push_back({-z_val, x_val, x_val});

        myReturn.push_back({-x_val, x_val, z_val});
        myReturn.push_back({-x_val, z_val, x_val});
        myReturn.push_back({ z_val,-x_val, x_val});

        myReturn.push_back({-x_val, x_val,-z_val});
        myReturn.push_back({-x_val,-z_val, x_val});
        myReturn.push_back({-z_val,-x_val, x_val});

        myReturn.push_back({ x_val,-x_val, z_val});
        myReturn.push_back({ x_val, z_val,-x_val});
        myReturn.push_back({ z_val, x_val,-x_val});

        myReturn.push_back({ x_val,-x_val,-z_val});
        myReturn.push_back({ x_val,-z_val,-x_val});
        myReturn.push_back({-z_val, x_val,-x_val});

        myReturn.push_back({-x_val,-x_val, z_val});
        myReturn.push_back({-x_val, z_val,-x_val});
        myReturn.push_back({ z_val,-x_val,-x_val});

        myReturn.push_back({-x_val,-x_val,-z_val});
        myReturn.push_back({-x_val,-z_val,-x_val});
        myReturn.push_back({-z_val,-x_val,-x_val});
    }
    else if (y_val == z_val) // x != y == z; (24)
    {
        myReturn.push_back({ x_val, y_val, y_val});
        myReturn.push_back({ y_val, x_val, y_val});
        myReturn.push_back({ y_val, y_val, x_val});

        myReturn.push_back({ x_val, y_val,-y_val});
        myReturn.push_back({ y_val, x_val,-y_val});
        myReturn.push_back({ y_val,-y_val, x_val});

        myReturn.push_back({ x_val,-y_val, y_val});
        myReturn.push_back({-y_val, x_val, y_val});
        myReturn.push_back({-y_val, y_val, x_val});

        myReturn.push_back({ x_val,-y_val,-y_val});
        myReturn.push_back({-y_val, x_val,-y_val});
        myReturn.push_back({-y_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, y_val});
        myReturn.push_back({ y_val,-x_val, y_val});
        myReturn.push_back({ y_val, y_val,-x_val});

        myReturn.push_back({-x_val, y_val,-y_val});
        myReturn.push_back({ y_val,-x_val,-y_val});
        myReturn.push_back({ y_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-y_val, y_val});
        myReturn.push_back({-y_val,-x_val, y_val});
        myReturn.push_back({-y_val, y_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-y_val});
        myReturn.push_back({-y_val,-x_val,-y_val});
        myReturn.push_back({-y_val,-y_val,-x_val});
    }
    else if (x_val == z_val) // x == z; x != y; (24)
    {
        myReturn.push_back({ x_val, y_val, x_val});
        myReturn.push_back({ y_val, x_val, x_val});
        myReturn.push_back({ x_val, x_val, y_val});

        myReturn.push_back({ x_val, y_val,-x_val});
        myReturn.push_back({ y_val,-x_val, x_val});
        myReturn.push_back({-x_val, x_val, y_val});

        myReturn.push_back({-x_val, y_val, x_val});
        myReturn.push_back({-y_val, x_val, x_val});
        myReturn.push_back({ x_val,-x_val, y_val});

        myReturn.push_back({-x_val, y_val,-x_val});
        myReturn.push_back({-y_val,-x_val, x_val});
        myReturn.push_back({-x_val,-x_val, y_val});

        myReturn.push_back({ x_val,-y_val, x_val});
        myReturn.push_back({ y_val, x_val,-x_val});
        myReturn.push_back({ x_val, x_val,-y_val});

        myReturn.push_back({ x_val,-y_val,-x_val});
        myReturn.push_back({ y_val,-x_val,-x_val});
        myReturn.push_back({-x_val, x_val,-y_val});

        myReturn.push_back({-x_val,-y_val, x_val});
        myReturn.push_back({-y_val, x_val,-x_val});
        myReturn.push_back({ x_val,-x_val,-y_val});

        myReturn.push_back({-x_val,-y_val,-x_val});
        myReturn.push_back({-y_val,-x_val,-x_val});
        myReturn.push_back({-x_val,-x_val,-y_val});
    }
    else // x != y != z != x; (48)
    {
        myReturn.push_back({ x_val, y_val, z_val});
        myReturn.push_back({ y_val, x_val, z_val});
        myReturn.push_back({ y_val, z_val, x_val});

        myReturn.push_back({ x_val, y_val,-z_val});
        myReturn.push_back({ y_val, x_val,-z_val});
        myReturn.push_back({ y_val,-z_val, x_val});

        myReturn.push_back({ x_val,-y_val, z_val});
        myReturn.push_back({-y_val, x_val, z_val});
        myReturn.push_back({-y_val, z_val, x_val});

        myReturn.push_back({ x_val,-y_val,-z_val});
        myReturn.push_back({-y_val, x_val,-z_val});
        myReturn.push_back({-y_val,-z_val, x_val});

        myReturn.push_back({ x_val, z_val, y_val});
        myReturn.push_back({ z_val, x_val, y_val});
        myReturn.push_back({ z_val, y_val, x_val});

        myReturn.push_back({ x_val, z_val,-y_val});
        myReturn.push_back({ z_val, x_val,-y_val});
        myReturn.push_back({ z_val,-y_val, x_val});

        myReturn.push_back({ x_val,-z_val, y_val});
        myReturn.push_back({-z_val, x_val, y_val});
        myReturn.push_back({-z_val, y_val, x_val});

        myReturn.push_back({ x_val,-z_val,-y_val});
        myReturn.push_back({-z_val, x_val,-y_val});
        myReturn.push_back({-z_val,-y_val, x_val});

        myReturn.push_back({-x_val, y_val, z_val});
        myReturn.push_back({ y_val,-x_val, z_val});
        myReturn.push_back({ y_val, z_val,-x_val});

        myReturn.push_back({-x_val, y_val,-z_val});
        myReturn.push_back({ y_val,-x_val,-z_val});
        myReturn.push_back({ y_val,-z_val,-x_val});

        myReturn.push_back({-x_val,-y_val, z_val});
        myReturn.push_back({-y_val,-x_val, z_val});
        myReturn.push_back({-y_val, z_val,-x_val});

        myReturn.push_back({-x_val,-y_val,-z_val});
        myReturn.push_back({-y_val,-x_val,-z_val});
        myReturn.push_back({-y_val,-z_val,-x_val});

        myReturn.push_back({-x_val, z_val, y_val});
        myReturn.push_back({ z_val,-x_val, y_val});
        myReturn.push_back({ z_val, y_val,-x_val});

        myReturn.push_back({-x_val, z_val,-y_val});
        myReturn.push_back({ z_val,-x_val,-y_val});
        myReturn.push_back({ z_val,-y_val,-x_val});

        myReturn.push_back({-x_val,-z_val, y_val});
        myReturn.push_back({-z_val,-x_val, y_val});
        myReturn.push_back({-z_val, y_val,-x_val});

        myReturn.push_back({-x_val,-z_val,-y_val});
        myReturn.push_back({-z_val,-x_val,-y_val});
        myReturn.push_back({-z_val,-y_val,-x_val});
    }
    return myReturn;
}

////////// extract Rotation from total deformation //////////
/// Implemented from https://dl.acm.org/doi/abs/10.1145/1028523.1028541
void Tools::jacobiRotate(dMatrix3x3 &A, dMatrix3x3 &R, int p, int q)
{
    if (A(p, q) == 0.0)
    {
        return;
    }
    double d = (A(p, p) - A(q, q)) / (2.0*A(p, q));
    double t = 1.0 / (fabs(d) + sqrt(d*d + 1.0));
    if (d < 0.0)
    {
        t = -t;
    }
    double c = 1.0 / sqrt(t*t + 1);
    double s = t*c;
    A(p, p) += t*A(p, q);
    A(q, q) -= t*A(p, q);
    A(p, q) = A(q, p) = 0.0;
    int k;
    for (k = 0; k < 3; k++)
    {
        if (k != p && k != q)
        {
            double Akp = c*A(k, p) + s*A(k, q);
            double Akq = -s*A(k, p) + c*A(k, q);
            A(k, p) = A(p, k) = Akp;
            A(k, q) = A(q, k) = Akq;
        }
    }
    for (k = 0; k < 3; k++)
    {
        double Rkp = c*R(k, p) + s*R(k, q);
        double Rkq = -s*R(k, p) + c*R(k, q);
        R(k, p) = Rkp;
        R(k, q) = Rkq;
    }
}
void Tools::eigenDecomposition(const dMatrix3x3 &A, dMatrix3x3 &EigenVectors, dVector3 &EigenValues)
{
    const int numJacobiIterations = 10;
    dMatrix3x3 D = A;
    EigenVectors.set_to_unity();
    int iter = 0;
    while (iter < numJacobiIterations)
    {
        int p, q;
        double a, max;
        max = fabs(D(0, 1));
        p = 0;
        q = 1;
        a = fabs(D(0, 2));
        if (a > max)
        {
            p = 0;
            q = 2;
            max = a;
        }
        a = fabs(D(1, 2));
        if (a > max)
        {
            p = 1;
            q = 2;
            max = a;
        }
        if (max < DBL_EPSILON)
        {
            break;
        }
        jacobiRotate(D, EigenVectors, p, q);
        iter++;
    }
    EigenValues[0] = D(0, 0);
    EigenValues[1] = D(1, 1);
    EigenValues[2] = D(2, 2);
}
void Tools::rotationMatrixIrving(const dMatrix3x3 &A, dMatrix3x3 &R)
{
    dMatrix3x3 AT_A, V;
    AT_A = A.transposed() * A;
    dVector3 S;
    eigenDecomposition(AT_A, V, S);
    const double detV = V.determinant();
    if (detV < 0.0)
    {
        double minLambda = DBL_MAX;
        unsigned char pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (S[l] < minLambda)
            {
                pos = l;
                minLambda = S[l];
            }
        }
        V(0, pos) = -V(0, pos);
        V(1, pos) = -V(1, pos);
        V(2, pos) = -V(2, pos);
    }
    if (S[0] < 0.0f)
        S[0] = 0.0f;
    if (S[1] < 0.0f)
        S[1] = 0.0f;
    if (S[2] < 0.0f)
        S[2] = 0.0f;

    dVector3 sigma;
    sigma[0] = sqrt(S[0]);
    sigma[1] = sqrt(S[1]);
    sigma[2] = sqrt(S[2]);
    unsigned char chk = 0;
    unsigned char pos = 0;
    dMatrix3x3 U;
    for (unsigned char l = 0; l < 3; l++)
    {
        if (fabs(sigma[l]) < 1.0e-4)
        {
            pos = l;
            chk++;
        }
    }
    if (chk > 0)
    {
        if (chk > 1)
        {
            U.set_to_unity();
        }
        else
        {
            U = A * V;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    for (unsigned char m = 0; m < 3; m++)
                    {
                        U(m, l) *= 1.0f / sigma[l];
                    }
                }
            }
            dVector3 v[2];
            unsigned char index = 0;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    v[index++] = dVector3({U(0, l), U(1, l), U(2, l)});
                }
            }
            dVector3 vec = v[0].cross(v[1]);
            vec.normalized();
            U(0, pos) = vec[0];
            U(1, pos) = vec[1];
            U(2, pos) = vec[2];
        }
    }
    else
    {
        dVector3 sigmaInv;
        sigmaInv = {1.0 / sigma[0], 1.0 / sigma[1], 1.0 / sigma[2]};
        U = A * V;
        for (unsigned char l = 0; l < 3; l++)
        {
            for (unsigned char m = 0; m < 3; m++)
            {
                U(m, l) *= sigmaInv[l];
            }
        }
    }
    const double detU = U.determinant();
    if (detU < 0.0)
    {
        double minLambda = DBL_MAX;
        unsigned char pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (sigma[l] < minLambda)
            {
                pos = l;
                minLambda = sigma[l];
            }
        }
        sigma[pos] = -sigma[pos];
        U(0, pos) = -U(0, pos);
        U(1, pos) = -U(1, pos);
        U(2, pos) = -U(2, pos);
    }
    R = U * V.transpose();
}
void Tools::GetAxisAngleFromRotationMatrix(const dMatrix3x3 RotMatrix, double& Angle, dVector3& Axis)
{
    ///  Implemented from https://www.euclideanspace.com/maths/geometry/rotations/index.htm
    double epsilon = 0.01; // margin to allow for rounding errors
    double epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees
    if ((fabs(RotMatrix(0,1)-RotMatrix(1,0)) < epsilon)
     && (fabs(RotMatrix(0,2)-RotMatrix(2,0)) < epsilon)
     && (fabs(RotMatrix(1,2)-RotMatrix(2,1)) < epsilon))
    {
        // singularity found
        // first check for identity matrix which must have +1 for all terms
        // in leading diagonal and zero in other terms
        if ((fabs(RotMatrix(0,1)+RotMatrix(1,0)) < epsilon2)
         && (fabs(RotMatrix(0,2)+RotMatrix(2,0)) < epsilon2)
         && (fabs(RotMatrix(1,2)+RotMatrix(2,1)) < epsilon2)
         && (fabs(RotMatrix(0,0)+RotMatrix(1,1)+RotMatrix(2,2)-3.0) < epsilon2))
        {
            // this singularity is identity matrix so angle = 0
            Angle = 0.0;
            Axis = {1,0,0};
            // zero angle, arbitrary axis
        }
        // otherwise this singularity is angle = 180
        Angle = Pi;
        double xx = (RotMatrix(0,0)+1)/2;
        double yy = (RotMatrix(1,1)+1)/2;
        double zz = (RotMatrix(2,2)+1)/2;
        double xy = (RotMatrix(0,1)+RotMatrix(1,0))/4;
        double xz = (RotMatrix(0,2)+RotMatrix(2,0))/4;
        double yz = (RotMatrix(1,2)+RotMatrix(2,1))/4;
        if ((xx > yy) && (xx > zz))
        {
            // RotMatrix(0,0) is the largest diagonal term
            if (xx< epsilon)
            {
                Axis = {0,0.7071,0.7071};
            }
            else
            {
                double x = sqrt(xx);
                Axis = {x,xy/x,xz/x};
            }
        }
        else if(yy > zz)
        {
            // m[1][1] is the largest diagonal term
            if (yy< epsilon)
            {
                Axis = {0.7071,0,0.7071};
            }
            else
            {
                double y = sqrt(yy);
                Axis = {xy/y,y,yz/y};
            }
        }
        else
        {
            // m[2][2] is the largest diagonal term so base result on this
            if (zz< epsilon)
            {
                Axis = {0.7071,0.7071,0};
            }
            else
            {
                double z = sqrt(zz);
                Axis = {xz/z,yz/z,z};
            }
        }
    }
    // as we have reached here there are no singularities so we can handle normally
    // First, normalise the rotation matrix
    double Norm = pow((RotMatrix(1,2) - RotMatrix(2,1)),2) + pow((RotMatrix(2,0) - RotMatrix(0,2)),2) + pow((RotMatrix(0,1) - RotMatrix(1,0)),2);
    Norm =  sqrt(Norm);

    if (fabs(Norm) < 0.001)
    {
        Norm = 1.0;
        // prevent divide by zero, should not happen if matrix is orthogonal and should be
        // caught by singularity test above, but I've left it in just in case
    }
    Angle = acos(0.5*(RotMatrix.trace() - 1));
    Axis[0] = ( RotMatrix(2,1) - RotMatrix(1,2) ) / Norm;
    Axis[1] = ( RotMatrix(0,2) - RotMatrix(2,0) ) / Norm;
    Axis[2] = ( RotMatrix(1,0) - RotMatrix(0,1) ) / Norm;
}
void Tools::getAxisAngle(const dMatrix3x3& TransformationMatrix, dVector3& Axis, double &Angle )
{
    dMatrix3x3 RotationMatrix;
    rotationMatrixIrving(TransformationMatrix,RotationMatrix);
    GetAxisAngleFromRotationMatrix(RotationMatrix,Angle,Axis);
}

double Tools::getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b
    return acos(((RotMatB*RotMatA.transposed()).trace()-1.0)/2.0);
}

double Tools::getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b

    // Only for cubic system. Needs to be generalized!

    double misorientation = 10.0;
    for (int i = 0; i < 24; i++)
    {
        double misorientationloc = acos(((CR.SymmetriesCubic[i]*RotMatB*
                RotMatA.transposed()).trace()-1.0)/2.0);
        if (abs(misorientationloc) <= abs(misorientation))
        {
            misorientation = misorientationloc;
        }
    }
    return misorientation;
}

double Tools::getDisorientationCubic(const Quaternion OrientationA, const Quaternion OrientationB)
{
    // Taken from Grimmer 1974
    // Only for cubic system. Needs to be generalized!

    double Angle[24];
    double misorientation = 0.0;
    Quaternion deltaQ;
    deltaQ  = OrientationA* OrientationB.inverted();

    Angle[0]  = deltaQ[0];
    Angle[1]  = deltaQ[1];
    Angle[2]  = deltaQ[2];
    Angle[3]  = deltaQ[3];

    Angle[4]  = (deltaQ[0] + deltaQ[1])*0.707;
    Angle[5]  = (deltaQ[0] - deltaQ[1])*0.707;
    Angle[6]  = (deltaQ[2] + deltaQ[3])*0.707;
    Angle[7]  = (deltaQ[2] - deltaQ[3])*0.707;

    Angle[8]  = (deltaQ[0] + deltaQ[2])*0.707;
    Angle[9]  = (deltaQ[0] - deltaQ[2])*0.707;
    Angle[10] = (deltaQ[1] + deltaQ[3])*0.707;
    Angle[11] = (deltaQ[1] - deltaQ[3])*0.707;

    Angle[12] = (deltaQ[0] + deltaQ[3])*0.707;
    Angle[13] = (deltaQ[0] - deltaQ[3])*0.707;
    Angle[14] = (deltaQ[1] + deltaQ[2])*0.707;
    Angle[15] = (deltaQ[1] - deltaQ[2])*0.707;

    Angle[16] = (deltaQ[0] + deltaQ[1] + deltaQ[2] + deltaQ[3])*0.5;
    Angle[17] = (deltaQ[0] + deltaQ[1] - deltaQ[2] - deltaQ[3])*0.5;
    Angle[18] = (deltaQ[0] - deltaQ[1] + deltaQ[2] - deltaQ[3])*0.5;
    Angle[19] = (deltaQ[0] - deltaQ[1] - deltaQ[2] + deltaQ[3])*0.5;

    Angle[20] = (deltaQ[0] + deltaQ[1] + deltaQ[2] - deltaQ[3])*0.5;
    Angle[21] = (deltaQ[0] + deltaQ[1] - deltaQ[2] + deltaQ[3])*0.5;
    Angle[22] = (deltaQ[0] - deltaQ[1] + deltaQ[2] + deltaQ[3])*0.5;
    Angle[23] = (deltaQ[0] - deltaQ[1] - deltaQ[2] - deltaQ[3])*0.5;

    for (int a=0; a<24; a++)
    {
        Angle[a] = abs(Angle[a]);
    }
    double DisAngle = 0;
    for (int a=0; a<24; a++)
    {
        if (DisAngle < Angle[a])
            DisAngle =  Angle[a];
    }
    misorientation = 2*acos(DisAngle);
    return misorientation;
}

EulerAngles Tools::RotationToEuler(const dMatrix3x3& Rot, const EulerConvention EConvention )
{
    EulerAngles Euler;
    Euler.Convention = EConvention;
    dMatrix3x3 RT = Rot.transposed();

    if(EConvention == ZXZ)
    {
        // Euler ZXZ passive (following formulation by Martin Boeff)
        double squvw = sqrt(RT(0,0)*RT(0,0)
            + RT(1,0)*RT(1,0) + RT(2,0)*RT(2,0));
        double sqhk = sqrt(RT(0,2)*RT(0,2)
            + RT(1,2)*RT(1,2));
        double sqhkl = sqrt(RT(0,2)*RT(0,2)
            + RT(1,2)*RT(1,2) + RT(2,2)*RT(2,2));
        double tempval = RT(2,2)/sqhkl;

        if(tempval >  1.0) {tempval =  1.0;}
        if(tempval < -1.0) {tempval = -1.0;}
        Euler.Q[1] = acos(tempval);

        if(Euler.Q[1] < 1.0e-8)
        {
            // calculate phi2
            Euler.Q[2] = 0.0;
            // calculate phi1
            tempval = RT(0,0)/squvw;
            if(tempval >  1.0) {tempval =  1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[0] = acos(tempval);
            if(RT(1,0) > 0.0) {Euler.Q[0] = 2.0*Pi - Euler.Q[0];}
        }
        else
        {
            // calculate phi2
            tempval = RT(1,2)/sqhk;
            if(tempval >  1.0) {tempval =  1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[2] = acos(tempval);
            if(RT(0,2) < 0.0) {Euler.Q[2] = 2.0*Pi - Euler.Q[2];}
            // calculate phi1
            tempval = - RT(2,1)/Euler.SinQ[1];
            if(tempval >  1.0) {tempval = 1.0;}
            if(tempval < -1.0) {tempval = -1.0;}

            Euler.Q[0] = acos(tempval);
            if(RT(2,0) < 0.0) {Euler.Q[0] = 2.0*Pi - Euler.Q[0];}
        }
        Euler.setTrigonometricFunctions();
    }
    else
    {
        std::stringstream message;
        message<< "Wrong/Unknown/None Euler convention used: " << EConvention ;

        ConsoleOutput::WriteExit(message.str(), "Tools", "RotationToEulerAngles()");
        exit(13);
    }
    return Euler;
}

dVector3 Tools::MillerConversion(const std::vector<double>& hkil, dVector3 hkl, bool PlaneNormal)
{
    if(PlaneNormal)
    {
        hkl[0] = hkil[0];
        hkl[1] = hkil[1];
        hkl[2] = hkil[3];
    }
    else
    {
        hkl[0] = hkil[0]-hkil[2];
        hkl[1] = hkil[1]-hkil[2];
        hkl[2] = hkil[3];
    }

    // Titanium c/a ratio ( to be generalized depends on the element )
    double a = 0.2950;
    double c = 0.4683;

    dMatrix3x3 locMat {a, -a/2.0, 0,
                       0, a*sqrt(3)/2.0, 0,
                       0, 0, c};
    hkl = locMat * hkl;

    hkl.normalize();
    return hkl;
}

void Tools::SetRandomGrainOrientations(PhaseField& Phase, const int seed)
{
    // uniform sampling of quaternion parameters using the method presented in (Shoemake, 1992)
    // taken from :: http://planning.cs.uiuc.edu/node198.html#eqn:shoemake

    default_random_engine generator(seed);

    uniform_real_distribution <double> Q1Distribution(0, 1.0);
    uniform_real_distribution <double> Q2Distribution(0, 1.0);
    uniform_real_distribution <double> Q3Distribution(0, 1.0);

    uniform_real_distribution <double> A1Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A2Distribution(0.0, 2.0*Pi);
    uniform_real_distribution <double> A3Distribution(0.0, 2.0*Pi);

    for(size_t alpha = 0; alpha < Phase.FieldsProperties.size(); alpha++)
    if(Phase.FieldsProperties[alpha].Exist)
    {
        Quaternion tempQuat;
        switch(Phase.Grid.Active())
        {
            case 1: // 1D
            {
                Phase.FieldsProperties[alpha].Orientation.set(1.0, 0.0, 0.0, 0.0);
                break;
            }
            case 2: // 2D
            {
                double a1 = 0.0;
                double a2 = 0.0;
                double a3 = 0.0;
                if(Phase.Grid.dNx == 0) a1 = A1Distribution(generator);
                if(Phase.Grid.dNy == 0) a2 = A2Distribution(generator);
                if(Phase.Grid.dNz == 0) a3 = A3Distribution(generator);
                EulerAngles ph1({a1,a2,a3},XYZ);
                Phase.FieldsProperties[alpha].Orientation = ph1.getQuaternion().normalized();
                break;
            }
            case 3: // 3D Full Rotation
            {
                double u1 = Q1Distribution(generator);
                double u2 = Q2Distribution(generator);
                double u3 = Q3Distribution(generator);

                tempQuat.set(sqrt(1.0-u1)*sin(2.0*Pi*u2),
                             sqrt(1.0-u1)*cos(2.0*Pi*u2),
                             sqrt(u1)*sin(2.0*Pi*u3),
                             sqrt(u1)*cos(2.0*Pi*u3));

                Phase.FieldsProperties[alpha].Orientation = tempQuat.normalized();
                break;
            }
        }
    }
}

/*
 * A Robust Method to Extract the Rotational Part of Deformations
 * https://animation.rwth-aachen.de/media/papers/2016-MIG-StableRotation.pdf
 */

/* The quat can be :
 *
 * 1- The solution of the previous step of an iterative solve.
 * 2- If such a solution is not available, start with q = Quat(Mat)/|Quat(Mat)|.
 * 3- q=(1,0,0,0).
 *
 *
 * No of iterations: tested, maxIter: 20 for best results!
 *
 *
*/
void Tools::ExtractRotation(dMatrix3x3 &Mat, Quaternion& Quat, const size_t maxIter)
{
    for (size_t iter = 0; iter < maxIter; iter++)
    {
        dMatrix3x3 Rot = Quat.getRotationMatrix();

        std::vector<dVector3> Rot_col = Col(Rot);
        std::vector<dVector3> Mat_col = Col(Mat);

        double b = 1.0 / fabs(Rot_col[0]*Mat_col[0] +
                              Rot_col[1]*Mat_col[1] +
                              Rot_col[2]*Mat_col[2]) + 1.0e-9;

        dVector3 omega = (Rot_col[0].cross(Mat_col[0]) +
                          Rot_col[1].cross(Mat_col[1]) +
                          Rot_col[2].cross(Mat_col[2]))*b;

        double Angle = omega.length();

        if (Angle < 1.0e-9) break;

        dVector3 Axis = omega * (1.0 / Angle);

        Quaternion tempQuat;
        tempQuat.set(Axis, Angle);
        Quat = tempQuat*Quat;
        Quat.normalize();
    }
}

double Tools::getMemoryUsageMB()
{
    std::ifstream file("/proc/self/status");
    std::string line;
    while (getline(file, line))
    {
        if (line.substr(0, 6) == "VmRSS:")
        {
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stol(value) / 1024.0;  // Convert kB to MB
        }
    }
    return -1.0;
}

}// namespace openphase
