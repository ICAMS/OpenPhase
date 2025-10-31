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
 *                         Philipp Engels; Muhammad Adil Ali; Hesham Salama
 *
 */

#ifndef MATRIXVECTOROPERATORS_H
#define MATRIXVECTOROPERATORS_H

#include "Globals.h"
#include "dMatrix3x3.h"
#include "dMatrix6x6.h"
#include "dVector3.h"
#include "iVector3.h"
#include "vStrain.h"
#include "vStress.h"
#include "vStress.h"

namespace openphase
{

extern inline dVector3 operator*(const dMatrix3x3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
        tmp[i] += lhs(i,j)*rhs[j];
    }
    return tmp;
}

extern inline dVector3 operator*(const dMatrix3x3& lhs, const iVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
        tmp[i] += lhs(i,j)*rhs[j];
    }
    return tmp;
}

extern inline double operator*(const iVector3& lhs, const dVector3& rhs)
{
    double tmp = 0.0;
    for(int i = 0; i < 3; i++)
    {
        tmp += lhs[i]*rhs[i];
    }
    return tmp;
}

extern inline double operator*(const dVector3& lhs, const iVector3& rhs)
{
    double tmp = 0.0;
    for(int i = 0; i < 3; i++)
    {
        tmp += lhs[i]*rhs[i];
    }
    return tmp;
}

extern inline dVector3 operator+(const iVector3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]+rhs[i];
    }
    return tmp;
}

extern inline dVector3 operator+(const dVector3& lhs, const iVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]+rhs[i];
    }
    return tmp;
}

extern inline dVector3 operator-(const iVector3& lhs, const dVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]-rhs[i];
    }
    return tmp;
}

extern inline dVector3 operator-(const dVector3& lhs, const iVector3& rhs)
{
    dVector3 tmp;
    for(int i = 0; i < 3; i++)
    {
        tmp[i] = lhs[i]-rhs[i];
    }
    return tmp;
}

extern inline dMatrix6x6 outer(const dMatrix3x3& lhs, const dMatrix3x3& rhs)
{
    dMatrix6x6 tmp;

    tmp(0,0) += lhs(0,0)*rhs(0,0);
    tmp(0,5) += lhs(0,0)*rhs(0,1);
    tmp(0,4) += lhs(0,0)*rhs(0,2);
    tmp(0,5) += lhs(0,0)*rhs(1,0);
    tmp(0,1) += lhs(0,0)*rhs(1,1);
    tmp(0,3) += lhs(0,0)*rhs(1,2);
    tmp(0,4) += lhs(0,0)*rhs(2,0);
    tmp(0,3) += lhs(0,0)*rhs(2,1);
    tmp(0,2) += lhs(0,0)*rhs(2,2);
    tmp(5,0) += lhs(0,1)*rhs(0,0);
    tmp(5,5) += lhs(0,1)*rhs(0,1);
    tmp(5,4) += lhs(0,1)*rhs(0,2);
    tmp(5,5) += lhs(0,1)*rhs(1,0);
    tmp(5,1) += lhs(0,1)*rhs(1,1);
    tmp(5,3) += lhs(0,1)*rhs(1,2);
    tmp(5,4) += lhs(0,1)*rhs(2,0);
    tmp(5,3) += lhs(0,1)*rhs(2,1);
    tmp(5,2) += lhs(0,1)*rhs(2,2);
    tmp(4,0) += lhs(0,2)*rhs(0,0);
    tmp(4,5) += lhs(0,2)*rhs(0,1);
    tmp(4,4) += lhs(0,2)*rhs(0,2);
    tmp(4,5) += lhs(0,2)*rhs(1,0);
    tmp(4,1) += lhs(0,2)*rhs(1,1);
    tmp(4,3) += lhs(0,2)*rhs(1,2);
    tmp(4,4) += lhs(0,2)*rhs(2,0);
    tmp(4,3) += lhs(0,2)*rhs(2,1);
    tmp(4,2) += lhs(0,2)*rhs(2,2);
    tmp(5,0) += lhs(1,0)*rhs(0,0);
    tmp(5,5) += lhs(1,0)*rhs(0,1);
    tmp(5,4) += lhs(1,0)*rhs(0,2);
    tmp(5,5) += lhs(1,0)*rhs(1,0);
    tmp(5,1) += lhs(1,0)*rhs(1,1);
    tmp(5,3) += lhs(1,0)*rhs(1,2);
    tmp(5,4) += lhs(1,0)*rhs(2,0);
    tmp(5,3) += lhs(1,0)*rhs(2,1);
    tmp(5,2) += lhs(1,0)*rhs(2,2);
    tmp(1,0) += lhs(1,1)*rhs(0,0);
    tmp(1,5) += lhs(1,1)*rhs(0,1);
    tmp(1,4) += lhs(1,1)*rhs(0,2);
    tmp(1,5) += lhs(1,1)*rhs(1,0);
    tmp(1,1) += lhs(1,1)*rhs(1,1);
    tmp(1,3) += lhs(1,1)*rhs(1,2);
    tmp(1,4) += lhs(1,1)*rhs(2,0);
    tmp(1,3) += lhs(1,1)*rhs(2,1);
    tmp(1,2) += lhs(1,1)*rhs(2,2);
    tmp(3,0) += lhs(1,2)*rhs(0,0);
    tmp(3,5) += lhs(1,2)*rhs(0,1);
    tmp(3,4) += lhs(1,2)*rhs(0,2);
    tmp(3,5) += lhs(1,2)*rhs(1,0);
    tmp(3,1) += lhs(1,2)*rhs(1,1);
    tmp(3,3) += lhs(1,2)*rhs(1,2);
    tmp(3,4) += lhs(1,2)*rhs(2,0);
    tmp(3,3) += lhs(1,2)*rhs(2,1);
    tmp(3,2) += lhs(1,2)*rhs(2,2);
    tmp(4,0) += lhs(2,0)*rhs(0,0);
    tmp(4,5) += lhs(2,0)*rhs(0,1);
    tmp(4,4) += lhs(2,0)*rhs(0,2);
    tmp(4,5) += lhs(2,0)*rhs(1,0);
    tmp(4,1) += lhs(2,0)*rhs(1,1);
    tmp(4,3) += lhs(2,0)*rhs(1,2);
    tmp(4,4) += lhs(2,0)*rhs(2,0);
    tmp(4,3) += lhs(2,0)*rhs(2,1);
    tmp(4,2) += lhs(2,0)*rhs(2,2);
    tmp(3,0) += lhs(2,1)*rhs(0,0);
    tmp(3,5) += lhs(2,1)*rhs(0,1);
    tmp(3,4) += lhs(2,1)*rhs(0,2);
    tmp(3,5) += lhs(2,1)*rhs(1,0);
    tmp(3,1) += lhs(2,1)*rhs(1,1);
    tmp(3,3) += lhs(2,1)*rhs(1,2);
    tmp(3,4) += lhs(2,1)*rhs(2,0);
    tmp(3,3) += lhs(2,1)*rhs(2,1);
    tmp(3,2) += lhs(2,1)*rhs(2,2);
    tmp(2,0) += lhs(2,2)*rhs(0,0);
    tmp(2,5) += lhs(2,2)*rhs(0,1);
    tmp(2,4) += lhs(2,2)*rhs(0,2);
    tmp(2,5) += lhs(2,2)*rhs(1,0);
    tmp(2,1) += lhs(2,2)*rhs(1,1);
    tmp(2,3) += lhs(2,2)*rhs(1,2);
    tmp(2,4) += lhs(2,2)*rhs(2,0);
    tmp(2,3) += lhs(2,2)*rhs(2,1),
    tmp(2,2) += lhs(2,2)*rhs(2,2);

    return tmp;
}

extern inline vStrain VoigtStrain(const dMatrix3x3& locStrainTensor)
{
    vStrain locStrain;
    locStrain[0] = locStrainTensor(0,0);
    locStrain[1] = locStrainTensor(1,1);
    locStrain[2] = locStrainTensor(2,2);
    locStrain[3] = locStrainTensor(1,2)*2.0;
    locStrain[4] = locStrainTensor(0,2)*2.0;
    locStrain[5] = locStrainTensor(0,1)*2.0;
    return locStrain;
}

extern inline vStress VoigtStress(const dMatrix3x3& locStressTensor)
{
    vStress locStress;
    locStress[0] = locStressTensor(0,0);
    locStress[1] = locStressTensor(1,1);
    locStress[2] = locStressTensor(2,2);
    locStress[3] = locStressTensor(1,2);
    locStress[4] = locStressTensor(0,2);
    locStress[5] = locStressTensor(0,1);
    return locStress;
}

// Return columns of a matrix
extern inline std::vector<dVector3> Col(const dMatrix3x3& Mat)
{
    std::vector<dVector3> Col;
    for(int i = 0; i < 3; i++)
    {
        dVector3 tmp;
        for(int j = 0; j < 3; j++)
        {
            tmp[j] = Mat(j,i);
        }
        Col.push_back(tmp);
    }
    return Col;
}

extern inline vStrain operator*(const dMatrix6x6& locCompliance, const vStress& locStress)
{
    vStrain locStrain;
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        locStrain[i] += locCompliance(i,j)*locStress[j];
    }
    return locStrain;
}

extern inline vStress operator*(const dMatrix6x6& locStiffness, const vStrain& locStrain)
{
    vStress locStress;
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        locStress[i] += locStiffness(i,j)*locStrain[j];
    }
    return locStress;
}

extern inline double operator*(const vStrain& locStrain, const vStress& locStress)
{
    double locEnergy = locStrain[0] * locStress[0]
                     + locStrain[1] * locStress[1]
                     + locStrain[2] * locStress[2]
                     + locStrain[3] * locStress[3]
                     + locStrain[4] * locStress[4]
                     + locStrain[5] * locStress[5];
    return locEnergy;
}

extern inline double operator*(const vStress& locStress, const vStrain& locStrain)
{
    double locEnergy = locStrain[0] * locStress[0]
                     + locStrain[1] * locStress[1]
                     + locStrain[2] * locStress[2]
                     + locStrain[3] * locStress[3]
                     + locStrain[4] * locStress[4]
                     + locStrain[5] * locStress[5];
    return locEnergy;
}

extern inline void Eigensystem(dMatrix3x3 M, dMatrix3x3& Eigenvectors,
                                             dMatrix3x3& Eigenvalues,
                                             double accuracy = FLT_EPSILON)     ///< Calculates eigenvectors and eigenvalues (output as a diagonal matrix for convenience) of a real symmetric matrix M using Jacobi method
{
    // M MUST BE _SYMMETRIC_
    assert((std::fabs(M(0,1) - M(1,0)) <= std::numeric_limits<double>::epsilon() and
            std::fabs(M(0,2) - M(2,0)) <= std::numeric_limits<double>::epsilon() and
            std::fabs(M(1,2) - M(2,1)) <= std::numeric_limits<double>::epsilon()) and
            "Matrix is not symmetric in MatrixVectorOperators::Eigensystem()");

    // Set elements which are already below specified accuracy to zero:
    for(int x = 0; x < 3; x++)
    for(int y = 0; y < 3; y++)
    {
        if(fabs(M(x,y)) < accuracy)
        {
            M(x,y) = 0.0;
        }
    }

    Eigenvectors.set_to_unity();

    do
    {
        // Maximum off-diagonal element value:
        double maxValue = 0.0;
        // Givens rotation angle:
        double theta = 0.0;
        // Indices of the maximum off-diagonal element:
        int indI = -1;
        int indJ = -1;

        // Find maximum off-diagonal element:
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        if((i != j))
        {
            if(fabs(M(i,j)) > fabs(maxValue))
            {
                maxValue = M(i,j);
                indI = i;
                indJ = j;
            }
        }

        // Perform iteration if there are non-zero off-diagonal elements
        if(indI != -1 and indJ != -1)
        {
            // Set rotation angle for Givens rotation:
            if(M(indI,indI) == M(indJ,indJ))
            {
                if(M(indI,indJ) > 0) theta =  Pi/4.0;
                else                 theta = -Pi/4.0;
            }
            else
            {
                theta = 0.5 * atan((2 * maxValue) / (M(indI,indI) - M(indJ,indJ)));
            }

            // Set Givens rotation matrix:
            dMatrix3x3 S;
            S.set_to_unity();
            S(indI,indI) =  cos(theta);
            S(indJ,indJ) =  cos(theta);
            S(indI,indJ) = -sin(theta);
            S(indJ,indI) =  sin(theta);

            // Compute update of the eigenvalue matrix:
            M = S.transposed() * M * S;

            // Compute update of the eigenvector matrix:
            Eigenvectors = Eigenvectors * S;

            // Set elements which have reached specified accuracy to zero:
            for(int x = 0; x < 3; x++)
            for(int y = 0; y < 3; y++)
            {
                if(fabs(M(x,y)) < accuracy)
                {
                    M(x,y) = 0.0;
                }
            }
        }
    }
    while (M(0,1) != 0.0 or M(0,2) != 0.0 or M(1,2) != 0.0 or
           M(1,0) != 0.0 or M(2,0) != 0.0 or M(2,1) != 0.0);

    // Assign eigenvalues:
    Eigenvalues = M;
}

extern inline dMatrix3x3 sqrtM3x3(const dMatrix3x3& M,
                                  double accuracy = FLT_EPSILON)                ///< Calculates sqrt(M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    Eigensystem(M,Eigenvectors,Eigenvalues,accuracy);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = sqrt(Eigenvalues(i,i));
    }
    // sqrt(M) = Eigenvectors_column . sqrt(Eigenvalues) . Eigenvectors_row
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();
}

extern inline dMatrix3x3 qurtM3x3(const dMatrix3x3& M,
                                  double accuracy = FLT_EPSILON)                ///< Calculates fourth order (quartic) root of M
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    Eigensystem(M,Eigenvectors,Eigenvalues,accuracy);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = sqrt(sqrt(Eigenvalues(i,i)));
    }
    // qurt(M) = Eigenvectors_column . sqrt(Eigenvalues) . Eigenvectors_row
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();
}

extern inline dMatrix3x3 logM3x3(const dMatrix3x3& M,
                                 double accuracy = FLT_EPSILON)                 ///< Calculates log(M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    Eigensystem(M,Eigenvectors,Eigenvalues,accuracy);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = log(Eigenvalues(i,i));
    }
    // log(M) = Eigenvectors_column . log(Eigenvalues) . Eigenvectors_row
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();
}

extern inline dMatrix3x3 expM3x3(const dMatrix3x3& M,
                                 double accuracy = FLT_EPSILON)                 ///< Calculates exp(M)
{
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

    if(M.norm() < DBL_EPSILON)
    {
        return dMatrix3x3::UnitTensor();
    }

    Eigensystem(M,Eigenvectors,Eigenvalues,accuracy);

    for(int i = 0; i < 3; i++)
    {
        Eigenvalues(i,i) = exp(Eigenvalues(i,i));
    }
    // exp(M) = Eigenvectors_column . exp(Eigenvalues) . Eigenvectors_row
    return Eigenvectors * Eigenvalues * Eigenvectors.transposed();
}

extern inline void PolarDecomposition(const dMatrix3x3& D, dMatrix3x3& R, dMatrix3x3& U,
                                      double accuracy = FLT_EPSILON)            ///< Calculates polar decomposition of a positive semi-definite matrix D = RU (R - rotation, U - stretch)
{
    dMatrix3x3 M = D.transposed() * D;
    U = sqrtM3x3(M,accuracy);
    R = D*U.inverted();
}

extern inline dMatrix3x3 AlignBaseAxes(dMatrix3x3 M)                            /// Aligns X, Y and Z axes of tensor M with the base vectors.
{
    /*
     * This function uses Givens rotation algorithm to zero (1,0), (2,0) and
     * (2,1) elements of M.
    */
    double c = 0.0;
    double s = 0.0;
    double den = 0.0;
    dMatrix3x3 G;

    auto denom = [](double x, double y) {return sqrt(x*x + y*y);};

    // Zeroing element (1,0)
    if(M(1,0) > DBL_EPSILON or M(1,0) < -DBL_EPSILON)
    {
        den = denom (M(1,0),M(0,0));
        if(den > DBL_EPSILON)
        {
            c = M(0,0) / den;
            s = M(1,0) / den;
        }
        G(0,0) = c;  G(0,1) = s; G(0,2) = 0;
        G(1,0) = -s; G(1,1) = c; G(1,2) = 0;
        G(2,0) = 0;  G(2,1) = 0; G(2,2) = 1;
        M = G*M;
    }

    // Zeroing element (2,0)
    if(M(2,0) > DBL_EPSILON or M(2,0) < -DBL_EPSILON)
    {
        den = denom (M(2,0),M(0,0));
        if(den > DBL_EPSILON)
        {
            c = M(0,0) / den;
            s = M(2,0) / den;
        }
        G(0,0) = c;  G(0,1) = 0; G(0,2) = s;
        G(1,0) = 0;  G(1,1) = 1; G(1,2) = 0;
        G(2,0) = -s; G(2,1) = 0; G(2,2) = c;
        M = G*M;
    }

    // Zeroing element (2,1)
    if(M(2,1) > DBL_EPSILON or M(2,1) < -DBL_EPSILON)
    {
        den = denom (M(2,1),M(1,1));
        if(den > DBL_EPSILON)
        {
            c = M(1,1) / den;
            s = M(2,1) / den;
        }
        G(0,0) = 1; G(0,1) = 0;  G(0,2) = 0;
        G(1,0) = 0; G(1,1) = c;  G(1,2) = s;
        G(2,0) = 0; G(2,1) = -s; G(2,2) = c;
        M = G*M;
    }

    return M;
}

}// namespace openphase
#endif
