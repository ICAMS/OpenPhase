/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2022
 *
 *   Authors:    Raphael Schiedung
 *
 *   Copyright (c) 2009-2022 Interdisciplinary Centre for Advanced Materials
 *                 Simulation (ICAMS). Ruhr-Universitaet Bochum. Germany
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LOCALLBM_H
#define LOCALLBM_H

#include "BoundaryConditions.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "ConsoleOutput.h"
#include "Initializations.h"
#include "Macros.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Velocities.h"
#include <complex>

using namespace openphase;

class LocalLBM : public FlowSolverLBM
{
public:
    LocalLBM(){};
    LocalLBM(Settings& locSettings, double in_dt):
        FlowSolverLBM(locSettings, in_dt){};

    double CalculateContactAngle(const double sigmaSL, const double sigmaSV, const double sigmaLV) const;
    double CalculateMass(void) const;                                          ///< Calculates the total mass
    double CalculateRadius(const Settings& locSettings,
            const int j, const int k) const;                                   ///< Calculates the diameter of the capillary bridge
    double CalculateSigmaX(const int j, const int k) const;                    ///< Calculates the surface tension
    double CalculateSigmaY(const int j0, const int jmax, const int k) const;   ///< Calculates the surface tension
    std::array<double,2> CalculateFluidVolumes(void) const;                    ///< Calculates the volume of liquid

    void EnforceVolume(BoundaryConditions& BC, const double& lbVCurrent, const double& lbVGoal, const double precision = 1E-8);
    void EnforceLiquidDiameter(BoundaryConditions& BC, const double& RCurrent, const double& RGoal, const double precision = 1E-8);
    void EnforceCurvature(double& Wetting, const double& kappa, const double& kappaGoal, const double precision = 1E-8) const;

    dVector3 CalculateNormal(const int i, const int j, const int k) const;
    std::array<double,2> PrincipleCurvatures(const int i, const int j, const int k) const;
    double AverageCurvature(const int j0, const int jmax, const int k0,
            const int kmax) const;

};
double LocalLBM::CalculateRadius(const Settings& locSettings, const int j, const int k) const
{
    double result = 0.0;
    const double lbRhoLG = VaporDensity[0]/dRho + (LiquidDensity[0]/dRho - VaporDensity[0]/dRho)/2;
    // Calculate diameter of liquid phase across the x-axis
    for (long int i =  0; i < Grid.Nx; ++i)
    {
        const double lbRho0 = DensityWetting(i  ,j,k,{0})/dRho;
        const double lbRho1 = DensityWetting(i+1,j,k,{0})/dRho;
        if (lbRho0 > lbRhoLG and lbRho1 < lbRhoLG)
        {
            // solve lbRhoLG = aa * xx + bb
            const double bb = lbRho0;
            const double aa = (lbRho1 - lbRho0);
            const double xx = i  + (lbRhoLG - bb)/aa;
            result = (locSettings.Grid.OffsetX + xx)*Grid.dx;
            break;
        }
    }
    #ifdef MPI_PARALLEL
    double tmp = result;
    OP_MPI_Allreduce(&tmp, &(result), 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    #endif
    return result;
}
double LocalLBM::CalculateSigmaX(const int j, const int k) const
{
    // Calculate difference between stress tensor components
    double sigma = 0;
    for (int i = 0; i < Grid.Nx; ++i)
    {
        const dMatrix3x3 locP = PressureTensor(i,j,k);
        sigma += (locP(0,0) - locP(1,1));
    }
    #ifdef MPI_PARALLEL
    double tmp = sigma;
    OP_MPI_Allreduce(&tmp, &(sigma), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif
    sigma *= Grid.dx;
    return sigma;
}
double LocalLBM::CalculateSigmaY(const int j0, const int jmax, const int k) const
{
    // Calculate difference between stress tensor components
    double sigma = 0;
    for (int j = j0; j < std::min(jmax,Grid.Ny-1); ++j)
    {
        const dMatrix3x3 locP = PressureTensor(0,j,k);
        sigma += (locP(1,1) - locP(0,0));
    }
    sigma *= Grid.dx;
    return sigma;
}
double LocalLBM::CalculateContactAngle(const double sigmaSL, const double sigmaSV, const double sigmaLV) const
{
    if (sigmaLV > DBL_EPSILON)
    {
        double theta = 0;
        const double cosTheta = (sigmaSV - sigmaSL)/sigmaLV;
        if      (cosTheta < - 1) theta = 180;
        else if (cosTheta >   1) theta =   0;
        else                     theta = 180/Pi*acos(cosTheta);
        return theta;

    }
    else return 0;
}
std::array<double,2> LocalLBM::CalculateFluidVolumes(void) const
{
    double VolumeLiquid = 0;
    double VolumeVapor  = 0;

    STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0)
    if (not Obstacle(i,j,k))
    {
        double DistLiquid = LiquidDensity[0]/dRho - DensityWetting(i,j,k,{0})/dRho;
        double DistVapor  = DensityWetting(i,j,k,{0})/dRho - VaporDensity[0]/dRho;

        if (DistLiquid < 0) DistLiquid = 0;
        if (DistVapor  < 0)  DistVapor = 0;

        if (DistLiquid < DistVapor)
        {
            VolumeLiquid += (1.0 - DistLiquid/LiquidDensity[0]*dRho);
        }
        else
        {
            VolumeVapor  += (1.0 - DistVapor/VaporDensity[0]*dRho);
        }
    }
    STORAGE_LOOP_END

    return std::array<double,2>({VolumeLiquid, VolumeVapor});
}
double LocalLBM::CalculateMass(void) const
{
    double lbMass = 0;
    STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0)
    if (not Obstacle(i,j,k))
    {
        lbMass += DensityWetting(i,j,k,{0})/dRho;
    }
    STORAGE_LOOP_END
    return lbMass;
}
void LocalLBM::EnforceVolume(BoundaryConditions& BC, const double& lbVCurrent, const double& lbVGoal, const double precision)
{
    //if (std::abs((RCurrent-RGoal)/dx) > precision)
    {
        const double lbDeltaV       = lbVGoal-lbVCurrent;
        const double lbDeltaM       = lbDeltaV*(LiquidDensity[0]/dRho-VaporDensity[0]/dRho);
        const double lbDeltaDensity = lbDeltaM/CountFluidNodes();

        const double prefactor = 10;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
        if (not Obstacle(i,j,k))
        {
            lbPopulations(i,j,k,{0})(0,0,0) += prefactor*lbDeltaDensity;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        SetBoundaryConditions(BC);
    }
}
void LocalLBM::EnforceLiquidDiameter(BoundaryConditions& BC, const double& RCurrent, const double& RGoal, const double precision)
{
    const double lbH            = 2*RGoal/Grid.dx;
    const double lbV            = 0.25*Pi*(RGoal*RGoal-RCurrent*RCurrent)/Grid.dx/Grid.dx*lbH;
    const double lbDeltaM       = lbV*(LiquidDensity[0]/dRho-VaporDensity[0]/dRho);
    const double lbDeltaDensity = lbDeltaM/CountFluidNodes();

    const double prefactor = 10/RGoal/RGoal*Grid.dx*Grid.dx;

    // Fix density
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DensityWetting,0,)
    if (not Obstacle(i,j,k))
    {
        lbPopulations(i,j,k,{0})(0,0,0) += prefactor*lbDeltaDensity;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}
void LocalLBM::EnforceCurvature(double& Wetting, const double& kappa, const double& kappaGoal, const double precision) const
{
    const double prefactor = 20*kappaGoal*kappaGoal/Grid.dx/Grid.dx;

    //if (std::abs((kappa-kappaGoal)*dx) > precision)
    {
        Wetting -= prefactor*(kappaGoal-kappa)*Grid.dx*dRho;
    }
}
dVector3 LocalLBM::CalculateNormal(const int i, const int j, const int k) const
{
    dVector3 normal;
    dVector3 grad;
    normal.set_to_zero();
    grad.set_to_zero();

    const double GWeights[3] = {-0.5/Grid.dx, 0.0, 0.5/Grid.dx};
    for (int ii = -Grid.dNx; ii <= Grid.dNx; ii += 2) grad[0] += DensityWetting(i+ii,j,k,{0})/dRho*GWeights[ii+1];
    for (int ii = -Grid.dNy; ii <= Grid.dNy; ii += 2) grad[1] += DensityWetting(i,j+ii,k,{0})/dRho*GWeights[ii+1];
    for (int ii = -Grid.dNz; ii <= Grid.dNz; ii += 2) grad[2] += DensityWetting(i,j,k+ii,{0})/dRho*GWeights[ii+1];

    const double NormGrad = sqrt(grad*grad);

    if (NormGrad > 1E-10) normal = grad /NormGrad;

    return normal;
}
std::array<double,2> LocalLBM::PrincipleCurvatures(const int i, const int j, const int k) const
{
    using namespace std;

    double kappa1 = 0;
    double kappa2 = 0;

    const double GWeights[3] = {-0.5/Grid.dx, 0.0, 0.5/Grid.dx};

    // Check if neighbour cell are also in the interface
    // Calculate gradients of phase normal fields (Jacobian matrix)
    Tensor<double, 2 > Mat;
    Mat.Allocate({3,3});
    Mat.set_to_zero();

    // Calculate gradients of normals
    for (int ii = -1; ii <= +1; ii += 2)
    {
        int jj = ii;
        int kk = ii;

        if      (ii > Grid.dNx) ii =  Grid.dNx;
        else if (ii < Grid.dNx) ii = -Grid.dNx;

        if      (jj > Grid.dNy) jj =  Grid.dNy;
        else if (jj < Grid.dNy) jj = -Grid.dNy;

        if      (kk > Grid.dNz) kk =  Grid.dNz;
        else if (kk < Grid.dNz) kk = -Grid.dNz;

        const double GWeight = GWeights[ii+1];

        const dVector3 locNormalX = CalculateNormal(i+ii,j,k);
        const dVector3 locNormalY = CalculateNormal(i,j+jj,k);
        const dVector3 locNormalZ = CalculateNormal(i,j,k+kk);

        Mat({0,0}) += GWeight*(locNormalX[0]);
        Mat({1,0}) += GWeight*(locNormalX[1]);
        Mat({2,0}) += GWeight*(locNormalX[2]);
        Mat({0,1}) += GWeight*(locNormalY[0]);
        Mat({1,1}) += GWeight*(locNormalY[1]);
        Mat({2,1}) += GWeight*(locNormalY[2]);
        Mat({0,2}) += GWeight*(locNormalZ[0]);
        Mat({1,2}) += GWeight*(locNormalZ[1]);
        Mat({2,2}) += GWeight*(locNormalZ[2]);
    }

    // Calculate basis of tangent space at (i,j,k)
    // the basis of a spherical coordinate system at (i,j,k) will be
    // used here. The phase normal vector is in this case equal to the
    // radial normal vector of the spherical coordinate and the
    // remaining basis vectors will be the basis of the tangent space

    dVector3 locN = CalculateNormal(i,j,k);
    const double lbRho  = DensityWetting(i,j,k,{0})/dRho;
    const double lbDRho = std::abs(LiquidDensity[0]/dRho - VaporDensity[0]/dRho);
    const double lbRho0 = VaporDensity[0]/dRho + lbDRho/2;

    if ( std::abs(lbRho-lbRho0) < 0.35 * lbDRho)
    {
        const double phi   = atan2(locN[1],locN[0]);
        const double theta = acos(locN[2]);

        // Calculate basis vectors of tangent space
        dVector3 e_n;
        e_n.set_to_zero();
        e_n[0] = locN[0];
        e_n[1] = locN[1];
        e_n[2] = locN[2];

        dVector3 e_theta;
        e_theta.set_to_zero();
        e_theta[0] = cos(theta) * cos(phi);
        e_theta[1] = cos(theta) * sin(phi);
        e_theta[2] =            - sin(theta);

        dVector3 e_phi;
        e_phi.set_to_zero();
        e_phi[0] = - sin(phi);
        e_phi[1] =   cos(phi);

        // Calculate projection and inclusion matrices
        double Projection [2][3];
        double Inclusion  [3][2];
        for (int m = 0; m < 3; ++m)
        {
            Projection[1][m] = e_phi  [m];
            Projection[0][m] = e_theta[m];

            Inclusion[m][1] = Projection[1][m];
            Inclusion[m][0] = Projection[0][m];
        }

        // Calculate local Weingarten map W
        std::complex<double> W[2][2];
        for (unsigned int l = 0; l < 2; ++l)
        for (unsigned int m = 0; m < 2; ++m)
        {
            W[l][m] = 0.0;
            for (unsigned int n = 0; n < 3; ++n)
            for (unsigned int o = 0; o < 3; ++o)
            {
                W[l][m] -= Projection[l][n] * Mat({n,o})
                    * Inclusion[o][m];
            }
        }

        // Calculate eigenvalues of local Weingarten map
        const double part1 = real((W[0][0] + W[1][1])/2.);
        const std::complex<double> b1 = (W[0][0] - W[1][1]);
        const std::complex<double> b2 = b1*b1; /// NOTE !!! do NOT use pow(b1,2)!!
        const double part2 = std::real(0.5 + sqrt(b2 + 4.*W[0][1] * W[1][0]));;

        kappa1 = real(part1 - part2);
        kappa2 = real(part1 + part2);
    }

    return std::array<double,2>{{kappa1, kappa2}};
}
double LocalLBM::AverageCurvature(const int j0, const int jmax, const int k0, const int kmax)  const
{
    size_t ii = 0;
    double kappa = 0;

    for (int k = k0; k < kmax; ++k)
    for (int j = j0; j < jmax; ++j)
    for (int i = 0;  i < Grid.Nx; ++i)
    {
        std::array<double,2> pkappa = PrincipleCurvatures(i,j,k);
        const double locKappa = (pkappa[0] + pkappa[1]);
        if (locKappa != 0)
        {
            kappa += locKappa;
            ii ++;
        }
    }
    #ifdef MPI_PARALLEL
    double tmpkappa = kappa;
    size_t tmpii    = ii;
    OP_MPI_Allreduce(&tmpkappa, &(kappa), 1, OP_MPI_DOUBLE,        OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(&tmpii,    &(ii),    1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    return (ii > 0) ? kappa/ii : 0;
}
#endif
