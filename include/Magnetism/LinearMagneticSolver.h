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

 *   File created :   2019
 *   Main contributors :   Marvin Tegeler; Raphael Schiedung
 *
 */

#ifndef MAGNETICSOLVER_H
#define MAGNETICSOLVER_H

#include "Includes.h"
#include "BoundaryConditions.h"

namespace openphase
{
class ElasticProperties;
class PhaseField;
class Settings;
class dVector3;

class OP_EXPORTS LinearMagneticSolver : public OPObject
{
public:
    using OPObject::ReadInput;

    LinearMagneticSolver(){};
    LinearMagneticSolver(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;
    void ReadInput(const std::string InputFileName) override;
    void ReadInput(std::stringstream& inp) override;
    void SetEffectiveSusceptibility(const PhaseField& Phase, const BoundaryConditions& BC);
    void Solve(const BoundaryConditions& BC, const double MaxResidual);
    void CalcForceDensity(ElasticProperties& EP, const BoundaryConditions& BC) const;
    void WriteVTK(const int tStep, const Settings& locSettings);

    Storage3D<double, 0> chi;                                                   ///< Magnetic susceptibility
    Storage3D<double, 0> phi;                                                   ///< Magnetic potential field

    double H0x;                                                                 ///< x-component of external H-Field
    double H0y;                                                                 ///< y-component of external H-Field
    double H0z;                                                                 ///< z-component of external H-Field
    double dH0x_dx;                                                             ///< x-component of external H-Field gradient in x-direction
    double dH0y_dx;                                                             ///< y-component of external H-Field gradient in x-direction
    double dH0z_dx;                                                             ///< z-component of external H-Field gradient in x-direction
    double dH0x_dy;                                                             ///< x-component of external H-Field gradient in y-direction
    double dH0y_dy;                                                             ///< y-component of external H-Field gradient in y-direction
    double dH0z_dy;                                                             ///< z-component of external H-Field gradient in y-direction
    double dH0x_dz;                                                             ///< x-component of external H-Field gradient in z-direction
    double dH0y_dz;                                                             ///< y-component of external H-Field gradient in z-direction
    double dH0z_dz;                                                             ///< z-component of external H-Field gradient in z-direction

    size_t Nphases;

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t MaxIterations;                                                       ///< Maximal residual value for iteration to stop

    std::vector<double> PhaseChi;
    std::vector<double> b,x;
    SparseMatrixCSR A;

    inline size_t idx(const int i, const int j, const int k) const
    {
        return i+Grid.Nx*(j+Grid.Ny*k);
    }

    /// Magnetic permeability
    inline double mu (const int i, const int j, const int k) const
    {
        return (chi(i,j,k) + 1)*PhysicalConstants::mu_0;
    }

    /// Magnetic susceptibility
    //inline double chi (const int i, const int j, const int k) const
    //{
    //    return (mu(i,j,k)/PhysicalConstants::mu_0 - 1);
    //}

    /// Magnetic field strength or H-Field
    inline dVector3 H0 (const int i, const int j, const int k) const
    {
        const double Hx = H0x + (dH0x_dx*i + dH0x_dy*j + dH0x_dz*k)*Grid.dx;
        const double Hy = H0y + (dH0y_dx*i + dH0y_dy*j + dH0y_dz*k)*Grid.dx;
        const double Hz = H0z + (dH0z_dx*i + dH0z_dy*j + dH0z_dz*k)*Grid.dx;
        return dVector3({Hx,Hy,Hz});
    }

    /// Magnetic field strength or H-Field
    inline dVector3 H (const int i, const int j, const int k) const
    {
        const double dchi_dx = (Grid.dNx) ? (phi(i+1,j,k)-phi(i-1,j,k))/2.0/Grid.dx : 0.0;
        const double dchi_dy = (Grid.dNy) ? (phi(i,j+1,k)-phi(i,j-1,k))/2.0/Grid.dx : 0.0;
        const double dchi_dz = (Grid.dNz) ? (phi(i,j,k+1)-phi(i,j,k-1))/2.0/Grid.dx : 0.0;

        const double Hx = H0x + (dH0x_dx*i + dH0x_dy*j + dH0x_dz*k)*Grid.dx - dchi_dx;
        const double Hy = H0y + (dH0y_dx*i + dH0y_dy*j + dH0y_dz*k)*Grid.dx - dchi_dy;
        const double Hz = H0z + (dH0z_dx*i + dH0z_dy*j + dH0z_dz*k)*Grid.dx - dchi_dz;
        return dVector3({Hx,Hy,Hz});
    }

    /// Magnetic flux density or B-Field
    inline dVector3 B (const int i, const int j, const int k) const
    {
        return H(i,j,k)*mu(i,j,k);
    }

    /// Material magnetization
    inline dVector3 M (const int i, const int j, const int k) const
    {
        return H(i,j,k)*chi(i,j,k);
    }

    /// Magnetic Force Density Field
    inline dVector3 ForceDensity (const int i, const int j, const int k) const
    {
        dMatrix3x3 gradB;
        gradB(0,0) = (Grid.dNx) ? (B(i+1,j,k)[0]-B(i-1,j,k)[0])/(2.0*Grid.dx) : 0.0;
        gradB(0,1) = (Grid.dNy) ? (B(i,j+1,k)[0]-B(i,j-1,k)[0])/(2.0*Grid.dx) : 0.0;
        gradB(0,2) = (Grid.dNz) ? (B(i,j,k+1)[0]-B(i,j,k-1)[0])/(2.0*Grid.dx) : 0.0;
        gradB(1,0) = (Grid.dNx) ? (B(i+1,j,k)[1]-B(i-1,j,k)[1])/(2.0*Grid.dx) : 0.0;
        gradB(1,1) = (Grid.dNy) ? (B(i,j+1,k)[1]-B(i,j-1,k)[1])/(2.0*Grid.dx) : 0.0;
        gradB(1,2) = (Grid.dNz) ? (B(i,j,k+1)[1]-B(i,j,k-1)[1])/(2.0*Grid.dx) : 0.0;
        gradB(2,0) = (Grid.dNx) ? (B(i+1,j,k)[2]-B(i-1,j,k)[2])/(2.0*Grid.dx) : 0.0;
        gradB(2,1) = (Grid.dNy) ? (B(i,j+1,k)[2]-B(i,j-1,k)[2])/(2.0*Grid.dx) : 0.0;
        gradB(2,2) = (Grid.dNz) ? (B(i,j,k+1)[2]-B(i,j,k-1)[2])/(2.0*Grid.dx) : 0.0;
        return gradB*M(i,j,k);
    }
};
}// namespace openphase
#endif
