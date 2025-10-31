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
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef INTERFACEDIFFUSION_H
#define INTERFACEDIFFUSION_H

#include "Includes.h"

namespace openphase
{
class InterfaceProperties;
class PhaseField;
class Settings;

class OP_EXPORTS InterfaceDiffusion : public OPObject                           ///< Calculates the interface diffusion inside the interface of multiple phases
{
 public:
    InterfaceDiffusion(){};                                                     ///< Empty constructor
    InterfaceDiffusion(Settings& locSetting,
                       const std::string InputFileName = DefaultInputFileName); ///< Constructor, initializes the class form input file
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes, just to indicate that module has been created.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads the initial parameters for the interface diffusion
	void ReadInput(std::stringstream& inp) override;                            ///< Reads the initial parameters for the interface diffusion
    void CalculatePhaseFieldIncrements(PhaseField& Phase,
            const InterfaceProperties& SigmaMu);                                ///< Calculates the movement of the interface due to the surface diffusion

    double PotentialDerivative(const double alpha, const double beta) const;
    void CalculateDiffusionPotential(PhaseField& Phase, const InterfaceProperties& IP);///< Calculates diffusion potential to the interface curvature
    void CalculateDiffusionPotentialGradients();                                ///< Calculates the gradient/diffusion-flux
    void CalculateDiffusionPotentialLaplacian(PhaseField& Phase);               ///< Laplacian of the diffusion potential

    Matrix<double>        Coefficients;                                         ///< Interface diffusion coefficient of each interface
    Storage3D< NodeAB<double,double>,  0 > DiffusionPotential;                  ///< Interface diffusion-potential of the tangential diffusion-flux

    GridParameters Grid;                                                        ///< Simulation grid parameters

    double                DoubleObstacleSmoothnessRange;                        ///< Defines the smoothing of the double obstacle potential
    unsigned int          Nphases;                                              ///< Number of thermodynamic phases

    LaplacianStencil LStencil;                                                  ///< Laplacian stencil. Uses user specified stencil as the basis
 protected:
 private:
};

class InterfaceDiffusionAnisotropic : public InterfaceDiffusion                 ///< Calculates the interface diffusion inside the interface of multiple phases
{
 public:
    InterfaceDiffusionAnisotropic(){};                                          ///< Empty constructor
    InterfaceDiffusionAnisotropic(Settings& locSetting,
            const std::string InputFileName = DefaultInputFileName);            ///< Constructor that initialises the class form input file
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;                            ///< Initializes, just to indicate that module has been created.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads the initial parameters for the interface diffusion

    void CalculatePhaseFieldIncrements(PhaseField& Phase,
            const InterfaceProperties& IP);                                     ///< Calculates the movement of the interface due to the surface diffusion
 protected:
    double PotentialDerivative(const double alpha, const double beta) const;
    void CalculateDiffusionFlux(PhaseField& Phase);                             ///< Calculates the interface diffusion flux
    void CalculateDiffusionFluxDivergence(PhaseField& Phase);                   ///< Calculates the divergence of the diffusion potential
    void CalculateDiffusionPotentialGradients();                                ///< Calculates the gradient/diffusion-flux

    Storage3D< NodeAB<dVector3,dVector3>, 0 > DiffusionFlux;                    ///< Interface diffusion-potential of the tangential diffusion-flux
 private:
};
}
#endif

