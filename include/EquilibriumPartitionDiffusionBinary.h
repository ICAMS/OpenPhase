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

 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef EQUILIBRIUMPARTITIONDIFFUSIONBINARY_H
#define EQUILIBRIUMPARTITIONDIFFUSIONBINARY_H

#include "Includes.h"

namespace openphase
{

class Settings;
class PhaseField;
class DrivingForce;
class Composition;
class Composition1Dextension;
class Temperature;
class InterfaceProperties;
class ElasticProperties;
class FractureField;
class BoundaryConditions;

enum class BoundaryConditionTypes;

enum class OP_EXPORTS BinaryDrivingForceModels
{
    Standard,
    LowerSlope,
    Weighted
};

class OP_EXPORTS EquilibriumPartitionDiffusionBinary : public OPObject                     ///<  Solver for diffusion using binary linearized phase diagrams
{
    friend ElasticProperties;
    friend FractureField;

 public:
    EquilibriumPartitionDiffusionBinary(){}
    EquilibriumPartitionDiffusionBinary(Settings& locSettings,
            const std::string InputFileName = DefaultInputFileName)             ///<  Initializes global settings
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings,
                    std::string ObjectNameSuffix = "") override;                ///<  Initializes global settings
    void ReadInput(const std::string InputFileName) override;                   ///<  Reads input parameters from a file
    void ReadInput(std::stringstream& inp) override;                            ///<  Reads input parameters from a file
    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///<  Remesh the storage while keeping the data

    void CalculateDrivingForce(PhaseField& Phase,
                               Composition& Cx,
                               Temperature& Tx,
                               DrivingForce& dGab);                             ///<  Calculates the driving force for each point
    void SolveDiffusion(PhaseField& Phase,
                        Composition& Cx,
                        Temperature& Tx,
                        BoundaryConditions& BC,
                        double dt);                                             ///<  Calculates the change of total concentrations in one time step taking into account cross terms
    void CalculateInterfaceMobility(PhaseField& Phase,
                                    Temperature& Tx,
                                    BoundaryConditions& BC,
                                    InterfaceProperties& IP);                   ///<  Calculates concentration-dependent mobility
    double ReportMaximumTimeStep(Temperature& Tx);                              ///<  Returns maximum time step for diffusion solver

 protected:
    Storage3D<double, 2> dMu;                                                   ///<  External contribution to the chemical potential.
    size_t Comp;                                                                ///<  Index of the chemical component for a given instance of the diffusion solver

 private:

    GridParameters Grid;                                                        ///<  Simulation grid parameters

    double R;                                                                   ///<  Universal gas constant

    BinaryDrivingForceModels DrivingForceModel;                                 ///<  Driving force models selector

    size_t RefComp;                                                             ///<  Index of the reference chemical component
    size_t Nphases;                                                             ///<  Number of thermodynamic phases
    size_t Ncomp;                                                               ///<  Number of chemical components
    double maxDC;                                                               ///<  Maximum diffusion coefficient in the simulation.

    std::vector<std::string> ElementNames;                                      ///<  Names of chemical elements

    std::vector<bool  > Stoichiometric;                                         ///<  Stoichiometry flags for each thermodynamic phase
    std::vector<double> IDC;                                                    ///<  Diffusivity enhancement coefficient in the interface.
    std::vector<double> DC;                                                     ///<  Solute temperature dependent diffusion coefficients for each thermodynamic phase
    std::vector<double> DC0;                                                    ///<  Solute diffusion coefficients for each thermodynamic phase
    std::vector<double> AE;                                                     ///<  Diffusion activation energies for each thermodynamic phase

    Tensor<double, 2> Slope;                                                    ///<  Pairwise liquidus slopes for a linear phase diagram
    Tensor<double, 2> Tcross;                                                   ///<  Pairwise temperatures of the liquidus-solidus or solvus lines intersections for a linear phase diagram
    Tensor<double, 2> Ccross;                                                   ///<  Pairwise concentrations of the liquidus-solidus or solvus lines intersections for a linear phase diagram

    std::vector<double> Entropy;                                                ///<  Entropies of different phases. Their differences are used in the driving force

    LaplacianStencils DiffusionStencil;                                         ///<  Diffusion stencil selector
    LaplacianStencil DStencil;                                                  ///<  Diffusion stencil. Uses Laplacian stencil as the basis

    bool EnableAntiTrapping;                                                    ///<  Enables antitrapping current calculation
    bool ConsiderChemicalPotential;                                             ///<  Enables chemical potential consideration

    void CalculateDiffusionCoefficients(Temperature& Tx);                       ///<  Calculates temperature dependent diffusion coefficients depending on the average temperature

    void CalculateLocalPhaseConcentrations(PhaseField& Phase,
                                           Composition& Cx,
                                           Temperature& Tx,
                                           int i, int j, int k);                ///<  Distributes the total concentrations into concentrations in each phase for a given point
    double EquilibriumComposition(size_t n, size_t m, double Temp);             ///<  Returns equilibrium composition for phase {n} in {n,m} phase pair
    double PartitioningCoefficient(size_t n, size_t m);                         ///<  Returns partitioning coefficient for phases {n,m}

    void CalculatePhaseConcentrations(PhaseField& Phase,
                                      Composition& Cx,
                                      Temperature& Tx);                         ///<  Distributes the total concentrations into concentrations in each phase
    void CalculateDiffusionIncrements(PhaseField& Phase,
                                      Composition& Cx);                         ///<  Calculates Fick's diffusion composition increments
    void CalculateChemicalPotentialContribution(PhaseField& Phase,
                                                Composition& Cx,
                                                Temperature& Tx);               ///<  Calculates diffusion due to externally supplied chemical potential contribution
    void ClearChemicalPotential(void);                                          ///<  Clears chemical potential storage
    void CalculateAntitrappingIncrements(PhaseField& Phase, Composition& Cx);   ///<  Calculates antitrapping concentration increments
    bool CalculateLimits(PhaseField& Phase, Composition& Elements, double dt);  ///<  Calculated limits for increments in order to keep concentrations within physical range
    void LimitDiffusionIncrements(PhaseField& Phase,
                                  Composition& Cx,
                                  Temperature& Tx);                             ///<  Limits Fick's diffusion concentration increments
    void ApplyIncrements(PhaseField& Phase, Composition& Cx, double dt);        ///<  Apply the limited increments
    void Solve1Dextension(Composition& Cx,
                          Temperature& Tx,
                          Composition1Dextension& CxExt,
                          BoundaryConditionTypes& lextBC,
                          double dt);                                           ///<  Solves diffusion in 1D composition extension
};

} // namespace openphase
#endif
