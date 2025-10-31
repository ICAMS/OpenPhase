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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#ifndef ELASTICPROPERTIES_H
#define ELASTICPROPERTIES_H

#include "Includes.h"
#include "SymmetryVariants.h"
#include "Tools.h"

namespace openphase
{

class AdvectionHR;
class BoundaryConditions;
class Composition;
class DrivingForce;
class EquilibriumPartitionDiffusionBinary;
class InterfaceProperties;
class PhaseField;
class Settings;
class Temperature;
class Velocities;

enum class ElasticityModels : int                                               ///< Elasticity homogenization and materials models
{
    Khachaturyan,                                                               ///< Khachaturyan's material model
    Steinbach,                                                                  ///< Steinbach's material model
    Reuss,                                                                      ///< Reuss/Sachs homogenization model
    Voigt,                                                                      ///< Voigt/Taylor homogenization model
    Rank1,                                                                      ///< Rank 1 homogenization model
    Rank1NL                                                                     ///< Rank 1 homogenization model with non-local contribution
};

enum class StrainModels : int                                                   ///< Models of mechanical strain
{
    Small                                                                       ///< Small strain model
};

inline vStrain StrainSmall(const dMatrix3x3& locDefGrad)                        ///< Small strain model
{
    return VoigtStrain((locDefGrad.transposed() + locDefGrad)*0.5 - dMatrix3x3::UnitTensor());
};

class OP_EXPORTS ElasticProperties : public OPObject                            ///< Module which stores and handles elastic properties
{
 public:
    ElasticProperties(){};                                                      ///< Constructor (does nothing)
    ElasticProperties(Settings& locSettings,
                      const std::string InputFileName = DefaultInputFileName);  ///< Constructs and initializes the elastic properties

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates storages, initializes internal parameters
    void InitializeLD(Settings& locSettings);                                   ///< Initialized parameters required by the large deformation framework
    void ReadInput(const std::string InputFileName) override;                   ///< Reads elastic properties from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads elastic properties from the input stream
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the elasticity storages

    void SetBoundaryConditions(const BoundaryConditions& BC)  override;         ///< Sets boundary conditions for strains, stresses and displacements
    void SetBoundaryConditionsPlastic(const BoundaryConditions& BC);            ///< Sets boundary conditions Plastic

    void SetEffectiveProperties(const PhaseField& Phase);                       ///< Sets effective elastic properties
    void SetEffectiveProperties(const PhaseField& Phase, const Composition& Cx);///< Sets effective elastic properties considering chemo-mechanical coupling
    void SetEffectiveProperties(const PhaseField& Phase, const Temperature& Tx);///< Sets effective elastic properties considering thermo-mechanical coupling
    void SetEffectiveProperties(const PhaseField& Phase, const Composition& Cx,
                                                         const Temperature& Tx);///< Sets effective elastic properties  considering chemo- and thermo-mechanical coupling
    void SetEffectiveProperties(const PhaseField& Phase, const InterfaceProperties& IP);///< Sets effective elastic properties considering interface stress contribution

    void CalculateDrivingForce(const PhaseField& Phase, DrivingForce& dGab);    ///< Calculates elastic driving force

    void CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceProperties& IP) const;///< Calculates interface energy contribution
    void CalculateChemicalPotentialContribution(const PhaseField& Phase,
                                                EquilibriumPartitionDiffusionBinary& DF);///< Calculates chemical potential contribution
    
    double Energy(const PhaseField& Phase) const;                               ///< Elastic energy in simulation domain [Joule]
    double EnergyDensity(const PhaseField& Phase, int i, int j, int k) const;   ///< Returns elastic energy density in a given point [Joule/m^3]
    double AverageEnergyDensity(const PhaseField& Phase) const;                 ///< Average elastic energy density in [Joule/<simulation cell>]

    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep) override;

    // Binary input/
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tSep) override;           ///< Reads raw deformation gradients from a file

    // Binary output
    bool Write(const Settings& locSettings, const int tSep) const override;     ///< Writes restart output in binary format

    // VTK output
    void WriteEffectiveElasticConstantsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes effective elastic constants in VTK format
    void WriteEigenStrainsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes effective eigenstrains using Voigt index notation, there is a factor 2 difference to the tensor notation
    void WriteDeformationGradientsTotalVTK(Settings& locSettings, const int tStep,
            const int precision=20) const;                                      ///< Writes EffectiveEigenstrains using Voigt index notation, there is a factor 2 difference to the tensor notation
    void WriteElasticStrainsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes elastic strains in VTK format
    void WriteEnergyDensityVTK(Settings& locSettings, const PhaseField& Phase,
            const int tStep, const int precision=16) const;                     ///< Writes local elastic energy in VTK format
    void WritePlasticStrainsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes plastic strains in VTK format
    void WriteStressIncrementsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes stress increments in VTK format
    void WriteStressesVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes stresses in VTK format
    void WriteTotalStrainsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes total strains in VTK format
    void WriteForceDensityVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes force densities in VTK format
    void WriteForceDensityDistortedVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes force densities in distorted configuration in VTK format
    void WriteVelocityGradientVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes velocity gradients in VTK format
    void WriteTotalRotationsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes local rotations (in VTK format) extracted from total deformation gradient as axis/angle pair
    void WriteDisplacementsVTK(Settings& locSettings, const int tStep,
            const int precision=16) const;                                      ///< Writes displacement vectors in VTK format
    void WriteDeformationJumpsVTK(Settings& locSettings, const int tStep,
            const int idxA, const int idxB,
            const int precision=16) const;                                      ///< Writes Hadamard jump vectors in VTK format
    void WriteLocalRotationsVTK(Settings& locSettings, const int tStep,
            const int precision = 16);                                          ///< Writes local rotations in VTK format
    // ASCII output
    void WriteStressStrainData(int time_step, double time,
                               std::string filename = "StessStrainCurve.dat") const;///< Writes stress-strain curve into an ASCII file
    void WriteDilatometerCurve(std::string filename,
                                 double temperature,
                                 double real_size = 1.0) const;                 ///< Write dilatometer curve

    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints various properties at a given point (x, y, z) to screen

    // Auxiliary methods
    void SetVelocityGradients(double dt);                                       ///< Sets velocity gradients
    void SetLocalRotations(double dt);                                          ///< Sets local rotations
    void SetStressesRX(PhaseField& Phase,
                       BoundaryConditions& BC, std::vector<int> targetPhases);
    ElasticProperties& operator= (const ElasticProperties& rhs);

    // Parameters and storages
    vStress AverageStress;                                                      ///< Average stress calculated by the spectral solver.
    vStrain AverageStrain;                                                      ///< Average strain calculated by the spectral solver.
    vStrain StrainToRemesh;                                                     ///< Accumulates the strain that has to be applied by remeshing
    vStrain RemeshedStrain;                                                     ///< Accumulates the strain that has been applied by remeshing
    vStress AppliedStress;                                                      ///< Applied stress tensor
    vStrain AppliedStrain;                                                      ///< Applied strain tensor
    vStrain AppliedStrainRate;                                                  ///< Applied strain rate tensor

    dVector6 AppliedStrainMask;                                                 ///< Applied strain markers for individual strain components
    dVector6 AppliedStrainRateMask;                                             ///< Applied strain rate markers for individual strain components
    dVector6 AppliedStressMask;                                                 ///< Applied stress markers for individual stress components

    Storage<dMatrix3x3> GrainTransformationStretches;                           ///< Reference transformation stretches for each phase field
    Storage<dMatrix6x6> GrainElasticConstants;                                  ///< Reference elastic constants for each phase field

    Storage<dMatrix3x3> PhaseTransformationStretches;                           ///< Reference transformation stretches gradients of each phase (relative to Cref), grain orientation not considered
    Storage<dMatrix6x6> PhaseElasticConstants;                                  ///< Reference elastic constants of each phase (relative to Cref), grain orientation not considered
    Storage<double>     PoissonRatio;                                           ///< Poisson Ratio of each phase

    Storage3D<dMatrix3x3,0> VelocityGradientsTotal;                             ///< Storage for deformation gradient rates
    Storage3D<dMatrix3x3,0> DeformationGradientsTotal;                          ///< Storage for total deformation gradients
    Storage3D<dMatrix3x3,0> DeformationGradientsTotalAdvectionDot;              ///< Temporary storage for total deformation gradients
    Storage3D<dMatrix3x3,0> VelocityGradientsPlastic;                           ///< Storage for advection of plastic deformation gradient rates
    Storage3D<dMatrix3x3,0> DeformationGradientsPlastic;                        ///< Storage for plastic deformation gradients
    Storage3D<dMatrix3x3,0> DeformationGradientsPlasticAdvectionDot;            ///< Temporary storage for advection of plastic deformation gradients
    Storage3D<dMatrix3x3,0> VelocityGradientsEigen;                             ///< Storage for transformation induced deformation gradient rates
    Storage3D<dMatrix3x3,0> DeformationGradientsEigen;                          ///< Storage for transformation induced deformation gradients
    Storage3D<dMatrix3x3,0> DeformationGradientsEigenAdvectionDot;              ///< Temporary storage for advection of transformation induced deformation gradients
    Storage3D<NodeAB<dVector3,dMatrix3x3>,0> DeformationJumps;                  ///< Storage for deformation jumps in the interfaces
    Storage3D<NodeAB<dVector3,dVector3>,0>   NonLocalDeformationJumpTerm;       ///< Storage for non-local deformation jump contribution to the driving force

    Storage3D<vStress,0>    Stresses;                                           ///< Storage for 2PK stresses in the intermediate configuration
    Storage3D<vStress,0>    StressesAdvectionDot;                               ///< Temporary storage for advection of stresses
    Storage3D<vStress,0>    StressIncrements;                                   ///< Storage for stress increments

    Storage3D<dMatrix6x6,0> EffectiveElasticConstants;                          ///< Storage for effective elastic constants

    Storage3D<NodeA<dMatrix6x6>,0> ElasticConstants;                            ///< Storage for elastic constants for each phase field in each grid point
    Storage3D<NodeA<dMatrix3x3>,0> TransformationStretches;                     ///< Storage for transformation stretches for each phase field in each grid point

    Storage3D<dVector3,0>   ForceDensity;                                       ///< Force density storage
    Storage3D<dVector3,0>   Displacements;                                      ///< Displacements storage

    Storage3D<Quaternion,0> LocalRotations;                                     ///< Storage for local rotations due to large deformation
    Storage3D<Quaternion,0> LocalRotationsAdvectionDot;                         ///< Temporary storage for advection of local rotations

    SymmetryVariants        Variants;                                           ///< Symmetry variants storage

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    bool KeepAspectRatio;                                                       ///< Restricts mechanical conditions to only volume relaxation
    bool PreventShear;                                                          ///< Prevents shear deformation of the simulation domain
    bool LargeDeformations;                                                     ///< Large deformation flag.
    bool AnyPlasticity;                                                         ///< True if any plasticity model is active
    bool ConsiderExternalForces;                                                ///< Enables external force density consideration

    dMatrix6x6 MAXElasticConstants(void) const;                                 ///< Calculates maximum elastic constants
    dMatrix6x6 AverageElasticConstants(void) const;                             ///< Calculates average elastic constants
    dMatrix3x3 AverageDeformationGradient(double accuracy = 1.0e-6) const;      ///< Calculates average deformation gradient
    vStrain    AveragePlasticStrain(void) const;                                ///< Calculates average plastic strain

    dMatrix3x3 DeformationGradientsElastic(int i, int j, int k) const;          ///< Returns elastic deformation gradient considering the strain model in use
    vStrain TotalStrains(int i, int j, int k) const;                            ///< Returns strain considering strain operation mode (small, finite)
    vStrain EigenStrains(int i, int j, int k) const;                            ///< Returns strain considering strain operation mode (small, finite)
    vStrain PlasticStrains(int i, int j, int k) const;                          ///< Returns strain considering strain operation mode (small, finite)
    vStrain StressFreeStrains(int i, int j, int k) const;                       ///< Returns strain considering strain operation mode (small, finite)
    vStrain ElasticStrains(int i, int j, int k) const;                          ///< Returns strain considering strain operation mode (small, finite)
    vStress ShearStresses(int i, int j, int k) const;                           ///< Returns shear stress
    vStress Stress2PK(int i, int j, int k, const dMatrix3x3& locDefGrad) const; ///< Returns second P-K stress for a given total deformation gradient in the intermediate configuration
    dMatrix3x3 Stress2PKpullBack(int i, int j, int k) const;                    ///< Returns second P-K stress in the reference configuration
    dMatrix3x3 Stress1PK(int i, int j, int k) const;                            ///< Returns first P-K stress in a given point

    vStrain EigenStrainDifference(
            const dMatrix3x3 locStretchesAlpha,
            const dMatrix3x3 locStretchesBeta,
            const int i, const int j, const int k) const;                       ///< Returns local eigenstrain difference based on the given transformation stretch input
            
    ElasticityModels ElasticityModel;                                           ///< Elasticity model selector
    StrainModels     StrainModel;                                               ///< Strain model selector

 protected:
 private:

    bool ThermoMechanicalCoupling;                                              ///< Thermo-mechanical coupling flag
    Storage<dMatrix3x3> GrainAlpha;                                             ///< Linear thermal expansion coefficients for each phase field
    Storage<dMatrix3x3> PhaseAlpha;                                             ///< Linear thermal expansion coefficients for each phase
    Storage<dMatrix6x6> GrainGamma;                                             ///< Linear temperature dependence coefficients for the elasticity parameters for each phase field
    Storage<dMatrix6x6> PhaseGamma;                                             ///< Linear temperature dependence coefficients for the elasticity parameters for each phase
    Storage<double>     Tref;                                                   ///< Reference temperature for temperature dependence parameters for each phase

    bool ChemoMechanicalCoupling;                                               ///< Chemo-mechanical coupling flag
    std::vector < std::string > ElementNames;                                   ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)
    /* template specialization integer number stands for
     * the number of extra dimensions
     * order:
     * Rank = 2: phase, component
     * Rank = 1: phase or component
     */
    Tensor<dMatrix6x6,2> GrainKappa;                                            ///< Linear composition dependence coefficients for the elasticity parameters for each phase field
    Tensor<dMatrix6x6,2> PhaseKappa;                                            ///< Linear composition dependence coefficients for the elasticity parameters for each phase
    Tensor<dMatrix3x3,2> GrainLambda;                                           ///< Linear composition dependence coefficients for the transformation stretches for each phase field
    Tensor<dMatrix3x3,2> PhaseLambda;                                           ///< Linear composition dependence coefficients for the transformation stretches for each phase
    Tensor<double,2>     Cref;                                                  ///< Reference composition for composition dependence parameters for each phase

    void SetGrainsProperties(const PhaseField& Phase);                          ///< Sets elastic properties for each grain according to its phase and orientation

    void SetBaseTransformationStretches(const PhaseField& Phase);               ///< Sets effective transformation stretches base values
    void AddThermalExpansion(const PhaseField& Phase, const Temperature& Tx);   ///< Calculates thermal expansion contribution to the transformation stretches
    void AddVegardsExpansion(const PhaseField& Phase, const Composition& Cx);   ///< Calculates Vegard's expansion contribution to the transformation stretches
    void AddInterfaceStressContribution(const PhaseField& Phase, const InterfaceProperties& IP);///< Calculates interface contribution to the transformation stretches
    void SetDeformationGradientsEigen(const PhaseField& Phase);                 ///< Sets effective transformation-induced deformations

    void SetBaseElasticConstants(const PhaseField& Phase);                      ///< Sets base elastic constants for each grain
    void AddThermalSoftening(const PhaseField& Phase, const Temperature& Tx);   ///< Calculates thermal softening contribution to the elastic constants
    void AddChemicalSoftening(const PhaseField& Phase, const Composition& Cx);  ///< Calculates chemical softening contribution to the elastic constants
    void SetEffectiveElasticConstants(const PhaseField& Phase);                 ///< Sets effective elastic constants

    void ApplyLocalRotationsToDeformationGradientsEigen(void);                  ///< Applies local rotations to transformation stretches base values
    void ApplyLocalRotationsToElasticConstants(void);                           ///< Applies local rotations to elastic constants

    void CalculateDeformationJumps(const PhaseField& Phase);                    ///< Calculates deformation jumps from the Hadamard jump condition for all phase-field pairs
    void CalculateNonlocalJumpContribution(const PhaseField& Phase);            ///< Calculates non-local deformation jump contribution to the driving force
};
}// namespace openphase
#endif
