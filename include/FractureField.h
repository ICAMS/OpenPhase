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

 *   File created :   2022
 *   Main contributors :   Oleg Shchyglo; Hossein Jafarzadeh, Muhammad Adil Ali
 *
 */

#ifndef FRACTUREFIELD_H
#define FRACTUREFIELD_H

#include "Includes.h"

namespace openphase
{

class ElasticProperties;
class BoundaryConditions;
class Velocities;
class InterfaceProperties;
class PhaseField;
class AdvectionHR;
class Settings;
class DoubleObstacle;
class Composition;
class Temperature;
class DrivingForce;
class Orientations;
class RunTimeControl;

enum FractureModels
{   
    Obstacle,               ///< Fracture model from  Schneider et al. http://dx.doi.org/10.1016/j.cma.2016.04.009
    Well,                   ///< Fracture model from  Houssein et al.  https://doi.org/10.1007/s10704-024-00762-x
};

/******************************************************************************/
class FractureField : public OPObject                                           ///< Fracture field class. It stores the crack phase field, performs basic operations on them.
{
 public:
    FractureField(){ };

    FractureField(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }

    ///< Initializes the storage and initial variables of the driving force class.
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;
    void ReadInput(const std::string InputFileName) override;                            ///< Reads input for the fracture field from the input file
    void ReadInput(std::stringstream& inp) override;                                     ///< Reads input for the fracture field from the input stream

    void CreateCrack(dVector3 &CrackStart,
                     dVector3 &CrackEnd,
                     Settings &OP,
                     BoundaryConditions &BC);                                            ///< Create a crack in the system From Houssein
    
    void Solve(PhaseField &Phase,
               DoubleObstacle &DO,
               InterfaceProperties &IP,
               BoundaryConditions &BC,
               ElasticProperties &EP,
               double dt);                                                               ///< Solves the fracture field equation

    void CalculateIncrementsWell(const PhaseField& Phase,
                            InterfaceProperties& IP,
                            ElasticProperties& EP,
                            DoubleObstacle& DO);                                         ///< Calculates interface curvature related driving force.

    void CalculateIncrementsObstacle(const PhaseField& Phase,
                            InterfaceProperties& IP,
                            ElasticProperties& EP,
                            DoubleObstacle& DO);                                         ///< Calculates interface curvature related driving force.
    
    void CorrectFractureProfile(BoundaryConditions &BC, double locMobility, size_t nSteps);

    void CreatePF(PhaseField& Phase, size_t PhaseIndex, BoundaryConditions& BC);
    void SetMobilityPF(PhaseField &Phase, InterfaceProperties &IP);
    void SetDrivingForcePF(const PhaseField &Phase, DrivingForce &dGab);

    void SetSurfaceEnergy(Composition& Cx, Temperature& Tx);

    void Advect(AdvectionHR& Adv,
                const Velocities& Vel,
                PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt,
                const double tStep) override;                                            ///< Advects phase-fields
    
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///< Changes the mesh size while keeping the data.
    void SetBoundaryConditions(const BoundaryConditions& BC) override;                   ///< Set boundary conditions

    bool Read(const Settings& locSettings,const BoundaryConditions& BC, 
                                            const int tStep) override;                   ///< Reads raw (binary) phase fields from the file
    bool Read(std::string FileName);                                                     ///< Reads raw (binary) phase fields from the file
    bool Write(const Settings& locSettings,const int tStep) const override;              ///< Writes raw (binary) phase fields to the file

    void WriteLaplacianVTK(const int tStep, const Settings& locSettings,
                                               const int precision = 16);                ///< Writes Laplacian of a fracture field to the file in VTK format
    void WriteVTK(const Settings& locSettings, const int tStep,
                                        const int precision = 16) const;                 ///< Writes fracture fields to the file in VTK format
    void WriteEnergyVTK(const Settings& locSettings, const int tStep,
                        const InterfaceProperties& IP, const int precision=16);          ///< Writes fracture energy in VTK format for a given timestep tStep
    void WriteDistortedVTK(const Settings& locSettings, const ElasticProperties& EP,
                           const int tStep,  const int precision = 16) const;            ///< Writes fracture fields on the distorted grid to the file in VTK format in MPI environment
    void ApplyAdvection(Settings &locSettings,
                        PhaseField& Phase,
                        ElasticProperties& EP,
                        Orientations& OR,
                        AdvectionHR &Adv,
                        Velocities &Vel,
                        BoundaryConditions& BC,
                        RunTimeControl &RTC,
                        double dt);                                                      ///< Applies advection to the phase fields
    
    void FixAdvectionFractureProfile(ElasticProperties &EP);

    void SaveDisplacements(ElasticProperties &EP);
    void SetDisplacementsIncrement(Velocities& Vel, ElasticProperties &EP, double dt);
    double Energy(const InterfaceProperties &IP);                                        ///< Provides total interfacial energy in the simulation domain.
    double AverageEnergyDensity(const InterfaceProperties& IP);                          ///< Provides the average interfacial energy density in the simulation domain.
    double PointEnergy(const InterfaceProperties& IP, int i, int j, int k);              ///< Provides the interfacial energy in a given point (i,j,k).
    void PrintVolumeFractions(void);                                                  ///< Calculates the volume of the fracture field

    double sWidth; ///< (crack) surface width in grid points

    Storage3D<double, 0> Fields;                                                         /// Storage for the fracture field values
    Storage3D<double, 0> Fields_dot;                                                     /// Storage for the rate of change of the fracture field values
    Storage3D<double, 0> Laplacian;                                                      /// Storage for the Laplacian of the fracture field
    Storage3D<double, 0> Flag;                                                           /// Storage for flags indicating the state of the fracture field
    Storage3D<dVector3, 0> DisplacementsOLD;                                             /// Storage for the Old displacements for the Advection
    Storage3D<double, 0> SurfaceEnergy;                                                     ///< Storage for the crack tip values
    Storage3D<double, 0> TestOutput;                                                     ///< Storage for the crack tip values
    Storage3D<double, 0> TestOutput2;                                                     ///< Storage for the crack tip values
    
    FractureField& operator= (const FractureField& rhs);

    double Eta;                                                                          ///< Interface width in physical units
    double Mobility;                                                                     ///< Crack mobility
    void Finalize(const BoundaryConditions& BC);                                         ///< Finalizes the phase fields calculations in this time step
    void SetFracturedElasticConstants(ElasticProperties& EP);                            ///< Applies degredation to ElasticProperties::EffectiveElasticConstants
    FractureModels FractureModel;                                                          ///< Fracture model type
    GridParameters Grid;
    double PFcutOff;                                                               ///< Cut-off value for the Fracture Field for creating a PhaseField
    void MergeIncrements(const BoundaryConditions& BC, const double dt);                 ///< Merges the increments into the phase fields
  protected:
  private:
    void Clear();                                                                        ///< Clears the phase field storage
    void SetFlags();                                                                     ///< Sets the flags which mark interfaces
    void CalculateLaplacians(void);                                                      ///< Calculates Laplacians and stores them in the phasefields nodes
    dVector3 Gradient(const int i, const int j, const int k) const;                      ///< Returns gradient of the phase field at a given position
    dVector3 Gradient_v2(const int i, const int j, const int k) const;                   ///< Returns gradient of the phase field at a given position
    bool Interface(const int x, const int y, const int z) const                          ///< Indicates the location of the interface.
    {
        return (Flag(x, y, z) > 1);
    };
    void CalculateCrackTip();

    double iWidth;                                                              ///< Interface width in grid points
    double RefVolume;                                                           ///< Reference volume for the nucleation
    double Xi;                                                                  ///< hydrogen degrading coefficient for fracture energy
    double delta_g_0;                                                           ///< Gibbs free energy difference between the decohering interface and the surrounding material
    double maxGradient;
    double InitialSurfaceEnergy;

    LaplacianStencils FractureFieldLaplacianStencil;                            ///< Laplacian stencil selector
    GradientStencils  FractureFieldGradientStencil;                             ///< Gradient stencil selector
    LaplacianStencil LStencil;                                                  ///< Laplacian stencil. Uses user specified stencil as the basis
    GradientStencil GStencil;                                                   ///< Gradient stencil. Uses user specified stencil as the basis

    bool ChemoMechanicalCoupling;
    std::vector<std::string> ElementNames;                                      ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)
    size_t Ncomp;                                                               ///< Number of components
    size_t Comp;

};

}// namespace openphase
#endif
