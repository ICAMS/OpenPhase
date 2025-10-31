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
 *   Main contributors :   Oleg Shchyglo; Dmitri Medvedev; Amol Subhedar;
 *                         Marvin Tegeler; Raphael Schiedung
 *
 */

#ifndef FLOWSOLVERLBM_H
#define FLOWSOLVERLBM_H

#include "Includes.h"
#include "FluidDynamics/D3Q27.h"

namespace openphase
{

class BoundaryConditions;
class Composition;
class D3Q27;
class DrivingForce;
class GrainsProperties;
class Node;
class PhaseField;
class Settings;
class Temperature;
class Velocities;

class OP_EXPORTS FlowSolverLBM : public OPObject                                ///< Calculation of fluid flow and advective solute transport
{
public:
    FlowSolverLBM(void){};                                                      ///< Empty constructor
    FlowSolverLBM(Settings& locSettings, double in_dt,
                  const std::string InputFileName = DefaultInputFileName);      ///< Read Settings and InputFile with constructor
    void Initialize(Settings& locSettings, std::stringstream& inp, double in_dt);                       ///< Allocates memory and initializes global settings
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from a file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from a stringstream
    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;

    static D3Q27 EquilibriumDistribution(double lbDensity,
            const double weights[3][3][3], dVector3 lbMomentum = {0.0,0.0,0.0});

    double DensityProfile(const double x, const size_t n=0) const;              ///< Simple tangent hyperbolic for the liquid-vapor interface

    static double OptimalParaKuper(const double ReducedTemperature,
            const double GasPrameter = 0.01);                                   ///< Calculates the optimal Parameter ParaKuper

    dMatrix3x3 PressureTensor(const int i, const int j, const int k) const;     ///< Calculates the pressure tensor at the point (i,j,k)
    std::vector<dVector3> CalculateFluidMomentum(const PhaseField& Phase) const;///< Calculates the momentum of the entire fluid in the system
    double BounceBack(const int i, const int j, const int k,
            const int ii, const int jj, const int kk, const size_t n,
            PhaseField& Phase, const BoundaryConditions& BC,
            double& lbDensityChange); 
    double SecondOrderBounceBack(const int i, const int j, const int k,
            const int ii, const int jj, const int kk, const size_t n,
            PhaseField& Phase, const BoundaryConditions& BC,
            double& lbDensityChange);                                            ///< BounceBack at Solid Interfaces
    double Pressure(const int i, const int j, const int k) const;               ///< Calculates the pressure at the point (i,j,k)
    //double PressureTensorXX(const double i0, const double j , const double k ) const; ///< Calculates the xx-component of the pressure tensor

    size_t CountFluidNodes(void) const;                                         ///< Counts the number of fluid nodes
    size_t CountObstacleNodes(void) const;                                      ///< Counts the number of obstacle nodes
    std::vector<double> CalculateFluidMass(void) const;                         ///< Calculates mass of fluid components
    double CalculateLiquidVolume(size_t Fluid_Comp);                            ///< Calculates volume of liquid phase
    double CalculateLiquidMass(size_t Fluid_Comp);                              ///< Calculates mass of liquid phase
    double CalculateVaporVolume(size_t Fluid_Comp);                             ///< Calculates volume of vapor phase
    double CalculateVapordMass(size_t Fluid_Comp);                              ///< Calculates mass of vapor phase
    virtual void SetObstacleNodes(const PhaseField& Phase,
            const Velocities& Vel);                                             ///< Sets density and momentum in obstacles to be used when solid moves
    void ApplyForces(PhaseField& Phase, const Velocities& Vel,
            const Composition& Cx);                                             ///< Applies forces
    void ApplyForces(PhaseField& Phase, const Velocities& Vel);                 ///< Applies force without composition
    void CalculateDensityAndMomentum(void);                                     ///< Calculates Density and momentum from lbDensity
    void CalculateFluidVelocities(Velocities& Vel, const PhaseField& Phase,
            const BoundaryConditions& BC) const;                                ///< Calculates Fluid velocities
    void CalculateForceDrag(const int i, const int j, const int k,
            PhaseField& Phase, const Velocities& Vel);                          ///< Force contribution of surface drag
    void CalculateForceGravity(PhaseField& Phase);                              ///< Force contribution of Gravitation
    void CalculateForceGravity(PhaseField& Phase, const Composition& Cx);       ///< Force contribution of Gravitation considering liquid composition
    void CalculateForceTwoPhase(const int i, const int j, const int k,
            PhaseField& Phase);                                                 ///< Force contribution according to Benzi
    double CalculateLocalIncomingSolidDensityChange(PhaseField& Phase, int i, int j, int k);
    void AcountForAndLimitPhaseTransformation(PhaseField& Phase, int tStep);    ///< Add/Removes fluid mass according to phase changes
    void Collision();                                                           ///< Processes the collisions. Adjusts center of mass properties of colliding particles
    void DetectObstaclesAdvection(const PhaseField& Phase,
            const Velocities &Vel, const BoundaryConditions& BC);               ///< Detects Obstacles (Detects change of obstacles for advection)
    void DetectObstacles(const PhaseField& Phase);                              ///< Detects Obstacles
    void EnforceMassConservation(void);                                         ///< Enforces fluid mass conservation for advectional problems
    void EnforceTotalMassConservation(PhaseField& Phase);                       ///< Enforces total mass conservation (solid and fluid)
    void EnforceSolidMomentum(PhaseField& Phase, const dVector3 value);         ///< Enforces solid momentum conservation

    void FixPopulations(void);                                                  ///< Ensures positive lbPopulations
    void Propagation(PhaseField& Phase, const BoundaryConditions& BC);          ///< Propagates Populations
    void PropagationSecondOrderBB(PhaseField& Phase, const BoundaryConditions& BC);          ///< Propagates Populations
    bool Read(const Settings& locSettings, const BoundaryConditions& BC,
              const int tStep = -1) override;                                   ///< Read raw (binary) fields from file
    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions for the particle distribution functions and flow velocities
    void SetFluidNodesNearObstacle();                                           ///< Sets lbPoblulations of vanishing Obstacles and of appearing Obstacles
    void SetUniformVelocity(const BoundaryConditions& BC, const dVector3 U0);   ///< Sets the initial values of particle distribution functions
    void Solve(PhaseField& Phase, Velocities& Vel,
            const BoundaryConditions& BC);                                      ///< Calculates one time step of the Navier-Stokes solver
    void Solve(PhaseField& Phase, const Composition& Cx, Velocities& Vel,
            const BoundaryConditions& BC);                                      ///< Calculates one time step of the Navier-Stokes solver
    void SolveTC3D(PhaseField& Phase, Velocities& Vel,
             Composition& Cx,  Temperature& Tx,
            const BoundaryConditions& BC);
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Writes raw (binary) fields to file
    void WriteVTK(const Settings& locSettings, const PhaseField& Phase,
                  const int tStep, const int precision = 16) const;             ///< Writes physical fields values into file with VTK format
    void lbWriteVTK(const Settings& locSettings, const PhaseField& Phase,
                    const int tStep, const int precision = 16) const;           ///< Writes physical fields values into file with VTK format
    bool SingleSolid(const int i, const int j, const int k,
             const PhaseField& Phase) const;                                    ///< Returns true if single solid is present at (i,j,k)
    double SolidFraction(const int i, const int j, const int k,
             const PhaseField& Phase) const;                                    ///< Returns if solid fraction at (i,j,k)
    bool LocalObstacle(const int i, const int j, const int k,
            const PhaseField& Phase) const;                                     ///< Returns true if obstacle has been detected at (i,j,k)
    void Update_dt(double _dt);                                                 ///< Calculates one time step of the Navier-Stokes solver
    void SetInitialPopulationsTC(BoundaryConditions& BC , Velocities& Vel);     ///< Initializing the populations when thermal compressibility considered
    void CollisionTC(Velocities& Vel);
    void CalculateHydrodynamicPressureAndMomentum(Velocities& Vel);
    void CalculateDensityTC(Temperature& Tx, Composition& Cx, const BoundaryConditions& BC);   ///<  Calculates Density from ideal gas law
    void SetOutFlow(int order);
    void SetPressureOutlet(Velocities& Vel, const PhaseField& Phase, double KMaxv);
    void SetDivUnearObstZero();
    void SetDivVelBCNXZero();
    void SetDivUnearBCsZeroXZPlane();
    void SetVelocityInletUniform(const BoundaryConditions& BC, const PhaseField& Phase);
    void SetCornerCorrection(const BoundaryConditions& BC, const PhaseField& Phase);
    void Set1DNonReflectingPressureOutletatBCXN(Velocities& Vel, const PhaseField& Phase);
    void CalculateDivergenceVelocity(Temperature& Tx, PhaseField& Phase,Velocities& Vel, BoundaryConditions& BC, double dt);
    void FinalizeInitialiation(const PhaseField& Phase, const Velocities& Vel,
            const BoundaryConditions& BC);                                      ///< Needs to be called after the microstructure has been initialized
    D3Q27 EquilibriumDistributionTC(dVector3 lbvel, double lbph, double lbnut,  double lbDivVel, 
                                               dVector3 lbFb, const double weights[3][3][3], double lbrho,
                                               dVector3 lbGradRho, size_t n); 

    double CalculateSutherlandViscosity(double p, double mw, double Temp);
    double CalculateIdealGasDensity(double p, double Rm, double Temp);
    void CalculateDensityGradient(PhaseField& Phase, Velocities& Vel);

    dVector3 FluidVelocity(const int i, const int j, const int k, const size_t comp) const; ///< Returns local fluid velocity of component comp
    dVector3 FluidVelocity(const int i, const int j, const int k) const ;       ///< Returns local fluid velocity summed over all components
    double FluidDensity(const int i, const int j, const int k, const size_t comp) const;  ///< Returns local fluid mass density of component comp
    double FluidDensity(const int i, const int j, const int k) const;          ///< Returns local fluid mass density summed over all components
    double Density(const PhaseField& Phase, const int i, const int j, const int k) const; ///< Returns local mass density including solids

    Storage3D< D3Q27,    1 > lbPopulations;                                     ///< Populations (discretized particle distribution functions) PDF)
    Storage3D< D3Q27,    1 > lbPopulationsTMP;                                  ///< Temporary array for Populations PDF propagation
    Storage3D< int,      0 > Obstacle;                                          ///< 1 if Node is solid
    Storage3D< int,      0 > ObstacleAppeared;                                  ///< True if obstacle node appeared
    Storage3D< int,      0 > ObstacleChangedDensity;                            ///< True if nearby obstacle changed local density
    Storage3D< int,      0 > ObstacleVanished;                                  ///< True if obstacle node vanished
    Storage3D< dVector3, 1 > ForceDensity;                                      ///< Force density
    Storage3D< dVector3, 1 > MomentumDensity;                                   ///< Momentum density
    Storage3D< double,   1 > DensityWetting;                                    ///< Fluid density / Solid wetting parameter
    Storage3D< double,   1 > nut;                                               ///< kinematic viscosity in lattice units when thermal compressibility considered
    Storage3D< double,   1 > HydroPressure;                                     ///< Hydrolic Pressure for each lattice when thermal compressibility considered
    Storage3D< double,   1 > DivVel;                                            ///< Divergence of velocity when thermal compressibility considered
    Storage3D< dVector3, 1 > GradRho;                                           ///< Gradient of density when thermal compressibility considered

    bool Do_Benzi;                                                              ///< Set to "true" if Benzi force should be calculated
    bool Do_BounceBack;                                                         ///< Set to "true" for the fluid bounce back at the interface (not energy conserving with mobile solids)
    bool Do_BounceBackElastic;                                                  ///< Set to "true" for the elastic fluid bounce back (energy conserving with mobile solids)
    bool Do_Drag;                                                               ///< Set to "true" if drag force at the interfaces should be calculated
    bool Do_EDForcing;                                                          ///< Set to "true" if exact difference method should be used to apply the calculated force
    bool Do_FixPopulations;                                                     ///< Set to "true" if negative Densities and Populations should be fixes
    bool Do_Gravity;                                                            ///< Set to "true" if gravitational force should be calculated
    bool Do_GuoForcing;
    bool Do_Kupershtokh;                                                        ///< Set to "true" if Kupershtokh's method for two phase fluid flow should be used
    bool Do_TwoPhase;                                                           ///< Set to "true" if two phase flow should be calculated
    bool Do_ThermalComp;                                                        ///< Set to "true" if thermal compressibility considered
    bool Do_DI;                                                        			///< Set to "true" if the interface considered diffuse
    bool ObstaclesChanged;                                                      ///< True if an obstacle changed

    bool GradRho_Upwind; 
    bool GradRho_Central; 
    bool GradRho_VanLeer; 

    double Pth;                                                                 ///< Thermodynamic Pressure [Pa]
    double PthOld;                                                              ///< Thermodynamic Pressure for previous time step [Pa]
    double Pth0;                                                                ///< Initial Thermodynamic Pressure [Pa]
    double Poutlet;                                                             ///< Outlet Hydrodynamic Pressure [Pa]

    double U0X;                                                                 ///< Inlet Velocity in x direction  [m/s]
    double U0Y;                                                                 ///< Inlet Velocity in y direction  [m/s]
    double U0Z;                                                                 ///< Inlet Velocity in z direction  [m/s]

    dVector3 GA;                                                                ///< Gravity acceleration in true units [m/s^2]

    double ParaKuper;                                                           ///< Kupershtokh parameter (to tune numerical accuracy!)
    double h_star;                                                              ///< dimensionless drag force parameter

    double dRho;                                                                ///< Coefficient to calculate physical density          from lattice density [Kg/m^3]
    double dP;                                                                  ///< Coefficient to calculate physical pressure         from lattice pressure [Pa]
    double dm;                                                                  ///< Coefficient to calculate physical momentum density from lattice momentum density [Kg/(m^2 s)]
    double df;                                                                  ///< Coefficient to calculate physical force density    from lattice force density [Kg/(m^2 s^2)]
    double dnu;                                                                 ///<  Coefficient to calculate physical kinematic viscosity from lattice kinematic viscosity [m^2/s]
    double dv;
    
    double dt;                                                                  ///< Coefficient to calculate physical time   from lattice time [s]
    double dM;                                                                  ///< Coefficient to calculate physical mass   from lattice mass [Kg]

    double lbWeights[3][3][3];                                                  ///< Lattice Boltzmann stencil weights

    size_t Ncomp;                                                               ///< Number of chemical components
    size_t N_Fluid_Comp;                                                        ///< Number of fluid components
    size_t Nphases;                                                             ///< Number of thermodynamic phases per grain

    GridParameters Grid;                                                        ///< Simulation grid parameters

    int Bcells;                                                                 ///< Number of boundary cells
    int FluidRedistributionRange;                                               ///< Fluid distribution range in case of 2 phase flow

    double cs2;                                                                 ///< Speed of sound squared [m/s]

    double InitialTotalMass;                                                    ///< Solid mass (used for mass conservation)
    std::vector<double> GasParameter;                                           ///< Gas parameter for Kupersthok (interface width)
    std::vector<double> InterfaceWidth;                                         ///< Equilibrium interface width [TODO]
    std::vector<double> SurfaceTension;                                         ///< Surface tension (physical units) [kg/s^2]
    std::vector<double> dDensity_dc;                                            ///< Density change according to the change of concentration [Kg/(mol/mol)]
    std::vector<double> CriticalDensity;                                        ///< Critical density of phase separation (Two phase flow)
    std::vector<double> lbCriticalPressure;                                     ///< Critical Pressure for Kupersthok (Van der Waal)
    std::vector<double> lbCriticalTemperature;                                  ///< Critical temperature of phase separation (Two phase flow)
    std::vector<double> FluidMass;                                              ///< Fluid mass (used for mass conservation)
    std::vector<double> LiquidDensity;                                          ///< Equilibrium density of the liquid phase
    std::vector<double> lbTemperature;                                          ///< Temperature (Two phase flow)
    std::vector<double> VaporDensity;                                           ///< Equilibrium density of the solid phase
    std::vector<double> nu;                                                     ///< Kinematic viscosity [m^2/s]
    std::vector<double> rho_0;                                                  ///< Reference density (Parameter for Benzi) [Kg/m^3]
    std::vector<double> lbtau;                                                  ///< Relaxation time for BGK collision operator
    std::vector<std::string> ElementNames;                                      ///< Names of chemical elements

    std::vector<std::vector<double>> Gb;                                        ///< Parameter for Benzi-force [Kg/(m^3)]
    std::vector<std::vector<double>> Wetting;                                   ///< Wetting Parameter [Kg/m^3]
    std::vector<std::vector<double>> lbGK;                                      ///< Parameter for Kupershtohk

protected:
};

inline double FlowSolverLBM::FluidDensity(const int i, const int j, const int k,
        const size_t comp) const
{
    if (Obstacle(i,j,k)) return 0.0;
    else                 return DensityWetting(i,j,k,{comp});
}

} //namespace openphase
#endif
