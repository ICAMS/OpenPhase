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

 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Raphael Schiedung;
 *                         Johannes Goerler; Marvin Tegeler
 *
 */

#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include "GrainsProperties.h"
#include "Includes.h"
#include "H5Interface.h"

namespace openphase
{

class ElasticProperties;
class BoundaryConditions;
class ChemicalProperties;
class Velocities;
class FlowSolverLBM;
class Settings;
class H5Interface;

class OP_EXPORTS PhaseField : public OPObject                                   ///< Phase field class. It stores the phase fields, performs basic operations on them.
{
 public:

    PhaseField(){};
    PhaseField(Settings& locSettings, std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override;///< Initializes the storage and initial variables of the driving force class.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input from the file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input from the stringstream
    void ReadJSON(const std::string InputFileName);

    void AllocateStorages(GridParameters& Grid);                                ///< Allocates internal storages

    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///< Changes the mesh size while keeping the data.

    void Clear(void);                                                           ///< Clears the phase field storage

    void Refine(void);                                                          ///< Calculates double resolution phase-fields from the single resolution ones

    double CalculateReferenceVolume(double radius);                             ///< Calculates reference volume (RefVolume) for nucleation

    size_t AddGrainInfo(size_t PhaseIndex);                                     ///< Adds new grain information for a phase "PhaseIndex", returns resulting phase field index
    size_t PlantGrainNucleus(size_t PhaseIndex, int x, int y, int z);           ///< Plants grain nucleus at position (x,y,z) of a phase "PhaseIndex", returns resulting phase field index

    dVector3 Normal(NodePF::citerator alpha, NodePF::citerator beta) const;     ///< Returns interface normal for a given pair of phase fields at a given position in double resolution
    //dMatrix3x3 dNormal_dGradientAlpha(NodePF::citerator alpha, NodePF::citerator beta) const;
    //dMatrix3x3 dNormal_dGradientBeta (NodePF::citerator alpha, NodePF::citerator beta) const {return dNormal_dGradientAlpha(alpha,beta)*(-1);};
    NodeAB<dVector3,dVector3> Normals  (const int i, const int j, const int k) const;///< Returns interface normals for all pairs of phase fields at a given position
    NodeAB<dVector3,dVector3> NormalsDR(const int i, const int j, const int k) const;///< Returns interface normals for all pairs of phase fields at a given position in double resolution
    NodeA<dVector3> NormalsPhase  (const int i, const int j, const int k) const;///< Returns interface normals for all pairs of phase fields at a given position
    NodeA<dVector3> NormalsPhaseDR(const int i, const int j, const int k) const;///< Returns interface normals for all pairs of phase fields at a given position

    template<class T>
    double LocalN(const T& locPF) const
    {
        if (!ConsiderNucleusVolume) return locPF.size();
        double norm = 0.0;
        for(auto alpha  = locPF.cbegin();
                 alpha != locPF.cend(); ++alpha)
        {
            norm += FieldsProperties[alpha->index].VolumeRatio;
        }
        return norm;
    };

    void ClearGrainsForcesAndAccelerations()                                    ///< Resets forces and accelerations for each grain in dynamic simulations. It has to be called at the beginning of each time-step when advection is used!
    {
        for (size_t pIdx = 0ul; pIdx < FieldsProperties.size(); pIdx++)
        if(FieldsProperties[pIdx].Exist)
        {
            FieldsProperties[pIdx].ClearForcesAndAccelerations();
        }
    }

    bool PhaseFieldPresent(const int i, const int j, const int k, const size_t Index) const;///< Returns true if phase field is present in a grid cell with coordinates (i,j,k)
    bool ThermodynamicPhasePresent(size_t alpha);
    bool ThermodynamicPhasePairPresent(size_t alpha, size_t beta);

    std::vector<int> GetPresentPhaseFields() const;                             ///< Returns an ordered list with indices of present phase fields
    std::vector<int> VicinityPhaseFields(const int i, const int j, const int k) const; ///< Returns list of phase fields at point (i,j,k) and neighboring nodes

    void Advect(AdvectionHR& Adv, const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dt, const double tStep,
            const bool finalize = true);                                        ///< Advects phase-field when no LB solver is used

    void AdvectALE(AdvectionHR& Adv, const Velocities& Vel,
                const BoundaryConditions& BC,
                const double dt, const double tStep,
                const bool finalize = true);                                    ///< Advects phase-field in the ALE scheme

    void Advect(AdvectionHR& Adv, const Velocities& Vel,
            const BoundaryConditions& BC,
            FlowSolverLBM& LBM, const double dt, const double tStep,
            const bool finalize = true);                                        ///< Advects phase-field when a LB solver is used

    //NodePF FieldsDerivatives  (const int i, const int j, const int k) const;    ///< Calculates local phase-field derivatives
    //NodePF FieldsDerivativesDR(const int i, const int j, const int k) const;    ///< Calculates local phase-field derivatives in double resolution

    //NodePF FieldsLaplacians  (const int i, const int j, const int k) const;     ///< Calculates local phase-field Laplacians
    //NodePF FieldsLaplaciansDR(const int i, const int j, const int k) const;     ///< Calculates local phase-field Laplacians in double resolution
    //NodePF FieldsGradients  (const int i, const int j, const int k) const;      ///< Calculates local phase-field gradients
    //NodePF FieldsGradientsDR(const int i, const int j, const int k) const;      ///< Calculates local phase-field gradients in double resolution

    dMatrix3x3 EffectiveOrientation(const int i, const int j, const int k) const;///< Returns averaged grain rotations, averaging is done via quaternions

    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Set boundary conditions

    void Finalize(const BoundaryConditions& BC, const bool finalize = true);    ///< Finalizing the phase fields calculations

    void FinalizeInitialization(const BoundaryConditions& BC);                  ///< Finalizing the phase fields after initialization

    void FixSpreading(const BoundaryConditions& BC, double cutoff);             ///< Removes phase-fields with values below the cutoff (used for reducing parasitic diffusion in advection)
    void KeepPhaseFieldsVolume(void);                                           ///< Keeps phase fields volume constant by allowing only grain shape change. Should be called before NormalizeIncrements()

    void NormalizeIncrements(const BoundaryConditions& BC, const double dt);    ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.

    void MergeIncrements(const BoundaryConditions& BC,
                         const double dt,
                         const bool finalize = true,
                         const bool clear = true);                              ///< Merges the increments into the phase fields

    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC) override;                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y and or z directions correspondingly.

    void ConsumePlane(const int dx, const int dy, const int dz,
                      const int x, const int y, const int z,
                      const BoundaryConditions& BC);                            ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.

    NodeA<double> Dot (const int i, const int j, const int k, const double dt) const;   ///< Calculates the upcoming derivatives of Fields with respect to time including finalization! (NOTE: Dot = Dot1 + Dot2)
    NodeA<double> Dot1(const int i, const int j, const int k, const double dt) const;   ///< (used for semi-implicit solver Grand Potential solver!)
    NodeA<double> Dot2(const int i, const int j, const int k, const double dt) const;   ///< (used for semi-implicit solver Grand Potential solver!)

    Storage3D<double,1> NewFractions(double dt);                                ///< Phase fractions at next time step
    Tensor<double,1> CalculateNewFractions(Tensor<double,1> oldFractions,
                                           NodeAB<double,double>& FieldsDot,
                                           double locdt);                       ///< Calculates new thermodynamic phase fractions as if Phi.MergeIncrements would have been applied.
    Tensor<double,2> CalculatePsi(NodeAB<double,double>& FieldsDot, double locdt);///< Calculates Phi.FieldsDot for thermodynamic phases.

    Tensor<double,2> FractionsPsi(const int x, const int y, const int z, double locdt) const;///< Calculates Phi.FieldsDot for thermodynamic phases.
    //Tensor<double,1> Fractions(const int x, const int y, const int z) const;    ///< Returns thermodynamic phase fractions
    Tensor<double,1> NewFractions(const int x, const int y, const int z, double locdt) const;///< Calculates new thermodynamic phase fractions as if Phi.MergeIncrements would have been applied.

    //double Fraction(const int x, const int y, const int z, const size_t idx) const;///< Returns fraction of thermodynamic phase idx
    double SolidPhaseFraction(const int i, const int j, const int k) const;     ///< Returns local solid phase fraction
    double SolidMassDensity(const int i, const int j, const int k) const;       ///< Returns local solid mass density
    double CalculateSolidMass(void) const;                                      ///< Calculates and returns SolidMass

    void SelectiveCombinePhaseFields(const BoundaryConditions& BC,              ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
                                     const size_t TargetPFIndex,
                                     const size_t SourcePFIndex);

    std::vector<int> GetMaxPhaseFieldOverlap(const size_t thPhase1,             ///< Returns a vector containing the phase indices thPhase1 and thPhase2 that have the biggest overlap, also the number of interface points between the phasefields. vector contains 1. phasefield of thPhase1, 2. phasefield of thePhase2, 3. number of overlap points. All values are -1 if no overlap is found.
                                             const size_t thPhase2);
    Matrix<int> GetPhaseFieldOverlap(const size_t thPhase1,                     ///< Returns a Matrix containing overlap of phase indices denoted in the columns and rows of the matrix. Only values above the diagonal are set.
                                     const size_t thPhase2);

    void PrintPFVolumes() const;                                                ///< Prints volume of all phase fields
    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints the locally present phase fields and their values to the screen.
    void PrintVolumeFractions(void);                                            ///< Prints volume fraction of all thermodynamic phases
    void WriteVolumePercentages(const std::string& filename, double time,
            char separator = ',') const;                                        ///< Writes volume fraction of all thermodynamic phases to the file FileName

    void WriteGrainsVolume(int time_step, double time, std::string filename);   ///< Writes grains' volume in grid cells over time into a file

    bool Read(std::string FileName);                                            ///< Read raw (binary) phase fields from the file named FileName
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep) override;          ///< Read raw (binary) phase fields from the file of specific time step
    bool ReadH5(const BoundaryConditions& BC, H5Interface& H5, const int tStep);///< Read raw data from HDF5 file
    
	bool WriteMPI(const std::string& Path); 
    bool Write(const std::string& FileName) const;                              ///< Write raw (binary) phase fields to the file FileName
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Write raw (binary) phase fields to the file PhaseField_tStep.dat
    void WriteH5(H5Interface& H5, const int tStep);                             ///< Writes output in HDF5 format

    void WriteAverageVolume(const int tStep, const size_t PhaseIndex) const;    ///< Writes to the file the average volume of a given phase
    void WriteVTK(Settings& locSettings,
                  const int tStep,
                  const bool CurvatureOutput = false,
                  const int precision = 16) const;                              ///< Writes phase fields to the file in VTK format
    void WriteDistortedVTK(Settings& locSettings,
                           const ElasticProperties& EP,
                           const int tStep,
                           const bool CurvatureOutput = false,
                           const int precision = 16) const;                     ///< Writes phase fields on distorted grid to the file in VTK format in MPI environment
    void WriteIndividualPhaseFieldValuesVTK(Settings& locSettings,
                            const int tStep,
                            const std::initializer_list<size_t> FieldIndices,
                            const int precision = 16) const;                    ///< Writes values of individual phase fields to file in VTK format
    void WriteLaplacianVTK(Settings& locSettings, const int tStep,
                           size_t PhiIndex, const int precision = 16);          ///< Writes Laplacian of a given phase field to the file in VTK format

    PhaseField& operator= (const PhaseField& rhs);                              ///< Copy operator for PhaseField class

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases

    double NucleusVolumeFactor;                                                 ///< Scales the critical nucleus volume

    InterfaceNormalModels InterfaceNormalModel;                                 ///< Model selector for the interface normal calculation
    LaplacianStencils PhaseFieldLaplacianStencil;                               ///< Phase-field Laplacian stencil selector
    GradientStencils  PhaseFieldGradientStencil;                                ///< Phase-field gradient stencil selector

    NodeAB<double,double> PairwiseGrowthFactors;                                ///< Pairwise growth reduction factors to deal with stoichiometric and other constraints

    std::vector<std::string> PhaseNames;                                        ///< Names of phases for subset selection
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases
    std::vector<double> PhaseEquilibriumDensities;                              ///< Equilibrium densities for all phases
    std::vector<double> FractionsTotal;                                         ///< Total phase fractions in the entire simulation domain

    bool NucleationPresent;                                                     ///< True if there are nuclei of any phase, false otherwise
    bool ConsiderNucleusVolume;

    LaplacianStencil LStencil;                                                  ///< Laplacian stencil. Uses user specified stencil as the basis
    GradientStencil  GStencil;                                                  ///< Gradient stencil. Uses user specified stencil as the basis

    GrainsProperties FieldsProperties;                                          ///< Phase fields properties. Contains information about location, velocity, orientation, volume, aggregate state etc. of each phase field

    Storage3D< double, 1 > Fractions;                                           ///< Phase fractions storage

    Storage3D< NodePF, 0 > Fields;                                              ///< Phase-field storage
    Storage3D< NodeAB<double,double>, 0 > FieldsDot;                            ///< Phase-field increments storage

    Storage3D< NodeA<double>,  0 > FieldsAdvectionDot;                          ///< Phase-field storage for advection increments
    Storage3D< NodePF, 0 > FieldsAdvectionBackup;                               ///< Phase-field backup storage for advection

    Storage3D< NodePF, 0 > FieldsDR;                                            ///< Phase-field storage for double resolution mode
    Storage3D< NodeAB<double,double>, 0 > FieldsDotDR;                          ///< Phase-field increments storage for double resolution mode
    
    // VTK output helper methods:
    double CurvaturePhase  (const int i, const int j, const int k, const size_t Index) const;///< Curvature for each thermodynamic phase VTK output
    double CurvaturePhaseDR(const int i, const int j, const int k, const size_t Index) const;///< Curvature for each thermodynamic phase VTK output
    double Interfaces  (const int i, const int j, const int k) const;           ///< Indicates interfaces for VTK output
    double InterfacesDR(const int i, const int j, const int k) const;           ///< Indicates interfaces for VTK output in double resolution
    double Variants  (const int i, const int j, const int k) const;             ///< Variants for VTK output
    double VariantsDR(const int i, const int j, const int k) const;             ///< Variants for VTK output in double resolution
    double ParentGrain  (const int i, const int j, const int k) const;          ///< Parent grain index for VTK output
    double ParentGrainDR(const int i, const int j, const int k) const;          ///< Parent grain index for VTK output in double resolution
    std::array<double,2> PrincipalCurvatures(const int i, const int j, const int k,
                                             const size_t phase) const;         ///< Principle Curvatures for plotting VTK
    std::pair<NodeA<double>,NodeA<double>> PrincipalCurvatures(const int i, const int j, const int k) const; ///< Returns principal curvatures of all phase fields at a given position
    // End VTK output helper methods
 protected:
 private:

    void SetStencils(GridParameters& Grid);                                     ///< Sets gradient and Laplacian stencils

    void SetFlagsSR();                                                          ///< Sets the flags which mark interfaces
    void SetFlagsDR();                                                          ///< Sets the flags which mark interfaces in double resolution case

    void CalculateDerivativesSR(void);                                          ///< Calculates local phase-field derivatives
    void CalculateDerivativesDR(void);                                          ///< Calculates local phase-field derivatives in double resolution

    void SetBoundaryConditionsSR(const BoundaryConditions& BC);                 ///< Set boundary conditions in single resolution case
    void SetBoundaryConditionsDR(const BoundaryConditions& BC);                 ///< Set boundary conditions in double resolution case

    void SetIncrementsBoundaryConditionsSR(const BoundaryConditions& BC);       ///< Set boundary conditions for phase field increments
    void SetIncrementsBoundaryConditionsDR(const BoundaryConditions& BC);       ///< Set boundary conditions for phase field increments in double resolution case

    void FinalizeSR(const BoundaryConditions& BC, const bool finalize = true);  ///< Finalizing the phase fields calculations in this time step
    void FinalizeDR(const BoundaryConditions& BC, const bool finalize = true);  ///< Finalizing the phase fields calculations in this time step in double resolutions

    void KeepPhaseVolumeSR(Tensor<bool,2> AllowedTransitions);                  ///< Keeps phases volume constant by allowing only grain shape change and grain transformations between the grains of the same phase (emulates coexistence of immiscible phases). Should be called before NormalizeIncrements()
    void KeepPhaseVolumeDR(Tensor<bool,2> AllowedTransitions);                  ///< Keeps phases volume constant by allowing only grain shape change and grain transformations between the grains of the same phase (emulates coexistence of immiscible phases). Should be called before NormalizeIncrements()

    void NormalizeIncrementsSR(const BoundaryConditions& BC, const double dt);  ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.
    void NormalizeIncrementsDR(const BoundaryConditions& BC, const double dt);  ///< Normalizes the interface fields such that resulting phase fields after merge do not escape interval [0,1] and sum up to 1.

    void MergeIncrementsSR(const BoundaryConditions& BC,
                           const double dt,
                           const bool finalize = true,
                           const bool clear = true);                            ///< Merges the increments into the phase fields in single resolution
    void MergeIncrementsDR(const BoundaryConditions& BC,
                           const double dt,
                           const bool finalize = true,
                           const bool clear = true);                            ///< Merges the increments into the phase fields in double resolution

    void MoveFrameSR(const int dx, const int dy, const int dz,
                     const BoundaryConditions& BC);                             ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y and or z directions correspondingly.
    void MoveFrameDR(const int dx, const int dy, const int dz,
                     const BoundaryConditions& BC);                             ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -2 or +2) in x, y and or z directions correspondingly.

    void ConsumePlaneSR(const int dx, const int dy, const int dz,
                        const int x, const int y, const int z,
                        const BoundaryConditions& BC);                          ///< Shifts the data in the single resolution storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.
    void ConsumePlaneDR(const int dx, const int dy, const int dz,
                        const int x, const int y, const int z,
                        const BoundaryConditions& BC);                          ///< Shifts the data in the double resolution storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.

    void SelectiveCombinePhaseFieldsSR(const BoundaryConditions& BC,            ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
                                       const size_t TargetPFIndex,
                                       const size_t SourcePFIndex);
    void SelectiveCombinePhaseFieldsDR(const BoundaryConditions& BC,            ///< Merge phase fields of sourceIndex to phase field with index targetPhaseIndex
                                       const size_t TargetPFIndex,
                                       const size_t SourcePFIndex);

    void Coarsen(void);                                                         ///< Calculates single resolution phase-fields from the double resolution ones
    void CoarsenDot(void);                                                      ///< Calculates single resolution phase-field increments from the double resolution ones

    void CombinePhaseFields(void);                                              ///< Merge phase fields of same phase to phase field with index PhaseIndex
    std::vector<bool> Combine;                                                  ///< Indicates which phase's phase fields should be combined

    void CalculateFractions(void);                                              ///< Calculates phase fractions from phase-fields, populates Fractions storage
    void CalculateGrainsVolume(void);                                           ///< Collects volume for each phase field.

    void Advect(AdvectionHR& Adv, const Velocities& Vel,
                PhaseField& Phi, const BoundaryConditions& BC,
                const double dt, const double tStep) override                   ///< Hidden standard advection method. Should not be used for PhaseField
    {
        (void) Adv; //unused
        (void) Vel; //unused
        (void) Phi; //unused
        (void) BC; //unused
        (void) dt; //unused
        (void) tStep; //unused
    };
};

}// namespace openphase
#endif

