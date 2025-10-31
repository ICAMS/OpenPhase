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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Raphael Schiedung
 *
 */

#ifndef INTERFACEPROPERTIES_H
#define INTERFACEPROPERTIES_H

#include "Includes.h"
#include "InterfaceProperties/InterfaceEnergyModel.h"
#include "InterfaceProperties/InterfaceMobilityModel.h"
#include "Thermodynamics/Element.h"
#include "Thermodynamics/ThermodynamicPhase.h"

#include <set>

namespace openphase
{
class ElasticProperties;
class PhaseField;
class Settings;
class Temperature;
class DoubleObstacle;
class DrivingForce;
enum class TripleJunctionModels                                                 ///< Models for triple junction term calculation
{
    Model,                                                                      ///< Factorised sum of pair interface energies
    Value,                                                                      ///< User specified triple junction energy value
    None                                                                        ///< Triple junction term is not considered
};

class OP_EXPORTS InterfaceProperties : public OPObject                          ///< Interface properties module.
{
    friend DoubleObstacle;
    friend DrivingForce;

 public:

    InterfaceProperties(){};
    InterfaceProperties(Settings& locSettings,
                        const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes all variables, allocates storages
    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remeshes/reallocates the storage
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input from a file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input from a stringstream
    void ReadJSON(const std::string InputFileName);

    void SetBoundaryConditions(const BoundaryConditions& BC) override;
    void SetBoundaryConditionsDR(const BoundaryConditions& BC);

    void Coarsen(const PhaseField& Phase);                                      ///< Populates single resolution storage using the double resolution values
    void Refine(const PhaseField& Phase);                                       ///< Populates double resolution storage using the single resolution values
    void Set(const PhaseField& Phase, const BoundaryConditions& BC);            ///< Sets both, interface energy and mobility

    void Set(const PhaseField& Phase,
             const Temperature& Tx,
             const BoundaryConditions& BC);                                     ///< Sets both, interface energy and mobility considering temperature effect
             
    void SetStoichiometry( std::vector<ThermodynamicPhase> Phase);

    void WriteFacetArea(const PhaseField& Phase, double& tStep,
                        double& DegreeTolerance, std::string OutFile);

    void WriteVTK(const Settings& locSettings,
                  const int tStep,
                  const int precision = 16);                                    ///< Writes effective interface properties in VTK format

    double ReportMaximumTimeStep(void);                                         ///< Returns maximum time step for phase-field equation

    Storage3D< NodeIP, 0 > Properties;                                          ///< 3D interface properties storage
    Storage3D< NodeIP, 0 > PropertiesDR;                                        ///< 3D interface properties storage
    Matrix<InterfaceEnergyModel> InterfaceEnergy;                               ///< Interface energy model storage for all pairs of phases
    Matrix<InterfaceMobilityModel> InterfaceMobility;                           ///< Interface mobility model storage for all pairs of phases

    Table< bool > RespectParentBoundaries;                                      ///< Prevents growth of product phase beyond parent grain boundaries

    //Table <bool > PhaseInteractions;                                            ///< Indicates which phase pairs are allowed to interact
    //std::vector <ActivationModes> PhaseActivationModes;                         ///< Phase activation modes. The phase can be enabled or disabled
    
    Matrix<double> maxEnergies;                                                 ///< Maximum interface energy for each phase pair in the simulation
    Matrix<double> maxMobilities;                                               ///< Maximum interface mobility for each phase pair in the simulation

  protected:

    dVector3 dEnergy_dGradientAlpha(const PhaseField& Phase,
                            NodePF::citerator alpha,
                            NodePF::citerator beta) const;                      ///< Partial derivative of interface energy/Set function with respect to the gradient of alpha

    TripleJunctionModels TripleJunctionModel;                                   ///< Triple junction model
    double TripleJunctionEnergy;                                                ///< Triple junction energy term
    double TripleJunctionFactor;                                                ///< Parameter in front of the triple junction energy term

    double CurvatureFactor;                                                     ///< Curvature reduction parameter. It should have value between 0 (no curvature) to 1 (full curvature -> default)
    double RegularizationFactor;                                                ///< Increases the interface stabilizing force. It should have value >= 0.0. This parameter is an extension of CurvatureFactor

    Matrix<double> NucleiMobilityFactor;                                        ///< Reduction of mobilities for growing nuclei

    Storage3D<NodeAB<dVector3,dVector3>, 0> InterfaceStiffnessTMP;              ///< Divergence of this will be the interface stiffness driving force

 private:

    size_t Nphases;

    GridParameters Grid;                                                        ///< Simulation grid parameters

    bool   FullAnisotropy;                                                      ///< True if full interface energy anisotropy model is used, false otherwise

    void SetSR(const PhaseField& Phase);                                        ///< Sets both, interface energy and mobility
    void SetDR(const PhaseField& Phase);                                        ///< Sets both, interface energy and mobility

    void SetMobilityThermalEffectSR(const PhaseField& Phase,
                                    const Temperature& Tx);                     ///< Sets interface mobility in the entire simulation domain using user selected interface mobility models and considering temperature effect
    void SetMobilityThermalEffectDR(const PhaseField& Phase,
                                    const Temperature& Tx);                     ///< Sets interface mobility in the entire simulation domain using user selected interface mobility models and considering temperature effect
    enum class ExtrapolationModes { None, ZerothOrder, FirstOrder};             ///< Curvature extrapolation mode
    ExtrapolationModes ExtrapolationMode;                                       ///< Interface energy extrapolation mode
    size_t MaxExtrapolationInterations;                                         ///< Maximum number of extrapolation iterations
    Storage3D< NodeIP, 0 > PropertiesExtrapolations;                            ///< Temporary storage for extrapolated properties
    std::set<size_t> localPhaseIdx( const int i, const int j, const int k, const PhaseField& Phase) const; ///< Returns the phase indices next to the local grid point
    void Extrapolate(const PhaseField& Phase, const BoundaryConditions& BC);    ///< Extrapolates the interface energy
    bool FirstOrderExtrapolation(const PhaseField& Phase);                      ///< Extrapolates the interface energy using first order approximation
};

}// namespace openphase
#endif
