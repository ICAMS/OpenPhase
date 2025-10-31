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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
 *
 */
#ifndef FLOWFEATURES_H_INCLUDED
#define FLOWFEATURES_H_INCLUDED
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"

#include "Tools/TimeInfo.h"
#include "RunTimeControl.h"

#include "SpeciesTransport.h"
#include "EnergyTransport.h"

using namespace std;

namespace openphase { class SpeciesTransport;
class EnergyTransport;  }

namespace openphase {
class FlowMixture
{
	public:
	FlowMixture(const std::string InputFileName)                                   ///< Constructor
    {
		ReadInput(InputFileName);
    }
	 
	void ReadInput(string InputFile);
	void Initialize(Settings& locSettings, PhaseField& Phase, FlowSolverLBM& FL, Velocities& Vel, BoundaryConditions& BC);
	void InitializingFlowProperties(PhaseField& Phase, FlowSolverLBM& FL,Velocities& Vel, BoundaryConditions& BC);
	void UpdatingInletVelocity(FlowSolverLBM& FL, double ConsumptionRate, double FuelMF);
	void DetectObstacles(FlowSolverLBM& FL, const PhaseField& Phase, bool DI);
	void SetVelocityInletUniformBCX0(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um);
	void SetInletVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um);
	void SetInletUniformVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um);
	void SetInletPoiseuilleVelocity(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um);

	void SetCornerAtInlet(FlowSolverLBM& FL, const PhaseField& Phase);
	void CalculatingMassFlowRate(FlowSolverLBM& FL, PhaseField& Phase, double &Sum, int i);
	void SetVelocityInletPoiseuilleFlow(FlowSolverLBM& FL, const BoundaryConditions& BC, const PhaseField& Phase, double Um);
	void SetInitialVelocity(PhaseField& Phase, BoundaryConditions& BC,FlowSolverLBM& FL, Velocities& Vel, double Um);
	void SettingDivergenceVelocityZero(FlowSolverLBM& FL);

	void Collision(FlowSolverLBM& FL, Velocities& Vel, BoundaryConditions& BC);
	void Propagation(FlowSolverLBM& FL, PhaseField& Phase, BoundaryConditions& BC);

	void LBMLimits(PhaseField& Phase, FlowSolverLBM& FL, double Umax, double lbnu, double dt);
	double CalculateInletArea(PhaseField& Phase, FlowSolverLBM& FL, int InletIndex);

	double CalculateSutherlandViscosity(double Temp);
	double CalculateIdealGasDensity(double p, double Rm, double Temp);
	void UpdateFluidProperties(FlowSolverLBM& FL,EnergyTransport& ET, double mw);
	void CalculateForceDragbyPF(PhaseField& Phase, FlowSolverLBM& FL, const Velocities& Vel, double hHyd);
	void CalculateForceDragbyPFNEW(PhaseField& Phase, FlowSolverLBM& FL, const Velocities& Vel, double hHyd);
	double CalculateRayleighNumber(PhaseField& Phase, FlowSolverLBM& FL, double Pr, double Th, double Tc, double Mw, BoundaryConditions& BC);
	void WriteVTKMixtureVelocity(Settings& locSettings, const int tStep, const int precision = 16);
	void UpdateMixtureVelocity(Velocities& Vel);

	void CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  EnergyTransport& ET, 
												BoundaryConditions& BC, double dt);
	void CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  EnergyTransport& ET,
								SpeciesTransport& ST, BoundaryConditions& BC, double dt);

    Storage3D <dVector3, 0 > V_Mixture;                                         ///< Convective Velocity used to advect the mixture

	bool DI;
	bool DO_YCORNER;
	bool DO_ZCORNER;
	bool SecOrdBB;

	double lbnu;
	double MaxU;
	double LengthScale;
	double InletDensity;
	double InletViscosity;
	double Uavg;
	bool Updating_Velocity;
	bool UniformVel;
	bool PoiseuilleVel;
	double Austar;
	bool Species;
	double InletMassFlowRate;

	double InletArea;
	double ModifiedInletArea;
	double Kp;

};
}
#endif

