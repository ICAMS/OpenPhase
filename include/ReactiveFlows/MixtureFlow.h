/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2025
 *  Main contributors :   Oleg Shchyglo; Reza Namdar
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

#include "Species.h"
#include "Energy.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { class Species;
class Energy;  }
enum class PressureOutletMode{Fixed, NonReflecting};
enum class ExtrapolateOrder{ZeroOrder, FirstOrder};

namespace openphase {
class OP_EXPORTS MixtureFlow
{
	public:
	MixtureFlow(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		ReadInput(InputFileName);
		Initialize(locSettings);
    }
	void Initialize(Settings& locSettings);
	void ReadInput(string InputFile);
	void SetInitial(PhaseField &Phase, FlowSolverLBM& FL, BoundaryConditions& BC);
	void SetFlowProperties(PhaseField& Phase, FlowSolverLBM& FL, BoundaryConditions& BC);
	void UpdateIncomingVelocityValue(FlowSolverLBM& FL, double ConsumptionRate, double FuelMF);
	void DetectObstacles(FlowSolverLBM& FL, const PhaseField& Phase, BoundaryConditions &BC, bool DI);
	void SetVanishedObstacleNodes(PhaseField &Phase, FlowSolverLBM& FL, SolidBody &SB, Energy& EN, Species& SP);		
	void SetAppeardObstacleNodes(FlowSolverLBM& FL, SolidBody &SB);
	void SetInitialPopulations(FlowSolverLBM& FL, BoundaryConditions& BC);

	void SetPressureBoundary(FlowSolverLBM& FL, double Pout, 
							 PressureOutletMode mode,   // Fixed or NonReflecting
    						 const dVector3& A,         // start point (INDEX space)
    						 const dVector3& B,         // end point   (INDEX space)
    						 int n );     

	void SetVelocityBoundary(FlowSolverLBM& FL, const dVector3& UinMean,
		    				 const dVector3& A,         // start point (GLOBAL index space)
    						 const dVector3& B,         // end point   (GLOBAL index space)
    						 int n,                     // +1 or -1 : direction INTO the domain
    						 bool usePoiseuille         // false = uniform, true = Poiseuille along A->B
							 );
	void SetOpenBoundary(FlowSolverLBM& FL, ExtrapolateOrder order, const dVector3& A, const dVector3& B, int n);
	void SetInitialVelocity(PhaseField& Phase, BoundaryConditions& BC,FlowSolverLBM& FL);

	void Collision(FlowSolverLBM& FL, BoundaryConditions &BC);
	void Propagation(FlowSolverLBM& FL, PhaseField& Phase, SolidBody& SB, BoundaryConditions& BC);

	void LBMLimits(PhaseField& Phase, FlowSolverLBM& FL, double Umax, double lbnu, double dt);


	void CalculateForceDragbyPF(PhaseField& Phase, FlowSolverLBM& FL, double hHyd);
	void WriteVTK(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision = 16);
	void WriteVTKDivergenceVelocity(Settings& locSettings, FlowSolverLBM& FL, const int tStep, const int precision = 16);
	double CalculateMassFlowRate(const FlowSolverLBM& FL, const dVector3& A, const dVector3& B, int nSign);
	void CalculateDensityGradient(FlowSolverLBM &FL, Energy &EN);
	void CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  Energy& EN, SolidBody &SB, double dt);
	void CalculateDivergenceVelocity(PhaseField& Phase, FlowSolverLBM& FL,  Energy& EN, 
									Species& SP, SolidBody &SB, double dt);
	void PropagationSecondOrderBB(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB, const BoundaryConditions& BC);
	void ApplyForces(PhaseField& Phase, FlowSolverLBM& FL);
	void CalculateForceGravity(PhaseField& Phase, FlowSolverLBM &FL);
	void CalculateVelocityAndPressure(FlowSolverLBM& FL);
	double SecondOrderBounceBack(const int i, const int j, const int k, const int ii, const int jj, const int kk, const size_t n,
            PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB, const BoundaryConditions& BC, double& lbDensityChange);    ///< BounceBack at Solid Interfaces

    Storage3D <dVector3, 0 > VelocityMix;                                         ///< Convective Velocity used to advect the mixture
	
    GridParameters Grid;                                                        ///< Simulation grid parameters

	bool DI;
	bool DO_YCORNER;
	bool DO_ZCORNER;
	bool SecOrdBB;

    Storage3D <double, 0 > DensityMix;                                        
    Storage3D <dVector3, 0 > GradDensity;                                        
    Storage3D <double, 0 > ViscosityMix;                                        


	double lbnu;
	double MaxU;
	double LengthScale;
	double ReynoldsNumber;
	double InletDensity;
	double InletViscosity;
	double Uavg;

	bool UPDATE_VELOCITY;
	bool POISEFLOW;

	double Austar;
	double InletMassFlowRate;

	double InletArea;
	double ModifiedInletArea;
	double Kp;

	double TrackingValue=0.0;

};
}
#endif

