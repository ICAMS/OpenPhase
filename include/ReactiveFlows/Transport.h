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

#ifndef TRANSPORT_H_INCLUDED
#define TRANSPORT_H_INCLUDED

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <array>
#include <random>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <tuple>

#include "Initializations.h"

#include "Energy.h"
#include "Species.h"
#include "MixtureFlow.h"
#include "SolidBody.h"

using namespace std;

namespace openphase { 
class MixtureFlow;
class Energy;
class Species;
class SolidBody;
 }
namespace openphase  {

class OP_EXPORTS Transport 
{
    public:
	Transport(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		Initialize(locSettings);
    }
	void Initialize(Settings& locSettings);
	
	void CalculatePhaseFieldCoupling(PhaseField& Phase, Energy &EN, MixtureFlow& MF, double hT, double dt);
	void CalculateSolidDiffusion(PhaseField& Phase, Energy &EN, FlowSolverLBM& FL, 
                                SolidBody& SB, double dt);
	void CalculateAdvectionDiffusion(PhaseField& Phase, Energy &EN, FlowSolverLBM& FL, MixtureFlow &MF,
                                    SolidBody& SB, double dt);
	void CalculateAdvectionDiffusion(PhaseField& Phase, Species &SP, Energy& EN, FlowSolverLBM& FL, MixtureFlow &MF,
                                    SolidBody& SB, double dt);
	void CalculateReaction(Species &SP, Energy& EN, MixtureFlow &MF, FlowSolverLBM& FL, 
                                    int nDim, double &FuelConsumptionRate, size_t GasFuelIndex, double dt);
	
    GridParameters Grid;                                                        ///< Simulation grid parameters

	double CalculatePhiVanLeer(double r);
};
}
#endif
