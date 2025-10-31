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

 *   File created :   2016
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#ifndef ELASTICITYSOLVERSPECTRAL_H
#define ELASTICITYSOLVERSPECTRAL_H

/******************************************************************************
*  This module is based on the iterative algorithm of [S.Y.Hu, L.Q.Chen,      *
*  Acta. Mater. 49 (2001) 1879 - 1890] with modifications allowing to use     *
*  inhomogeneous elasticity parameters not limited to composition dependence. *
*  The boundary condition allowing free volume expansion is also introduced   *
*  following the approach outlined in [I.Steinbach, M.Apel, Physica D 217     *
*  (2006) 153 - 160].                                                         *
*  Finite strain extension considers Green-Lagrangian strain within the       *
*  St. Venant-Kirchhoff hyper-elastic model. The resulting solver algorithm   *
*  is similar to [P. Eisenlohr et al., International Journal of Plasticity 46 *
*  (2013) 37–53]                                                              *
******************************************************************************/

#ifdef MPI_PARALLEL
#include "mpi_wrapper.h"
#else
#include "fftw3.h"
#endif
#include "Includes.h"

namespace openphase
{
class Settings;
class ElasticProperties;
class ElasticitySolverSpectralImpl;

class OP_EXPORTS ElasticitySolverSpectral : public OPObject                     ///< Elastic problem solver based on spectral method.
{
 public:
    ElasticitySolverSpectral();
    ElasticitySolverSpectral(Settings& locSettings,
                        const std::string InputFileName = DefaultInputFileName);///< Constructor calls Initialize and ReadInput
    ElasticitySolverSpectral(Settings& locSettings,
                        const BoundaryConditions& BC,
                        const std::string InputFileName = DefaultInputFileName);///< Constructor calls Initialize and ReadInput
    ~ElasticitySolverSpectral(void);
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Named constructor, allocates internal storages and initialized variables
    void Initialize(Settings& locSettings, const BoundaryConditions& BC, std::string ObjectNameSuffix = "");       ///< Named constructor (considering non periodic boundary conditions)
    void ReadInput(const std::string InputFileName) override;                   ///< Reads elastic properties from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads elastic properties from the input stream
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the storages, reinitilizes the solver

    void Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC);     ///< Needed if the system dimensions have been altered (considering non periodic boundary conditions)

    int  Solve(ElasticProperties& EP, BoundaryConditions& BC, double dt,
               std::function<bool()> EndSolve = [](){return false;});           ///< Solves mechanical equilibrium problem
 private:
    std::unique_ptr<ElasticitySolverSpectralImpl> impl_;
};
} // namespace openphase
#endif //SpectralElasticitySolver
