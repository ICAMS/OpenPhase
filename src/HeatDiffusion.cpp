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
 *   Main contributors :   Oleg Shchyglo; Philipp Engels; Raphael Schiedung;
 *                         Marvin Tegeler; Helge Schaar
 *
 */

#include "HeatDiffusion.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Composition.h"

namespace openphase
{
using namespace std;

inline pair<int, int> three_state_bounds_selector(int selector, int lower_bound, int upper_bound)
{
    pair<int, int> result;
    switch (selector)
    {
        case -1: // lower bound only
        {
            result = make_pair(lower_bound, lower_bound + 1);
            break;
        }
        case 0: // both bounds
        {
            result = make_pair(lower_bound, upper_bound);
            break;
        }
        case 1: // upper bound only
        {
            result = make_pair(upper_bound - 1, upper_bound);
            break;
        }
    }
    return result;
}

void HeatDiffusion::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "HeatDiffusion";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    Tolerance = 1.0e-9;
    MaxIterations = 10000;
    SolverCallsInterval = 1;
    SolverCallsCounter = 0;
    VerboseIterations = false;

    PhaseThermalConductivity.Allocate(Nphases);
    PhaseVolumetricHeatCapacity.Allocate(Nphases);

    size_t Bcells = Grid.Bcells;
    EffectiveThermalConductivity.Allocate(Grid, Bcells);
    EffectiveHeatCapacity.Allocate(Grid, Bcells);

    TxOld.Allocate(Grid, Bcells);
    dTx.Allocate(Grid, Bcells);
    Qdot.Allocate(Grid, Bcells);

    locSettings.AddForRemeshing(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void HeatDiffusion::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("HeatDiffusion input");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();
    ReadInput(inp);
}

void HeatDiffusion::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    for(size_t n = 0; n < Nphases; n++)
    {
        stringstream converter;
        converter << n;
        string counter = converter.str();

        ConsoleOutput::WriteBlankLine();

        PhaseThermalConductivity[n]    = FileInterface::ReadParameterD(inp, moduleLocation, string("ThermalConductivity_") + counter);
        PhaseVolumetricHeatCapacity[n] = FileInterface::ReadParameterD(inp, moduleLocation, string("VolumetricHeatCapacity_") + counter, false, 0.0);
    }

    Tolerance       = FileInterface::ReadParameterD(inp, moduleLocation, string("Tolerance"), false, Tolerance);
    MaxIterations   = FileInterface::ReadParameterI(inp, moduleLocation, string("MaxIterations"), false, MaxIterations);
    SolverCallsInterval = FileInterface::ReadParameterI(inp, moduleLocation, string("SolverCallsInterval"), false, SolverCallsInterval);
    VerboseIterations = FileInterface::ReadParameterB(inp, moduleLocation, string("VerboseIterations"), false, VerboseIterations);

    MaxThermalConductivity = 0.0;
    for(size_t n = 0; n < Nphases; n++)
    {
        MaxThermalConductivity = max(MaxThermalConductivity, PhaseThermalConductivity[n]);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void HeatDiffusion::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    EffectiveThermalConductivity.Reallocate(newNx, newNy, newNz);
    EffectiveHeatCapacity.Reallocate(newNx, newNy, newNz);
    TxOld.Reallocate(newNx, newNy, newNz);
    dTx.Reallocate(newNx, newNy, newNz);
    Qdot.Reallocate(newNx, newNy, newNz);

    Grid.SetDimensions(newNx, newNy, newNz);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void HeatDiffusion::SetLocalLatentHeat(const PhaseField& Phase, const Temperature& Tx)
{
    /** This function accounts for the latent heat release due to
        phase-transformations. */

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Qdot,0,)
    if(Phase.Fields(i,j,k).interface())
    {
        for(auto it  = Phase.FieldsDot(i,j,k).cbegin();
                 it != Phase.FieldsDot(i,j,k).cend(); ++it)
        if(fabs(it->value1) > DBL_EPSILON)
        {
            size_t PIdxA = Phase.FieldsProperties[it->indexA].Phase;
            size_t PIdxB = Phase.FieldsProperties[it->indexB].Phase;

            Qdot(i,j,k) += Tx.LatentHeat({PIdxA,PIdxB})*it->value1;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::SetLocalLatentHeatAndApply(const PhaseField& Phase, Temperature& Tx, double dt)
{
    /** This function accounts for the latent heat release due to
        phase-transformations. */

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Qdot,0,)
    if(Phase.Fields(i,j,k).interface())
    {
        for(auto it  = Phase.FieldsDot(i,j,k).cbegin();
                 it != Phase.FieldsDot(i,j,k).cend(); ++it)
        if(fabs(it->value1) > DBL_EPSILON)
        {
            size_t PIdxA = Phase.FieldsProperties[it->indexA].Phase;
            size_t PIdxB = Phase.FieldsProperties[it->indexB].Phase;

            Qdot(i,j,k) += Tx.LatentHeat({PIdxA,PIdxB})*it->value1;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Qdot,0,)
    if(Phase.Fields(i,j,k).interface())
    {
        Tx(i,j,k) += dt*Qdot(i,j,k)/EffectiveHeatCapacity(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Qdot.Clear();
}
void HeatDiffusion::SetEffectiveProperties(const PhaseField& Phase, const Temperature& Tx)
{
    if(Tx.LatentHeatMode == LatentHeatModes::Local)
    {
        SetLocalLatentHeat(Phase, Tx);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveThermalConductivity,EffectiveThermalConductivity.Bcells(),)
    {
        EffectiveHeatCapacity(i,j,k) = 0.0;
        EffectiveThermalConductivity(i,j,k) = 0.0;

        for(size_t alpha = 0; alpha < Phase.Nphases; ++alpha)
        if(Phase.Fractions(i,j,k,{alpha}) > 0.0)
        {
            EffectiveHeatCapacity(i,j,k)
            += PhaseVolumetricHeatCapacity[alpha]*Phase.Fractions(i,j,k,{alpha});

            EffectiveThermalConductivity(i,j,k)
            += PhaseThermalConductivity[alpha]*Phase.Fractions(i,j,k,{alpha});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::Perform1DIteration(Temperature& Tx,
                                       Temperature1Dextension& TxExt,
                                       BoundaryConditionTypes extBC,
                                       double& residual, double dt)
{
    /** Calculates one iteration of the heat diffusion using a Jacobi implicit
        algorithm for different 1D-Extensions:

        double dT = (RhoCp*TxOld[x]+Lambda*dt/(dx*dx)*(Tx[x-1]+Tx[x+1]))
                   /(RhoCp + 2.0*Lambda*dt/(dx*dx)) - Tx[x];
    */
    std::pair<int, int> Xbounds = three_state_bounds_selector(TxExt.Direction[0], 0, Tx.Grid.Nx);
    std::pair<int, int> Ybounds = three_state_bounds_selector(TxExt.Direction[1], 0, Tx.Grid.Ny);
    std::pair<int, int> Zbounds = three_state_bounds_selector(TxExt.Direction[2], 0, Tx.Grid.Nz);

    double boundaryValue = 0.0;
    double averageRhoCP  = 0.0;
    double averageLambda = 0.0;
    double area          = 0.0;

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        averageRhoCP  += EffectiveHeatCapacity(i,j,k);                          // Average volumetric heat capacity [J/(m^3 K)]
        averageLambda += EffectiveThermalConductivity(i,j,k);                   // Average thermal conductivity [J/(m s K)]
        boundaryValue += Tx(i,j,k);
        area++;
    }

#ifdef MPI_PARALLEL
    if(TxExt.Direction[0] == 0)
    {
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &averageRhoCP, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &averageLambda, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &boundaryValue, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &area, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
#endif

    averageRhoCP  /= area;
    averageLambda /= area;

    TxExt.Data(-1) = boundaryValue/area;

    double dt_dx2 = dt/(Grid.dx*Grid.dx);

    for (size_t i = 0; i < TxExt.size(); ++i)
    {
        double locQ = 0.0;
        if(i == TxExt.size()-1)
        {
            locQ = TxExt.Qdot;
        }

        double loc_deltaT = (averageRhoCP*TxExt.DataTMP(i) +
                             averageLambda*dt_dx2*(TxExt.Data(i+1) + TxExt.Data(i-1)) + locQ*dt)
                           /(averageRhoCP + 2.0*averageLambda*dt_dx2) - TxExt.Data(i);

        residual = max(residual,loc_deltaT*loc_deltaT);
        TxExt.Data(i) += loc_deltaT;
    }

    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        Tx(i+TxExt.Direction[0],j+TxExt.Direction[1],k+TxExt.Direction[2]) = TxExt.Data(0);
    }

    TxExt.setBC(extBC);
}

int HeatDiffusion::SolveImplicit(const PhaseField& Phase,
                                 const BoundaryConditions& BC,
                                 Temperature& Temp,
                                 const double dt_in)
{
    int iteration = 0;                                                          // Iterations counter
    double dt = dt_in;
    if(SolverCallsCounter != 0) dt = SolverCallsInterval*dt_in;

    if((SolverCallsCounter % SolverCallsInterval) == 0)
    {
        /** This function calculates the heat conductivity in the simulation domain
        using an implicit solver method. */

        /* Populate "old" temperature for the iterations. */

        if(Temp.ExtensionsActive)
        {
            if(Temp.ExtensionX0.isActive()) {Temp.ExtensionX0.store_temporary();}
            if(Temp.ExtensionXN.isActive()) {Temp.ExtensionXN.store_temporary();}
            if(Temp.ExtensionY0.isActive()) {Temp.ExtensionY0.store_temporary();}
            if(Temp.ExtensionYN.isActive()) {Temp.ExtensionYN.store_temporary();}
            if(Temp.ExtensionZ0.isActive()) {Temp.ExtensionZ0.store_temporary();}
            if(Temp.ExtensionZN.isActive()) {Temp.ExtensionZN.store_temporary();}
        }
        Temp.CalculateMinMaxAvg();
        TxOld = Temp.Tx;

        double residual = 0.0;                                                      // Residual for calculating the accuracy of the converged solution
        const double dt_dx2 = dt/(Grid.dx*Grid.dx);
        const double dimension = 2.0*double(Grid.Active());                       // Laplacian stencil dimension parameter

        do
        {
            iteration++;
            residual = 0.0;

            if(Temp.ExtensionsActive)
            {
                /* Calculation of heat diffusion in 1D extension using Gauss–Seidel implicit method. */
    #ifdef MPI_PARALLEL
                if(MPI_RANK == 0)
    #endif
                if (Temp.ExtensionX0.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionX0, BC.BC0X, residual, dt);
                }
    #ifdef MPI_PARALLEL
                if(MPI_RANK == MPI_SIZE - 1)
    #endif
                if (Temp.ExtensionXN.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionXN, BC.BCNX, residual, dt);
                }
                if (Temp.ExtensionY0.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionY0, BC.BC0Y, residual, dt);
                }
                if (Temp.ExtensionYN.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionYN, BC.BCNY, residual, dt);
                }
                if (Temp.ExtensionZ0.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionZ0, BC.BC0Z, residual, dt);
                }
                if (Temp.ExtensionZN.isActive())
                {
                    Perform1DIteration(Temp, Temp.ExtensionZN, BC.BCNZ, residual, dt);
                }
            }

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,reduction(max:residual))
            {
                /* Calculation of heat diffusion using Jacobi implicit method. */

                double RhoCp      = EffectiveHeatCapacity(i,j,k);                   // Volumetric heat capacity [J/(m^3 K)]
                double Lambda     = EffectiveThermalConductivity(i,j,k);            // Thermal conductivity [J/(m s K)]
                double locQdot    = Qdot(i,j,k);                                    // Heat source [J/(m^3 s)]
                double locStencil = 0.0;                                            // Temperature stencil from the updated field [K]

                if (Grid.dNx > 0) locStencil += (Temp(i+1,j,k) + Temp(i-1,j,k));
                if (Grid.dNy > 0) locStencil += (Temp(i,j+1,k) + Temp(i,j-1,k));
                if (Grid.dNz > 0) locStencil += (Temp(i,j,k+1) + Temp(i,j,k-1));

                dTx(i,j,k) = (RhoCp*TxOld(i,j,k) + Lambda*locStencil*dt_dx2 + locQdot*dt)
                            /(RhoCp + dimension*Lambda*dt_dx2) - Temp(i,j,k);

                residual = max(residual,dTx(i,j,k)*dTx(i,j,k));
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            /* Updating the temperature with the calculated increment.*/

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,)
            {
                Temp(i,j,k) += dTx(i,j,k);
                dTx(i,j,k) = 0.0;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            /* Assigning values for the boundary cells. */
            Temp.SetBoundaryConditions(BC);

            /* Warning output, if number of iterations exceeds MaxIterations. */

            residual = sqrt(residual)/(Temp.Tavg + Tolerance); // Added Tolerance in denominator in case if Tavg == 0

    #ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &residual, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    #endif

            if(VerboseIterations)
            {
                std::string message  = ConsoleOutput::GetStandard("Solver call",
                                                 std::to_string(SolverCallsCounter));
                            message += ConsoleOutput::GetStandard("Iteration",
                                                 std::to_string(iteration));
                            message += ConsoleOutput::GetStandard("Residual",
                                       ConsoleOutput::to_string_with_precision(residual));
                            message += ConsoleOutput::GetStandard("Tolerance",
                                       ConsoleOutput::to_string_with_precision(Tolerance));
                ConsoleOutput::WriteWithinMethod(message, thisclassname, "SolveImplicit()");
            }

            if(iteration == MaxIterations)
            {
                stringstream value;
                value.precision(16);
                value << residual;

                string tmp = "Maximum number of " + to_string(MaxIterations) + " iterations reached.\n"
                           + "Solver failed to converge to the requested tolerance of " + to_string(Tolerance) + ".\n"
                           + "Current convergence residual is " + value.str() + ".\n"
                           + "Simulation will continue!\n";

                ConsoleOutput::WriteWarning(tmp,thisclassname,"SolveImplicit()");
                break;
            }
        }
        while (residual >= Tolerance);

        /* Evaluate temperature statistics. */
        Temp.CalculateMinMaxAvg();

        if(Temp.ExtensionsActive)
        {
            if(Temp.ExtensionX0.isActive()) {Temp.ExtensionX0.Qdot = 0.0;}
            if(Temp.ExtensionXN.isActive()) {Temp.ExtensionXN.Qdot = 0.0;}
            if(Temp.ExtensionY0.isActive()) {Temp.ExtensionY0.Qdot = 0.0;}
            if(Temp.ExtensionYN.isActive()) {Temp.ExtensionYN.Qdot = 0.0;}
            if(Temp.ExtensionZ0.isActive()) {Temp.ExtensionZ0.Qdot = 0.0;}
            if(Temp.ExtensionZN.isActive()) {Temp.ExtensionZN.Qdot = 0.0;}
        }
        Qdot.Clear();
    }
    SolverCallsCounter++;
    return iteration;
}
void HeatDiffusion::SolveSemiImplicit(const PhaseField& Phase,
                                      const BoundaryConditions& BC,
                                      Temperature& Temp,
                                      const double dt)
{
    /** Solves the heat diffusion equation using a semi-implicit scheme with the
        Thomas Algorithm (Tridiagonal Matrix Algorithm, TDMA).

        The equation is discretized as:
            RhoCp * (T^{n+1} - T^n) / dt = Lambda * Laplacian(T^{n+1}) + Qdot

        In each directional sweep the operator along that axis is treated
        implicitly while contributions from the remaining axes are taken
        explicitly from the current iterate.  Three sweeps (X → Y → Z)
        advance the solution one full time step.

        For a sweep along X at row (j, k) the per-node equations are:
            a_i * T_{i-1} + b_i * T_i + c_i * T_{i+1} = rhs_i

        where:
            r      = Lambda * dt / (RhoCp * dx^2)
            a_i    = -r
            b_i    =  1 + 2*r          (1 + 2*r for interior, adjusted at boundaries)
            c_i    = -r
            rhs_i  = T_old_i + (explicit Y/Z contributions) * r + Qdot * dt / RhoCp
    */

    if ((SolverCallsCounter % SolverCallsInterval) != 0)
    {
        SolverCallsCounter++;
        return;
    }

    const double dx2 = Grid.dx * Grid.dx;

    // --- X-direction sweep: implicit along X, explicit along Y and Z ---
    if (Grid.dNx > 0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,)
        {
            TxOld(i,j,k) = Temp(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int j = 0; j < Grid.Ny; ++j)
        for (int k = 0; k < Grid.Nz; ++k)
        {
            const int N = Grid.Nx;
            std::vector<double> a(N), b(N), c(N), rhs(N);
            std::vector<double> cp(N), dp(N);

            for (int i = 0; i < N; ++i)
            {
                double RhoCp  = EffectiveHeatCapacity(i,j,k);
                double Lambda = EffectiveThermalConductivity(i,j,k);
                double r      = Lambda * dt / (RhoCp * dx2);

                double explicitContrib = 0.0;
                if (Grid.dNy > 0) explicitContrib += Lambda * (TxOld(i,j+1,k) + TxOld(i,j-1,k)) / dx2;
                if (Grid.dNz > 0) explicitContrib += Lambda * (TxOld(i,j,k+1) + TxOld(i,j,k-1)) / dx2;

                a[i]   = -r;
                b[i]   =  1.0 + 2.0 * r;
                c[i]   = -r;
                rhs[i] = TxOld(i,j,k) + (explicitContrib + Qdot(i,j,k)) * dt / RhoCp;
            }

            // Dirichlet-like boundary adjustment at the first and last node
            // (zero-flux: ghost values equal boundary node → coefficient folds into diagonal)
            a[0]     = 0.0;
            c[N - 1] = 0.0;

            // Forward sweep
            cp[0] = c[0] / b[0];
            dp[0] = rhs[0] / b[0];
            for (int i = 1; i < N; ++i)
            {
                double denom = b[i] - a[i] * cp[i - 1];
                cp[i] = c[i] / denom;
                dp[i] = (rhs[i] - a[i] * dp[i - 1]) / denom;
            }

            // Back substitution
            Temp(N - 1, j, k) = dp[N - 1];
            for (int i = N - 2; i >= 0; --i)
            {
                Temp(i, j, k) = dp[i] - cp[i] * Temp(i + 1, j, k);
            }
        }
        Temp.SetBoundaryConditions(BC);
    }

    // --- Y-direction sweep: implicit along Y, explicit along X and Z ---
    if (Grid.dNy > 0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,)
        {
            TxOld(i,j,k) = Temp(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) collapse(2)
        #endif
        for (int i = 0; i < Grid.Nx; ++i)
        for (int k = 0; k < Grid.Nz; ++k)
        {
            const int N = Grid.Ny;
            std::vector<double> a(N), b(N), c(N), rhs(N);
            std::vector<double> cp(N), dp(N);

            for (int j = 0; j < N; ++j)
            {
                double RhoCp  = EffectiveHeatCapacity(i,j,k);
                double Lambda = EffectiveThermalConductivity(i,j,k);
                double r      = Lambda * dt / (RhoCp * dx2);

                double explicitContrib = 0.0;
                if (Grid.dNx > 0) explicitContrib += Lambda * (TxOld(i+1,j,k) + TxOld(i-1,j,k)) / dx2;
                if (Grid.dNz > 0) explicitContrib += Lambda * (TxOld(i,j,k+1) + TxOld(i,j,k-1)) / dx2;

                a[j]   = -r;
                b[j]   =  1.0 + 2.0 * r;
                c[j]   = -r;
                rhs[j] = TxOld(i,j,k) + (explicitContrib + Qdot(i,j,k)) * dt / RhoCp;
            }

            a[0]     = 0.0;
            c[N - 1] = 0.0;

            cp[0] = c[0] / b[0];
            dp[0] = rhs[0] / b[0];
            for (int j = 1; j < N; ++j)
            {
                double denom = b[j] - a[j] * cp[j - 1];
                cp[j] = c[j] / denom;
                dp[j] = (rhs[j] - a[j] * dp[j - 1]) / denom;
            }

            Temp(i, N - 1, k) = dp[N - 1];
            for (int j = N - 2; j >= 0; --j)
            {
                Temp(i, j, k) = dp[j] - cp[j] * Temp(i, j + 1, k);
            }
        }
        Temp.SetBoundaryConditions(BC);
    }

    // --- Z-direction sweep: implicit along Z, explicit along X and Y ---
    if (Grid.dNz > 0)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temp.Tx,0,)
        {
            TxOld(i,j,k) = Temp(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) collapse(2)
        #endif
        for (int i = 0; i < Grid.Nx; ++i)
        for (int j = 0; j < Grid.Ny; ++j)
        {
            const int N = Grid.Nz;
            std::vector<double> a(N), b(N), c(N), rhs(N);
            std::vector<double> cp(N), dp(N);

            for (int k = 0; k < N; ++k)
            {
                double RhoCp  = EffectiveHeatCapacity(i,j,k);
                double Lambda = EffectiveThermalConductivity(i,j,k);
                double r      = Lambda * dt / (RhoCp * dx2);

                double explicitContrib = 0.0;
                if (Grid.dNx > 0) explicitContrib += Lambda * (TxOld(i+1,j,k) + TxOld(i-1,j,k)) / dx2;
                if (Grid.dNy > 0) explicitContrib += Lambda * (TxOld(i,j+1,k) + TxOld(i,j-1,k)) / dx2;

                a[k]   = -r;
                b[k]   =  1.0 + 2.0 * r;
                c[k]   = -r;
                rhs[k] = TxOld(i,j,k) + (explicitContrib + Qdot(i,j,k)) * dt / RhoCp;
            }

            a[0]     = 0.0;
            c[N - 1] = 0.0;

            cp[0] = c[0] / b[0];
            dp[0] = rhs[0] / b[0];
            for (int k = 1; k < N; ++k)
            {
                double denom = b[k] - a[k] * cp[k - 1];
                cp[k] = c[k] / denom;
                dp[k] = (rhs[k] - a[k] * dp[k - 1]) / denom;
            }

            Temp(i, j, N - 1) = dp[N - 1];
            for (int k = N - 2; k >= 0; --k)
            {
                Temp(i, j, k) = dp[k] - cp[k] * Temp(i, j, k + 1);
            }
        }
        Temp.SetBoundaryConditions(BC);
    }

    Qdot.Clear();
    Temp.CalculateMinMaxAvg();
    SolverCallsCounter++;
}

}// namespace openphase
