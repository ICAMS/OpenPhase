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

 *   File created : 2025
 *   Main contributors : Efim Borukhovich; Oleg Shchyglo; Raphael Schiedung;
 *
 */

#include "Settings.h"
#include "Electrics/ElectricProperties.h"
#include "Electrics/ElectricSolverSpectral.h"
#include <complex>

namespace openphase
{

ElectricSolverSpectral::ElectricSolverSpectral(Settings& locSettings,
        const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void ElectricSolverSpectral::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ElectricSolverSpectral";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;
    Nz2 = locSettings.Grid.Nz/2+1;
    rlSize = Grid.Nx*Grid.Ny*Grid.Nz;
    ftSize = Grid.Nx*Grid.Ny*Nz2;
    freq = new fftw_complex [ftSize];
    rhs = new double [rlSize];

    ForwardPlan  = fftw_plan_dft_r2c_3d(Grid.Nx, Grid.Ny, Grid.Nz, rhs, freq, FFTW_ESTIMATE);
    BackwardPlan = fftw_plan_dft_c2r_3d(Grid.Nx, Grid.Ny, Grid.Nz, freq, rhs, FFTW_ESTIMATE);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

ElectricSolverSpectral::~ElectricSolverSpectral()
{
    delete [] freq;
    delete [] rhs;

    fftw_destroy_plan(ForwardPlan);
    fftw_destroy_plan(BackwardPlan);
}

void ElectricSolverSpectral::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Electric properties input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream inp;
    inp << inpF.rdbuf();

    ReadInput(inp);

    inpF.close();
}

void ElectricSolverSpectral::ReadInput(std::stringstream& inp)
{
}

void  ElectricSolverSpectral::Solve(ElectricProperties& EP, BoundaryConditions& BC)
{
    const double Lx = Grid.TotalNx*Grid.dx;
    const double Ly = Grid.TotalNy*Grid.dx;
    const double Lz = Grid.TotalNz*Grid.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.ChargeDensity,0,)
    {
        rhs[rlidx(i,j,k)] = - EP.ChargeDensity(i,j,k)/PhysicalConstants::epsilon_0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    fftw_execute(ForwardPlan);

    for (long int i = 0; i < Grid.Nx; i++)
    for (long int j = 0; j < Grid.Ny; j++)
    for (long int k = 0; k < Nz2; k++)
    {
        const long int ii = Grid.OffsetX + i;
        const long int jj = Grid.OffsetY + j;
        const long int kk = Grid.OffsetZ + k;

        assert(ii < Grid.TotalNx);
        assert(jj < Grid.TotalNy);
        assert(kk < Grid.TotalNz);

        const double ni = (ii <= Grid.TotalNx/2) ? ii : ii - Grid.TotalNx;
        const double nj = (jj <= Grid.TotalNy/2) ? jj : jj - Grid.TotalNy;
        const double nk = kk;

        const double Qx = 2.0*Pi*ni/Lx;
        const double Qy = 2.0*Pi*nj/Ly;
        const double Qz = 2.0*Pi*nk/Lz;

        const double QQ = Qx*Qx + Qy*Qy + Qz*Qz;

        const size_t XYZ = k + Nz2*(j + Grid.Ny*i);
        assert( XYZ < ftSize ); 
        if (QQ == 0.0)
        {
            freq[XYZ][0] = 0.0;
            freq[XYZ][1] = 0.0;
        }
        else
        {
            freq[XYZ][0] /= -QQ;
            freq[XYZ][1] /= -QQ;
        }
    }

    fftw_execute(BackwardPlan);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.ElectricPotential,0,)
    {
        EP.ElectricPotential(i,j,k) = rhs[rlidx(i,j,k)]/Grid.LocalNumberOfCells();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    BC.SetX(EP.ElectricPotential);
    BC.SetY(EP.ElectricPotential);
    BC.SetZ(EP.ElectricPotential);
}

}
