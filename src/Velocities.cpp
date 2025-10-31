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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitri Medvedev
 *
 */

#include "Velocities.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include "ElasticProperties.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;
/*************************************************************************/

Velocities::Velocities(Settings& locSettings)
{
    Initialize(locSettings);
}

void Velocities::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Velocities";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases    = locSettings.Nphases;
    PhaseNames = locSettings.PhaseNames;

    size_t Bcells = Grid.Bcells;
    Phase.Allocate(Grid, {Nphases}, Bcells);
    Average.Allocate(Grid, Bcells);

    Clear();

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Velocities::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(Average);
    BC.SetYVector(Average);
    BC.SetZVector(Average);

    BC.SetXVector(Phase);
    BC.SetYVector(Phase);
    BC.SetZVector(Phase);
}

void Velocities::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();
        for(size_t alpha = 0; alpha < Nphases; ++alpha)
        {
            Phase(i, j, k, {alpha}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void Velocities::SetAverage(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::SetAverage(ElasticProperties& EP, double dt)
{
    const double dx_dt = Grid.dx / dt;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.Displacements,EP.Displacements.Bcells(),)
    {
        Average(i,j,k) = EP.Displacements(i,j,k) * dx_dt;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::SetAllPhases(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(size_t alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k,{alpha}) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::PrescribePhaseVelocities(PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(size_t alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k,{alpha}) = Average(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::CalculateAverage(const PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();

        if(Phi.Fields(i,j,k).interface())
        {
            double locDensity = 0.;
            dVector3 locMomentumDensity = {0.,0.,0.};
            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t locPindex = Phi.FieldsProperties[alpha->index].Phase;
                locDensity += Phi.FieldsProperties[alpha->index].Density*alpha->value;
                locMomentumDensity += Phase(i,j,k,{locPindex})*Phi.FieldsProperties[alpha->index].Density*alpha->value;
            }
            Average(i,j,k) += locMomentumDensity/locDensity;
        }
        else
        {
            size_t pInd = Phi.FieldsProperties[Phi.Fields(i,j,k).front().index].Phase;
            Average(i,j,k) = Phase(i,j,k,{pInd});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::WriteVTK(Settings& locSettings, int tStep) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"VelocityAverage_",  [this](int i,int j,int k){return Average(i,j,k);}});

    for (size_t n = 0; n < Nphases; n++)
    {
        ListOfFields.push_back((VTK::Field_t){"VelocityPhase_" + std::to_string(n), [n,this](int i,int j,int k){return Phase(i,j,k,{n});}});
    }
    VTK::Write(Filename, locSettings, ListOfFields);
}

void Velocities::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);

    Grid.SetDimensions(newNx, newNy, newNz);

    Phase.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    Average.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    SetBoundaryConditions(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void Velocities::MoveFrame(int dx, int dy, int dz, const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    int xEnd = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    int yEnd = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    int zEnd = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Phase(i, j, k) = Phase(i + dx, j + dy, k + dz);
        Average(i, j, k) = Average(i + dx, j + dy, k + dz);
    }

    SetBoundaryConditions(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Frame moved");

}

void Velocities::PrintPointStatistics(int x, int y, int z)
{
    ConsoleOutput::WriteStandard("Point", iVector3{x,y,z});
    for (size_t PhaseIdx = 0; PhaseIdx < Nphases;  ++PhaseIdx)
    {
        std::string vname = "Phase Velocity "+PhaseNames[PhaseIdx];
        ConsoleOutput::WriteStandard(vname, Phase(x,y,z,{PhaseIdx}),3);
    }
    ConsoleOutput::WriteStandard("Average Velocity", Average(x,y,z),3);
    ConsoleOutput::WriteBlankLine();
}

double Velocities::GetMaxVelocity()
{
    double MAXVelocity = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,0, reduction(max: MAXVelocity))
    {
        MAXVelocity = max(MAXVelocity, Average(i,j,k).abs());
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &MAXVelocity, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif

    return MAXVelocity;
}

void Velocities::WriteStatistics(RunTimeControl& RTC)
{
    double MAXVelocity = GetMaxVelocity();

    ofstream output_file;
    if (RTC.TimeStep == 0)
    {
        output_file.open("VelocityStatistics.txt", ios::out);
        output_file << "tStep" << "\t\t\t" << "sim_time" << "\t\t\t"
                    << "max velocity" << endl;
        output_file.close();
    }

    output_file.open("VelocityStatistics.txt", ios::app);
    output_file << RTC.TimeStep << "\t\t\t" << RTC.SimulationTime << "\t\t\t" << MAXVelocity << endl;
    output_file.close();
}

Velocities& Velocities::operator= (const Velocities& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Velocities")
    {
        thisclassname = rhs.thisclassname;
        Grid = rhs.Grid;

        Nphases = rhs.Nphases;

        if (Phase.IsNotAllocated())
        {
            Phase.Allocate(Grid, {Nphases}, rhs.Phase.Bcells());
            Average.Allocate(Grid, rhs.Average.Bcells());
        }
        else if (not Phase.IsSize(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz))
        {
            Phase.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Average.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
        {
            Average(i,j,k) = rhs.Average(i,j,k);
            for(size_t alpha = 0; alpha < Nphases; ++alpha)
            {
                Phase(i,j,k, {alpha}) = rhs.Phase(i,j,k, {alpha});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void Velocities::Advect(BoundaryConditions& BC, double dt, AdvectionSchemes scheme)
{
    if(AverageDot.IsNotAllocated())
    {
        AverageDot.Allocate(Grid, 1);
    }
    if(not AverageDot.IsSize(Grid.Nx, Grid.Ny, Grid.Nz))
    {
        AverageDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case AdvectionSchemes::Upwind:
        {
            const double dx = Grid.dx;
            const double dx2 = 0.5/dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
            for (auto dir = 0; dir < 3; dir++)
            {
                AverageDot(i,j,k)[dir] = dx2 *
                     ((fabs(Average(i-1,j,k)[0]) + Average(i-1,j,k)[0])*Average(i-1,j,k)[dir] +
                      (fabs(Average(i,j-1,k)[1]) + Average(i,j-1,k)[1])*Average(i,j-1,k)[dir] +
                      (fabs(Average(i,j,k-1)[2]) + Average(i,j,k-1)[2])*Average(i,j,k-1)[dir] +
                      (fabs(Average(i+1,j,k)[0]) - Average(i+1,j,k)[0])*Average(i+1,j,k)[dir] +
                      (fabs(Average(i,j+1,k)[1]) - Average(i,j+1,k)[1])*Average(i,j+1,k)[dir] +
                      (fabs(Average(i,j,k+1)[2]) - Average(i,j,k+1)[2])*Average(i,j,k+1)[dir]) -
                      (fabs(Average(i,j,k)[0]) +
                       fabs(Average(i,j,k)[1]) +
                       fabs(Average(i,j,k)[2])) * Average(i, j, k)[dir]/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        default:
        {
            ConsoleOutput::WriteExit("Wrong advection scheme selected",
                             thisclassname, "Velocities::Advect(BC, dt, scheme)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
    for (int dir = 0; dir < 3; ++dir)
    {
        Average(i, j, k)[dir] += AverageDot(i, j, k)[dir]*dt;
        AverageDot(i, j, k)[dir] = 0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Advected");
}

bool Velocities::Write(const Settings& OPSettings, const int tStep) const
{
    std::string FileName = FileInterface::MakeFileName(OPSettings.RawDataDir,"Velocity_",tStep,".dat");

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::stringstream message;
        message << "File \"" << FileName << "\" could not be created! Terminating!!!" << std::endl;
        throw std::runtime_error(message.str());
        return false;
    };

    int Nx = Grid.Nx;
    int Ny = Grid.Ny;
    int Nz = Grid.Nz;

    out.write(reinterpret_cast<const char*>(&Nx), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Ny), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Nz), sizeof(int));

    STORAGE_LOOP_BEGIN(i,j,k,Average,0)
    {
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Average(i,j,k)[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
    for(size_t n = 0; n < Nphases; n++)
    {
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k,{n})[0]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k,{n})[1]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&Phase(i,j,k,{n})[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    out.close();
    return true;
}

bool Velocities::Read(const Settings& OPSettings, const BoundaryConditions& BC, const int tStep)
{
    std::string FileName = FileInterface::MakeFileName(OPSettings.RawDataDir,"Velocity_", tStep, ".dat");

    std::fstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        std::stringstream message;
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
        return false;
    };

    int locNx = 0;
    int locNy = 0;
    int locNz = 0;

    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));

    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz)
    {
        std::stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: ("
                << locNx << ", "
                << locNy << ", "
                << locNz << ") grid points.\n"
                << "Required data dimensions: ("
                << Grid.Nx << ", "
                << Grid.Ny << ", "
                << Grid.Nz << ") grid points.\n";
        throw std::runtime_error(message.str());
    }

    STORAGE_LOOP_BEGIN(i,j,k,Average,0)
    {
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[0]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[1]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Average(i,j,k)[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
    for(size_t n = 0; n < Nphases; n++)
    {
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k,{n})[0]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k,{n})[1]), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Phase(i,j,k,{n})[2]), sizeof(double));
    }
    STORAGE_LOOP_END

    inp.close();

    SetBoundaryConditions(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

}// namespace openphase
