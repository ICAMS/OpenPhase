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
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "Includes.h"
#include "DrivingForce.h"
#include "Noise.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"
#include "fftw3.h"

namespace openphase
{
using namespace std;

Noise::~Noise(void)
{
    if (initialized)
    {
        delete[] RandomFourier;
        delete[] RandomReal;
        fftw_destroy_plan(FFTBackward);
//#ifdef MPI_PARALLEL // TODO implement MPI
//        op_fftw_mpi_cleanup();
//#endif
#ifdef _OPENMP
        fftw_cleanup_threads();
#endif
    }
}

void Noise::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname  = "Noise";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases   = locSettings.Nphases;

    // Allocate three-dimensional Fourier space noise
    RandomFourier = new complex<double>[Grid.LocalNumberOfCells()]();
    RandomReal    = new complex<double>[Grid.LocalNumberOfCells()]();

    // Allocate raw Noise Storage
    Raw.Allocate(Grid, 0);

    // Assign storages to fftw pointers
    fftw_In  = reinterpret_cast<fftw_complex*>(RandomFourier);
    fftw_Out = reinterpret_cast<fftw_complex*>(RandomReal);

    // Create FFT-Plan, which determines how the FFT will be executed
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    FFTBackward = fftw_plan_dft_3d(Grid.Nx, Grid.Ny, Grid.Nz, fftw_In,fftw_Out,FFTW_BACKWARD,FFTW_ESTIMATE);

    // Initialize random number distribution
    distribution = normal_distribution<double>(0.0,0.5);

    size_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator = std::default_random_engine(seed);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Noise::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Noise input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);
    inp.close();
}

void Noise::ReadInput(std::stringstream& inp)
{
    // Read Parameters
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    gamma            = FileInterface::ReadParameterD(inp, moduleLocation, string("dgamma"),false,0);
    TimeSteps        = FileInterface::ReadParameterI(inp, moduleLocation, string("iTime"));
    CutOffWaveLength = FileInterface::ReadParameterI(inp, moduleLocation, string("dWaveLength"));
    FixAmpl          = FileInterface::ReadParameterI(inp, moduleLocation, string("dFixAmpl"));

    // Check if the input TimeSteps makes sense
    if (TimeSteps < 1)
    {
        std::string message = "A input value TimeSteps < 1 does not make sense! "
            "Set TimeSteps to 1!";
        ConsoleOutput::WriteWarning( message, thisclassname, "ReadInput");
        TimeSteps = 0;
    }

    // Check if the input CutOffWaveLength makes sense
    if (CutOffWaveLength < 0)
    {
        std::string message = "A input value CutOffWaveLength < 0 does not make sense! "
            "Set CutOffWaveLength to 0!";
        ConsoleOutput::WriteWarning( message, thisclassname, "ReadInput");
        CutOffWaveLength = 0;
    }
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Noise::Generate(void)
{
    double Lx = 2.0 * Pi / (double(Grid.Nx) * Grid.dx);
    double Ly = 2.0 * Pi / (double(Grid.Ny) * Grid.dx);
    double Lz = 2.0 * Pi / (double(Grid.Nz) * Grid.dx);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        dVector3 WaveVector;
        WaveVector[0] = Lx*(i*(i <= Grid.Nx/2) - (Grid.Nx-i)*(i > Grid.Nx/2));
        WaveVector[1] = Ly*(j*(j <= Grid.Ny/2) - (Grid.Ny-j)*(j > Grid.Ny/2));
        WaveVector[2] = Lz*(k*(k <= Grid.Nz/2) - (Grid.Nz-k)*(k > Grid.Nz/2));

        if(WaveVector.abs() < 2.0 * Pi/CutOffWaveLength)
        {
            // Generate Gaussian white noise coefficients
            complex<double> cijk{distribution(generator),distribution(generator)};
            //cijk /= abs(cijk);

            // Store computed coefficient
            RandomFourier[k + Grid.Nz*(j + Grid.Ny*i)] = cijk;
        }
        else
        {
            RandomFourier[k + Grid.Nz*(j + Grid.Ny*i)] = complex<double>(0.0,0.0);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    // Do Fourier transform
    fftw_execute(FFTBackward);

    // Calculate global maximum of RandomReal
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {
        Max = max(Max, abs(real(RandomReal[k + Grid.Nz*(j + Grid.Ny*i)])));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        RandomReal[k + Grid.Nz*(j + Grid.Ny*i)] /= Max;
        RandomReal[k + Grid.Nz*(j + Grid.Ny*i)] -= Raw(i,j,k);
        RandomReal[k + Grid.Nz*(j + Grid.Ny*i)] /= TimeSteps;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::UpdateRaw(int tStep)
{
    if (tStep == 0)
    {
        // Generate initial noise distribution
        Generate();

        // Save initial noise distribution and add it to the driving force

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Set initial noise
            Raw(i,j,k) = TimeSteps * real(RandomReal[k + Grid.Nz*(j + Grid.Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        // Generate new noise distribution in order to able to interpolate
        // between the old distribution and the new one.
        Generate();
    }
    else
    {
        // First check if a new noise Distribution has to be generated
        if (!(tStep%TimeSteps)) Generate();

        // Interpolate between old and new noise distribution by adding RandomReal
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Compute noise for this time step
            Raw(i,j,k) += real(RandomReal[k + Grid.Nz*(j + Grid.Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    // Calculate the global maximum of Raw
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {
        Max = max(Max, abs(Raw(i,j,k)));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        Raw(i,j,k) /= Max;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::AddToDrivingForce(int tStep, Temperature& Temp, DrivingForce& DF)
{
    UpdateRaw(tStep);

    // Interpolate between old and new noise distribution by adding RandomReal
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        // Add noise to driving force
        for (auto it = DF.Force(i,j,k).cbegin(); it < DF.Force(i,j,k).cend(); ++it)
        {
            double Ampl = sqrt( 2 * PhysicalConstants::k_Boltzmann * Temp(i,j,k)/gamma);
            DF.Force(i,j,k).add_raw(it->indexA, it->indexB, Ampl * Raw(i,j,k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::AddToDrivingForce(int tStep, DrivingForce& DF)
{
    UpdateRaw(tStep);

    // Interpolate between old and new noise distribution by adding RandomReal
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        // Add noise to driving force
        for (auto it = DF.Force(i,j,k).cbegin(); it < DF.Force(i,j,k).cend(); ++it)
        {
            DF.Force(i,j,k).add_raw(it->indexA, it->indexB, FixAmpl * Raw(i,j,k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::WriteVTK(const Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t){"Noise",  [this](int i,int j,int k){return Raw(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

}// end of name space
