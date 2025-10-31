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

#ifdef MPI_PARALLEL
#include "mpi_wrapper.h"
#else
#include "fftw3.h"
#endif

#include "Settings.h"
#include "BoundaryConditions.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"

namespace openphase
{
using namespace std;

class OP_EXPORTS ElasticitySolverSpectralImpl : public OPObject                     ///< Elastic problem solver based on spectral method.
{
 public:
    ElasticitySolverSpectralImpl();
    ElasticitySolverSpectralImpl(Settings& locSettings,
                        const std::string InputFileName = DefaultInputFileName);///< Constructor calls Initialize and ReadInput
    ElasticitySolverSpectralImpl(Settings& locSettings,
                        const BoundaryConditions& BC,
                        const std::string InputFileName = DefaultInputFileName);///< Constructor calls Initialize and ReadInput
    ~ElasticitySolverSpectralImpl(void);
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Named constructor, allocates internal storages and initialized variables
    void Initialize(Settings& locSettings, const BoundaryConditions& BC, std::string ObjectNameSuffix = ""); ///< Named constructor (considering non periodic boundary conditions)
    void ReadInput(const std::string InputFileName) override;                   ///< Reads elastic properties from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads elastic properties from the input stream
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override;///< Remeshes the storages, reinitilizes the solver

    void Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC);     ///< Needed if the system dimensions have been altered (considering non periodic boundary conditions)

    int  Solve(ElasticProperties& EP, BoundaryConditions& BC, double dt,
               std::function<bool()> EndSolve = [](){return false;});           ///< Solves mechanical equilibrium problem
    
    void WriteDeformationIncrementsVTK(const int tStep,
                                         const Settings& locSettings,
                                         const int precision=16) const;         ///< Writes deformation increments in VTK format
 private:
    
    GridParameters Grid;                                                        ///< Simulation grid parameters

    double IncrementScaling;                                                    ///< Improves solver convergence for extreme deformations by reducing deformation increment in the current iteration. Recommended values 1 (default) for normal cases, 0.5 ~ 1.0 for extreme deformations.

    double StrainAccuracy;                                                      ///< Convergence parameter for the elasticity solver
    double StressAccuracy;                                                      ///< Convergence parameter for the elasticity solver
    int    MAXIterations;                                                       ///< Maximum iterations per time step for the elasticity solver
    bool   VerboseIterations;                                                   ///< If true enables iterations statistics output to console
    bool   DiscreteDerivatives;                                                 ///< If true, uses sin(Q) instead of Q (wave vector) for reciprocal space derivatives calculations

    double *  UandForce[3];                                                     ///< Force and displacements in real and reciprocal space
    double *  RHSandDefGrad[9];                                                 ///< RHS and deformation gradient in real and reciprocal space

    fftw_plan ForwardPlanRHS[9];                                                ///< Forward FFT plans for the RHSide
    fftw_plan ForwardPlanForce[3];                                              ///< Forward FFT plans for the force
    fftw_plan BackwardPlanDefGrad[9];                                           ///< Backward FFT plans for the deformation gradients
    fftw_plan BackwardPlanU[3];                                                 ///< Backward FFT plans for the displacements

    void CopyForceDensity(ElasticProperties& EP);                               ///< Copies force density into the internal storage of the solver
    void CalculateRHS(const ElasticProperties& EP, const dMatrix6x6& Cij);      ///< Calculates right hand side of the mechanical equilibrium equation
    void ExecuteForwardFFT(void);                                               ///< Executes forward FFT of the RHS
    void CalculateFourierSolution(const ElasticProperties& EP, const dMatrix6x6& Cij);///< Calculates solution in Fourier space
    void ExecuteBackwardFFT(void);                                              ///< Executes backward FFT of the deformation gradients
    void ExecuteBackwardFFTdisplacements(void);                                 ///< Execute forward FFT of the displacements
    void ExecuteForwardFFTforces(void);                                         ///< Execute forward FFT of the external forces
    void CalculateTargetStress(ElasticProperties& EP,
                               vStress& TargetStress);                          ///< Calculates average stress
    void CopyDisplacements(ElasticProperties& EP);                              ///< Copies displacements to the ElasticProperties
    void SetElasticProperties(ElasticProperties& EP,
                              double& MAXStrainDifference,
                              double& MAXStressDifference);                     ///< Sets stresses and deformation gradients

    void ApplyMechanicalBC(ElasticProperties& EP, vStress TargetStress);        ///< Applies mechanical boundary conditions
    bool OPMECH;
};

ElasticitySolverSpectral::ElasticitySolverSpectral() : impl_(new ElasticitySolverSpectralImpl) {}

ElasticitySolverSpectral::ElasticitySolverSpectral(Settings& locSettings, const std::string InputFileName) : impl_(new ElasticitySolverSpectralImpl) 
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}
ElasticitySolverSpectral::ElasticitySolverSpectral(Settings& locSettings, const BoundaryConditions& BC,  const std::string InputFileName) : impl_(new ElasticitySolverSpectralImpl) 
{
    Initialize(locSettings, BC);
    ReadInput(InputFileName);
}

ElasticitySolverSpectral::~ElasticitySolverSpectral() = default;
            
void ElasticitySolverSpectral::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    impl_->Initialize(locSettings, ObjectNameSuffix);
}                         
void ElasticitySolverSpectral::Initialize(Settings& locSettings, const BoundaryConditions& BC, std::string ObjectNameSuffix)
{
    impl_->Initialize(locSettings, BC, ObjectNameSuffix);
}       
void ElasticitySolverSpectral::Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC)     ///< Needed if the system dimensions have been altered (considering non periodic boundary conditions)
{
    impl_->Reinitialize(EP, BC);
}

void ElasticitySolverSpectral::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    impl_->Remesh(newNx, newNy, newNz, BC);
}

void ElasticitySolverSpectral::ReadInput(const std::string InputFileName) 
{
    impl_->ReadInput(InputFileName) ;
}
void ElasticitySolverSpectral::ReadInput(std::stringstream& inp) 
{
    impl_->ReadInput(inp);
}
int ElasticitySolverSpectral::Solve(ElasticProperties& EP, BoundaryConditions& BC, double dt, std::function<bool()> EndSolve)
{
    return impl_->Solve(EP,BC,dt, EndSolve);
}

ElasticitySolverSpectralImpl::ElasticitySolverSpectralImpl()
{
    IncrementScaling = 1.0;
    StrainAccuracy   = 1.0e-6;
    StressAccuracy   = 0.0;
    MAXIterations    = 100;

    DiscreteDerivatives = false;
    VerboseIterations = false;
    OPMECH = false;
};

ElasticitySolverSpectralImpl::ElasticitySolverSpectralImpl(Settings& locSettings,
                                                const std::string InputFileName)
{
    OPMECH = false;
    Initialize(locSettings);
    ReadInput(InputFileName);
}

ElasticitySolverSpectralImpl::ElasticitySolverSpectralImpl(Settings& locSettings,
                         const BoundaryConditions& BC, const std::string InputFileName)
{
    OPMECH = false;
    Initialize(locSettings, BC);
    ReadInput(InputFileName);
} 

void ElasticitySolverSpectralImpl::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ElasticitySolverSpectral";
    thisobjectname = thisclassname + ObjectNameSuffix;

    IncrementScaling = 1.0;
    StressAccuracy   = 0.0;
    StrainAccuracy   = 1.0e-6;
    MAXIterations    = 100;

    VerboseIterations = false;
    DiscreteDerivatives = false;
    //Initialize FFTW OpenMP threads
#ifdef _OPENMP
    fftw_init_threads();
#endif

    // Set system dimensions
    Grid = locSettings.Grid;

    // Initialize FFTW MPI parallelism
#ifdef MPI_PARALLEL
    op_fftw_mpi_init();

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*op_fftw_mpi_local_size_3d(Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, OP_MPI_COMM_WORLD,
                                              &local_n0, &local_0_start);
#else
    size_t SIZE = Grid.Nx*Grid.Ny*Grid.Nz2*2;
#endif

    // Allocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Set number of FFTW OpenMP threads
#ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        OP_MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = op_fftw_mpi_plan_dft_c2r_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         OP_MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = op_fftw_mpi_plan_dft_c2r_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                         reinterpret_cast<fftw_complex*> (UandForce[n]),
                         UandForce[n],
                         OP_MPI_COMM_WORLD,
                         FFTW_ESTIMATE);

        ForwardPlanForce[n] = op_fftw_mpi_plan_dft_r2c_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            OP_MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz,
                         reinterpret_cast<fftw_complex*> (UandForce[n]),
                         UandForce[n],
                         FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz, UandForce[n],
                         reinterpret_cast<fftw_complex*> (UandForce[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    locSettings.AddForRemeshing(*this);
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectralImpl::Initialize(Settings& locSettings, const BoundaryConditions& BC, std::string ObjectNameSuffix)
{
    thisclassname = "ElasticitySolverSpectral";
    thisobjectname = thisclassname + ObjectNameSuffix;

    IncrementScaling = 1.0;
    StressAccuracy   = 0.0;
    StrainAccuracy   = 1.0e-6;
    MAXIterations    = 100;

    VerboseIterations = false;
    DiscreteDerivatives = false;

    // Set system dimensions
    Grid = locSettings.Grid;

    if(BC.BC0X != BoundaryConditionTypes::Periodic or
       BC.BCNX != BoundaryConditionTypes::Periodic)
    {
#ifdef MPI_PARALLEL
        std::cerr << "ElasticitySolverSpectralImpl::Initialize()\n"
                  << "\tNonperiodic boundary conditions in X-direction are not permitted in MPI parallel mode!" << std::endl;
        OP_Exit(EXIT_FAILURE);
#else
        Grid.TotalNx *= 2;
        Grid.Nx *= 2;
#endif
    }
    if(BC.BC0Y != BoundaryConditionTypes::Periodic or
       BC.BCNY != BoundaryConditionTypes::Periodic)
    {
        Grid.TotalNy *= 2;
        Grid.Ny *= 2;
    }
    if(BC.BC0Z != BoundaryConditionTypes::Periodic or
       BC.BCNZ != BoundaryConditionTypes::Periodic)
    {
        Grid.TotalNz *= 2;
        Grid.Nz *= 2;
    }

    Grid.Nz2 = Grid.Nz/2+1;

    // Initialize FFTW MPI parallelism
#ifdef MPI_PARALLEL
    op_fftw_mpi_init();

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*op_fftw_mpi_local_size_3d(Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, OP_MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Grid.Nx*Grid.Ny*Grid.Nz2*2;
#endif

    // Initialize FFTW OpenMP threads
#ifdef _OPENMP
   fftw_init_threads();
#endif

   // Set number of FFTW OpenMP threads
#ifdef _OPENMP
   fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Allocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Create FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        OP_MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = op_fftw_mpi_plan_dft_c2r_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         OP_MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = op_fftw_mpi_plan_dft_c2r_3d
                            (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             OP_MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            OP_MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectralImpl::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    // Destroy old FFT plans and free allocated memory
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        fftw_free(RHSandDefGrad[n]);
    }
    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);
        fftw_destroy_plan(ForwardPlanForce[n]);

        fftw_free(UandForce[n]);
    }

#ifdef MPI_PARALLEL

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*op_fftw_mpi_local_size_3d(Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, OP_MPI_COMM_WORLD,
                                  &local_n0, &local_0_start);
#else
    size_t SIZE = Grid.Nx*Grid.Ny*Grid.Nz2*2;
#endif

    // Reallocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Create new FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        OP_MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = op_fftw_mpi_plan_dft_c2r_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         OP_MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = op_fftw_mpi_plan_dft_c2r_3d
                            (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             OP_MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            OP_MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    ConsoleOutput::Write(thisclassname, "Remeshed");
}

void ElasticitySolverSpectralImpl::Reinitialize(ElasticProperties& EP, const BoundaryConditions& BC)
{
    // Set new system dimensions
    Grid = EP.Grid;

    if(BC.BC0X != BoundaryConditionTypes::Periodic or
       BC.BCNX != BoundaryConditionTypes::Periodic)
    {
#ifdef MPI_PARALLEL
        std::cerr << "ElasticitySolverSpectralImpl::Initialize()\n"
                  << "\tNonperiodic boundary conditions in X-direction are not permitted in MPI parallel mode!" << std::endl;
        OP_Exit(EXIT_FAILURE);
#else
        Grid.TotalNx *= 2;
        Grid.Nx *= 2;
#endif
    }
    if(BC.BC0Y != BoundaryConditionTypes::Periodic or
       BC.BCNY != BoundaryConditionTypes::Periodic)
    {
        Grid.TotalNy *= 2;
        Grid.Ny *= 2;
    }
    if(BC.BC0Z != BoundaryConditionTypes::Periodic or
       BC.BCNZ != BoundaryConditionTypes::Periodic)
    {
        Grid.TotalNz *= 2;
        Grid.Nz *= 2;
    }

    Grid.Nz2 = Grid.Nz/2 + 1;

    // Destroy old FFT plans and free allocated memory
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        fftw_free(RHSandDefGrad[n]);
    }
    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);
        fftw_destroy_plan(ForwardPlanForce[n]);

        fftw_free(UandForce[n]);
    }

#ifdef _OPENMP
    fftw_cleanup_threads();
#else
    fftw_cleanup();
#endif

#ifdef MPI_PARALLEL
    op_fftw_mpi_cleanup();
#endif

    // Initialize FFTW MPI parallelism
#ifdef MPI_PARALLEL
    op_fftw_mpi_init();

    ptrdiff_t local_n0      = 0;
    ptrdiff_t local_0_start = 0;

    size_t SIZE = 2*op_fftw_mpi_local_size_3d(Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, OP_MPI_COMM_WORLD,
                                              &local_n0, &local_0_start);
#else
    size_t SIZE = Grid.Nx*Grid.Ny*Grid.Nz2*2;
#endif

    //Initialize FFTW OpenMP threads
#ifdef _OPENMP
    fftw_init_threads();
#endif

    // Set number of FFTW OpenMP threads
#ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Reallocate arrays
    for(int n = 0; n < 9; n++)
    {
        RHSandDefGrad[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }
    for(int n = 0; n < 3; n++)
    {
        UandForce[n] = (double *)fftw_malloc(sizeof(double)*SIZE);
    }

    // Create new FFT plans
    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        ForwardPlanRHS[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, RHSandDefGrad[n],
                        reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                        OP_MPI_COMM_WORLD,
                        FFTW_ESTIMATE);
        #else
        ForwardPlanRHS[n] = fftw_plan_dft_r2c_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz, RHSandDefGrad[n],
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 9; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanDefGrad[n] = op_fftw_mpi_plan_dft_c2r_3d
                        (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         OP_MPI_COMM_WORLD,
                         FFTW_ESTIMATE);
        #else
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Grid.Nx, Grid.Ny, Grid.Nz,
                         reinterpret_cast<fftw_complex*> (RHSandDefGrad[n]),
                         RHSandDefGrad[n],
                         FFTW_ESTIMATE);
        #endif
    }

    for(int n = 0; n < 3; n++)
    {
        #ifdef MPI_PARALLEL
        BackwardPlanU[n] = op_fftw_mpi_plan_dft_c2r_3d
                            (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             OP_MPI_COMM_WORLD,
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = op_fftw_mpi_plan_dft_r2c_3d (Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, UandForce[n],
                            reinterpret_cast<fftw_complex*> (UandForce[n]),
                            OP_MPI_COMM_WORLD,
                            FFTW_ESTIMATE);
        #else
        BackwardPlanU[n] = fftw_plan_dft_c2r_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz,
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             UandForce[n],
                             FFTW_ESTIMATE);

        ForwardPlanForce[n] = fftw_plan_dft_r2c_3d
                            (Grid.Nx, Grid.Ny, Grid.Nz, UandForce[n],
                             reinterpret_cast<fftw_complex*> (UandForce[n]),
                             FFTW_ESTIMATE);
        #endif
    }
    ConsoleOutput::WriteStandard(thisclassname, "Reinitialized");
}

void ElasticitySolverSpectralImpl::ReadInput(const string InputFileName)
{
    ConsoleOutput::Write("Source", InputFileName);
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + InputFileName + " could not be opened",thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ReadInput(data);
}

void ElasticitySolverSpectralImpl::ReadInput(stringstream& inp)
{
    ConsoleOutput::WriteLineInsert("ElasticitySolverSpectralImpl input");

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    StrainAccuracy      = FileInterface::ReadParameterD(inp, moduleLocation, "StrainAccuracy",      false, StrainAccuracy);
    StressAccuracy      = FileInterface::ReadParameterD(inp, moduleLocation, "StressAccuracy",      false, StressAccuracy);
    MAXIterations       = FileInterface::ReadParameterI(inp, moduleLocation, "MAXIterations",       false, MAXIterations);
    IncrementScaling    = FileInterface::ReadParameterD(inp, moduleLocation, "IncrementScaling",    false, IncrementScaling);
    VerboseIterations   = FileInterface::ReadParameterB(inp, moduleLocation, "VerboseIterations",   false, VerboseIterations);
    DiscreteDerivatives = FileInterface::ReadParameterB(inp, moduleLocation, "DiscreteDerivatives", false, DiscreteDerivatives);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

//++++++++++++ Destructor +++++++++++++++++++++++++

ElasticitySolverSpectralImpl::~ElasticitySolverSpectralImpl(void)
{
    if(initialized)
    {
        for(int n = 0; n < 9; n++)
        {
            fftw_destroy_plan(ForwardPlanRHS[n]);
            fftw_destroy_plan(BackwardPlanDefGrad[n]);

            fftw_free(RHSandDefGrad[n]);
        }

        for(int n = 0; n < 3; n++)
        {
            fftw_destroy_plan(BackwardPlanU[n]);
            fftw_destroy_plan(ForwardPlanForce[n]);

            fftw_free(UandForce[n]);
        }

#ifdef MPI_PARALLEL
        op_fftw_mpi_cleanup();
#endif

#ifdef _OPENMP
        fftw_cleanup_threads();
#else
        fftw_cleanup();
#endif
    }
}

int ElasticitySolverSpectralImpl::Solve(ElasticProperties& EP,
                                    BoundaryConditions& BC, double dt,
                                    std::function<bool()> EndSolve)
{
    dMatrix6x6 Cij = EP.MAXElasticConstants();

    vStress TargetStress;
    vStrain oldAverageStrain;

    int    IterationsCount = 0;
    double MAXLocalStressDifference = 0.0;
    double MAXLocalStrainDifference = 0.0;
    double MAXAverageStrainDifference = 0.0;

    if (!EP.LargeDeformations)
    {
        for(int n = 0; n < 6; n++)
        if(EP.AppliedStrainRateMask[n])
        {
            EP.AppliedStrain[n] += EP.AppliedStrainRate[n]*dt;
        }
    }

    do // Iteration loop begin
    {
        MAXLocalStressDifference = 0.0;
        MAXLocalStrainDifference = 0.0;
        MAXAverageStrainDifference = 0.0;
        if(EP.ConsiderExternalForces)
        {
            CopyForceDensity(EP);
            ExecuteForwardFFTforces();
        }
        CalculateRHS(EP, Cij);
        ExecuteForwardFFT();
        CalculateFourierSolution(EP, Cij);
        ExecuteBackwardFFT();
        CalculateTargetStress(EP, TargetStress);

        oldAverageStrain = EP.AverageStrain;

        ApplyMechanicalBC(EP, TargetStress);
        SetElasticProperties(EP, MAXLocalStrainDifference, MAXLocalStressDifference);
        MAXAverageStrainDifference = max(MAXAverageStrainDifference,(EP.AverageStrain - oldAverageStrain).norm());

#ifdef MPI_PARALLEL
        OP_MPI_Allreduce(OP_MPI_IN_PLACE,&MAXAverageStrainDifference,1,OP_MPI_DOUBLE,OP_MPI_MAX,OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE,&MAXLocalStrainDifference,1,OP_MPI_DOUBLE,OP_MPI_MAX,OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE,&MAXLocalStressDifference,1,OP_MPI_DOUBLE,OP_MPI_MAX,OP_MPI_COMM_WORLD);
#endif

        if(VerboseIterations)
        {
            cout << EP.AverageStrain.print() << endl;
            cout << EP.AverageDeformationGradient(StrainAccuracy*0.1).print() << endl;

            std::string message  = "Iteration " + std::to_string(IterationsCount) + ":\n";
                        message += ConsoleOutput::GetStandard("Average strains converged to",
                                   ConsoleOutput::to_string_with_precision(MAXAverageStrainDifference));
                        message += ConsoleOutput::GetStandard("Local strains converged to",
                                   ConsoleOutput::to_string_with_precision(MAXLocalStrainDifference));
                        message += ConsoleOutput::GetStandard("Local stresses converged to",
                                   ConsoleOutput::to_string_with_precision(MAXLocalStressDifference));
            ConsoleOutput::WriteWithinMethod(message, thisclassname, "Solve()");
        }

        IterationsCount++;

        if(IterationsCount > MAXIterations)
        {
            std::string message  = "Maximum number of iterations (" + std::to_string(MAXIterations) + ") reached\n";
                        message += ConsoleOutput::GetStandard("Average strains converged to",
                                   ConsoleOutput::to_string_with_precision(MAXAverageStrainDifference));
                        message += ConsoleOutput::GetStandard("Local strains converged to",
                                   ConsoleOutput::to_string_with_precision(MAXLocalStrainDifference));
                        message += ConsoleOutput::GetStandard("Local stresses converged to",
                                   ConsoleOutput::to_string_with_precision(MAXLocalStressDifference));
            ConsoleOutput::WriteWarning(message, thisclassname, "Solve()");
            break;
        }
    } // Iteration loop end
    while(not EndSolve() and ((MAXLocalStressDifference > StressAccuracy and
                               MAXLocalStrainDifference > StrainAccuracy) or
                               MAXAverageStrainDifference > StrainAccuracy));

    ExecuteBackwardFFTdisplacements();
    CopyDisplacements(EP);

    if(EP.LargeDeformations)
    {
        EP.SetVelocityGradients(dt);
    }

    EP.SetBoundaryConditions(BC);

    return IterationsCount;
}

void ElasticitySolverSpectralImpl::CalculateRHS(const ElasticProperties& EP, const dMatrix6x6& Cij)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
    {
        dMatrix3x3 locRHStensor = (Cij*StrainSmall(EP.DeformationGradientsTotal(i,j,k))).tensor() - EP.Stress1PK(i,j,k);

        const long int xyz = k + 2*Grid.Nz2*(j + Grid.Ny*i);

        RHSandDefGrad[0][xyz] = locRHStensor(0,0);
        RHSandDefGrad[1][xyz] = locRHStensor(0,1);
        RHSandDefGrad[2][xyz] = locRHStensor(0,2);

        RHSandDefGrad[3][xyz] = locRHStensor(1,0);
        RHSandDefGrad[4][xyz] = locRHStensor(1,1);
        RHSandDefGrad[5][xyz] = locRHStensor(1,2);

        RHSandDefGrad[6][xyz] = locRHStensor(2,0);
        RHSandDefGrad[7][xyz] = locRHStensor(2,1);
        RHSandDefGrad[8][xyz] = locRHStensor(2,2);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    if(Grid.Ny > EP.Grid.Ny)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int y = Grid.Ny - j - 1;
            RHSandDefGrad[n][k + 2*Grid.Nz2*(y + Grid.Ny*i)] = RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Grid.Nz > EP.Grid.Nz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int z = Grid.Nz - k - 1;
            RHSandDefGrad[n][z + 2*Grid.Nz2*(j + Grid.Ny*i)] = RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#else
    if(Grid.Nx > EP.Grid.Nx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int x = Grid.Nx - i - 1;
            RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*x)] = RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Grid.Ny > EP.Grid.Ny)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int y = Grid.Ny - j - 1;
            RHSandDefGrad[n][k + 2*Grid.Nz2*(y + Grid.Ny*i)] = RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Grid.Nz > EP.Grid.Nz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.DeformationGradientsTotal, 0,)
        for(int n = 0; n < 9; n++)
        {
            int z = Grid.Nz - k - 1;
            RHSandDefGrad[n][z + 2*Grid.Nz2*(j + Grid.Ny*i)] = RHSandDefGrad[n][k + 2*Grid.Nz2*(j + Grid.Ny*i)];
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
#endif
}

void ElasticitySolverSpectralImpl::CalculateFourierSolution(const ElasticProperties& EP, const dMatrix6x6& Cij)
{
    double Norm = 1.0/double(Grid.TotalNumberOfCells());

    double DPi_Nx = 2.0*Pi/double(Grid.TotalNx);
    double DPi_Ny = 2.0*Pi/double(Grid.TotalNy);
    double DPi_Nz = 2.0*Pi/double(Grid.TotalNz);

    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Grid.Nx ; i++)
    for(int j = 0; j < Grid.Ny ; j++)
    for(int k = 0; k < Grid.Nz2; k++)
    {
        long int XYZ = k + Grid.Nz2*(j + Grid.Ny*i);
#ifdef MPI_PARALLEL
        long int ii = i + Grid.OffsetX;
        double Qx = DPi_Nx*(ii*(ii <= Grid.TotalNx/2) - (Grid.TotalNx-ii)*(ii > Grid.TotalNx/2));
        long int jj = j + Grid.OffsetY;
        double Qy = DPi_Ny*(jj*(jj <= Grid.TotalNy/2) - (Grid.TotalNy-jj)*(jj > Grid.TotalNy/2));
        long int kk = k + Grid.OffsetZ;
        double Qz = DPi_Nz*(kk*(kk <= Grid.TotalNz/2) - (Grid.TotalNz-kk)*(kk > Grid.TotalNz/2));
#else
        double Qx = DPi_Nx*(i*(i <= Grid.Nx/2) - (Grid.Nx-i)*(i > Grid.Nx/2));
        double Qy = DPi_Ny*(j*(j <= Grid.Ny/2) - (Grid.Ny-j)*(j > Grid.Ny/2));
        double Qz = DPi_Nz*(k*(k <= Grid.Nz/2) - (Grid.Nz-k)*(k > Grid.Nz/2));

        if(DiscreteDerivatives)
        {
            Qx = sin(Qx);
            Qy = sin(Qy);
            Qz = sin(Qz);
        }
#endif

        complex<double> rhsX = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[0])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[1])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[2])[XYZ]);
        complex<double> rhsY = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[3])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[4])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[5])[XYZ]);
        complex<double> rhsZ = -I*(Qx*reinterpret_cast<complex<double>*>(RHSandDefGrad[6])[XYZ] +
                                   Qy*reinterpret_cast<complex<double>*>(RHSandDefGrad[7])[XYZ] +
                                   Qz*reinterpret_cast<complex<double>*>(RHSandDefGrad[8])[XYZ]);

        if(EP.ConsiderExternalForces)
        {
            rhsX += reinterpret_cast<complex<double>*>(UandForce[0])[XYZ];
            rhsY += reinterpret_cast<complex<double>*>(UandForce[1])[XYZ];
            rhsZ += reinterpret_cast<complex<double>*>(UandForce[2])[XYZ];
        }

        double a11 = (Cij(0,0)*Qx*Qx + 2.0*Cij(0,5)*Qx*Qy + Cij(5,5)*Qy*Qy +
                  2.0*Cij(0,4)*Qx*Qz + 2.0*Cij(4,5)*Qy*Qz + Cij(4,4)*Qz*Qz);

        double a21 = (Cij(0,5)*Qx*Qx + Cij(0,1)*Qx*Qy + Cij(5,5)*Qx*Qy +
                      Cij(1,5)*Qy*Qy + Cij(0,3)*Qx*Qz + Cij(4,5)*Qx*Qz +
                      Cij(1,4)*Qy*Qz + Cij(3,5)*Qy*Qz + Cij(3,4)*Qz*Qz);

        double a31 = (Cij(0,4)*Qx*Qx + Cij(0,3)*Qx*Qy + Cij(4,5)*Qx*Qy +
                      Cij(3,5)*Qy*Qy + Cij(0,2)*Qx*Qz + Cij(4,4)*Qx*Qz +
                      Cij(2,5)*Qy*Qz + Cij(3,4)*Qy*Qz + Cij(2,4)*Qz*Qz);

        double a12 = a21;

        double a22 = (Cij(5,5)*Qx*Qx + 2.0*Cij(1,5)*Qx*Qy + Cij(1,1)*Qy*Qy +
                  2.0*Cij(3,5)*Qx*Qz + 2.0*Cij(1,3)*Qy*Qz + Cij(3,3)*Qz*Qz);

        double a32 = (Cij(4,5)*Qx*Qx + Cij(1,4)*Qx*Qy + Cij(3,5)*Qx*Qy +
                      Cij(1,3)*Qy*Qy + Cij(2,5)*Qx*Qz + Cij(3,4)*Qx*Qz +
                      Cij(1,2)*Qy*Qz + Cij(3,3)*Qy*Qz + Cij(2,3)*Qz*Qz);

        double a13 = a31;

        double a23 = a32;

        double a33 = (Cij(4,4)*Qx*Qx + 2.0*Cij(3,4)*Qx*Qy + Cij(3,3)*Qy*Qy +
                  2.0*Cij(2,4)*Qx*Qz + 2.0*Cij(2,3)*Qy*Qz + Cij(2,2)*Qz*Qz);

        double denominator = (-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 -
                               a11*a23*a32 - a12*a21*a33 + a11*a22*a33);

        if(std::abs(denominator) > DBL_EPSILON and std::abs(denominator) < DBL_MAX)
        {
            denominator = 1.0/denominator;
        }
        else
        {
            denominator = 0.0;
        }

        complex<double> locUrcX = (-a23*a32*rhsX + a22*a33*rhsX + a13*a32*rhsY -
                                    a12*a33*rhsY - a13*a22*rhsZ + a12*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcY = ( a23*a31*rhsX - a21*a33*rhsX - a13*a31*rhsY +
                                    a11*a33*rhsY + a13*a21*rhsZ - a11*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcZ = (-a22*a31*rhsX + a21*a32*rhsX + a12*a31*rhsY -
                                    a11*a32*rhsY - a12*a21*rhsZ + a11*a22*rhsZ)*denominator*Norm;

        reinterpret_cast<complex<double>*>(UandForce[0])[XYZ] = locUrcX;
        reinterpret_cast<complex<double>*>(UandForce[1])[XYZ] = locUrcY;
        reinterpret_cast<complex<double>*>(UandForce[2])[XYZ] = locUrcZ;

        reinterpret_cast<complex<double>*>(RHSandDefGrad[0])[XYZ] = I*(Qx*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[1])[XYZ] = I*(Qy*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[2])[XYZ] = I*(Qz*locUrcX);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[3])[XYZ] = I*(Qx*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[4])[XYZ] = I*(Qy*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[5])[XYZ] = I*(Qz*locUrcY);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[6])[XYZ] = I*(Qx*locUrcZ);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[7])[XYZ] = I*(Qy*locUrcZ);
        reinterpret_cast<complex<double>*>(RHSandDefGrad[8])[XYZ] = I*(Qz*locUrcZ);
    }
}

void ElasticitySolverSpectralImpl::CalculateTargetStress(ElasticProperties& EP,
                                                    vStress& TargetStress)
{
    const dMatrix3x3 avgDefGrad = EP.AverageDeformationGradient(StrainAccuracy*0.1);
    vStress locAverageStress;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.DeformationGradientsTotal,0,reduction(vStressSUM:locAverageStress))
    {
        const long int xyz = k + 2*Grid.Nz2*(j + Grid.Ny*i);

        dMatrix3x3 locDefGrad = avgDefGrad;
        locDefGrad(0,0) += RHSandDefGrad[0][xyz];
        locDefGrad(0,1) += RHSandDefGrad[1][xyz];
        locDefGrad(0,2) += RHSandDefGrad[2][xyz];

        locDefGrad(1,0) += RHSandDefGrad[3][xyz];
        locDefGrad(1,1) += RHSandDefGrad[4][xyz];
        locDefGrad(1,2) += RHSandDefGrad[5][xyz];

        locDefGrad(2,0) += RHSandDefGrad[6][xyz];
        locDefGrad(2,1) += RHSandDefGrad[7][xyz];
        locDefGrad(2,2) += RHSandDefGrad[8][xyz];

        vStress locStress = EP.Stress2PK(i,j,k,locDefGrad);
        locAverageStress += locStress;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, locAverageStress.data(),6, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif
    locAverageStress /= double(Grid.TotalNumberOfCells());

    EP.AverageStress = locAverageStress;

    TargetStress[0] = -EP.AverageStress[0]*(1.0 - EP.AppliedStressMask[0]) + (EP.AppliedStress[0]-EP.AverageStress[0])*EP.AppliedStressMask[0];
    TargetStress[1] = -EP.AverageStress[1]*(1.0 - EP.AppliedStressMask[1]) + (EP.AppliedStress[1]-EP.AverageStress[1])*EP.AppliedStressMask[1];
    TargetStress[2] = -EP.AverageStress[2]*(1.0 - EP.AppliedStressMask[2]) + (EP.AppliedStress[2]-EP.AverageStress[2])*EP.AppliedStressMask[2];
    TargetStress[3] = -EP.AverageStress[3]*(1.0 - EP.AppliedStressMask[3]) + (EP.AppliedStress[3]-EP.AverageStress[3])*EP.AppliedStressMask[3];
    TargetStress[4] = -EP.AverageStress[4]*(1.0 - EP.AppliedStressMask[4]) + (EP.AppliedStress[4]-EP.AverageStress[4])*EP.AppliedStressMask[4];
    TargetStress[5] = -EP.AverageStress[5]*(1.0 - EP.AppliedStressMask[5]) + (EP.AppliedStress[5]-EP.AverageStress[5])*EP.AppliedStressMask[5];
}

void ElasticitySolverSpectralImpl::SetElasticProperties(ElasticProperties& EP,
                                                    double& MAXStrainDifference,
                                                    double& MAXStressDifference)
{
    const dMatrix3x3 avgDefGrad = EP.AverageDeformationGradient(StrainAccuracy*0.1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.DeformationGradientsTotal,0,reduction(max:MAXStrainDifference) reduction(max:MAXStressDifference) )
    {
        const long int xyz = k + 2*Grid.Nz2*(j + Grid.Ny*i);

        dMatrix3x3 locDefGrad = avgDefGrad;
        locDefGrad(0,0) += RHSandDefGrad[0][xyz];
        locDefGrad(0,1) += RHSandDefGrad[1][xyz];
        locDefGrad(0,2) += RHSandDefGrad[2][xyz];

        locDefGrad(1,0) += RHSandDefGrad[3][xyz];
        locDefGrad(1,1) += RHSandDefGrad[4][xyz];
        locDefGrad(1,2) += RHSandDefGrad[5][xyz];

        locDefGrad(2,0) += RHSandDefGrad[6][xyz];
        locDefGrad(2,1) += RHSandDefGrad[7][xyz];
        locDefGrad(2,2) += RHSandDefGrad[8][xyz];

        dMatrix3x3 locDefGradDelta = (locDefGrad - EP.DeformationGradientsTotal(i,j,k));
        MAXStrainDifference = max(MAXStrainDifference, locDefGradDelta.norm());

        if(IncrementScaling != 1.0)
        {
            locDefGradDelta *= IncrementScaling;
        }

        EP.DeformationGradientsTotal(i,j,k) += locDefGradDelta;

        vStress locStress = EP.EffectiveElasticConstants(i,j,k)*EP.ElasticStrains(i,j,k);
        MAXStressDifference = max(MAXStressDifference, (EP.Stresses(i,j,k) - locStress).norm());

        if(EP.LargeDeformations)
        {
            EP.StressIncrements(i,j,k) = locStress;
        }
        else
        {
            EP.Stresses(i,j,k) = locStress;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectralImpl::ExecuteForwardFFT(void)
{
    fftw_execute(ForwardPlanRHS[0]);
    fftw_execute(ForwardPlanRHS[1]);
    fftw_execute(ForwardPlanRHS[2]);
    fftw_execute(ForwardPlanRHS[3]);
    fftw_execute(ForwardPlanRHS[4]);
    fftw_execute(ForwardPlanRHS[5]);
    fftw_execute(ForwardPlanRHS[6]);
    fftw_execute(ForwardPlanRHS[7]);
    fftw_execute(ForwardPlanRHS[8]);
}

void ElasticitySolverSpectralImpl::ExecuteBackwardFFT(void)
{
    fftw_execute(BackwardPlanDefGrad[0]);
    fftw_execute(BackwardPlanDefGrad[1]);
    fftw_execute(BackwardPlanDefGrad[2]);
    fftw_execute(BackwardPlanDefGrad[3]);
    fftw_execute(BackwardPlanDefGrad[4]);
    fftw_execute(BackwardPlanDefGrad[5]);
    fftw_execute(BackwardPlanDefGrad[6]);
    fftw_execute(BackwardPlanDefGrad[7]);
    fftw_execute(BackwardPlanDefGrad[8]);
}

void ElasticitySolverSpectralImpl::ExecuteBackwardFFTdisplacements(void)
{
    fftw_execute(BackwardPlanU[0]);
    fftw_execute(BackwardPlanU[1]);
    fftw_execute(BackwardPlanU[2]);
}

void ElasticitySolverSpectralImpl::ExecuteForwardFFTforces(void)
{
    fftw_execute(ForwardPlanForce[0]);
    fftw_execute(ForwardPlanForce[1]);
    fftw_execute(ForwardPlanForce[2]);
}

void ElasticitySolverSpectralImpl::CopyForceDensity(ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.ForceDensity, 0,)
    {
        const long int xyz = k + 2*Grid.Nz2*(j + Grid.Ny*i);

        UandForce[0][xyz] = EP.ForceDensity(i,j,k)[0];
        UandForce[1][xyz] = EP.ForceDensity(i,j,k)[1];
        UandForce[2][xyz] = EP.ForceDensity(i,j,k)[2];
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectralImpl::CopyDisplacements(ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.Displacements,0,)
    {
        const long int xyz = k + 2*Grid.Nz2*(j + Grid.Ny*i);

        EP.Displacements(i,j,k)[0] = UandForce[0][xyz];
        EP.Displacements(i,j,k)[1] = UandForce[1][xyz];
        EP.Displacements(i,j,k)[2] = UandForce[2][xyz];
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectralImpl::ApplyMechanicalBC(ElasticProperties& EP, vStress TargetStress)
{
    dMatrix6x6 AverageElasticConstants = EP.AverageElasticConstants();
    size_t size = 6;
    for(int n = 0; n < 6; n++)
    if(EP.AppliedStrainMask[n] or EP.AppliedStrainRateMask[n])
    {
        size--;
    }

    dVectorN locRHS(size);
    dMatrixNxN locLHS(size);

    int x = 0;
    for(int n = 0; n < 6; n++)
    if(!(EP.AppliedStrainMask[n] or EP.AppliedStrainRateMask[n]))
    {
        locRHS[x] = TargetStress[n];
        int y = 0;
        for(int m = 0; m < 6; m++)
        if(!(EP.AppliedStrainMask[m] or EP.AppliedStrainRateMask[m]))
        {
            locLHS(x,y) = AverageElasticConstants(n,m);
            y++;
        }
        x++;
    }

    dVectorN Solution = locLHS.inverted()*locRHS;

    x = 0;
    for(int n = 0; n < 6; n++)
    if(!(EP.AppliedStrainMask[n] or EP.AppliedStrainRateMask[n]))
    {
        EP.AverageStrain[n] += Solution[x];
        x++;
    }
    else
    {
        EP.AverageStrain[n] = EP.AppliedStrain[n];
    }

    if(EP.KeepAspectRatio)
    {
        double trace = (1.0/(Grid.Active()))*(EP.AverageStrain[0]*Grid.dNx + EP.AverageStrain[1]*Grid.dNy + EP.AverageStrain[2]*Grid.dNz);
        EP.AverageStrain[0] = trace * Grid.dNx;
        EP.AverageStrain[1] = trace * Grid.dNy;
        EP.AverageStrain[2] = trace * Grid.dNz;
    }

    if(EP.KeepAspectRatio or EP.PreventShear)
    {
        EP.AverageStrain[3] = 0.0;
        EP.AverageStrain[4] = 0.0;
        EP.AverageStrain[5] = 0.0;
    }
}
} // namespace openphase
