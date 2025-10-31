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

 *   File created :   2022
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *                         Raphael Schiedung
 *
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <complex>
#include <iostream>

#ifndef _OPENMP
    // overloading several OpenMP methods used in OpenPhase for serial operation
    inline int omp_get_num_threads()
    {
        return 1;
    }
    inline int omp_get_max_threads()
    {
        return 1;
    }
    inline int omp_get_thread_num()
    {
        return 0;
    }
    inline void omp_set_num_threads(int)
    {
        // Do nothing
    };
#endif

#ifdef _WIN32
#    define NOMINMAX
#    define and &&
#    define or ||
#    define not !
#endif

#ifdef _WIN32
#    define OP_EXPORTS __declspec(dllexport)
#    define OP_IMPORTS __declspec(dllimport)
#    ifdef OPENPHASE_EXPORTS
#        define OPENPHASE_API __declspec(dllexport)
#    else
#        define OPENPHASE_API __declspec(dllimport)
#    endif
#else
#    define OP_EXPORTS
#    define OP_IMPORTS
#    ifdef OPENPHASE_EXPORTS
#        define OPENPHASE_API
#    else
#        define OPENPHASE_API
#    endif
#endif

#ifdef MPI_PARALLEL
    extern int MPI_RANK;                                                        ///< Local MPI RANK in MPI parallel mode
    extern int MPI_SIZE;                                                        ///< Total number of MPI RANKs in MPI parallel mode
    extern int MPI_CART_RANK[3];                                                ///< Cartesian coordinates of the current process if using MPI 3D domain decomposition
    extern int MPI_CART_SIZE[3];                                                ///< Number of processes used in each direction if using MPI 3D domain decomposition
    extern bool MPI_3D_DECOMPOSITION;                                           ///< "true" if MPI should decompose in 3 dimensions
#endif

namespace openphase
{

inline void simulation_end()
{
#ifdef _WIN32
    std::getchar();                                                             ///< Prevents closing of the terminal window at the end of the simulation
#endif
}

inline void OP_Exit(int exit_code)
{
#if defined(_WIN32)
//Who put that here??
//    std::cout << "Press enter to exit program!\n";
//    std::cin.get();
#endif

#if defined(NDEBUG) // Release mode (Normal)
    std::exit(exit_code);
#else 
    if (exit_code == EXIT_SUCCESS) std::exit(EXIT_SUCCESS);
    else std::abort();
#endif
}
static constexpr int EXIT_H5_ERROR = 5;                                         ///< OpenPhase not compiled for HDF5 support
static constexpr int EXIT_NUCMODE_ERROR = 16;
static constexpr int EXIT_CONTROLMODE_ERROR = 10;
static constexpr int EXIT_LATENTHEATMODE_ERROR = 13;

#ifdef M_PI
#    undef M_PI
#    undef M_PIl
#    undef M_PIf
#endif

static constexpr double Pi = 3.14159265358979323846;                            ///< Pi constant value

const std::complex< double > I(0.0, 1.0);                                       ///< sqrt(-1) declaration

// Special keywords
enum class Resolutions : int                                                    ///< Possible simulation resolutions
{
    Single,                                                                     ///< Single resolution for all fields
    Dual                                                                        ///< Double resolution for phase fields
};

enum class AdvectionSchemes : int                                               ///< Available advection schemes
{
    Upwind,
    Minmod,
//    VanLeer,
    Superbee,
    //LaxWendroff,
    MonotonizedCentral
};

enum class AggregateStates : int                                                ///< Aggregate state markers
{
    Solid,
    Liquid,
    Gas
};

enum class EventTriggers : int                                                  ///< Events triggers/conditions
{
    User,                                                                       ///< The event is manually controlled by the user
    Tmax,                                                                       ///< Trigger at a given maximum temperature approaching from below [K]
    Tmin,                                                                       ///< Trigger at a given minimum temperature approaching from above [K]
    Time,                                                                       ///< Trigger at a particular moment in time [s]
    TimeStep,                                                                   ///< Trigger at a particular time step
    PhaseFractionMax,                                                           ///< Trigger on particular phase fraction value approaching from below
    PhaseFractionMin,                                                           ///< Trigger on particular phase fraction value approaching from above
    Stress,                                                                     ///< Trigger on a given strain value
    Strain                                                                      ///< Trigger on a given stress value
};

enum class ClearingModes : int                                                  ///< Convenience feature to control major Clear() methods calls in various classes
{
    Automatic,                                                                  ///< Relevant storages are cleared automatically at the end of the time step
    Manual                                                                      ///< To clear the relevant storages Clear() method should be called explicitly
};

enum class ExtrapolationModes : int                                             ///< Property extrapolation modes. Used in Thermodynamics to accelerate the simulations
{
    None,                                                                       ///< No extrapolation, slow
    Local,                                                                      ///< Local extrapolation in every grid cell, fast
    Global                                                                      ///< Global extrapolation for entire simulation domain, fastest
};

enum class ActivationModes : int                                                ///< Activation modes for different entities in a simulation
{
    Enabled,                                                                    ///< The entity is enabled and interacts with all other entities.
    Disabled                                                                    ///< The entity is disabled in the simulation and does not require input.
};

enum class InterfaceNormalModels : int
{
    AverageGradient,                                                            ///< Use normalized average phase-field gradient for a pair of grains as their interface normal
    WeightedGradient                                                            ///< Use normalized weighted average phase-field gradient for a pair of grains as their interface normal
};

enum class LatticeSystems2D : int                                               ///< Lattice systems according to crystallography (in 2D)
{
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Hexagonal
};

enum class LatticeSystems : int                                                 ///< Lattice systems according to crystallography (in 3D)
{
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Rhombohedral,
    Trigonal = Rhombohedral,                                                    ///< Trigonal crystal system and Rhombohedral lattice system overlap but are not equal. Here they are made equal for simplicity
    Hexagonal,
    Cubic
};

enum class LatticeCentering : int                                               ///< Unit cell centering for lattice systems
{
    Simple,
    Base,
    Body,
    Face
};

enum class BravaisLattices : int                                                ///< Bravais lattices plus HCP (hexagonal close-packed) lattice
{
    Triclinic,
    Monoclinic,
    MonoclinicBaseCentered,
    Orthorhombic,
    OrthorhombicBaseCentered,
    OrthorhombicBodyCentered,
    OrthorhombicFaceCentered,
    Tetragonal,
    TetragonalBodyCentered,
    BCT = TetragonalBodyCentered,
    Rhombohedral,
    Hexagonal,
    HEX = Hexagonal,
    HexagonalClosePacked,
    HCP = HexagonalClosePacked,
    Cubic,
    CubicBodyCentered,
    BCC = CubicBodyCentered,
    CubicFaceCentered,
    FCC = CubicFaceCentered
};

enum class BravaisLattices2D : int
{
    Monoclinic,
    Orthorhombic,
    OrthorhombicCentered,
    Tetragonal,
    Hexagonal
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
inline void ignore_result(T /* unused result */) {}                             ///< Allows to suppress warnings on unused function return values

// Standard OpenPhase files and directories
#ifdef _WIN32
    const std::string dirSeparator = "\\";                                      ///< Windows style directory separator
#else
    const std::string dirSeparator = "/";                                       ///< Unix/Linux style directory separator
#endif

const std::string DefaultInputFileName = "ProjectInput.opi";
const std::string DefaultVTKDir        = "VTK" + dirSeparator;
const std::string DefaultRawDataDir    = "RawData" + dirSeparator;
const std::string DefaultTextDir       = "TextData" + dirSeparator;

}// namespace openphase
#endif
