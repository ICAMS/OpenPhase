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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Philipp Engels;
 *                         Raphael Schiedung
 *
 */

#ifndef INITIALIZATIONS_H
#define INITIALIZATIONS_H

#include "BoundaryConditions.h"
#include "Includes.h"
#include "Tools.h"
namespace openphase
{

class BoundaryConditions;
class ElasticProperties;
class Orientations;
class PhaseField;
class Settings;

class OP_EXPORTS Initializations// : public OPObject
{
 public:
    static size_t Single(PhaseField& Phase, size_t PhaseIndex,
            const BoundaryConditions& BC);                                      ///< Initializes a single space filling phase-field
    static std::vector<iVector3> QuasiRandomNuclei(PhaseField& Phase,
            size_t phaseIndex,
            iVector3 offset, iVector3 spacing, iVector3 deviation,
            double threshold = 0.05, int seed = 1);                             ///< Plants grain nuclei on a regular grid with the given "offset" from the origin, "spacing" and local random "deviation" from the ideal grid position for each seed
    static std::vector<iVector3> QuasiRandomSpheres(PhaseField& Phase,
            BoundaryConditions& BC, int phaseIndex1,
            int phaseIndex2, int dist, double radius1, double radius2,
            double probabilityPhase1, int offset = 0, int seed = 1);            ///< Populates spherical grains of two different phases on a regular grid with grid spacing "distance" and local random "offset" of each sphere
    static void RandomNuclei(PhaseField& Phase, const Settings& locSettings,
            size_t phaseIndex, size_t Nparticles, bool randomOrientation = true,
            bool randomVariants = false, std::string onPlane = "No", int seed = 1);
    static size_t Ellipsoid(PhaseField& Phase, size_t PhaseIndex,
            double RadiusX, double RadiusY, double RadiusZ,
            double x0, double y0, double z0, BoundaryConditions& BC);                     ///< Initializes a spherical grain
    static size_t SectionalPlane(PhaseField& Phase, const size_t PhaseIndex,
            const dVector3 Point, const dVector3 Orientation,
            const BoundaryConditions& BC,
            const bool NewGrain = true, const size_t FieldIndex = 0);           ///< Initializes a new phase on the negative side of a sectional plane
    static size_t Layer(PhaseField& Phase, size_t PhaseIndex,
            const dVector3 position, const dVector3 orientation,
            double thickness, BoundaryConditions& BC);   ///< Initializes a layer with thickness of PhaseIndex with orientation at position
    static size_t Sphere(PhaseField& Phase, const size_t PhaseIndex,
            const double Radius, double x0, double y0, double z0,
            const BoundaryConditions& BC, const bool Finalize = true);                                        ///< Initializes a spherical grain
    static std::vector<size_t> Fractional(PhaseField& Phase,
            size_t MajorityPhaseIndex, size_t MinorityPhaseIndex,
            double MinorityPhaseLayerThickness, BoundaryConditions& BC);                                             ///< ?? TODO
    static std::vector<size_t> ThreeFractionals(PhaseField& Phase,
            size_t MajorityPhaseIndex, double MajorityPhaseLayerThickness,
            size_t MinorityPhaseIndex1, double MinorityPhaseLayerThickness1,
            size_t MinorityPhaseIndex2, BoundaryConditions& BC);                                             ///< ?? TODO
    static std::vector<size_t> TwoDifferentWalls(PhaseField& Phase,
            size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
            double WallsThickness, BoundaryConditions& BC);
    static std::vector<size_t> TwoWalls(PhaseField& Phase,
            size_t ChannelPhaseIndex, size_t WallsPhaseIndex,
            double WallsThickness, BoundaryConditions& BC);
    static size_t Rectangular(PhaseField& Phase, const size_t PhaseIndex,
            const double Lx, const double Ly, const double Lz,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC, const bool Finalize = true);                                        ///< Initializes a rectangular grain
    static size_t Cylinder(PhaseField& Phase, const size_t PhaseIndex,
            const double Radius, const double length, const int Axis,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC);         ///< Initializes a cylindrical grain
    // Special initialization functions:
    static size_t  FillGrainWithSpheres(PhaseField& Phase,
            size_t ParrentPhaseFieldIndex, size_t SpheresPhaseIndex,
            double MinR, double MaxR, BoundaryConditions& BC,
            size_t Nspheres = 0, double MinDistance = -1);                                           ///< Fills grain "loosely" with spheres spaced MinDistance apart
    /// Fills rectangular domain densely with spherical grains of normal distributed size.
    static std::vector<size_t> FillRectangularWithSpheres(
            PhaseField& Phase,
            const BoundaryConditions& BC, const Settings& locSettings,
            std::function<size_t(long int, long int, long int)> PhaseIndex,     ///< function that determines the phases of the grain which may depend on the position (i,j,k).
            double MeanRadius,                                                  ///< Mean radius of grains
            double StdRadius,                                                   ///< Standard deviation grain radius
            long int x_min,                                                     ///< First x-index of rectangular
            long int x_max,                                                     ///< Last x-index of rectangular
            long int y_min,                                                     ///< First y-index of rectangular
            long int y_max,                                                     ///< Last y-index of rectangular
            long int z_min,                                                     ///< First a-index of rectangular
            long int z_max,                                                     ///< Last z-index of rectangular
            double MinDistance,                                                 ///< Minimal distance of grains which can be negative for even denser microstructure!
            size_t RandomSeedInput,                                             ///< Random number seed for initialization
            bool Verbose = true
            );                                                                  ///< Fills rectangular domain densely with spherical grains of normal distributed size.
    /// Fills rectangular domain densely with spherical grains of normal distributed size with specified relative density
    static std::vector<size_t> RectangularGreenBody(
            PhaseField& Phase,
            const BoundaryConditions& BC,const Settings& locSettings,
            std::function<size_t(long, long, long)> SolidPhaseIndex,            ///< function that determines the phases of the grain which may depend on the position (i,j,k).
            size_t FluidPhaseIndex,                                             ///< Index of fluid phase
            long SizeX,                                                         ///< First x-index of rectangular
            long SizeY,                                                         ///< First y-index of rectangular
            long SizeZ,                                                         ///< First a-index of rectangular
            double MeanRadius,                                                  ///< Mean radius of grains
            double StdRadius,                                                   ///< Standard deviation grain radius
            size_t RandomSeed,                                                  ///< Random number seed for initialization
            double RelativeDensity,                                             ///< Relative Density of initialized green body
            double accuracy = 0.001,                                            ///< Accuracy of relative Density of initialized green body
            size_t MaxItrations = 100                                           ///< Maximum number of attempt to initialize green body of specified density
            );
    static size_t SphereInGrain(PhaseField& Phase,
            const size_t ParrentPhaseFieldIndex,
            const size_t PhaseIndex, const double Radius,
            const double x0, const double y0, const double z0,
            const BoundaryConditions& BC,
            const bool Finalize = true);
    static int TwoDimEBSD(std::string filename, std::vector<int> columns,
            PhaseField& Phase, BoundaryConditions& BC);  ///< Reads microstructure from EBSD result files
    static int TwoDimEBSDWithOrientations(std::string filename,
            std::vector<int> columns, std::string anglerepresentation,
            PhaseField& Phase, ElasticProperties& EP, Orientations& OR,
            BoundaryConditions& BC);                     ///< Reads microstructure from EBSD result files
    static void VoronoiTessellation(PhaseField& Phase, BoundaryConditions& BC,
            const size_t seedpoints, const size_t matrixphase);                 ///< Creates Voronoi grain structures using the algorithm based on coordination shells.
    static void TripleJunction(PhaseField& Phase, size_t PhaseIndex,
                BoundaryConditions& BC);
    static std::vector<size_t> Young3(PhaseField& Phase, size_t alpha,
            size_t beta, size_t gamma, size_t delta, BoundaryConditions& BC);
    static std::vector<size_t> Young4(PhaseField& Phase, size_t alpha,
            size_t beta, size_t gamma, size_t delta,
            BoundaryConditions& BC);
    static void Young4Periodic(PhaseField& Phase, size_t PhaseIndex,
            BoundaryConditions& BC);
    static std::vector<size_t> ThermalGrooving(PhaseField& Phase,
            size_t PhaseIndex1, size_t PhaseIndex2, BoundaryConditions& BC);

    static void Read(PhaseField& Phase, std::string FileInputName,
            BoundaryConditions& BC);
    static void ReadCSV(PhaseField& Phase, BoundaryConditions& BC,
            std::filesystem::path FilePath, char Separator = ',');

    /// Calculates global coordinates from local mpi coordinates
    template<typename T>
    static inline void ApplyPeriodicBoundaryConditionsOnCoordinates
        (T& i, T& j, T& k, const BoundaryConditions& BC, const GridParameters& Dimensions)
    {
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicX or BC.BCNX == BoundaryConditionTypes::Periodic) while (i >= Dimensions.TotalNx) i -= Dimensions.TotalNx;
        if (BC.MPIperiodicY or BC.BCNY == BoundaryConditionTypes::Periodic) while (j >= Dimensions.TotalNy) j -= Dimensions.TotalNy;
        if (BC.MPIperiodicZ or BC.BCNZ == BoundaryConditionTypes::Periodic) while (k >= Dimensions.TotalNz) k -= Dimensions.TotalNz;
        if (BC.MPIperiodicX or BC.BC0X == BoundaryConditionTypes::Periodic) while (i <  0) i += Dimensions.TotalNx;
        if (BC.MPIperiodicY or BC.BC0Y == BoundaryConditionTypes::Periodic) while (j <  0) j += Dimensions.TotalNy;
        if (BC.MPIperiodicZ or BC.BC0Z == BoundaryConditionTypes::Periodic) while (k <  0) k += Dimensions.TotalNz;
#else
        if (BC.BCNX == BoundaryConditionTypes::Periodic) while (i >= Dimensions.Nx) i -= Dimensions.Nx;
        if (BC.BCNY == BoundaryConditionTypes::Periodic) while (j >= Dimensions.Ny) j -= Dimensions.Ny;
        if (BC.BCNZ == BoundaryConditionTypes::Periodic) while (k >= Dimensions.Nz) k -= Dimensions.Nz;
        if (BC.BC0X == BoundaryConditionTypes::Periodic) while (i <  0) i += Dimensions.Nx;
        if (BC.BC0Y == BoundaryConditionTypes::Periodic) while (j <  0) j += Dimensions.Ny;
        if (BC.BC0Z == BoundaryConditionTypes::Periodic) while (k <  0) k += Dimensions.Nz;
#endif
    };

    /// Checks if (i,j,k) is the
    template<typename T, typename field_t>
    static inline bool CoordinatesInBoundaries(T& i, T& j, T& k, const field_t& field, const GridParameters& Dimensions)
    {
#ifdef MPI_PARALLEL
        if ((i < Dimensions.OffsetX) or (i >= Dimensions.OffsetX + field.sizeX())) return false;
        if ((j < Dimensions.OffsetY) or (j >= Dimensions.OffsetY + field.sizeY())) return false;
        if ((k < Dimensions.OffsetZ) or (k >= Dimensions.OffsetZ + field.sizeZ())) return false;
#else
        if ((i < 0) or (i >= field.sizeX())) return false;
        if ((j < 0) or (j >= field.sizeY())) return false;
        if ((k < 0) or (k >= field.sizeZ())) return false;
#endif
        return true;
    }

    /// Loops over all points (i,j,k) of a sphere with radius at (x0,y0,z0) and executes func(i,j,k,radius)
    template <class T, size_t Rank>
    static void loop_sphere(Storage3D<T,Rank>& field,
            const std::function<bool(long int, long int, long int, double)>& func,
            const double x0, const double y0, const double z0,
            const double radius, GridParameters& Dimensions, const BoundaryConditions& BC)
    {
        // <<< It is safe to put sphere with its center coordinates outside
        // <<< of the simulation domain and its radius bigger than system
        // <<< dimensions thus not need for assertion.
        //assert(x0 >= 0);
        //assert(y0 >= 0);
        //assert(z0 >= 0);
        //assert(x0 < locSettings.TotalNx);
        //assert(y0 < locSettings.TotalNy);
        //assert(z0 < locSettings.TotalNz);
        //assert(radius < std::max(locSettings.TotalNx,std::max(locSettings.TotalNy,locSettings.TotalNz)));

        dVector3 pos0 ({double(x0),double(y0),double(z0)});
        const long int iRadius = std::ceil(radius);
        const long int ii_min = (field.dNx() != 0 and field.sizeX() > 1) ? -std::min(iRadius,field.sizeX()) : 0;
        const long int ii_max = (field.dNx() != 0 and field.sizeX() > 1) ?  std::min(iRadius,field.sizeX()) : 0;
        const long int jj_min = (field.dNy() != 0 and field.sizeY() > 1) ? -std::min(iRadius,field.sizeY()) : 0;
        const long int jj_max = (field.dNy() != 0 and field.sizeY() > 1) ?  std::min(iRadius,field.sizeY()) : 0;
        const long int kk_min = (field.dNz() != 0 and field.sizeZ() > 1) ? -std::min(iRadius,field.sizeZ()) : 0;
        const long int kk_max = (field.dNz() != 0 and field.sizeZ() > 1) ?  std::min(iRadius,field.sizeZ()) : 0;

        std::vector<std::array<long int,3>> list;

        #pragma omp parallel for collapse(3) //NOTE has crashes MPI calculation in older version of loop!
        for (long int ii = ii_min; ii <= ii_max; ++ii)
        for (long int jj = jj_min; jj <= jj_max; ++jj)
        for (long int kk = kk_min; kk <= kk_max; ++kk)
        {
            // Calculate global coordinates (i,j,k)
            long int i = std::round(x0 + ii);
            long int j = std::round(y0 + jj);
            long int k = std::round(z0 + kk);

            ApplyPeriodicBoundaryConditionsOnCoordinates(i,j,k,BC,Dimensions);
            if(not CoordinatesInBoundaries(i,j,k,field,Dimensions)) continue;

            // Calculate distance between (i,j,k) and (x0,y0,z0)
            dVector3 pos ({double(i),double(j),double(k)});
            double rad = Tools::Distance<dVector3>(pos, pos0, Dimensions.TotalNx, Dimensions.TotalNy, Dimensions.TotalNz, BC).abs();
            if (rad > radius) continue;

            // Convert global coordinates to local block coordinates
#ifdef MPI_PARALLEL
            i -= Dimensions.OffsetX;
            j -= Dimensions.OffsetY;
            k -= Dimensions.OffsetZ;
#endif
            func(i,j,k,rad);
        }
    }

    /// Loops over all points (i,j,k) of a sphere with radius at (x0,y0,z0) and executes func(i,j,k,radius)
    /// Exits loop when func(i,j,k,radius) return true
    template <class T, size_t Rank>
    static bool loop_sphere_with_exit(Storage3D<T,Rank>& field,
            const std::function<bool(long int, long int, long int, double)>& func,
            const double x0, const double y0, const double z0,
            const double radius, const GridParameters& Dimensions, const BoundaryConditions& BC)
    {
        assert(x0 >= 0);
        assert(y0 >= 0);
        assert(z0 >= 0);
        assert(x0 < Dimensions.TotalNx);
        assert(y0 < Dimensions.TotalNy);
        assert(z0 < Dimensions.TotalNz);
        assert(radius < std::max(Dimensions.TotalNx,std::max(Dimensions.TotalNy,Dimensions.TotalNz)));

        dVector3 pos0 ({double(x0),double(y0),double(z0)});
        const long int iRadius = std::ceil(radius);
        const long int ii_min = (field.dNx() != 0) ? -iRadius : 0;
        const long int ii_max = (field.dNx() != 0) ?  iRadius : 0;
        const long int jj_min = (field.dNy() != 0) ? -iRadius : 0;
        const long int jj_max = (field.dNy() != 0) ?  iRadius : 0;
        const long int kk_min = (field.dNz() != 0) ? -iRadius : 0;
        const long int kk_max = (field.dNz() != 0) ?  iRadius : 0;

        int loc_exit_loop = 0; // is used to exit loop if a certain condition is met e.g. the presence of an existing phase
        #pragma omp parallel for shared(loc_exit_loop) collapse(3) //NOTE has crashes MPI calculation in older version of loop!
        for (long int ii = ii_min; ii <= ii_max; ++ii)
        for (long int jj = jj_min; jj <= jj_max; ++jj)
        for (long int kk = kk_min; kk <= kk_max; ++kk)
        {
            if (loc_exit_loop) continue;
            // Calculate global coordinates (i,j,k)
            long int i = std::round(x0 + ii);
            long int j = std::round(y0 + jj);
            long int k = std::round(z0 + kk);

            ApplyPeriodicBoundaryConditionsOnCoordinates(i,j,k,BC,Dimensions);
            if(not CoordinatesInBoundaries(i,j,k,field,Dimensions)) continue;

            // Calculate distance between (i,j,k) and (x0,y0,z0)
            dVector3 pos ({double(i),double(j),double(k)});
            double rad = Tools::Distance<dVector3>(pos, pos0, Dimensions.TotalNx, Dimensions.TotalNy, Dimensions.TotalNz, BC).abs();
            if (rad > radius) continue;

            // Convert global coordinates to local block coordinates
#ifdef MPI_PARALLEL
            i -= Dimensions.OffsetX;
            j -= Dimensions.OffsetY;
            k -= Dimensions.OffsetZ;
#endif
            if (func(i,j,k,rad)) loc_exit_loop++;
        }
        int exit_loop = 0;
#ifdef MPI_PARALLEL
        OP_MPI_Allreduce(&loc_exit_loop, &exit_loop, 1, OP_MPI_INT, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#else
        exit_loop += loc_exit_loop;
#endif

        if (exit_loop) return true;
        else return false;
    }

    static constexpr auto thisclassname = "Initializations";

 private:
    static void SphereFixedIdx(PhaseField& Phase, size_t PhaseIndex,
            double Radius, double x0, double y0, double z0,
            BoundaryConditions& BC);
   };
}// namespace openphase
#endif
