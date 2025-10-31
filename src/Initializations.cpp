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

#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "Includes.h"
#include "Initializations.h"
#include "GrainsProperties.h"
#include "InterfaceProperties.h"
#include "ElasticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Tools.h"
#include "Tools/AnalysisSintering.h"
#include "NumericalMethods/RootFindingAlgorithms.h"
#include "UserDrivingForce.h"
#include <cwctype>
#include <set>

/***************************************************************/

namespace openphase
{
using namespace std;

vector<iVector3> Initializations::QuasiRandomNuclei(PhaseField& Phase,
                                                    size_t phaseIndex,
                                                    iVector3 offset,
                                                    iVector3 spacing,
                                                    iVector3 deviation,
                                                    double threshold,
                                                    int seed)
{
    vector<iVector3> result;

    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
#ifdef MPI_PARALLEL
        srand(pow(seed + MPI_RANK,(MPI_RANK + 1)));
#else
        srand(seed);
#endif
    }

    int Nx = Phase.Grid.TotalNx;
    int Ny = Phase.Grid.TotalNy;
    int Nz = Phase.Grid.TotalNz;

    int distx = (offset[0] < Nx and spacing[0] < Nx) ? offset[0] : 0;
    int disty = (offset[1] < Ny and spacing[1] < Ny) ? offset[1] : 0;
    int distz = (offset[2] < Nz and spacing[2] < Nz) ? offset[2] : 0;

    for (int i = distx; i < Nx; i += spacing[0])
    for (int j = disty; j < Ny; j += spacing[1])
    for (int k = distz; k < Nz; k += spacing[2])
    {
        double chance = double(rand()) / double(RAND_MAX);
#ifdef MPI_PARALLEL
        OP_MPI_Bcast(&(chance), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
#endif
        if (chance > threshold)
        {
            int di = (spacing[0] >= Nx or deviation[0] == 0) ? 0 : (rand() % (2 * deviation[0]) - deviation[0]);
            int dj = (spacing[1] >= Ny or deviation[1] == 0) ? 0 : (rand() % (2 * deviation[1]) - deviation[1]);
            int dk = (spacing[2] >= Nz or deviation[2] == 0) ? 0 : (rand() % (2 * deviation[2]) - deviation[2]);
#ifdef MPI_PARALLEL
            OP_MPI_Bcast(&(di), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
            OP_MPI_Bcast(&(dj), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
            OP_MPI_Bcast(&(dk), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
#endif
            if (i + di >= 0 and i + di < Nx and
                j + dj >= 0 and j + dj < Ny and
                k + dk >= 0 and k + dk < Nz)
            {
                Phase.PlantGrainNucleus(phaseIndex, i + di, j + dj, k + dk);

                iVector3 temp;
                temp[0] = i + di;
                temp[1] = j + dj;
                temp[2] = k + dk;
                result.push_back(temp);
            }
        }
    }
    return result;
}

vector<iVector3> Initializations::QuasiRandomSpheres(PhaseField& Phase,
                                                     BoundaryConditions& BC,
                                                     int phaseIndex1,
                                                     int phaseIndex2,
                                                     int dist,
                                                     double radius1,
                                                     double radius2,
                                                     double probabilityPhase1,
                                                     int offset,
                                                     int seed)
{
    vector<iVector3> result;

    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
#ifdef MPI_PARALLEL
        srand(pow(seed + MPI_RANK,(MPI_RANK + 1)));
#else
        srand(seed);
#endif
    }

    int Nx = Phase.Grid.TotalNx;
    int Ny = Phase.Grid.TotalNy;
    int Nz = Phase.Grid.TotalNz;

    int distx = (dist < Nx) ? dist : 0;
    int disty = (dist < Ny) ? dist : 0;
    int distz = (dist < Nz) ? dist : 0;
    double threshold = 0.05;

    for (int i = distx; i < Nx; i += 2 * dist)
    for (int j = disty; j < Ny; j += 2 * dist)
    for (int k = distz; k < Nz; k += 2 * dist)
    {
        double chance = double(rand()) / double(RAND_MAX);
#ifdef MPI_PARALLEL
        OP_MPI_Bcast(&(chance), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
#endif
        if (chance > threshold)
        {
            int di = (distx == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dj = (disty == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);
            int dk = (distz == 0 or offset == 0) ? 0 : (rand() % (2 * offset) - offset);

#ifdef MPI_PARALLEL
            OP_MPI_Bcast(&(di), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
            OP_MPI_Bcast(&(dj), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
            OP_MPI_Bcast(&(dk), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
#endif
            if (i + di >= 0 and i + di < Nx and
                j + dj >= 0 and j + dj < Ny and
                k + dk >= 0 and k + dk < Nz)
            {
                double randPhase = double(rand()) / double(RAND_MAX);
                if (randPhase < probabilityPhase1)
                {
                    Initializations::Sphere(Phase, phaseIndex1, radius1, i+di, j+dj, k+dk, BC);
                }
                else
                {
                    Initializations::Sphere(Phase, phaseIndex2, radius2, i+di, j+dj, k+dk, BC);
                }

                iVector3 temp;
                temp[0] = i + di;
                temp[1] = j + dj;
                temp[2] = k + dk;

                result.push_back(temp);
            }
        }
    }
    return result;
}

void Initializations::RandomNuclei(PhaseField& Phase,
                                   const Settings& locSettings,
                                   size_t phaseIndex,
                                   size_t Nparticles,
                                   bool randomOrientation,
                                   bool randomVariants,
                                   string onPlane,
                                   int seed)
{
    if (seed == -1)
    {
        srand(time(NULL));
    }
    else
    {
#ifdef MPI_PARALLEL
        srand(pow(seed + MPI_RANK,(MPI_RANK + 1)));
#else
        srand(seed);
#endif
    }

    default_random_engine OrientationsGenerator(2*seed);
    default_random_engine VariantsGenerator(3*seed);

    for (size_t n = 0; n < Nparticles; n++)
    {
        bool planted = false;
        int iterations = 0;
        while(!planted and iterations < 1000)
        {
            iterations++;

            int di = 0;
            if     (onPlane == "Xbottom") {di = 0;}
            else if(onPlane == "Xtop")    {di = Phase.Grid.Nx - 1;}
            else                          {di = rand() % Phase.Grid.Nx;};

            int dj = 0;
            if     (onPlane == "Ybottom") {dj = 0;}
            else if(onPlane == "Ytop")    {dj = Phase.Grid.Ny - 1;}
            else                          {dj = rand() % Phase.Grid.Ny;};

            int dk = 0;
            if     (onPlane == "Zbottom") {dk = 0;}
            else if(onPlane == "Ztop")    {dk = Phase.Grid.Nz - 1;}
            else                          {dk = rand() % Phase.Grid.Nz;};

            bool freecell = true;

            for (int ii = - Phase.Grid.iWidth; ii <= Phase.Grid.iWidth; ii++)
            for (int jj = - Phase.Grid.iWidth; jj <= Phase.Grid.iWidth; jj++)
            for (int kk = - Phase.Grid.iWidth; kk <= Phase.Grid.iWidth; kk++)
            if ((di+ii > -Phase.Fields.Bcells() and di+ii < Phase.Grid.Nx + Phase.Fields.Bcells() and
                 dj+jj > -Phase.Fields.Bcells() and dj+jj < Phase.Grid.Ny + Phase.Fields.Bcells() and
                 dk+kk > -Phase.Fields.Bcells() and dk+kk < Phase.Grid.Nz + Phase.Fields.Bcells()) and
                Phase.Fields(di+ii,dj+jj,dk+kk).flag)
            {
                freecell = false;
                break;
            }
            if(freecell)
            {
                int locIndex = Phase.PlantGrainNucleus(phaseIndex, di, dj, dk);
                if(randomOrientation)
                {
                    uniform_real_distribution<double> loc_orientation(0, 180);
                    int Q1 = loc_orientation(VariantsGenerator);
                    int Q2 = loc_orientation(VariantsGenerator);
                    int Q3 = loc_orientation(VariantsGenerator);
                    EulerAngles locAngles({Q1*Pi/180, Q2*Pi/180, Q3*Pi/180}, XYZ);

                    Phase.FieldsProperties[locIndex].Orientation = locAngles.getQuaternion();
                }
                if(randomVariants)
                {
                    uniform_int_distribution<int> loc_variant(0, locSettings.Nvariants[phaseIndex]-1);
                    Phase.FieldsProperties[locIndex].Variant = loc_variant(VariantsGenerator);
                }
                planted = true;
            }
        }
    }
}

/**
 * Initializes a new phase on the negative side of a plane
 * @param Phase       Phase field
 * @param PhaseIndex  Phase index
 * @param Point       Point on plane
 * @param Orientation Plane normal
 * @param BC          Boundary conditions
 * @param locSettings Project settings
 */
size_t Initializations::SectionalPlane(PhaseField& Phase,
                                       const size_t PhaseIndex,
                                       const dVector3 Point,
                                       const dVector3 Orientation,
                                       const BoundaryConditions& BC,
                                       const bool NewGrain,
                                       const size_t FieldIndex)
{
    const dVector3 n      = Orientation.normalized();
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    int locIndex;
    if (NewGrain)
    {
         locIndex = Phase.AddGrainInfo(PhaseIndex);
    }
    else
    {
         locIndex = FieldIndex;
    }

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        const dVector3 x = {double(i),double(j),double(k)}; // position in space
        const double distance = (x-Point)*n;   // distance from plane

        if (distance < - 0.5*iWidth)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(locIndex, 1.0);
        }
        else if (distance <= 0.5*iWidth)
        {
            const double IntProf = 0.5 - 0.5*std::sin(Pi*(distance)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    return locIndex;
}

size_t Initializations::Sphere(PhaseField& Phase,
                               const size_t PhaseIndex,
                               const double Radius,
                               double x0, double y0, double z0,
                               const BoundaryConditions& BC,
                               const bool Finalize)
{
    // <<< It is safe to put sphere with its center coordinates outside
    // <<< of the simulation domain and its radius bigger than system
    // <<< dimensions thus not need for assertion.    //assert(x0 >= 0);
    //assert(y0 >= 0);
    //assert(z0 >= 0);
    //assert(x0 < locSettings.TotalNx);
    //assert(y0 < locSettings.TotalNy);
    //assert(z0 < locSettings.TotalNz);
    //assert(Radius < std::max(locSettings.TotalNx,std::max(locSettings.TotalNy,locSettings.TotalNz)));

    if (Phase.Grid.dNx == 0) x0 = 0.0;
    if (Phase.Grid.dNy == 0) y0 = 0.0;
    if (Phase.Grid.dNz == 0) z0 = 0.0;

    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    //assert(Phase.FieldsProperties[locIndex].is_solid());
    auto set_sphere = [&iWidth,&Radius,&locIndex,&Phase](long int i,long int j,long int k, double rad)
    {
        if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set_value(locIndex, 1.0);
        }
        else if (rad <= Radius + iWidth*0.5)
        {
            const double IntProf = 0.5 - 0.5*std::sin(Pi*(rad - Radius)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, Phase.Grid, BC);

    if (Finalize) Phase.FinalizeInitialization(BC);
    return locIndex;
}

void Initializations::Young4Periodic(PhaseField& Phase,
                                     size_t PhaseIndex,
                                     BoundaryConditions& BC)
{
    int Nx = Phase.Grid.Nx;
    int Ny = Phase.Grid.Ny;
    int Nz = Phase.Grid.Nz;

    size_t locIndex0 = Phase.AddGrainInfo(PhaseIndex);
    size_t locIndex1 = Phase.AddGrainInfo(PhaseIndex+1);
    size_t locIndex2 = Phase.AddGrainInfo(PhaseIndex+2);
    size_t locIndex3 = Phase.AddGrainInfo(PhaseIndex+3);
    size_t locIndex4 = Phase.AddGrainInfo(PhaseIndex+4);
    size_t locIndex5 = Phase.AddGrainInfo(PhaseIndex+5);
    size_t locIndex6 = Phase.AddGrainInfo(PhaseIndex+6);
    size_t locIndex7 = Phase.AddGrainInfo(PhaseIndex+7);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i,j,k).flag = 2;

        if (j <= Ny/2)
        {
            if      (i <= Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex1,1);
            else if (i >  Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex2,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k <= Nz/4)
                Phase.Fields(i, j, k).set_value(locIndex3,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k >  Nz*3./4)
                Phase.Fields(i, j, k).set_value(locIndex3,1);
            else Phase.Fields(i,j,k).set_value(locIndex0,1);
        }
        else
        {
            if      (i <= Nx/8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex5,1);
            else if (i >  Nx/8 and i <= Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex4,1);
            else if (i >  Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex5,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k <= Nz/8)
                Phase.Fields(i, j, k).set_value(locIndex6,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k >  Nz*5./8)
                Phase.Fields(i, j, k).set_value(locIndex6,1);
            else Phase.Fields(i,j,k).set_value(locIndex7,1);
        }
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
}

void Initializations::TripleJunction(PhaseField& Phase,
                                     size_t PhaseIndex,
                                     BoundaryConditions& BC)
{
    int Nx = Phase.Grid.Nx;
    int Nz = Phase.Grid.Nz;

    Rectangular(Phase, PhaseIndex  , 2.0*Nx/3.0, 0, Nz/2.0, 1.0*Nx/3.0, 0,     Nz/4.0,  BC, false);
    Rectangular(Phase, PhaseIndex+1, 2.0*Nx/3.0, 0, Nz/2.0, 1.0*Nx/3.0, 0, 3.0*Nz/4.0,  BC, false);
    Rectangular(Phase, PhaseIndex+2,     Nx/3.0, 0, Nz    , 5.0*Nx/6.0, 0,     Nz/2.0,  BC, false);

    Phase.FinalizeInitialization(BC);
}

vector<size_t> Initializations::Young4(PhaseField& Phase,
                                       size_t alpha, size_t beta,
                                       size_t gamma, size_t delta,
                                       BoundaryConditions& BC)
{
    vector<size_t> result;

    int Nx = Phase.Grid.Nx;
    int Ny = Phase.Grid.Ny;
    int Nz = Phase.Grid.Nz;

    result.push_back(Single     (Phase,alpha,                                                                                                                                BC));
    result.push_back(Rectangular(Phase,beta , 3.0*Nx/4.0+Phase.Grid.iWidth,    Ny/3.0+Phase.Grid.iWidth,Nz+2*Phase.Grid.iWidth,5.0*Nx/8.0,    Ny/6.0-Phase.Grid.iWidth,Nz/2, BC, false));
    result.push_back(Rectangular(Phase,gamma, 3.0*Nx/4.0+Phase.Grid.iWidth,2.0*Ny/3.0+Phase.Grid.iWidth,Nz/2+Phase.Grid.iWidth,5.0*Nx/8.0,4.0*Ny/6.0,                  Nz/4, BC, false));
    result.push_back(Rectangular(Phase,delta, 3.0*Nx/4.0+Phase.Grid.iWidth,2.0*Ny/3.0+Phase.Grid.iWidth,Nz/2+Phase.Grid.iWidth,5.0*Nx/8.0,4.0*Ny/6.0,              3.0*Nz/4, BC, false));

    Phase.FinalizeInitialization(BC);
    return result;
}

vector<size_t> Initializations::Young3(PhaseField& Phase,
                                       size_t alpha, size_t beta,
                                       size_t gamma, size_t delta,
                                       BoundaryConditions& BC)
{
    vector<size_t> result;

    int Nx = Phase.Grid.Nx;
    int Nz = Phase.Grid.Nz;

    double x = Nx/2.0;
    double y = 0;
    double z = Nz/2.0;

    result.push_back(Single     (Phase, alpha,                   BC));
    result.push_back(Rectangular(Phase,  beta, x,y,z, 0.5*x,0,z, BC, false));
    result.push_back(Rectangular(Phase, gamma, x,y,z,     x,0,0, BC, false));
    result.push_back(Rectangular(Phase, delta, x,y,z, 1.5*x,0,z, BC, false));

    Phase.FinalizeInitialization(BC);

    return result;
}

size_t Initializations::Single(PhaseField& Phase,
                               size_t PhaseIndex,
                               const BoundaryConditions& BC)
{
    size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set_value(locIndex, 1.0);
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);

    return locIndex;
}

size_t Initializations::Layer(PhaseField& Phase,
                              size_t PhaseIndex,
                              const dVector3 position,
                              const dVector3 orientation,
                              double thickness,
                              BoundaryConditions& BC)
{
    const dVector3 n = orientation.normalized();
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    size_t idx = Phase.AddGrainInfo(PhaseIndex);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        const dVector3 x = {double(i),double(j),double(k)};
        const double distance = std::abs((x-position)*n);

        if (distance < 0.5*thickness - 0.5*iWidth)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(idx, 1.0);
        }
        else if (distance < 0.5*thickness + 0.5*iWidth)
        {
            const double IntProf = 0.5 - 0.5*std::sin(Pi*(distance - 0.5*thickness)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set_value(idx, IntProf);
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);

    return idx;
}

vector<size_t> Initializations::Fractional(PhaseField& Phase,
                                           size_t MajorityPhaseIndex,
                                           size_t MinorityPhaseIndex,
                                           double MinorityPhaseLayerThickness,
                                           BoundaryConditions& BC)
{
    if( BC.BC0Z != BoundaryConditionTypes::NoFlux || BC.BCNZ != BoundaryConditionTypes::NoFlux )
    {
        ConsoleOutput::WriteWarning(
                "Consider NoFlux Boundary condition at least along \n"
                "Z-direction for Fractional initialization", thisclassname, "Fractional");
    }

    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    size_t index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(MinorityPhaseIndex);

    int offset = MinorityPhaseLayerThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);

    return vector<size_t>({index1,index2});
}

vector<size_t> Initializations::ThreeFractionals(PhaseField& Phase,
                                 size_t MajorityPhaseIndex,
                                 double MajorityPhaseLayerThickness,
                                 size_t MinorityPhaseIndex1,
                                 double MinorityPhaseLayerThickness1,
                                 size_t MinorityPhaseIndex2,
                                 BoundaryConditions& BC)
{

    if( BC.BC0Z != BoundaryConditionTypes::NoFlux or
        BC.BCNZ != BoundaryConditionTypes::NoFlux )
    {
        ConsoleOutput::WriteWarning(
                "Consider NoFlux Boundary condition atleast along \n"
                "Z-direction for Fractional initialization", thisclassname, "Fractional");
    }

    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    size_t index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(MinorityPhaseIndex1);
    size_t index3 = Phase.AddGrainInfo(MinorityPhaseIndex2);

    int offset1 = MajorityPhaseLayerThickness;
    int offset2 = MinorityPhaseLayerThickness1;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k < offset1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if ((k > offset1 + iWidth/2) and (k < offset1+offset2 - iWidth/2))
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k > offset1+offset2 + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index3, 1);
        }
        else if ((k <= offset1 + iWidth/2) and (k >= offset1 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset1)/iWidth);
            Phase.Fields(i, j, k).set_value(index2, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if ((k <= offset1+offset2 + iWidth/2) and
                 (k >= offset1+offset2 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - (offset1+offset2))/iWidth);
            Phase.Fields(i, j, k).set_value(index3, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);

    return vector<size_t>({index1,index2,index3});
}

vector<size_t> Initializations::TwoWalls(PhaseField& Phase,
                                         size_t ChannelPhaseIndex,
                                         size_t WallsPhaseIndex,
                                         double WallsThickness,
                                         BoundaryConditions& BC)
{
    int Nz = Phase.Grid.TotalNz;
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    size_t index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(WallsPhaseIndex);

    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k + Phase.Grid.OffsetZ > offset + iWidth/2 && k + Phase.Grid.OffsetZ < Nz - offset - 1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k+ Phase.Grid.OffsetZ < offset - iWidth/2 || k+ Phase.Grid.OffsetZ > Nz - offset - 1 + iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k+ Phase.Grid.OffsetZ >= offset - iWidth/2 && k+ Phase.Grid.OffsetZ <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*((k+ Phase.Grid.OffsetZ) - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);

            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k+ Phase.Grid.OffsetZ >= Nz - offset - 1 - iWidth/2 && k+ Phase.Grid.OffsetZ <= Nz - offset - 1 + iWidth/2)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*(Nz - (k + Phase.Grid.OffsetZ)- offset - 1)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i, j, k).set_value(index2, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    return vector<size_t>({index1,index2});
}

vector<size_t> Initializations::TwoDifferentWalls(PhaseField& Phase,
                                                  size_t ChannelPhaseIndex,
                                                  size_t WallsPhaseIndex,
                                                  double WallsThickness,
                                                  BoundaryConditions& BC)
{
    int Nz = Phase.Grid.Nz;
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;

    size_t index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    size_t index2 = Phase.AddGrainInfo(WallsPhaseIndex);
    size_t index3 = Phase.AddGrainInfo(WallsPhaseIndex);

    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2 && k < (Nz) - offset - iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set_value(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set_value(index2, 1);
        }
        else if (k > (Nz) - offset + iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set_value(index3, 1);
        }
        else if (k >= offset - iWidth/2 && k <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set_value(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k >= (Nz) - offset - iWidth/2 - 1 && k <= (Nz) - offset + iWidth/2 - 1)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*((Nz) - k - offset - 1)/iWidth);
            Phase.Fields(i, j, k).set_value(index1, IntProf);
            Phase.Fields(i, j, k).set_value(index3, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    return vector<size_t>({index1,index2,index3});
}

double PinE(const double x, const double y, const double z,
            const double a, const double b, const double c)                     /// Auxiliary function, check whether a point (x,y,z) lies within ellipsoid with axes a, b, c
{
  return (x*x/a/a + y*y/b/b + z*z/c/c - 1.0);
}

size_t Initializations::Ellipsoid(PhaseField& Phase,
                                  size_t PhaseIndex,
                                  double RadiusX, double RadiusY, double RadiusZ,
                                  double x0, double y0, double z0,
                                  BoundaryConditions& BC)
{
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    size_t index = Phase.AddGrainInfo(PhaseIndex);

    ///< Inner ellipsoid
    double Rxi = RadiusX - iWidth*0.5;
    double Ryi = RadiusY - iWidth*0.5;
    double Rzi = RadiusZ - iWidth*0.5;
    ///< Outer ellipsoid
    double Rxo = RadiusX + iWidth*0.5;
    double Ryo = RadiusY + iWidth*0.5;
    double Rzo = RadiusZ + iWidth*0.5;

    double Rxi2 = Rxi*Rxi;
    double Ryi2 = Ryi*Ryi;
    double Rzi2 = Rzi*Rzi;

    double Rxo2 = Rxo*Rxo;
    double Ryo2 = Ryo*Ryo;
    double Rzo2 = Rzo*Rzo;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double radX = (i+Phase.Grid.OffsetX-x0);
        double radY = (j+Phase.Grid.OffsetY-y0);
        double radZ = (k+Phase.Grid.OffsetZ-z0);

        if ((radX*radX/Rxi2 + radY*radY/Ryi2 + radZ*radZ/Rzi2) <= 1.0)  ///< Point in the pure phase
        {
            Phase.Fields(i, j, k).clear();     ///< remove all phase fields that may already be present
            Phase.Fields(i, j, k).set_value(index, 1.0);
        }
        else if ((radX*radX/Rxo2 + radY*radY/Ryo2 + radZ*radZ/Rzo2) <= 1.0)  ///< Point in the interface
        {
            // Find coordinates inside the interface
            double r1 = -iWidth*0.5;
            double r2 =  iWidth*0.5;
            double rr = (r1 + r2)*0.5;
            double tolerance = 1e-6;
            while (fabs(r2 - r1) > tolerance)
            {
                if(PinE(radX, radY, radZ, RadiusX+r1, RadiusY+r1, RadiusZ+r1) > 0.0)
                   r1 += 0.25*(r2 - r1);
                else r1 -= 0.25*(r2 - r1);
                if(PinE(radX, radY, radZ, RadiusX+r2, RadiusY+r2, RadiusZ+r2) < 0.0)
                    r2 -= 0.25*(r2 - r1);
                else r2 += 0.25*(r2 - r1);

                rr = (r1 + r2)*0.5;
            }

            double IntProf = 0.5 - 0.5*sin(Pi*rr/iWidth);

            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); ++alpha)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set_value(index, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    return index;
}

void Initializations::Read(PhaseField& Phase, string FileName,
                           BoundaryConditions& BC)
{
    Phase.Read(FileName);

    if(Phase.Grid.Resolution == Resolutions::Single)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
        {
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha != Phase.Fields(i,j,k).end(); alpha++)
            {
               if(alpha->index >= Phase.FieldsProperties.size())
               {
                   Phase.FieldsProperties.Reallocate(alpha->index+1);
               }
               Phase.FieldsProperties[alpha->index].Exist = true;
           }
        }
        STORAGE_LOOP_END
    }

    if(Phase.Grid.Resolution == Resolutions::Dual)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,Phase.FieldsDR.Bcells())
        {
            for(auto alpha = Phase.FieldsDR(i,j,k).begin();
                     alpha != Phase.FieldsDR(i,j,k).end(); alpha++)
            {
                if(alpha->index >= Phase.FieldsProperties.size())
                {
                    Phase.FieldsProperties.Reallocate(alpha->index+1);
                }
                Phase.FieldsProperties[alpha->index].Exist = true;
            }
        }
        STORAGE_LOOP_END
    }
    Phase.Finalize(BC);
}

size_t Initializations::Rectangular(PhaseField& Phase, size_t PhaseIndex,
                                    double Lx, double Ly, double Lz,
                                    double x0, double y0, double z0,
                                    const BoundaryConditions& BC,
                                    const bool Finalize)
{
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    // Determine the area where the phase-field is one
    const double SizeX = (Lx < iWidth) ? 1: (Lx - iWidth)/2.0;
    const double SizeY = (Ly < iWidth) ? 1: (Ly - iWidth)/2.0;
    const double SizeZ = (Lz < iWidth) ? 1: (Lz - iWidth)/2.0;

    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        // Determine relative coordinates
        const double x = double(i+Phase.Grid.OffsetX);
        const double y = double(j+Phase.Grid.OffsetY);
        const double z = double(k+Phase.Grid.OffsetZ);
        const dVector3 pos  ({x , y, z});
        const dVector3 pos0 ({x0,y0,z0});
        dVector3 dist = Tools::Distance<dVector3>(pos, pos0, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC);
        dist[0] = std::abs(dist[0]) - SizeX;
        dist[1] = std::abs(dist[1]) - SizeY;
        dist[2] = std::abs(dist[2]) - SizeZ;

        // Set cubic phase field
        if (dist[0] <= 0 and dist[1] <= 0 and dist[2] <= 0)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(locIndex, 1);
        }
        else if (dist[0] <= iWidth and dist[1] <= iWidth and dist[2] <= iWidth)
        {
            const double xx = std::max(dist[0], std::max(dist[1], dist[2]));
            const double IntProf = 0.5 + 0.5 * std::cos(Pi * xx / iWidth);
            Phase.Fields(i,j,k).flag = 2;

            double tmp = IntProf;
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha != Phase.Fields(i,j,k).end(); alpha++)
            if (alpha->value > tmp)
            {
                alpha->value -= tmp;
                tmp = 0;
            }
            else
            {
                tmp -= alpha->value;
                alpha->value = 0;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END

    if (Finalize)
    {
        Phase.Finalize(BC);

        if(Phase.Grid.Resolution == Resolutions::Dual)
        {
            Phase.Refine();
        }
    }
    return locIndex;
}

size_t Initializations::Cylinder(PhaseField& Phase,
                                       const size_t PhaseIndex,
                                       const double Radius,
                                       const double length,
                                       const int Axis,
                                       const double x0,
                                       const double y0,
                                       const double z0,
                                       const BoundaryConditions& BC)
{
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        dVector3 pos__plane ({double(i+Phase.Grid.OffsetX),double(j),double(k)});
        dVector3 pos0_plane ({double(x0),double(y0),double(z0)});
        dVector3 pos__axis  ({0,0,0});
        dVector3 pos0_axis  ({0,0,0});
        pos__axis [Axis] = pos__plane[Axis];
        pos0_axis [Axis] = pos0_plane[Axis];
        pos__plane[Axis] = 0;
        pos0_plane[Axis] = 0;

        dVector3 dist_plane = Tools::Distance<dVector3>(pos__plane, pos0_plane, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC);
        dVector3 dist_axis  = Tools::Distance<dVector3>(pos__axis,  pos0_axis , Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC);

        const double rad = dist_plane.abs();
        const double z   = dist_axis.abs();

        if ((rad < Radius - iWidth*0.5) and (z < length/2 - iWidth*0.5))
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set_value(locIndex, 1.0);
        }
        else if ((rad <= Radius + iWidth*0.5) and (z <= length/2 + iWidth*0.5))
        {
            double IntProf = 0.0;
            if (z >= length/2 - iWidth*0.5 and ((rad - Radius) <= (z - length/2)))
            {
                IntProf = 0.5 - 0.5*std::sin(Pi*(z - length/2)/iWidth);
            }
            if (rad >= Radius - iWidth*0.5 and ((rad - Radius) > (z - length/2)))
            {
                IntProf = 0.5 - 0.5*std::sin(Pi*(rad - Radius)/iWidth);
            }

            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i,j,k).set_value(locIndex, IntProf);
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    return locIndex;
}

void Initializations::SphereFixedIdx(PhaseField& Phase,
                                     size_t PhaseIndex,
                                     double Radius,
                                     double x0, double y0, double z0,
                                     BoundaryConditions& BC)
{
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    Phase.FieldsProperties[PhaseIndex].Exist = true;
    Phase.FieldsProperties[PhaseIndex].Stage = GrainStages::Stable;

    auto set_sphere = [&iWidth,&Radius,&PhaseIndex,&Phase,&z0](long int i,long int j,long int k, double rad)
    {
       if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set_value(PhaseIndex, 1.0);
        }
        else if (rad < Radius + iWidth*0.5)
        {
            const double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
            if(IntProf > Phase.Fields(i, j, z0).get_value(PhaseIndex))
            {
                for(auto alpha = Phase.Fields(i,j,k).begin();
                         alpha < Phase.Fields(i,j,k).end(); alpha++)
                {
                    alpha->value *= 1.0 - IntProf;
                }
                Phase.Fields(i, j, k).set_value(PhaseIndex, IntProf);
                Phase.Fields(i,j,k).flag = 2;
            }
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, Phase.Grid, BC);

    Phase.FinalizeInitialization(BC);
}

int Initializations::TwoDimEBSD(std::string filename, std::vector<int> columns,
        PhaseField& Phase, BoundaryConditions& BC)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    if(columns.size() != 3)
    {
        columns.clear(); columns = {1, 2, 3};
        ConsoleOutput::WriteWarning("Column list incomplete. Set to default {1, 2, 3}.", thisclassname, "TwoDimEBSD");
    }

    if(Phase.Grid.Nz != 1)
    {
        ConsoleOutput::WriteExit("Only 2d supported", thisclassname, "TwoDimEBSD");
        OP_Exit(EXIT_FAILURE);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + filename + " could not be opened",
                thisclassname, "TwoDimEBSD()");
        OP_Exit(EXIT_FAILURE);
    };

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("EBSD reader");
    ConsoleOutput::WriteStandard("Source", filename);
    ConsoleOutput::WriteStandard("Read x coordinates from column", columns[0]);
    ConsoleOutput::WriteStandard("Read y coordinates from column", columns[1]);
    ConsoleOutput::WriteStandard("Read phase indices from column", columns[2]);

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    ConsoleOutput::WriteStandard("Number of lines", nol);

    // Read input file

    vector<size_t> phaseindex;
    vector<size_t> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }
            if (str.size() > 0) spos++;
            if( !iline ) break;
        }
    }

    size_t individualgrainnum = individualphaseindices.size();
    ConsoleOutput::WriteStandard("Number of individual phases", individualgrainnum);
    size_t maxphaseindex = *max_element(individualphaseindices.begin(), individualphaseindices.end());
    ConsoleOutput::WriteStandard("Largest phase-field index", maxphaseindex);
    size_t pfsize = phaseindex.size();
    ConsoleOutput::WriteStandard("Phasefield storage size", pfsize);
    size_t NxEBSD = individualxcoord.size();
    ConsoleOutput::WriteStandard("Number of (individual) x coordinates", NxEBSD);
    double maxX = xcoord.back();
    ConsoleOutput::WriteStandard("Biggest X Coord", maxX);
    double distX = maxX/double(NxEBSD);
    ConsoleOutput::WriteStandard("Rastering X", distX);
    size_t NyEBSD = individualycoord.size();
    ConsoleOutput::WriteStandard("Number of (individual) y coordinates", NyEBSD);
    double maxY = ycoord.back();
    ConsoleOutput::WriteStandard("Biggest Y Coord", maxY);
    double distY = maxY/double(NyEBSD);
    ConsoleOutput::WriteStandard("Rastering Y", distY);

//    for(size_t i = 0 ; i < phaseindex.size(); ++i) cout << phaseindex[i] << " ";
//    cout << endl;
//
//    for(size_t i = 0 ; i < xcoord.size(); ++i) cout << xcoord[i] << " ";
//    cout << endl;
//
    stringstream outPhaseIndeces;
    for(size_t i = 0 ; i < individualphaseindices.size(); ++i) outPhaseIndeces << individualphaseindices[i] << " ";
    ConsoleOutput::WriteSimple(outPhaseIndeces.str());
    ConsoleOutput::WriteBlankLine();

    const int Nx = Phase.Grid.Nx;
    const int Ny = Phase.Grid.Ny;

    for (size_t idx = 0; idx < maxphaseindex+1; idx++)
    {
        if(std::find(individualphaseindices.begin(), individualphaseindices.end(), idx) != individualphaseindices.end())
        {
            size_t locIndex = Phase.AddGrainInfo(0) + idx%Phase.Nphases; // eigene zv

            ConsoleOutput::WriteStandard("Add grain/phase", locIndex);

            int incx = 1; // Skip points in x-direction
            int incy = 1; // Skip points in y-direction

            //int itery = 1;
            int colstep = 0;
            for(int j = 0; j < Ny; ++j)
            {
                int iterx = 1;

                for(int i = 0; i < Nx; ++i)
                {
                    if (idx + 1 == phaseindex[-1 + iterx + colstep])
    //                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                    {
    //                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                        Phase.Fields(i, j, 0).clear();
                        Phase.Fields(i, j, 0).flag = 2;
                        Phase.Fields(i, j, 0).set_value(locIndex, 1.0);
                    }
                    iterx += incx;
                }
                if(incx == 1 and incy == 1)
                {
                    colstep += Nx;
                }
                else
                {
                // Odd columns
                    if(j%2 == 1) colstep += Nx;
                    else colstep += (Nx-1);
                }
                //itery += incy;
            }
        }
    }

    Phase.FinalizeInitialization(BC);
    ConsoleOutput::WriteStandard("EBSD", "done");
    ConsoleOutput::WriteLine();

    return 0;
}

int Initializations::TwoDimEBSDWithOrientations(std::string filename, std::vector<int> columns, std::string anglerepresentation,
        PhaseField& Phase, ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    // Euler angle 1 -> columns[3]
    // Euler angle 2 -> columns[4]
    // Euler angle 3 -> columns[5]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    const int Nx = Phase.Grid.Nx;
    const int Ny = Phase.Grid.Ny;
    const int Nz = Phase.Grid.Nz;

    if(columns.size() != 6)
    {
        columns.clear(); columns = {1, 2, 3, 4, 5, 6};
        ConsoleOutput::WriteWarning("Column list incomplete. Set to default {1, 2, 3, 4, 5, 6}.",
                thisclassname, "TwoDimEBSD");
    }

    if(Phase.Grid.Nz != 1)
    {
        ConsoleOutput::WriteExit("Only 2d supported", thisclassname, "TwoDimEBSD");
        OP_Exit(EXIT_FAILURE);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + filename + " could not be opened",
                thisclassname, "TwoDimEBSD()");
        OP_Exit(EXIT_FAILURE);
    };

    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("EBSD reader");
    ConsoleOutput::WriteStandard("Source", filename);
    ConsoleOutput::WriteStandard("Read x coordinates from column", columns[0]);
    ConsoleOutput::WriteStandard("Read y coordinates from column", columns[1]);
    ConsoleOutput::WriteStandard("Read phase indices from column", columns[2]);
    ConsoleOutput::WriteStandard("Read Euler angle 1 from column", columns[3]);
    ConsoleOutput::WriteStandard("Read Euler angle 2 from column", columns[4]);
    ConsoleOutput::WriteStandard("Read Euler angle 3 from column", columns[5]);

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    ConsoleOutput::WriteStandard("Number of lines", nol);

    // Read input file

    vector<size_t> phaseindex;
    vector<size_t> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;
    vector<EulerAngles> Eang;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        EulerAngles tempEang;
        tempEang.set_to_zero();
        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }

            double anglefactor = 180.0/Pi;
            if(!anglerepresentation.compare("degree") or !anglerepresentation.compare("deg")
                    or !anglerepresentation.compare("Degree") or !anglerepresentation.compare("Deg"))
            {
                anglefactor = 1.0;
            }

            if(spos == columns[3] and str.size() > 0) // Read Euler angle 1
            {
                tempEang.Q[0] = std::stod(str)*anglefactor;
            }

            if(spos == columns[4] and str.size() > 0) // Read Euler angle 2
            {
                tempEang.Q[1] = std::stod(str)*anglefactor;
            }

            if(spos == columns[5] and str.size() > 0) // Read Euler angle 3
            {
                tempEang.Q[2] = std::stod(str)*anglefactor;
            }

            if (str.size() > 0) spos++;
            if( !iline ) break;
        }

        tempEang.set_convention(ZXZ);

        ConsoleOutput::WriteWarning("Using ZXZ angle convention", thisclassname, "TwoDimEBSDWithOrientations()");

        tempEang.setTrigonometricFunctions();
        Eang.push_back(tempEang);
    }

    size_t individualgrainnum = individualphaseindices.size();
    ConsoleOutput::WriteStandard("Number of individual phases", individualgrainnum);
    size_t pfsize = phaseindex.size();
    ConsoleOutput::WriteStandard("Phasefield storage size", pfsize);
    size_t NxEBSD = individualxcoord.size();
    ConsoleOutput::WriteStandard("Number of (individual) x coordinates", NxEBSD);
    double maxX = 0;
    if(xcoord.size()) maxX = xcoord.back();
    ConsoleOutput::WriteStandard("Biggest X Coord", maxX);
    double distX = maxX/double(NxEBSD);
    ConsoleOutput::WriteStandard("Rastering X", distX);
    size_t NyEBSD = individualycoord.size();
    ConsoleOutput::WriteStandard("Number of (individual) y coordinates", NyEBSD);
    double maxY = 0;
    if(ycoord.size()) maxY = ycoord.back();
    ConsoleOutput::WriteStandard("Biggest Y Coord", maxY);
    double distY = maxY/double(NyEBSD);
    ConsoleOutput::WriteStandard("Rastering Y", distY);

    Storage3D <EulerAngles, 0> StorageEulerAngles;
    StorageEulerAngles.Allocate(Nx, Ny, 1, 1,1,0, 0);

    for (size_t idx = 0; idx < individualgrainnum; idx++)
    {
        size_t locIndex = Phase.AddGrainInfo(0) + idx%Phase.Nphases;

        ConsoleOutput::WriteStandard("Add grain/phase", locIndex);

        int incx = 1; // Skip points in x-direction
        int incy = 1; // Skip points in y-direction

        //int itery = 1;
        int colstep = 0;
        for(int j = 0; j < Ny; ++j)
        {
            int iterx = 1;

            for(int i = 0; i < Nx; ++i)
            {
                // SetInitialOrientations
                EulerAngles EangLocal = Eang[-1 + iterx + colstep];
                OR.Quaternions(i,j,0) = EangLocal.getQuaternion();
                StorageEulerAngles(i,j,0) = EangLocal;

                // Set phase fields
                if (idx + 1 == phaseindex[-1 + iterx + colstep])
//                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                {
//                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                    Phase.Fields(i, j, 0).clear();
                    Phase.Fields(i, j, 0).flag = 2;
                    Phase.Fields(i, j, 0).set_value(locIndex, 1.0);
                }
                iterx += incx;
            }
            if(incx == 1 and incy == 1)
            {
                colstep += Nx;
            }
            else
            {
            // Odd columns
                if(j%2 == 1) colstep += Nx;
                else colstep += (Nx-1);
            }
            //itery += incy;
        }
    }

    // Create output to VTK

    stringstream outbuffer;

    outbuffer << "# vtk DataFile Version 3.0\n";
    outbuffer << "InitialEulerAngles\n";
    outbuffer << "ASCII\n";
    outbuffer << "DATASET STRUCTURED_GRID\n";
    outbuffer << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << "\n";
    outbuffer << "POINTS " <<  Nx*Ny << " int\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbuffer << i << " " << j << " " << k << "\n";
    }
    outbuffer << " \n";
    outbuffer << "POINT_DATA " << Nx*Ny << " \n";

    outbuffer << "SCALARS Eang_" << 1 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[0] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 2 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[1] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 3 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[2] << " ";
    }

    string FileName = DefaultVTKDir + "InitializedEulerAngles.vtk";

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbuffer.rdbuf();
    vtk_file.close();

    // Write table of orientations

    outbuffer.str("");

    for(size_t i = 0 ; i < Eang.size(); ++i)
    {
        outbuffer << i << "  " << Eang[i].Q[0] << "  "
                               << Eang[i].Q[1] << "  "
                               << Eang[i].Q[2] << endl;
    }
    FileName = "EBSDorientations.dat";

    ofstream orientations_file(FileName.c_str());
    orientations_file << outbuffer.rdbuf();
    orientations_file.close();

    // End write Euler angles

    Phase.FinalizeInitialization(BC);
    ConsoleOutput::WriteStandard("EBSD", "done");
    ConsoleOutput::WriteLine();

    return 0;
}

size_t Initializations::SphereInGrain(PhaseField& Phase,
                                      const size_t ParrentPhaseFieldIndex,
                                      const size_t PhaseIndex,
                                      const double Radius,
                                      const double x0, const double y0, const double z0,
                                      const BoundaryConditions& BC,
                                      const bool Finalize)
{
    const double iWidth = (Phase.Grid.Resolution == Resolutions::Dual) ? 0.5*Phase.Grid.iWidth : Phase.Grid.iWidth;
    const size_t locIndex = Phase.AddGrainInfo(PhaseIndex);
    auto set_sphere = [&iWidth,&Radius,&locIndex,&Phase,&ParrentPhaseFieldIndex](long int i,long int j,long int k, double rad)
    {
        const double ParrentPhaseAmount = Phase.Fields(i, j, k).get_value(ParrentPhaseFieldIndex);
        if (rad < (Radius - iWidth*0.5))
        {
            Phase.Fields(i,j,k).set_value(ParrentPhaseFieldIndex, 0.0);
            Phase.Fields(i,j,k).set_value(locIndex, ParrentPhaseAmount);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (rad < (Radius + iWidth*0.5))
        {
            const double Profile = (0.5 - 0.5*sin(Pi*(rad- Radius)/iWidth))*ParrentPhaseAmount;
            Phase.Fields(i, j, k).set_value(ParrentPhaseFieldIndex, ParrentPhaseAmount-Profile);
            Phase.Fields(i, j, k).set_value(locIndex, Profile);
            Phase.Fields(i,j,k).flag = 2;
        }
        return false;
    };
    const double loop_radius = Radius + 1.1*iWidth;
    loop_sphere(Phase.Fields, set_sphere, x0, y0, z0, loop_radius, Phase.Grid, BC);

    if (Finalize) Phase.FinalizeInitialization(BC);
    return locIndex;
}

size_t  Initializations::FillGrainWithSpheres(PhaseField& Phase,
                                              size_t ParrentPhaseFieldIndex,
                                              size_t SpheresPhaseIndex,
                                              double MinR, double MaxR,
                                              BoundaryConditions& BC,
                                              size_t Nspheres,
                                              double MinDistance)
{
    if (MinDistance < 0.0) MinDistance = MinR;

    vector<iVector3> spheres;
    vector<double> Radii;

    int counter = 0;
    int globalCounter = 0;
    const int MaxNumberOfTries = 10.0*(Phase.Grid.TotalNumberOfCells());

    //default_random_engine generator;
    std::random_device generator;
    //std::mt19937_64 generator(rd());
    uniform_int_distribution<int> distributionX(0, (Phase.Grid.TotalNx-1));
    uniform_int_distribution<int> distributionY(0, (Phase.Grid.TotalNy-1));
    uniform_int_distribution<int> distributionZ(0, (Phase.Grid.TotalNz-1));
    uniform_real_distribution<double> distributionR(MinR, MaxR);

    while(counter < MaxNumberOfTries)
    {
        int i = distributionX(generator);
        int j = distributionY(generator);
        int k = distributionZ(generator);
        double Radius  = distributionR(generator);
        #ifdef MPI_PARALLEL
        OP_MPI_Bcast(&(i), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&(j), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&(k), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&(Radius), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
        #endif

        bool overlapping = false;
        for(size_t n = 0; n < spheres.size(); n++)
        {
            double distance = Tools::Distance<iVector3>({i,j,k}, spheres[n], Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC).abs();
            if(distance < Radius + Radii[n] + MinDistance) overlapping = true;
        }

        #ifdef MPI_PARALLEL
        double locParentPhaseFraction = 0;
        double ParentPhaseFraction = 0;
        if (i > Phase.Grid.OffsetX and i < Phase.Grid.OffsetX + Phase.Grid.Nx)
        {
            locParentPhaseFraction = Phase.Fields(i - Phase.Grid.OffsetX, j, k).get_value(ParrentPhaseFieldIndex);
        }
        OP_MPI_Allreduce(&(locParentPhaseFraction), &(ParentPhaseFraction), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #else
            double ParentPhaseFraction = Phase.Fields(i,j,k).get_value(ParrentPhaseFieldIndex);
        #endif

        if(overlapping or ParentPhaseFraction < 1.0)
        {
            counter++;
        }
        else
        {
            size_t sphereIndex =
            SphereInGrain(Phase, ParrentPhaseFieldIndex, SpheresPhaseIndex,
                          Radius, i, j, k, BC, false);
            spheres.push_back({i,j,k});
            Radii.push_back(Radius);
            stringstream message;
            message << "*********************************************************\n";
            message << "after " << counter << " tries sphere nr: "
                    << spheres.size() << " at (" << i << "," << j << "," << k <<")"
                    << " with R = " << Radius<< " and phase field index = "
                    << sphereIndex << " set." << "\n";
            //Phase.PrintPointStatistics(i,j,k); //NOTE: not mpi-parallel
            message << "*********************************************************\n";
            ConsoleOutput::WriteSimple(message.str());
            globalCounter += counter;
            counter = 0;

            if(Nspheres != 0 and spheres.size() == Nspheres)
            {
                break;
            }
        }
    }
    cout << "after " << globalCounter << " tries " << spheres.size()
         << " spheres could be initialized!\n";

    Phase.FinalizeInitialization(BC);

    return spheres.size();
}

std::vector<size_t> Initializations::FillRectangularWithSpheres(
                PhaseField& Phase,
                const BoundaryConditions& BC,
                const Settings& locSettings,
                std::function<size_t(long int, long int, long int)> PhaseIndex,
                double MeanRadius, double StdRadius,
                long int x_min, long int x_max,
                long int y_min, long int y_max,
                long int z_min, long int z_max,
                double MinDistance,
                size_t RandomSeed,
                bool Verbose)
{
    const int Nx = Phase.Grid.TotalNx;
    const int Ny = Phase.Grid.TotalNy;
    const int Nz = Phase.Grid.TotalNz;
    const double MinRadius = Phase.Grid.iWidth/2;
    const double MaxRadius = std::max(Nx/2, std::max(Ny/2, Nz/2));

    if (StdRadius < 0.0)
    {
        ConsoleOutput::WriteExit("Standard deviation needs to be positive!", thisclassname, "FillWithSpheres");
        OP_Exit(EXIT_FAILURE);
    }
    if (MeanRadius < MinRadius)
    {
        ConsoleOutput::WriteExit("The mean radius needs to be at least half the interface width!", thisclassname, "FillWithSpheres");
        OP_Exit(EXIT_FAILURE);
    }
    if (MeanRadius > MaxRadius)
    {
        ConsoleOutput::WriteExit("The mean radius needs to be smaller than the simulation box!", thisclassname, "FillWithSpheres");
        OP_Exit(EXIT_FAILURE);
    }

    std::vector<size_t> SpheresIdx;
    std::mt19937_64 generator(RandomSeed);
    std::normal_distribution<double> radius_gen(MeanRadius,StdRadius);
    std::uniform_real_distribution<double> dist_01(0.0,1.0);
    std::uniform_real_distribution<double> phi_gen(0,2.0*Pi);
    std::uniform_real_distribution<double> theta_gen(0,Pi);

    auto rand01 = std::bind(dist_01,generator);

    Storage3D<int,0> Shielded(Phase.Grid, 0);
    struct Planted_t {int x; int y; int z; double radius;};
    std::vector<Planted_t> Planted;
    {
        // Die solid radius
        double radius = radius_gen(generator);
        #ifdef MPI_PARALLEL
        OP_MPI_Bcast(&(radius), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
        #endif

        while(radius < MinRadius or radius > MaxRadius) radius = radius_gen(generator);

        if (Verbose) ConsoleOutput::Write("Plant initial Grain");
        assert(radius < MaxRadius);

        // Select random position inside boundaries
        Planted_t PInit;
        PInit.x = (Nx > 1 ) ? x_min + rand01()*(x_max-x_min) : 0;
        PInit.y = (Ny > 1 ) ? y_min + rand01()*(y_max-y_min) : 0;
        PInit.z = (Nz > 1 ) ? z_min + rand01()*(z_max-z_min) : 0;
        PInit.radius = radius;

        // Check if simulation box is 2D
        while (PInit.x >= Nx) PInit.x -= Nx;
        while (PInit.y >= Ny) PInit.y -= Ny;
        while (PInit.z >= Nz) PInit.z -= Nz;
        while (PInit.x <   0) PInit.x += Nx;
        while (PInit.y <   0) PInit.y += Ny;
        while (PInit.z <   0) PInit.z += Nz;

        // Initialize phase-field
        #ifdef MPI_PARALLEL
        OP_MPI_Bcast(&(PInit.x), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&(PInit.y), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
        OP_MPI_Bcast(&(PInit.z), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
        #endif
        Planted.push_back(PInit);

        const size_t phaseIdx =  PhaseIndex(PInit.x, PInit.y, PInit.z);
        const size_t idx = Initializations::Sphere(Phase, phaseIdx, radius, PInit.x, PInit.y, PInit.z, BC, false);
        SpheresIdx.push_back(idx);

        // Shield surrounding against places of other solid particles
        auto set_shield = [&Shielded](long int i,long int j,long int k, double rad){Shielded(i,j,k) = true; return false;};
        Initializations::loop_sphere(Shielded, set_shield, PInit.x, PInit.y, PInit.z, radius + MinDistance/2.0, Phase.Grid, BC);

        if (Verbose) ConsoleOutput::Write("Planted number ", Planted.size());

        size_t i = 0;
        while(i < Planted.size())
        {
            size_t N = 1000;//4.0/3.0*M_PI*(std::pow(2*radius,3) - std::pow(radius,3));
            for (size_t j = 0; j < N; j++)
            {
                radius = radius_gen(generator);
                #ifdef MPI_PARALLEL
                //MPI_Bcast(&(i),      1, MPI_INT, 0, MPI_COMM_WORLD);
                //MPI_Bcast(&(j),      1, MPI_INT, 0, MPI_COMM_WORLD);
                OP_MPI_Bcast(&(radius), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
                #endif
                if (radius < MinRadius) continue;

                // Place new solid next to the old one
                double phi   = phi_gen(generator);
                double theta = theta_gen(generator);
                #ifdef MPI_PARALLEL
                OP_MPI_Bcast(&(phi),   1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
                OP_MPI_Bcast(&(theta), 1, OP_MPI_DOUBLE, 0, OP_MPI_COMM_WORLD);
                #endif

                Planted_t PNew;
                PNew.x = (Nx > 1) ? Planted[i].x + (Planted[i].radius+MinDistance+radius+3)*cos(phi)   : 0;//*sin(theta);
                PNew.y = (Ny > 1) ? Planted[i].y + (Planted[i].radius+MinDistance+radius+3)*sin(phi)   : 0;//*sin(theta);
                PNew.z = (Nz > 1) ? Planted[i].z + (Planted[i].radius+MinDistance+radius+3)*cos(theta) : 0;//TODO z is too large
                PNew.radius = radius;
                // Check if new solid is in simulation box
                if (PNew.x >= Nx) continue;
                if (PNew.y >= Ny) continue;
                if (PNew.z >= Nz) continue;
                if (PNew.x <   0) continue;
                if (PNew.y <   0) continue;
                if (PNew.z <   0) continue;

                // Check if floor and ceiling are not touched
                if (x_min != 0  and PNew.x-PNew.radius < x_min) continue;
                if (x_max != Nx and PNew.x+PNew.radius > x_max) continue;
                if (y_min != 0  and PNew.y-PNew.radius < y_min) continue;
                if (y_max != Ny and PNew.y+PNew.radius > y_max) continue;
                if (z_min != 0  and PNew.z-PNew.radius < z_min) continue;
                if (z_max != Nz and PNew.z+PNew.radius > z_max) continue;

                // Check if surrounding shielded against planting
                auto check_shield = [&Shielded](long int ii,long int jj,long int kk, double rad){return (Shielded(ii,jj,kk) == true) ? true : false;};
                int  IsShielded = Initializations::loop_sphere_with_exit(Shielded, check_shield, PNew.x, PNew.y, PNew.z, radius + MinDistance/2.0, Phase.Grid, BC);
                #ifdef MPI_PARALLEL
                int tmpIsShielded = IsShielded;
                OP_MPI_Allreduce(&tmpIsShielded, &(IsShielded), 1, OP_MPI_INT, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                #endif
                if (IsShielded) continue;

                // Shield surrounding against places of other solid particles
                Initializations::loop_sphere(Shielded, set_shield, PNew.x, PNew.y, PNew.z, radius + MinDistance/2.0, Phase.Grid, BC);

                // Initialize spherical phase field
                Planted.push_back(PNew);
                int phaseIdx2 = PhaseIndex(PNew.x, PNew.y, PNew.z);
                #ifdef MPI_PARALLEL
                OP_MPI_Bcast(&(phaseIdx2), 1, OP_MPI_INT, 0, OP_MPI_COMM_WORLD);
                #endif
                const size_t idx = Initializations::Sphere(Phase, phaseIdx2, radius, PNew.x, PNew.y, PNew.z, BC, false);
                SpheresIdx.push_back(idx);

                if (Verbose) ConsoleOutput::Write("Planted number ", Planted.size());
            }
            ++i;
        }
    }

    uniform_int_distribution <int> Q1Distribution(0, 360);
    uniform_int_distribution <int> Q2Distribution(0, 180);
    uniform_int_distribution <int> Q3Distribution(0, 360);
    for(size_t n = 0; n < Phase.FieldsProperties.size(); n++)
    {
        double Q1 = Q1Distribution(generator) * Pi/180.0;
        double Q2 = Q2Distribution(generator) * Pi/180.0;
        double Q3 = Q3Distribution(generator) * Pi/180.0;
        EulerAngles locAngles({Q1, Q2, Q3}, XYZ);
        Phase.FieldsProperties[n].Orientation = locAngles.getQuaternion();
    }

#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        const char sep = ',';
        ConsoleOutput::Write("Spheres planted", Planted.size());
        std::ofstream file(locSettings.TextDir + "InitialParticlePositions.csv");
        file << "x" << sep << "y" << sep << "z" << sep << "radius\n";
        for (auto sphere = Planted.cbegin(); sphere != Planted.end(); sphere++)
        {
            file << sphere->x << sep;
            file << sphere->y << sep;
            file << sphere->z << sep;
            file << sphere->radius << "\n";
        }
    }

    // Finalise initialization
    Phase.FinalizeInitialization(BC);
    return SpheresIdx;
}

std::vector<size_t> Initializations::RectangularGreenBody(
        PhaseField& Phase,
        const BoundaryConditions& BC,const Settings& locSettings,
        std::function<size_t(long, long, long)> SolidPhaseIndex,
        size_t FluidPhaseIndex,
        long SizeX,
        long SizeY,
        long SizeZ,
        double MeanRadius,
        double StdRadius,
        size_t RandomSeed,
        double RelativeDensity,
        double accuracy,
        size_t MaxItrations)
{
    // Clamp input parameters 
    SizeX = std::clamp(SizeX,0l,long(Phase.Grid.TotalNx));
    SizeY = std::clamp(SizeY,0l,long(Phase.Grid.TotalNy));
    SizeZ = std::clamp(SizeZ,0l,long(Phase.Grid.TotalNz));
    MeanRadius = std::clamp(MeanRadius,0.,std::numeric_limits<double>::max());
    StdRadius = std::clamp(StdRadius,0.,std::numeric_limits<double>::max());
    RelativeDensity = std::clamp(RelativeDensity,0.,1.0);
    accuracy = std::clamp(accuracy,0.,1.0);

    // Calculate rectangular coordinates
    long x_min = std::max(int((Phase.Grid.TotalNx - SizeX)/2), 0);
    long x_max = std::min(int((Phase.Grid.TotalNx - SizeX)/2 + SizeX), Phase.Grid.TotalNx);
    long y_min = std::max(int((Phase.Grid.TotalNy - SizeY)/2), 0);
    long y_max = std::min(int((Phase.Grid.TotalNy - SizeY)/2 + SizeY), Phase.Grid.TotalNx);
    long z_min = std::max(int((Phase.Grid.TotalNz - SizeZ)/2), 0);
    long z_max = std::min(int((Phase.Grid.TotalNz - SizeZ)/2 + SizeZ), Phase.Grid.TotalNx);

    ConsoleOutput::Write("Try to initialize micro structure with relative density", RelativeDensity);
    std::vector<size_t> SpheresIdx;
    auto SolidPhaseFraction = [&Phase](long i, long j, long k){return Phase.SolidPhaseFraction(i,j,k);};
    auto residual = [&] (double& MinDist)
    {
        Phase.Clear();
        Phase.FieldsProperties.Reallocate(0);
        Single(Phase,FluidPhaseIndex,BC);
        SpheresIdx = FillRectangularWithSpheres(Phase, BC,locSettings,
                SolidPhaseIndex, MeanRadius, StdRadius,
                x_min,x_max,y_min,y_max,z_min,z_max, MinDist, RandomSeed, false);
        double RelRho = AnalysisSintering::CalculateDensity_LineIntersection(Phase,SolidPhaseFraction);
        ConsoleOutput::Write("Density of micro structure: "+std::to_string(RelRho)+" MinDist:"+std::to_string(MinDist));
        return RelRho - RelativeDensity;
    };

    double MinDist = 0.0;
    try
    {
        double& x0 = MinDist;
        double  x1 = MinDist+0.1;
        RootFindingAlgorithms::Secant(residual, x0, x1, accuracy, MaxItrations);
    }
    catch (std::runtime_error& ecep)
    {
        ConsoleOutput::WriteWarning(ecep.what(),"Initializations","FillWithSpheresWithDensity");
    }

    // Finalise initialization
    if(Phase.Grid.Resolution == Resolutions::Dual)
    {
        Phase.Refine();
    }
    Phase.Finalize(BC);
    return SpheresIdx;
}

void RemoveParrentGrain(PhaseField& Phase, int ParrentPhaseFieldIndex,
                        BoundaryConditions& BC, Settings& locSettings)
{
    /*
     * Removes phase number ParrentPhaseFieldIndex and increases phase number 0
     * by the removed amount.
     */
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double ParrentPhaseAmount = Phase.Fields(i, j, k).get_value(ParrentPhaseFieldIndex);

        if(ParrentPhaseAmount != 0.0)
        {
            Phase.Fields(i, j, k).add_value(0, ParrentPhaseAmount);
            Phase.Fields(i, j, k).set_value(ParrentPhaseFieldIndex, 0);
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
}

std::vector<size_t> Initializations::ThermalGrooving(PhaseField& Phase,
                                                     size_t alpha, size_t beta,
                                                     BoundaryConditions& BC)
{
    size_t Nx = Phase.Grid.Nx;
    size_t Nz = Phase.Grid.Nz;

    int x = Nx/2.0;
    int y = 0;
    int z = Nz+1;
    int dx = 3.0*Nx/4.0+1;
    int dy = 0;
    int dz = Nz/2.0;

    vector<size_t> locIndex = Fractional(Phase,beta,beta,dz,BC);
    locIndex.push_back(Rectangular(Phase,alpha,x,y,z,dx,dy,dz,BC));

    Phase.SetBoundaryConditions(BC);

    return locIndex;
}

void Initializations::VoronoiTessellation(PhaseField& Phase,
                                         BoundaryConditions& BC,
                                         const size_t Ngrains,
                                         const size_t GrainsPhase)
{
    ConsoleOutput::WriteLineInsert("Voronoi creator");
    ConsoleOutput::WriteStandard("Number of grains", Ngrains);
    ConsoleOutput::WriteStandard("Thermodynamic phase", GrainsPhase);
    ConsoleOutput::WriteSimple("Generating grain seeds... ");

    Phase.Clear();

    for(size_t n = 0; n < Ngrains; n++)
    {
        Phase.AddGrainInfo(GrainsPhase);
    }

    int TotalNx = Phase.Grid.TotalNx;
    int TotalNy = Phase.Grid.TotalNy;
    int TotalNz = Phase.Grid.TotalNz;

    int Nx = Phase.Grid.Nx;
    int Ny = Phase.Grid.Ny;
    int Nz = Phase.Grid.Nz;

    size_t SeedX = 253;
    size_t SeedY = 4958;
    size_t SeedZ = 54861;

    mt19937_64 xPosGenerator(SeedX);
    mt19937_64 yPosGenerator(SeedY);
    mt19937_64 zPosGenerator(SeedZ);

    uniform_int_distribution <int> xPosDistribution(0, TotalNx - 1);
    uniform_int_distribution <int> yPosDistribution(0, TotalNy - 1);
    uniform_int_distribution <int> zPosDistribution(0, TotalNz - 1);

    for (size_t n = 0; n < Phase.FieldsProperties.size(); ++n)
    {
        Phase.FieldsProperties[n].Rcm[0] = xPosDistribution(xPosGenerator);
        Phase.FieldsProperties[n].Rcm[1] = yPosDistribution(yPosGenerator);
        Phase.FieldsProperties[n].Rcm[2] = zPosDistribution(zPosGenerator);
    }
    ConsoleOutput::WriteSimple("Done!");

    ConsoleOutput::WriteSimple("Creating coordination shells...");

    size_t NpointsTotal = TotalNx*TotalNy*TotalNz;
    size_t NpointsRun = 0;

    vector<iVector3> CoordinationShells;

    int MAXsize = std::max(std::max(Phase.Grid.TotalNx, Phase.Grid.TotalNy), Phase.Grid.TotalNz);
    for(int x = 0; x < MAXsize; x++)
    for(int y = x; y < MAXsize; y++)
    for(int z = y; z < MAXsize; z++)
    {
        CoordinationShells.push_back(iVector3({x,y,z}));
    }
    std::sort(CoordinationShells.begin(),CoordinationShells.end());

    ConsoleOutput::WriteSimple("Done!");

    ConsoleOutput::StartProgressIndicator("Creating Voronoi grain structure");
    int position = 0;

    for(size_t s = 0; s < CoordinationShells.size(); s++)
    {
        vector <iVector3> locShell = Tools::findIndexPermutations(CoordinationShells[s]);

        for (size_t n = 0; n < Phase.FieldsProperties.size(); ++n)
        {
            int Rx = Phase.FieldsProperties[n].Rcm[0];
            int Ry = Phase.FieldsProperties[n].Rcm[1];
            int Rz = Phase.FieldsProperties[n].Rcm[2];

            for(size_t m = 0; m < locShell.size(); m++)
            {
                int x = Rx + locShell[m][0];
                int y = Ry + locShell[m][1];
                int z = Rz + locShell[m][2];
#ifdef MPI_PARALLEL
                if(BC.BC0X == BoundaryConditionTypes::Periodic or BC.MPIperiodicX) x = (x + TotalNx) % TotalNx;
                if(BC.BC0Y == BoundaryConditionTypes::Periodic or BC.MPIperiodicY) y = (y + TotalNy) % TotalNy;
                if(BC.BC0Z == BoundaryConditionTypes::Periodic or BC.MPIperiodicZ) z = (z + TotalNz) % TotalNz;

                x -= Phase.Grid.OffsetX;
                y -= Phase.Grid.OffsetY;
                z -= Phase.Grid.OffsetZ;
#else
                if(BC.BC0X == BoundaryConditionTypes::Periodic) x = (x + Nx) % Nx;
                if(BC.BC0Y == BoundaryConditionTypes::Periodic) y = (y + Ny) % Ny;
                if(BC.BC0Z == BoundaryConditionTypes::Periodic) z = (z + Nz) % Nz;
#endif
                if(x >= 0 and x < Nx and
                   y >= 0 and y < Ny and
                   z >= 0 and z < Nz and
                   Phase.Fields(x,y,z).size() == 0)
                {
                    Phase.Fields(x,y,z).set_value(n, 1.0);
                    NpointsRun ++;
                }
            }
#ifdef MPI_PARALLEL
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &NpointsRun, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif
            ConsoleOutput::AdvanceProgressIndicator(0,NpointsTotal,NpointsRun, position);

            if(NpointsRun == NpointsTotal) break;
        }
        if(NpointsRun == NpointsTotal) break;
    }
    ConsoleOutput::EndProgressIndicator();

    Phase.SetBoundaryConditions(BC);

    // Spread interfaces by overlapping the grains near the interface
    ConsoleOutput::WriteSimple("Creating grain boundaries...");
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    if(Phase.Fields(i,j,k).size() == 1)
    {
        size_t alpha_index = Phase.Fields(i,j,k).begin()->index;

        for(int di = -1; di <= 1; di++)
        for(int dj = -1; dj <= 1; dj++)
        for(int dk = -1; dk <= 1; dk++)
        if(i+di >= 0 and i+di < Nx and
           j+dj >= 0 and j+dj < Ny and
           k+dk >= 0 and k+dk < Nz and
           (di != 0 or dj != 0 or dk != 0))
        {
            bool grain_not_present = true;
            for(auto beta  = Phase.Fields(i+di,j+dj,k+dk).begin();
                     beta != Phase.Fields(i+di,j+dj,k+dk).end(); ++beta)
            if(beta->index == alpha_index)
            {
                grain_not_present = false;
                break;
            }
            else
            {
                break;
            }
            if(grain_not_present)
            {
                for(auto beta  = Phase.Fields(i+di,j+dj,k+dk).begin();
                         beta != Phase.Fields(i+di,j+dj,k+dk).end(); ++beta)
                {
                    Phase.Fields(i,j,k).add_value(beta->index, 1.0);
                }
                Phase.Fields(i+di,j+dj,k+dk).add_value(alpha_index, 1.0);
            }
        }
    }
    STORAGE_LOOP_END

    Phase.FinalizeInitialization(BC);
    ConsoleOutput::WriteSimple("Done!");

    ConsoleOutput::WriteSimple("Generating random orientations...");

    mt19937_64 OrientGenerator1(45);
    mt19937_64 OrientGenerator2(697);
    mt19937_64 OrientGenerator3(255);

    uniform_int_distribution <int> Q1Distribution(0, 360);
    uniform_int_distribution <int> Q2Distribution(0, 180);
    uniform_int_distribution <int> Q3Distribution(0, 360);

    for(size_t n = 0; n < Phase.FieldsProperties.size(); n++)
    {
        double Q1 = Phase.Grid.dNy*Phase.Grid.dNz*Q1Distribution(OrientGenerator1) * Pi/180.0;
        double Q2 = Phase.Grid.dNx*Phase.Grid.dNz*Q2Distribution(OrientGenerator2) * Pi/180.0;
        double Q3 = Phase.Grid.dNx*Phase.Grid.dNy*Q3Distribution(OrientGenerator3) * Pi/180.0;

        EulerAngles locAngles({Q1, Q2, Q3}, XYZ);
        Phase.FieldsProperties[n].Orientation = locAngles.getQuaternion();
    }

    ConsoleOutput::WriteSimple("Done!");
}

void Initializations::ReadCSV(PhaseField& Phase, BoundaryConditions& BC, std::filesystem::path FilePath, char Separator)
{
    ConsoleOutput::WriteLineInsert("Reading CSV file");
    ConsoleOutput::WriteStandard("Filename", FilePath.string());

    // Open file
    std::ifstream file(FilePath);
    if (!file.is_open())
    {
        ConsoleOutput::WriteExit("Could not open file", FilePath.string(), "ReadCSV");
        OP_Exit(EXIT_FAILURE);
    }

    // Read data
    std::set<size_t> PhaseIndices;
    std::vector<std::vector<long int>> data;
    {
        std::string line;
        while (std::getline(file, line))
        {
            std::vector<long int> row;
            std::istringstream lineStream(line);
            std::string cell;

            while (std::getline(lineStream, cell, Separator))
            {
                try
                {
                    long int value = std::stol(cell);
                    row.push_back(value);
                }
                catch (std::invalid_argument& e)
                {
                    ConsoleOutput::WriteExit("Invalid argument: "+cell, thisclassname, "ReadCSV");
                    OP_Exit(EXIT_FAILURE);
                }
                catch (std::out_of_range& e)
                {
                    ConsoleOutput::WriteExit("Out of range: "+cell, thisclassname, "ReadCSV");
                    OP_Exit(EXIT_FAILURE);
                }
            }
            if (row.size() != 5u)
            {
                ConsoleOutput::WriteExit("Wrong data format. Expected 4 columns i,j,k,phaseIdx,grainIdx", thisclassname, "ReadCSV");
                OP_Exit(EXIT_FAILURE);
            }
            long int phaseIdx = row[3];
            if (phaseIdx < 0)
            {
                ConsoleOutput::WriteExit("Phase index must be positive", thisclassname, "ReadCSV");
                OP_Exit(EXIT_FAILURE);
            }
            if (row[4] < 0)
            {
                ConsoleOutput::WriteExit("Grain index must be positive", thisclassname, "ReadCSV");
                OP_Exit(EXIT_FAILURE);
            }
            PhaseIndices.insert(phaseIdx);
            data.push_back(row);
        }
        file.close();
    }
    if (data.size() != (size_t) Phase.Grid.TotalNumberOfCells())
    {
        std::stringstream message;
        message << "Wrong number of rows, " << data.size() << ", expected: Nx*Ny*Nz = " << Phase.Grid.TotalNumberOfCells();
        ConsoleOutput::WriteExit(message.str(), thisclassname, "ReadCSV");
        OP_Exit(EXIT_FAILURE);
    }

    // Set phase-field values
    std::map<size_t, size_t> grainIdxMap;
    for (auto& line : data) 
    {
        long int i = line[0] - Phase.Grid.OffsetX;
        long int j = line[1] - Phase.Grid.OffsetY;
        long int k = line[2] - Phase.Grid.OffsetZ;
        if (i < 0 or i >= Phase.Grid.Nx) continue;
        if (j < 0 or j >= Phase.Grid.Ny) continue;
        if (k < 0 or k >= Phase.Grid.Nz) continue;
        size_t phaseIdx = line[3];
        size_t grainIdxFile = line[4];
        size_t grainIdx = 0;
        if (grainIdxMap.find(grainIdxFile) == grainIdxMap.end())
        {
            grainIdxMap[grainIdxFile] = Phase.AddGrainInfo(phaseIdx);
            grainIdx = grainIdxMap[grainIdxFile];
        }
        else
        {
            grainIdx = grainIdxMap[grainIdxFile];
        }
        Phase.Fields(i,j,k).set_value(grainIdx, 1.0);
    }
    Phase.FinalizeInitialization(BC);
    Phase.SetBoundaryConditions(BC);

    // Smooth phase-field
    Storage3D<NodePF,0> tmp(Phase.Fields);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if (Phase.Grid.dNx)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i+1,j,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i+1,j,k).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i-1,j,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i-1,j,k).front().index, 1);
        }
        if (Phase.Grid.dNy)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j+1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j+1,k).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j-1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j-1,k).front().index, 1);
        }
        if (Phase.Grid.dNz)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j,k+1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j,k+1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j,k-1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j,k-1).front().index, 1);
        }
        if (Phase.Grid.dNx && Phase.Grid.dNy)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i-1,j+1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i-1,j+1,k).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i-1,j-1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i-1,j-1,k).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i+1,j+1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i+1,j+1,k).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i+1,j-1,k).front().index) tmp(i,j,k).add_value(Phase.Fields(i+1,j-1,k).front().index, 1);
        }
        if (Phase.Grid.dNx && Phase.Grid.dNz)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i-1,j,k+1).front().index) tmp(i,j,k).add_value(Phase.Fields(i-1,j,k+1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i-1,j,k-1).front().index) tmp(i,j,k).add_value(Phase.Fields(i-1,j,k-1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i+1,j,k+1).front().index) tmp(i,j,k).add_value(Phase.Fields(i+1,j,k+1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i+1,j,k-1).front().index) tmp(i,j,k).add_value(Phase.Fields(i+1,j,k-1).front().index, 1);
        }
        if (Phase.Grid.dNy && Phase.Grid.dNz)
        {
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j-1,k+1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j-1,k+1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j-1,k-1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j-1,k-1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j+1,k+1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j+1,k+1).front().index, 1);
            if (Phase.Fields(i,j,k).front().index != Phase.Fields(i,j+1,k-1).front().index) tmp(i,j,k).add_value(Phase.Fields(i,j+1,k-1).front().index, 1);
        }
    }
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        Phase.Fields(i,j,k) = tmp(i,j,k);
    }
    STORAGE_LOOP_END
    Phase.FinalizeInitialization(BC);
    Phase.SetBoundaryConditions(BC);
}

}// namespace openphase

