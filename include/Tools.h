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

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Philipp Engels;
 *                         Muhammad Adil Ali; Hesham Salama
 *
 */

#ifndef TOOLS_H
#define TOOLS_H

#include "Includes.h"
#include "BoundaryConditions.h"

namespace openphase
{

class Crystallography;

class OP_EXPORTS Tools
{
 public:
    static std::vector<iVector3> findIndexPermutations(iVector3& locIndices);   ///< Finds all index permutations of locIndices

    /////////// Extract Rotation Matrix from Total Deformation ////////
    static void jacobiRotate(dMatrix3x3 &A, dMatrix3x3 &R, int p, int q);
    static void eigenDecomposition(const dMatrix3x3 &A, dMatrix3x3 &eigenVecs, dVector3 &eigenVals);
    static void rotationMatrixIrving(const dMatrix3x3 &A, dMatrix3x3 &R);
    static void GetAxisAngleFromRotationMatrix(const dMatrix3x3 RotMatrix, double& Angle, dVector3& Axis);
    static void getAxisAngle(const dMatrix3x3& TransformationMatrix, dVector3& Axis, double &Angle);

    static double getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB); ///< Calculates missorientation between two matrices without consideration of symmetries
    static double getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR); ///< Calculates missorientation between two matrices considering cubic symmetry
    static double getDisorientationCubic(const Quaternion OrientationA, const Quaternion OrientationB); ///< Calculates Disorientation between two quaternions considering cubic symmetry
    static EulerAngles RotationToEuler(const dMatrix3x3& Rot, const EulerConvention locConvention);
    static dVector3 MillerConversion(const std::vector<double>& hkil, dVector3 hkl, bool PlaneNormal = true); ///< Hexagonal Miller plane Normal indices to Miller-Bravais // flase for directions --

    static void SetRandomGrainOrientations(PhaseField& Phase, const int seed = 0);///< Sets random grains orientations

    static void ExtractRotation(dMatrix3x3 &Mat,  Quaternion& Quat, const size_t maxIter = 20);

    static double getMemoryUsageMB();                                           ///< Returns memory used by the simulation (useful for memory leaks detection)

    template<typename vector>
    static vector Position(const vector loc_pos, int OffsetX, int OffsetY, int OffsetZ)///< Gives global coordinates in MPI mode.
    {
        #ifdef MPI_PARALLEL
        return vector({loc_pos[0] + OffsetX, loc_pos[1] + OffsetY, loc_pos[2] + OffsetZ});
        #else
        return loc_pos;
        #endif
    }

    template<typename vector>
    static vector Distance( const vector A, const vector B, int TotalNx, int TotalNy, int TotalNz, const BoundaryConditions& BC)
    {
        vector dist = A-B;

        if (TotalNx > 1)
        {
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicX or BC.BC0X == BoundaryConditionTypes::Periodic)
#else
            if (BC.BC0X == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[0] - TotalNx;
                if (std::abs(val) < std::abs(dist[0])) dist[0] = val;
            }
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicX or BC.BCNX == BoundaryConditionTypes::Periodic)
#else
            if (BC.BCNX == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[0] + TotalNx;
                if (std::abs(val) < std::abs(dist[0])) dist[0] = val;
            }
        }

        if (TotalNy > 1)
        {
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicY or BC.BC0Y == BoundaryConditionTypes::Periodic)
#else
            if (BC.BC0Y == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[1] - TotalNy;
                if (std::abs(val) < std::abs(dist[1])) dist[1] = val;
            }
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicY or BC.BCNY == BoundaryConditionTypes::Periodic)
#else
            if (BC.BCNY == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[1] + TotalNy;
                if (std::abs(val) < std::abs(dist[1])) dist[1] = val;
            }
        }

        if (TotalNz > 1)
        {
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicZ or BC.BC0Z == BoundaryConditionTypes::Periodic)
#else
            if (BC.BC0Z == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[2] - TotalNz;
                if (std::abs(val) < std::abs(dist[2])) dist[2] = val;
            }
#ifdef MPI_PARALLEL
            if (BC.MPIperiodicZ or BC.BCNZ == BoundaryConditionTypes::Periodic)
#else
            if (BC.BCNZ == BoundaryConditionTypes::Periodic)
#endif
            {
                double val = dist[2] + TotalNz;
                if (std::abs(val) < std::abs(dist[2])) dist[2] = val;
            }
        }
        return dist;
    }

    template<typename vector>
    static inline void MapIntoSimulationDomain(vector& A, int TotalNx, int TotalNy, int TotalNz, const BoundaryConditions& BC)
    {
#ifdef MPI_PARALLEL
        if (BC.MPIperiodicX or BC.BCNX == BoundaryConditionTypes::Periodic) while (A[0] >= TotalNx) A[0] -= TotalNx;
        if (BC.MPIperiodicY or BC.BCNY == BoundaryConditionTypes::Periodic) while (A[1] >= TotalNy) A[1] -= TotalNy;
        if (BC.MPIperiodicZ or BC.BCNZ == BoundaryConditionTypes::Periodic) while (A[2] >= TotalNz) A[2] -= TotalNz;
        if (BC.MPIperiodicX or BC.BC0X == BoundaryConditionTypes::Periodic) while (A[0] <  0) A[0] += TotalNx;
        if (BC.MPIperiodicY or BC.BC0Y == BoundaryConditionTypes::Periodic) while (A[1] <  0) A[1] += TotalNy;
        if (BC.MPIperiodicZ or BC.BC0Z == BoundaryConditionTypes::Periodic) while (A[2] <  0) A[2] += TotalNz;
#else
        if (BC.BCNX == BoundaryConditionTypes::Periodic) while (A[0] >= TotalNx) A[0] -= TotalNx;
        if (BC.BCNY == BoundaryConditionTypes::Periodic) while (A[1] >= TotalNy) A[1] -= TotalNy;
        if (BC.BCNZ == BoundaryConditionTypes::Periodic) while (A[2] >= TotalNz) A[2] -= TotalNz;
        if (BC.BC0X == BoundaryConditionTypes::Periodic) while (A[0] <  0) A[0] += TotalNx;
        if (BC.BC0Y == BoundaryConditionTypes::Periodic) while (A[1] <  0) A[1] += TotalNy;
        if (BC.BC0Z == BoundaryConditionTypes::Periodic) while (A[0] <  0) A[2] += TotalNz;
#endif
    };

    template <typename T>
    static void Smooth(Storage3D<T, 0> &Field, size_t SmoothIterations = 1)
    {
        int Nx = Field.sizeX();
        int Ny = Field.sizeY();
        int Nz = Field.sizeZ();
        int dNx = Field.dNx();
        int dNy = Field.dNy();
        int dNz = Field.dNz();
        double Stencil[3][3][3] = {{{1.0/64.0,   1.0/32.0, 1.0/64.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/64.0,   1.0/32.0, 1.0/64.0}},

                                  {{1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/16.0,   1.0/8.0, 1.0/16.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0}},

                                  {{1.0/64.0,   1.0/32.0, 1.0/64.0},
                                   {1.0/32.0,   1.0/16.0, 1.0/32.0},
                                   {1.0/64.0,   1.0/32.0, 1.0/64.0}}};

        std::function<bool(int,int,int)> condition = [Nx,Ny,Nz](int u,int v,int w)
        {
            return u < Nx and v < Ny and w < Nz and u >= 0 and v >= 0 and w >= 0;
        };
        Storage3D<T, 0> FieldTmp;
        FieldTmp.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 1);
        for(size_t start = 0; start < SmoothIterations; start++)
        {
            if(start)
            {
                FieldTmp.Reallocate(Nx, Ny, Nz);
            }
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Field,Field.Bcells(),)
            {
                        double tmpStencil[3][3][3] = {{{0.0,   0.0,  0.0},
                                                     {0.0,   0.0,  0.0},
                                                     {0.0,   0.0,  0.0}},

                                                     {{0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0}},

                                                     {{0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0},
                                                      {0.0,   0.0,  0.0}}};
                        double WS = 0;
                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if (condition(i+a,j+b,k+c))
                            {
                                WS += Stencil[a+1][b+1][c+1];
                            }
                        }
                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if (condition(i+a,j+b,k+c))
                            {
                                tmpStencil[a+1][b+1][c+1] = Stencil[a+1][b+1][c+1]/WS;
                            }
                            else
                            {
                                tmpStencil[a+1][b+1][c+1] = Stencil[a+1][b+1][c+1];
                            }
                        }

                        for(int a = -dNx; a <= +dNx; a++)
                        for(int b = -dNy; b <= +dNy; b++)
                        for(int c = -dNz; c <= +dNz; c++)
                        {
                            if(condition(i+a,j+b,k+c))
                            {
                                FieldTmp(i,j,k) += Field(i+a,j+b,k+c) * tmpStencil[a+1][b+1][c+1];
                            }
                        }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
            {
                Field(i,j,k) = FieldTmp(i,j,k);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
        }
    }
/*
    template <typename T>
    static T Average(Storage3D<T, 0> &Field)
    {
        int Nx = Field.sizeX();
        int Ny = Field.sizeY();
        int Nz = Field.sizeZ();

        size_t Nthreads = 1;

        #ifdef _OPENMP
        Nthreads = omp_get_max_threads();
        #endif

        std::vector<T> avgThreads;
        avgThreads.resize(Nthreads);
        T AverageField;

        double Norm = 1.0/double(Nx*Ny*Nz);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
        {
            size_t thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
            avgThreads[thread] += Field(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        for(size_t thread = 0; thread < Nthreads; thread++)
            AverageField += avgThreads[thread];

        #ifdef MPI_PARALLEL
        for(int m = 0; m < 6; m++)
        {
            double tmp = AverageField[m];
            MPI_Reduce(&tmp,&AverageField[m],1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
        }
        #endif
        return AverageField*Norm;
    }*/
 protected:
 private:
};
}
#endif
