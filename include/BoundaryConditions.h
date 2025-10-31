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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung;
 *                         Reza Namdar; Marvin Tegeler
 *
 */

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "Includes.h"

namespace openphase
{
class Settings;
/// Types of boundary conditions:

enum class BoundaryConditionTypes
{
    Periodic,                                                                   ///< Periodic boundary condition
    NoFlux,                                                                     ///< Adiabatic boundary condition, results in zero first order derivative at the boundary
    Free,                                                                       ///< Boundary condition resembling free surface, results in continuous gradient across the boundary
    Fixed,                                                                      ///< Fixed (Dirichlet or Neumann) boundary condition
    Mirror,                                                                     ///< Similar to NoFlux but places the domain boundary at the edge grid point
    MPIcomm,                                                                    ///< Adjacent boundaries communication between neighboring blocks in MPI parallel mode
};

class OP_EXPORTS BoundaryConditions : public OPObject
{
 public:

    BoundaryConditions(){};
    BoundaryConditions(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    template< class T, size_t Num >
    void Set(Storage3D<T, Num>& loc3Dstorage, int dir) const;                   ///< Set boundary conditions for standard values storage along dir boundaries

    template< class T, size_t Num >
    void SetX(Storage3D<T, Num>& loc3Dstorage) const;                           /// Set boundary conditions for standard values storage along X boundaries
    template< class T, size_t Num >
    void SetY(Storage3D<T, Num>& loc3Dstorage) const;                           /// Set boundary conditions for standard values storage along Y boundaries
    template< class T, size_t Num >
    void SetZ(Storage3D<T, Num>& loc3Dstorage) const;                           /// Set boundary conditions for standard values storage along Z boundaries

    template< class T>
    void SetXVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along X boundaries
    template< class T>
    void SetYVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Y boundaries
    template< class T>
    void SetZVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Z boundaries

/*    template< class T>
    void SetXVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// Set boundary conditions for vector values storage along X boundaries
    template< class T>
    void SetYVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// Set boundary conditions for vector values storage along Y boundaries
    template< class T>
    void SetZVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// Set boundary conditions for vector values storage along Z boundaries

    template< class T>
    void SetLBXVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// Set boundary conditions for LB vector values storage along X boundaries
    template< class T>
    void SetLBYVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// Set boundary conditions for LB vector values storage along Y boundaries
    template< class T>
    void SetLBZVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// Set boundary conditions for LB vector values storage along Z boundaries
*/
    template< class T>
    void SetXVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along X boundaries
    template< class T>
    void SetYVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Y boundaries
    template< class T>
    void SetZVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Z boundaries

    template< class T>
    void SetXVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along X boundaries
    template< class T>
    void SetYVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Y boundaries
    template< class T>
    void SetZVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Z boundaries

    template< class T>
    void SetXVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along X boundaries
    template< class T>
    void SetYVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Y boundaries
    template< class T>
    void SetZVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// Set boundary conditions for vector values storage along Z boundaries

    template< class T, size_t Num >
    void SetXFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// Set boundary conditions for storage of flags along X boundaries
    template< class T, size_t Num >
    void SetYFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// Set boundary conditions for storage of flags along Y boundaries
    template< class T, size_t Num >
    void SetZFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// Set boundary conditions for storage of flags along Z boundaries
    void Initialize(Settings& Settings, std::string ObjectNameSuffix = "") override;/// Initializes the class's variables
    void ReadInput(std::string InputFileName) override;                         /// Read boundary conditions
    void ReadInput(std::stringstream& inp) override;
	void ReadJSON(const std::string InputFileName);
    long int Index(const long int x, const long int y, const long int z,
                   const long int Nx, const long int Ny, const long int Nz,
                   const long int BcellsX,
                   const long int BcellsY,
                   const long int BcellsZ) const;                               /// Index of the (x, y, z) position using the boundary conditions.

    BoundaryConditionTypes BC0X;
    BoundaryConditionTypes BCNX;
    BoundaryConditionTypes BC0Y;
    BoundaryConditionTypes BCNY;
    BoundaryConditionTypes BC0Z;
    BoundaryConditionTypes BCNZ;

#ifdef MPI_PARALLEL
    bool Setup_MPI();                                                           /// Setting up 1D MPI domain decomposition along a single X dimension
    bool Setup_MPIX();                                                          /// Setting up 3D MPI domain decomposition along X dimension
    bool Setup_MPIY();                                                          /// Setting up 3D MPI domain decomposition along Y dimension
    bool Setup_MPIZ();                                                          /// Setting up 3D MPI domain decomposition along Z dimension

    bool MPIperiodicX;
    bool MPIperiodicY;
    bool MPIperiodicZ;

    template<typename A> void Communicate(A& storage) const;
    template<typename A> void CommunicateX(A& storage) const;
    template<typename A> void CommunicateY(A& storage) const;
    template<typename A> void CommunicateZ(A& storage) const;
#endif

    BoundaryConditionTypes TranslateBoundaryConditions(std::string Key);        ///< Translates input string into the valid boundary condition designation
 protected:
 private:

};

inline long int BoundaryConditions::Index(const long int x, const long int y, const long int z,
                                          const long int Nx, const long int Ny, const long int Nz,
                                          const long int BcellsX, const long int BcellsY, const long int BcellsZ) const
{
    long int xx = x;
    long int yy = y;
    long int zz = z;
    if(x < 0)
    {
        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                xx = (-x - 1);
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                xx = (-x);
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                xx = (x + Nx)%Nx;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                xx = 0;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                xx = std::max(xx,-BcellsX);
                break;
            }
        }
    }
    if(x > Nx - 1)
    {
        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                xx = 2*Nx - 1 - x;
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                xx = 2*Nx - 2 - x;
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                xx = x%Nx;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                xx = Nx - 1;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                xx = std::min(xx,Nx+BcellsX-1);
                break;
            }
        }
    }

    if(y < 0)
    {
        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                yy = (-y - 1);
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                yy = (-y);
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                yy = (y + Ny)%Ny;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                yy = 0;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                yy = std::max(yy,-BcellsY);
                break;
            }
        }
    }
    if(y > Ny - 1)
    {
        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                yy = 2*Ny - 1 - y;
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                yy = 2*Ny - 2 - y;
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                yy = y%Ny;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                yy = Ny - 1;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                yy = std::min(yy,Ny+BcellsY-1);
                break;
            }
        }
    }

    if(z < 0)
    {
        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                zz = (-z - 1);
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                zz = (-z);
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                zz = (z + Nz)%Nz;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                zz = 0;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                zz = std::max(zz,-BcellsZ);
                break;
            }
        }
    }
    if(z > Nz - 1)
    {
        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
                zz = 2*Nz - 1 - z;
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
                zz = 2*Nz - 2 - z;
                break;
            }
            case BoundaryConditionTypes::Periodic:
            {
                zz = z%Nz;
                break;
            }
            case BoundaryConditionTypes::Fixed:
            {
                zz = Nz - 1;
                break;
            }
            case BoundaryConditionTypes::Free: [[fallthrough]];
            default:
            {
                zz = std::min(zz,Nz+BcellsZ-1);
                break;
            }
        }
    }
    return ((Ny + 2*BcellsY)*(xx + BcellsX) + yy + BcellsY)*(Nz + 2*BcellsZ) + zz + BcellsZ;
}

// Flags
template< class T, size_t Num >
void BoundaryConditions::SetXFlags(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsX())
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field( 0, j, k) = std::max(Field( 0, j, k), Field(Field.sizeX(), j, k));
                Field(Field.sizeX() - 1, j, k) = std::max(Field(Field.sizeX() - 1, j, k), Field(-1, j, k));
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field( i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((- i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((- i - 1)%Field.sizeZ(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((- i)%Field.sizeZ(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - 1 - i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - 1 - i, j, k) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}
template< class T, size_t Num >
void BoundaryConditions::SetYFlags(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, 0, k) = std::max(Field(i, 0, k), Field(i, Field.sizeY(), k));
                Field(i, Field.sizeY() - 1, k) = std::max(Field(i, Field.sizeY() - 1, k), Field(i, -1, k));
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j -1, k) = Field(i, (- j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (- j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (- j)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - 1 - j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - 1 - j, k) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}
template< class T, size_t Num >
void BoundaryConditions::SetZFlags(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            {
                Field(i, j, 0) = std::max(Field(i, j, 0), Field(i, j, Field.sizeZ()));
                Field(i, j, Field.sizeZ() - 1) = std::max(Field(i, j, Field.sizeZ() - 1), Field(i, j, -1));
            }
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (- k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (- k - 1)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (- k)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - 1 - k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - 1 - k) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Free:  [[fallthrough]];
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

// Standard
template< class T, size_t Num >
void BoundaryConditions::SetX(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsX())
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((- i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((-i)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(0, j, k)
                        + (Field(1%Field.sizeX(), j, k) - Field(0, j, k)) * i;
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    const long int i0 = Field.sizeX() - 1;
                    Field(i0 - i, j, k) = Field(i0, j, k)
                        +(Field(i0, j, k) - Field((i0-1)%Field.sizeX(), j, k))*(-i);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}
template< class T, size_t Num >
void BoundaryConditions::SetY(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (-j)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, 0, k)
                              + (Field(i, 1%Field.sizeY(), k) - Field(i, 0, k)) * j;
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    const long int j0 = Field.sizeY() - 1;
                    Field(i, j0 - j, k) = Field(i, j0, k)
                        +(Field(i, j0, k) - Field(i, (j0-1)%Field.sizeY(), k))*(-j);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}
template< class T, size_t Num >
void BoundaryConditions::SetZ(Storage3D<T, Num> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, 0)
                        + (Field(i, j, 1%Field.sizeZ()) - Field(i, j, 0)) * k;
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    const long int k0 = Field.sizeZ() - 1;
                    Field(i, j, k0 - k) = Field(i, j, k0)
                        +(Field(i, j, k0) - Field(i, j, (k0-1)%Field.sizeZ()))*(-k);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

template< class T, size_t Num >
void BoundaryConditions::Set(Storage3D<T, Num>& loc3Dstorage, int dir) const
{
    switch(dir)
    {
        case 0:
        {
            SetX(loc3Dstorage);
            break;
        }
        case 1:
        {
            SetY(loc3Dstorage);
            break;
        }
        case 2:
        {
            SetZ(loc3Dstorage);
            break;
        }
        default:
        {
            ConsoleOutput::WriteExit("Wrong direction: \"dir\" has to be 0, 1 or 2!","BoundaryConditions","Set");
            OP_Exit(EXIT_FAILURE);
        }
    }
}

// Vectors
template< class T>
void BoundaryConditions::SetXVector(Storage3D<T, 0> &Field) const
{
    if(Field.BcellsX())
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((-i)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}

template< class T>
void BoundaryConditions::SetXVector(Storage3D<T, 1> &Field) const
{
    if(Field.BcellsX())
    {

        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n = 0; n < Field(0, j, k).size(0); n++)
            {
                Field(i, j, k)({n}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n});
                Field(Field.sizeX() - i - 1, j, k)({n}) = Field((-i - 1)%Field.sizeX(), j, k)({n});
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(0, j, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field((-i - 1)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(0, j, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field((-i)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }

            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(0, j, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field((-i - 1)%Field.sizeX(), j, k)({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(Field.sizeX() - 1, j, k).size(0); n++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(Field.sizeX() - 1, j, k).size(0); n++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n}) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(Field.sizeX() - 1, j, k).size(0); n++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}
template< class T>
void BoundaryConditions::SetXVector(Storage3D<T, 2> &Field) const
{
    if(Field.BcellsX())
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
            {
                Field(i, j, k)({n1, n2}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2});
                Field(Field.sizeX() - i - 1, j, k)({n1, n2}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2});
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field((-i)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2}) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}
template< class T>
void BoundaryConditions::SetXVector(Storage3D<T, 3> &Field) const
{
    if(Field.BcellsX())
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < 0; i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
            for(size_t n3 = 0; n3 < Field(0, j, k).size(2); n3++)
            {
                Field(i, j, k)({n1, n2, n3}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2, n3});
                Field(Field.sizeX() - i - 1, j, k)({n1, n2, n3}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2, n3});
            }
            return;
        }

        switch (BC0X)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(0, j, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(0, j, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field((-i)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(0, j, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(Field.sizeX() - 1, j, k).size(2); n3++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2, n3}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(Field.sizeX() - 1, j, k).size(2); n3++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2, n3}) = Field((Field.sizeX() + i - 1)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < 0; i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(Field.sizeX() - 1, j, k).size(2); n3++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2, n3}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateX(Field);
        }
        else
        {
            Communicate(Field);
        }
#endif
    }
}
template< class T>
void BoundaryConditions::SetYVector(Storage3D<T, 0> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() -j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (-j)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}

template< class T>
void BoundaryConditions::SetYVector(Storage3D<T, 1> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n = 0; n < Field(i, 0, k).size(0); n++)
            {
                Field(i, j, k)({n}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n});
                Field(i, Field.sizeY() - j - 1, k)({n}) = Field(i, (-j - 1)%Field.sizeY(), k)({n});
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, 0, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, (-j - 1)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, 0, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, (-j)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, 0, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, (-j - 1)%Field.sizeY(), k)({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, Field.sizeY() - 1, k).size(0); n++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, Field.sizeY() - 1, k).size(0); n++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n}) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n = 0; n < Field(i, Field.sizeY() - 1, k).size(0); n++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}
template< class T>
void BoundaryConditions::SetYVector(Storage3D<T, 2> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
            {
                Field(i, j, k)({n1, n2}) = Field(i, Field.sizeY() + j, k)({n1, n2});
                Field(i, Field.sizeY() - j - 1, k)({n1, n2}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2});
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, (-j)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2}) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.Bcellsz(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}
template< class T>
void BoundaryConditions::SetYVector(Storage3D<T, 3> &Field) const
{
    if(Field.BcellsY())
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < 0; j++)
            for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
            for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
            for(size_t n3 = 0; n3 < Field(i, 0, k).size(2); n3++)
            {
                Field(i, j, k)({n1, n2, n3}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2, n3});
                Field(i, Field.sizeY() - j - 1, k)({n1, n2, n3}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2, n3});
            }
            return;
        }

        switch (BC0Y)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, 0, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, 0, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, (-j)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, 0, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, Field.sizeY() - 1, k).size(2); n3++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2, n3}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, Field.sizeY() - 1, k).size(2); n3++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2, n3}) = Field(i, (Field.sizeY() + j - 1)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < 0; j++)
                for(long int k = -Field.BcellsZ(); k < Field.sizeZ() + Field.BcellsZ(); k++)
                for(size_t n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, Field.sizeY() - 1, k).size(2); n3++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2, n3}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateY(Field);
        }
#endif
    }
}
template< class T >
void BoundaryConditions::SetZVector(Storage3D<T, 0> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

template< class T >
void BoundaryConditions::SetZVector(Storage3D<T, 1> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            for(size_t n = 0; n < Field(i, j, 0).size(0); n++)
            {
                Field(i, j, k)({n}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n});
                Field(i, j, Field.sizeZ() - k - 1)({n}) = Field(i, j, (-k - 1)%Field.sizeZ())({n});
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, 0).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, j, (-k - 1)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, 0).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, j, (-k)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, 0).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, j, (-k - 1)%Field.sizeZ())({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, Field.sizeZ() - 1).size(0); n++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, Field.sizeZ() - 1).size(0); n++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n}) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n = 0; n < Field(i, j, 0).size(0); n++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

template< class T >
void BoundaryConditions::SetZVector(Storage3D<T, 2> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
            {
                Field(i, j, k)({n1, n2}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2});
                Field(i, j, Field.sizeZ() - k - 1)({n1, n2}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2});
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, j, (-k)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2}) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

template< class T >
void BoundaryConditions::SetZVector(Storage3D<T, 3> &Field) const
{
    if(Field.BcellsZ())
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
           BCNZ == BoundaryConditionTypes::Periodic)
        {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
            for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
            for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
            for(long int k = -Field.BcellsZ(); k < 0; k++)
            for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
            for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
            for(size_t n3 = 0; n3 < Field(i, j, 0).size(2); n3++)
            {
                Field(i, j, k)({n1, n2, n3}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2, n3});
                Field(i, j, Field.sizeZ() - k - 1)({n1, n2, n3}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2, n3});
            }
            return;
        }

        switch (BC0Z)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, 0).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, 0).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, j, (-k)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, 0).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case BoundaryConditionTypes::NoFlux:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, Field.sizeZ() - 1).size(2); n3++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2, n3}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Mirror:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, Field.sizeZ() - 1).size(2); n3++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2, n3}) = Field(i, j, (Field.sizeZ() + k - 1)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case BoundaryConditionTypes::Free:
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(OMP_SCHEDULING_TYPE,OMP_CHUNKSIZE)
#endif
                for(long int i = -Field.BcellsX(); i < Field.sizeX() + Field.BcellsX(); i++)
                for(long int j = -Field.BcellsY(); j < Field.sizeY() + Field.BcellsY(); j++)
                for(long int k = -Field.BcellsZ(); k < 0; k++)
                for(size_t n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(size_t n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                for(size_t n3 = 0; n3 < Field(i, j, Field.sizeZ() - 1).size(2); n3++)
                {
                    Field(i, j, Field.sizeZ() - k -1)({n1, n2, n3}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2, n3});
                }
                break;
            }
            case BoundaryConditionTypes::Fixed: [[fallthrough]];
            default:
            {
                break;
            }
        }
#ifdef MPI_PARALLEL
        if(MPI_3D_DECOMPOSITION)
        {
            CommunicateZ(Field);
        }
#endif
    }
}

#ifdef MPI_PARALLEL
template<typename A>
inline void BoundaryConditions::Communicate([[maybe_unused]] A& storage) const
{
    std::vector<double> sleft;
    std::vector<double> rleft;
    std::vector<double> sright;
    std::vector<double> rright;
    void* request_sleft = create_request();
    void* request_sleft_size = create_request();
    void* request_sright = create_request();
    void* request_sright_size = create_request();
    void* request_rleft = create_request();
    void* request_rleft_size = create_request();
    void* request_rright = create_request();
    void* request_rright_size = create_request();

    int sleft_size;
    int sright_size;
    const int RightProcess = (((MPI_RANK+1)%MPI_SIZE)+MPI_SIZE)%MPI_SIZE;
    const int LeftProcess  = (((MPI_RANK-1)%MPI_SIZE)+MPI_SIZE)%MPI_SIZE;

    const int LeftSizeTag  = 0; // used to identify data stream
    const int LeftDataTag  = 2; // used to identify data stream
    const int RightSizeTag = 4; // used to identify data stream
    const int RightDataTag = 8; // used to identify data stream

    if(MPI_RANK > 0 or MPIperiodicX)
    {
        std::vector<long int> window(6);
        window[0] = 0;
        window[1] = storage.BcellsX();
        window[2] = -storage.BcellsY();
        window[3] = storage.sizeY()+storage.BcellsY();
        window[4] = -storage.BcellsZ();
        window[5] = storage.sizeZ()+storage.BcellsZ();
        sleft = storage.pack(window);
        sleft_size = sleft.size();
        OP_MPI_Isend(&sleft_size , 1 , OP_MPI_INT , LeftProcess , RightSizeTag, OP_MPI_COMM_WORLD ,  request_sleft_size);
        OP_MPI_Isend(sleft.data() , sleft_size , OP_MPI_DOUBLE , LeftProcess , RightDataTag , OP_MPI_COMM_WORLD,  request_sleft);
    }
    if(MPI_RANK < MPI_SIZE-1 or MPIperiodicX)
    {
        std::vector<long int> window(6);
        window[0] = storage.sizeX()-storage.BcellsX();
        window[1] = storage.sizeX();
        window[2] = -storage.BcellsY();
        window[3] = storage.sizeY()+storage.BcellsY();
        window[4] = -storage.BcellsZ();
        window[5] = storage.sizeZ()+storage.BcellsZ();
        sright = storage.pack(window);
        sright_size = sright.size();
        OP_MPI_Isend(&sright_size , 1 , OP_MPI_INT , RightProcess , LeftSizeTag, OP_MPI_COMM_WORLD , request_sright_size);
        OP_MPI_Isend(sright.data() , sright_size , OP_MPI_DOUBLE , RightProcess , LeftDataTag, OP_MPI_COMM_WORLD , request_sright);
    }
    if(MPI_RANK > 0 or MPIperiodicX)
    {
        long int size = 0;
        OP_MPI_Irecv(&size, 1, OP_MPI_INT, LeftProcess, LeftSizeTag, OP_MPI_COMM_WORLD, request_rleft_size);
        OP_MPI_Wait(request_rleft_size, OP_MPI_STATUS_IGNORE);
        rleft.resize(size);
        OP_MPI_Irecv(rleft.data(), size, OP_MPI_DOUBLE, LeftProcess, LeftDataTag, OP_MPI_COMM_WORLD, request_rleft);
        OP_MPI_Wait(request_rleft, OP_MPI_STATUS_IGNORE);
        std::vector<long int> window(6);
        window[0] = -storage.BcellsX();
        window[1] = 0;
        window[2] = -storage.BcellsY();
        window[3] = storage.sizeY()+storage.BcellsY();
        window[4] = -storage.BcellsZ();
        window[5] = storage.sizeZ()+storage.BcellsZ();
        storage.unpack(rleft,window);
    }
    if (MPI_RANK < MPI_SIZE-1 or MPIperiodicX)
    {
        long int size = 0;
        OP_MPI_Irecv(&size, 1, OP_MPI_INT, RightProcess, RightSizeTag, OP_MPI_COMM_WORLD,  request_rright_size);
        OP_MPI_Wait(request_rright_size, OP_MPI_STATUS_IGNORE);
        rright.resize(size);
        OP_MPI_Irecv(rright.data(), size, OP_MPI_DOUBLE, RightProcess, RightDataTag, OP_MPI_COMM_WORLD, request_rright);
        OP_MPI_Wait(request_rright, OP_MPI_STATUS_IGNORE);
        std::vector<long int> window(6);
        window[0] = storage.sizeX();
        window[1] = storage.sizeX()+storage.BcellsX();
        window[2] = -storage.BcellsY();
        window[3] = storage.sizeY()+storage.BcellsY();
        window[4] = -storage.BcellsZ();
        window[5] = storage.sizeZ()+storage.BcellsZ();
        storage.unpack(rright,window);
    }
    if (MPI_RANK > 0 or MPIperiodicX)
    {
        OP_MPI_Wait(request_sleft, OP_MPI_STATUS_IGNORE);
        OP_MPI_Wait(request_sleft_size, OP_MPI_STATUS_IGNORE);
    }
    if (MPI_RANK < MPI_SIZE-1 or MPIperiodicX)
    {
        OP_MPI_Wait(request_sright, OP_MPI_STATUS_IGNORE);
        OP_MPI_Wait(request_sright_size, OP_MPI_STATUS_IGNORE);
    }
    free_request(request_sleft);
    free_request(request_sleft_size);
    free_request(request_sright);
    free_request(request_sright_size);
    free_request(request_rleft);
    free_request(request_rleft_size);
    free_request(request_rright);
    free_request(request_rright_size);
}

template<typename A>
inline void BoundaryConditions::CommunicateX([[maybe_unused]] A& storage) const
{
    if(MPI_CART_SIZE[0] > 1)
    {
        std::vector<double> sleft;
        std::vector<double> rleft;
        std::vector<double> sright;
        std::vector<double> rright;
        void* request_sleft = create_request();
        void* request_sleft_size = create_request();
        void* request_sright = create_request();
        void* request_sright_size = create_request();
        void* request_rleft = create_request();
        void* request_rleft_size = create_request();
        void* request_rright = create_request();
        void* request_rright_size = create_request();

        int sleft_size;
        int sright_size;
        int RightProcess = MPI_CART_RANK[2] + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0] + 1) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        int LeftProcess  = MPI_CART_RANK[2] + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0] - 1) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];

        if(MPI_CART_RANK[0] >= MPI_CART_SIZE[0]-1)
        {
            RightProcess = MPI_CART_RANK[2] + MPI_CART_RANK[1] * MPI_CART_SIZE[2];
        }
        if(MPI_CART_RANK[0] <= 0)
        {
            LeftProcess  = MPI_CART_RANK[2] + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_SIZE[0] - 1) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        }

        const int LeftSizeTag  = 0; // used to identify data stream
        const int LeftDataTag  = 2; // used to identify data stream
        const int RightSizeTag = 4; // used to identify data stream
        const int RightDataTag = 8; // used to identify data stream

        if (MPI_CART_RANK[0] > 0 or MPIperiodicX)
        {
            std::vector<long int> window(6);
            window[0] = 0;
            window[1] = storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            sleft = storage.pack(window);
            sleft_size = sleft.size();
            OP_MPI_Isend(&sleft_size , 1 , OP_MPI_INT , LeftProcess , RightSizeTag, OP_MPI_COMM_WORLD , request_sleft_size);
            OP_MPI_Isend(sleft.data() , sleft_size , OP_MPI_DOUBLE , LeftProcess , RightDataTag, OP_MPI_COMM_WORLD , request_sleft);
        }
        if (MPI_CART_RANK[0] < MPI_CART_SIZE[0]-1 or MPIperiodicX)
        {
            std::vector<long int> window(6);
            window[0] = storage.sizeX()-storage.BcellsX();
            window[1] = storage.sizeX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            sright = storage.pack(window);
            sright_size = sright.size();
            OP_MPI_Isend(&sright_size , 1 , OP_MPI_INT , RightProcess , LeftSizeTag, OP_MPI_COMM_WORLD , request_sright_size);
            OP_MPI_Isend(sright.data() , sright_size , OP_MPI_DOUBLE , RightProcess , LeftDataTag, OP_MPI_COMM_WORLD , request_sright);
        }
        if (MPI_CART_RANK[0] > 0 or MPIperiodicX)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, LeftProcess, LeftSizeTag, OP_MPI_COMM_WORLD, request_rleft_size);
            OP_MPI_Wait(request_rleft_size, OP_MPI_STATUS_IGNORE);
            rleft.resize(size);
            OP_MPI_Irecv(rleft.data(), size, OP_MPI_DOUBLE, LeftProcess, LeftDataTag, OP_MPI_COMM_WORLD, request_rleft);
            OP_MPI_Wait(request_rleft, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = 0;
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            storage.unpack(rleft,window);
        }
        if (MPI_CART_RANK[0] < MPI_CART_SIZE[0]-1 or MPIperiodicX)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, RightProcess, RightSizeTag, OP_MPI_COMM_WORLD, request_rright_size);
            OP_MPI_Wait(request_rright_size, OP_MPI_STATUS_IGNORE);
            rright.resize(size);
            OP_MPI_Irecv(rright.data(), size, OP_MPI_DOUBLE, RightProcess, RightDataTag, OP_MPI_COMM_WORLD, request_rright);
            OP_MPI_Wait(request_rright, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = storage.sizeX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            storage.unpack(rright,window);
        }
        if (MPI_CART_RANK[0] > 0 or MPIperiodicX)
        {
            OP_MPI_Wait(request_sleft, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sleft_size, OP_MPI_STATUS_IGNORE);
        }
        if (MPI_CART_RANK[0] < MPI_CART_SIZE[0]-1 or MPIperiodicX)
        {
            OP_MPI_Wait(request_sright, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sright_size, OP_MPI_STATUS_IGNORE);
        }
        free_request(request_sleft);
        free_request(request_sleft_size);
        free_request(request_sright);
        free_request(request_sright_size);
        free_request(request_rleft);
        free_request(request_rleft_size);
        free_request(request_rright);
        free_request(request_rright_size);
    }
}

template<typename A>
inline void BoundaryConditions::CommunicateY([[maybe_unused]] A& storage) const
{
    if(MPI_CART_SIZE[1]>1)
    {
        std::vector<double> sleft;
        std::vector<double> rleft;
        std::vector<double> sright;
        std::vector<double> rright;
        void* request_sleft = create_request();
        void* request_sleft_size = create_request();
        void* request_sright = create_request();
        void* request_sright_size = create_request();
        void* request_rleft = create_request();
        void* request_rleft_size = create_request();
        void* request_rright = create_request();
        void* request_rright_size = create_request();

        int sleft_size;
        int sright_size;

        int RightProcess = MPI_CART_RANK[2] + (MPI_CART_RANK[1] + 1) * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        int LeftProcess  = MPI_CART_RANK[2] + (MPI_CART_RANK[1] - 1) * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];

        if(MPI_CART_RANK[1] >= MPI_CART_SIZE[1]-1)
        {
            RightProcess = MPI_CART_RANK[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        }
        if(MPI_CART_RANK[1] <= 0 )
        {
            LeftProcess  = MPI_CART_RANK[2] + (MPI_CART_SIZE[1] - 1) * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        }

        const int LeftSizeTag  = 0; // used to identify data stream
        const int LeftDataTag  = 2; // used to identify data stream
        const int RightSizeTag = 4; // used to identify data stream
        const int RightDataTag = 8; // used to identify data stream

        if (MPI_CART_RANK[1] > 0 or MPIperiodicY)
        {
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = 0;
            window[3] = storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            sleft = storage.pack(window);
            sleft_size = sleft.size();
            OP_MPI_Isend(&sleft_size , 1 , OP_MPI_INT , LeftProcess , RightSizeTag, OP_MPI_COMM_WORLD , request_sleft_size);
            OP_MPI_Isend(sleft.data() , sleft_size , OP_MPI_DOUBLE , LeftProcess , RightDataTag , OP_MPI_COMM_WORLD, request_sleft);
        }
        if (MPI_CART_RANK[1] < MPI_CART_SIZE[1]-1 or MPIperiodicY)
        {
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = storage.sizeY()-storage.BcellsY();
            window[3] = storage.sizeY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            sright = storage.pack(window);
            sright_size = sright.size();
            OP_MPI_Isend(&sright_size , 1 , OP_MPI_INT , RightProcess , LeftSizeTag, OP_MPI_COMM_WORLD , request_sright_size);
            OP_MPI_Isend(sright.data() , sright_size , OP_MPI_DOUBLE , RightProcess , LeftDataTag, OP_MPI_COMM_WORLD , request_sright);
        }
        if (MPI_CART_RANK[1] > 0 or MPIperiodicY)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, LeftProcess, LeftSizeTag, OP_MPI_COMM_WORLD, request_rleft_size);
            OP_MPI_Wait(request_rleft_size, OP_MPI_STATUS_IGNORE);
            rleft.resize(size);
            OP_MPI_Irecv(rleft.data(), size, OP_MPI_DOUBLE, LeftProcess, LeftDataTag, OP_MPI_COMM_WORLD, request_rleft);
            OP_MPI_Wait(request_rleft, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = 0;
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            storage.unpack(rleft,window);
        }
        if (MPI_CART_RANK[1] < MPI_CART_SIZE[1]-1 or MPIperiodicY)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, RightProcess, RightSizeTag, OP_MPI_COMM_WORLD, request_rright_size);
            OP_MPI_Wait(request_rright_size, OP_MPI_STATUS_IGNORE);
            rright.resize(size);
            OP_MPI_Irecv(rright.data(), size, OP_MPI_DOUBLE, RightProcess, RightDataTag, OP_MPI_COMM_WORLD, request_rright);
            OP_MPI_Wait(request_rright, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = storage.sizeY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            storage.unpack(rright,window);
        }
        if (MPI_CART_RANK[1] > 0 or MPIperiodicY)
        {
            OP_MPI_Wait(request_sleft, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sleft_size, OP_MPI_STATUS_IGNORE);
        }
        if (MPI_CART_RANK[1] < MPI_CART_SIZE[1]-1 or MPIperiodicY)
        {
            OP_MPI_Wait(request_sright, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sright_size, OP_MPI_STATUS_IGNORE);
        }
        free_request(request_sleft);
        free_request(request_sleft_size);
        free_request(request_sright);
        free_request(request_sright_size);
        free_request(request_rleft);
        free_request(request_rleft_size);
        free_request(request_rright);
        free_request(request_rright_size);
    }
}

template<typename A>
inline void BoundaryConditions::CommunicateZ([[maybe_unused]] A& storage) const
{
    if(MPI_CART_SIZE[2] > 1)
    {
        std::vector<double> sleft;
        std::vector<double> rleft;
        std::vector<double> sright;
        std::vector<double> rright;
        void* request_sleft = create_request();
        void* request_sleft_size = create_request();
        void* request_sright = create_request();
        void* request_sright_size = create_request();
        void* request_rleft = create_request();
        void* request_rleft_size = create_request();
        void* request_rright = create_request();
        void* request_rright_size = create_request();

        int sleft_size;
        int sright_size;

        int RightProcess = MPI_CART_RANK[2] +1 + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        int LeftProcess  = MPI_CART_RANK[2] -1 + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];

        if(MPI_CART_RANK[2] >= MPI_CART_SIZE[2]-1)
        {
            RightProcess =  MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        }
        if(MPI_CART_RANK[2] <= 0)
        {
            LeftProcess = MPI_CART_SIZE[2] - 1 + MPI_CART_RANK[1] * MPI_CART_SIZE[2] + (MPI_CART_RANK[0]) * MPI_CART_SIZE[1] * MPI_CART_SIZE[2];
        }

        const int LeftSizeTag  = 0; // used to identify data stream
        const int LeftDataTag  = 2; // used to identify data stream
        const int RightSizeTag = 4; // used to identify data stream
        const int RightDataTag = 8; // used to identify data stream

        if(MPI_CART_RANK[2] > 0 or MPIperiodicZ)
        {
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = 0;
            window[5] = storage.BcellsZ();
            sleft = storage.pack(window);
            sleft_size = sleft.size();
            OP_MPI_Isend(&sleft_size , 1 , OP_MPI_INT , LeftProcess , RightSizeTag, OP_MPI_COMM_WORLD , request_sleft_size);
            OP_MPI_Isend(sleft.data() , sleft_size , OP_MPI_DOUBLE , LeftProcess , RightDataTag, OP_MPI_COMM_WORLD , request_sleft);
        }
        if (MPI_CART_RANK[2] < MPI_CART_SIZE[2]-1 or MPIperiodicZ)
        {
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = storage.sizeZ()-storage.BcellsZ();
            window[5] = storage.sizeZ();
            sright = storage.pack(window);
            sright_size = sright.size();
            OP_MPI_Isend(&sright_size , 1 , OP_MPI_INT , RightProcess , LeftSizeTag, OP_MPI_COMM_WORLD , request_sright_size);
            OP_MPI_Isend(sright.data() , sright_size , OP_MPI_DOUBLE , RightProcess , LeftDataTag, OP_MPI_COMM_WORLD , request_sright);
        }
        if (MPI_CART_RANK[2] > 0 or MPIperiodicZ)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, LeftProcess, LeftSizeTag, OP_MPI_COMM_WORLD, request_rleft_size);
            OP_MPI_Wait(request_rleft_size, OP_MPI_STATUS_IGNORE);
            rleft.resize(size);
            OP_MPI_Irecv(rleft.data(), size, OP_MPI_DOUBLE, LeftProcess, LeftDataTag, OP_MPI_COMM_WORLD, request_rleft);
            OP_MPI_Wait(request_rleft, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = -storage.BcellsZ();
            window[5] = 0;
            storage.unpack(rleft,window);
        }
        if (MPI_CART_RANK[2] < MPI_CART_SIZE[2]-1 or MPIperiodicZ)
        {
            long int size = 0;
            OP_MPI_Irecv(&size, 1, OP_MPI_INT, RightProcess, RightSizeTag, OP_MPI_COMM_WORLD, request_rright_size);
            OP_MPI_Wait(request_rright_size, OP_MPI_STATUS_IGNORE);
            rright.resize(size);
            OP_MPI_Irecv(rright.data(), size, OP_MPI_DOUBLE, RightProcess, RightDataTag, OP_MPI_COMM_WORLD, request_rright);
            OP_MPI_Wait(request_rright, OP_MPI_STATUS_IGNORE);
            std::vector<long int> window(6);
            window[0] = -storage.BcellsX();
            window[1] = storage.sizeX()+storage.BcellsX();
            window[2] = -storage.BcellsY();
            window[3] = storage.sizeY()+storage.BcellsY();
            window[4] = storage.sizeZ();
            window[5] = storage.sizeZ()+storage.BcellsZ();
            storage.unpack(rright,window);
        }
        if (MPI_CART_RANK[2] > 0 or MPIperiodicZ)
        {
            OP_MPI_Wait(request_sleft, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sleft_size, OP_MPI_STATUS_IGNORE);
        }
        if (MPI_CART_RANK[2] < MPI_CART_SIZE[2]-1 or MPIperiodicZ)
        {
            OP_MPI_Wait(request_sright, OP_MPI_STATUS_IGNORE);
            OP_MPI_Wait(request_sright_size, OP_MPI_STATUS_IGNORE);
        }
        free_request(request_sleft);
        free_request(request_sleft_size);
        free_request(request_sright);
        free_request(request_sright_size);
        free_request(request_rleft);
        free_request(request_rleft_size);
        free_request(request_rright);
        free_request(request_rright_size);
    }
}

#endif

}// namespace openphase
#endif
