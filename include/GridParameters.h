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

 *   File created :   2023
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#ifndef GRIDPARAMETERS_H
#define GRIDPARAMETERS_H

#include "Containers/iVector3.h"

namespace openphase
{

class Settings;

class OP_EXPORTS GridParameters
{
 public:
    std::string thisclassname;
    std::string thisobjectname;

    GridParameters(std::string ObjectNameSuffix = "");
    GridParameters(const GridParameters& rhs);
    GridParameters(int total_nx, int total_ny, int total_nz);

    GridParameters& operator=(const GridParameters& rhs);

    void ReadInput(const std::string InputFileName = DefaultInputFileName);     ///< Reads input from the specified input file
    void ReadInput(std::stringstream& data);                                    ///< Reads input from the specified input file
    void ReadJSON(const std::string InputFileName);

    void Read(const Settings& locSettings, int tStep);                          ///< Reads grid parameters history from the file
    void Write(const Settings& locSettings, int tStep);                         ///< Writes grid parameters history to the file

    GridParameters DoubleResolution(void) const;                                ///< Returns grid parameters in double resolution

    void SetDimensions(int total_nx, int total_ny, int total_nz);               ///< Sets grid parameters

    int Active(void) const                                                      ///< Indicates dimensionality, e.g. 1D, 2D or 3D
    {
        return dNx + dNy + dNz;
    }

    long int LocalNumberOfCells(void) const                                     ///< Number of grid cells in the local simulation domain
    {
        return Nx * Ny * Nz;
    }

    double LocalVolume(void) const                                              ///< Volume in 3D, area in 2D and length in 1D of the local simulation domain
    {
        return LocalNumberOfCells()*std::pow(dx, dNx + dNy + dNz);
    }

    long int TotalNumberOfCells(void)  const                                    ///< Number of grid cells in the entire simulation domain
    {
        return TotalNx * TotalNy * TotalNz;
    }

    double TotalVolume(void) const                                              ///< Volume in 3D, area in 2D and length in 1D of the entire simulation domain
    {
        return TotalNumberOfCells()*std::pow(dx, dNx + dNy + dNz);
    }

    double CellVolume(bool ignore_active_dimensions = false) const              ///< Volume in 3D, area in 2D and length in 1D of a single grid cell unless ignore_active_dimensions == true, in which case volume is returned regardless of active dimensions
    {
        if (ignore_active_dimensions) return std::pow(dx, 3);
        return std::pow(dx, dNx + dNy + dNz);
    }

    bool PositionInBounds(const int x, const int y, const int z) const          ///< Indicates if given local grid coordinates are within local domain bounds
    {
        return x >= 0 and x < Nx and
               y >= 0 and y < Ny and
               z >= 0 and z < Nz;
    }

    bool PositionInBounds(const iVector3& position) const                       ///< Indicates if given local grid coordinates are within local domain bounds
    {
        return PositionInBounds(position[0],position[1],position[2]);
    }

    bool PositionInLocalBounds(const int x, const int y, const int z) const     ///< Indicates if given global grid coordinates are within local domain bounds
    {
        return x - OffsetX >= 0 and x - OffsetX < Nx and
               y - OffsetY >= 0 and y - OffsetY < Ny and
               z - OffsetZ >= 0 and z - OffsetZ < Nz;
    }

    bool PositionInLocalBounds(const iVector3& position) const                  ///< Indicates if given global grid coordinates are within local domain bounds
    {
        return PositionInLocalBounds(position[0],position[1],position[2]);
    }

    bool PositionInTotalBounds(const int x, const int y, const int z) const     ///< Indicates if given global grid coordinates are within total domain bounds
    {
        return x >= 0 and x < TotalNx and
               y >= 0 and y < TotalNy and
               z >= 0 and z < TotalNz;
    }

    bool PositionInTotalBounds(const iVector3& position) const                  ///< Indicates if given global position coordinates are within total domain bounds
    {
        return PositionInTotalBounds(position[0],position[1],position[2]);
    }

    iVector3 ConvertToLocal(const iVector3& global_position) const              ///< Converts global grid coordinates from global to local domain coordinates (relevant for MPI mode)
    {
        return global_position - iVector3{OffsetX,OffsetY,OffsetZ};
    }

    iVector3 ConvertToGlobal(const iVector3& local_position) const              ///< Converts grid coordinates from local domain coordinates to global (relevant for MPI mode)
    {
        return local_position + iVector3{OffsetX,OffsetY,OffsetZ};
    }

    int Idx;                                                                    ///< Index of the current instance of the grid parameters (used to store grid history)

    int Nx;                                                                     ///< Number of grid cells in X direction
    int Ny;                                                                     ///< Number of grid cells in Y direction
    int Nz;                                                                     ///< Number of grid cells in Z direction
    int Nz2;                                                                    ///< Half of the number of grid cells in Z direction for FFT routines

    int dNx;                                                                    ///< Stencil step in X direction for reduced dimensions simulations
    int dNy;                                                                    ///< Stencil step in Y direction for reduced dimensions simulations
    int dNz;                                                                    ///< Stencil step in Z direction for reduced dimensions simulations

    int TotalNx;                                                                ///< Total number of grid cells in X direction (in case of MPI parallelism)
    int OffsetX;                                                                ///< X coordinate of the local domain coordinate system origin (in case of MPI parallelism)

    int TotalNy;                                                                ///< Total system size in Y direction (in case of MPI parallelism)
    int OffsetY;                                                                ///< Y coordinate of the local domain coordinate system origin (in case of MPI parallelism)

    int TotalNz;                                                                ///< Total system size in Z direction (in case of MPI parallelism)
    int OffsetZ;                                                                ///< Z coordinate of the local domain coordinate system origin (in case of MPI parallelism)

    int maxNx;                                                                  ///< Maximum number of grid cells in X direction during the simulation
    int maxNy;                                                                  ///< Maximum number of grid cells in Y direction during the simulation
    int maxNz;                                                                  ///< Maximum number of grid cells in Z direction during the simulation

    int maxTotalNx;                                                             ///< Maximum total number of grid cells in X direction during the simulation
    int maxTotalNy;                                                             ///< Maximum total number of grid cells in Y direction during the simulation
    int maxTotalNz;                                                             ///< Maximum total number of grid cells in Z direction during the simulation

    double dx;                                                                  ///< Grid spacing in true units

    int Bcells;                                                                 ///< Number of boundary cells for setting the boundary conditions
    double iWidth;                                                              ///< Interface width in grid cells
    double Eta;                                                                 ///< Interface width in true units

    Resolutions Resolution;                                                     ///< Phase field resolution (Single or Double)
};

}
#endif //GRIDPARAMETERS
