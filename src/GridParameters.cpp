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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "GridParameters.h"

#include "Includes.h"
#include "Settings.h"

namespace openphase
{

using namespace std;

GridParameters::GridParameters(std::string ObjectNameSuffix)
{
    thisclassname  = "GridParameters";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Idx = 0;

    Nx  = 1;
    Ny  = 1;
    Nz  = 1;
    Nz2 = 1;

    dNx = 1;
    dNy = 1;
    dNz = 1;

    TotalNx = 1;
    OffsetX = 0;

    TotalNy = 1;
    OffsetY = 0;

    TotalNz = 1;
    OffsetZ = 0;

    maxNx = 1;
    maxNy = 1;
    maxNz = 1;

    maxTotalNx = 1;
    maxTotalNy = 1;
    maxTotalNz = 1;

    dx = 1.0;

    Bcells = 1;
    iWidth = 1.0;
    Eta = 1.0;

    Resolution = Resolutions::Single;
}

GridParameters::GridParameters(const GridParameters& rhs)
{
    thisclassname  = rhs.thisclassname;
    thisobjectname = rhs.thisobjectname;

    Idx = rhs.Idx;

    Nx  = rhs.Nx;
    Ny  = rhs.Ny;
    Nz  = rhs.Nz;
    Nz2 = rhs.Nz2;

    dNx = rhs.dNx;
    dNy = rhs.dNy;
    dNz = rhs.dNz;

    TotalNx = rhs.TotalNx;
    OffsetX = rhs.OffsetX;

    TotalNy = rhs.TotalNy;
    OffsetY = rhs.OffsetY;

    TotalNz = rhs.TotalNz;
    OffsetZ = rhs.OffsetZ;

    maxNx = rhs.maxNx;
    maxNy = rhs.maxNy;
    maxNz = rhs.maxNz;

    maxTotalNx = rhs.maxTotalNx;
    maxTotalNy = rhs.maxTotalNy;
    maxTotalNz = rhs.maxTotalNz;

    dx = rhs.dx;

    Bcells = rhs.Bcells;
    iWidth = rhs.iWidth;
    Eta    = rhs.Eta;

    Resolution = rhs.Resolution;
}

GridParameters& GridParameters::operator=(const GridParameters& rhs)
{
    thisclassname  = rhs.thisclassname;
    thisobjectname = rhs.thisobjectname;

    Idx = rhs.Idx;

    Nx  = rhs.Nx;
    Ny  = rhs.Ny;
    Nz  = rhs.Nz;
    Nz2 = rhs.Nz2;

    dNx = rhs.dNx;
    dNy = rhs.dNy;
    dNz = rhs.dNz;

    TotalNx = rhs.TotalNx;
    OffsetX = rhs.OffsetX;

    TotalNy = rhs.TotalNy;
    OffsetY = rhs.OffsetY;

    TotalNz = rhs.TotalNz;
    OffsetZ = rhs.OffsetZ;

    maxNx = rhs.maxNx;
    maxNy = rhs.maxNy;
    maxNz = rhs.maxNz;

    maxTotalNx = rhs.maxTotalNx;
    maxTotalNy = rhs.maxTotalNy;
    maxTotalNz = rhs.maxTotalNz;

    dx = rhs.dx;

    Bcells = rhs.Bcells;
    iWidth = rhs.iWidth;
    Eta    = rhs.Eta;

    Resolution = rhs.Resolution;

    return *this;
}

GridParameters::GridParameters(int total_nx, int total_ny, int total_nz)
{
    SetDimensions(total_nx, total_ny, total_nz);
}

GridParameters GridParameters::DoubleResolution(void) const
{
    GridParameters DoubleGrid(*this);

    if(dNx)
    {
        DoubleGrid.Nx         *= 2;
        DoubleGrid.TotalNx    *= 2;
        DoubleGrid.OffsetX    *= 2;
        DoubleGrid.maxNx      *= 2;
        DoubleGrid.maxTotalNx *= 2;
    }
    if(dNy)
    {
        DoubleGrid.Ny         *= 2;
        DoubleGrid.TotalNy    *= 2;
        DoubleGrid.OffsetY    *= 2;
        DoubleGrid.maxNy      *= 2;
        DoubleGrid.maxTotalNy *= 2;
    }
    if(dNz)
    {
        DoubleGrid.Nz         *= 2;
        DoubleGrid.Nz2        *= 2;
        DoubleGrid.TotalNz    *= 2;
        DoubleGrid.OffsetZ    *= 2;
        DoubleGrid.maxNz      *= 2;
        DoubleGrid.maxTotalNz *= 2;
    }
    DoubleGrid.dx             *= 0.5;
    DoubleGrid.Eta            *= 0.5;

    return DoubleGrid;
}

void GridParameters::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ReadInput(data);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void GridParameters::ReadInput(std::stringstream& inp)
{
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Grid Parameters");

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    int tmpNx = FileInterface::ReadParameterI(inp, moduleLocation, string("Nx"));
    int tmpNy = FileInterface::ReadParameterI(inp, moduleLocation, string("Ny"));
    int tmpNz = FileInterface::ReadParameterI(inp, moduleLocation, string("Nz"));
    dx        = FileInterface::ReadParameterD(inp, moduleLocation, string("dx"));

    Bcells    = FileInterface::ReadParameterI(inp, moduleLocation, string("Bcells"), false, 1);
    iWidth    = FileInterface::ReadParameterD(inp, moduleLocation, {string("IWidth"), string("InterfaceWidth")});

    if (tmpNx + tmpNy + tmpNz == 0)
    {
        ConsoleOutput::WriteExit("All dimensions are suppressed: Nx=Ny=Nz=0!\n OpenPhase is not intended for simulations of zero-dimensional systems", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }

    string locResolution = FileInterface::ReadParameterK(inp, moduleLocation, string("Resolution"), false, "SINGLE");
    if(locResolution == "DOUBLE" or locResolution == "DUAL")
    {
        Resolution = Resolutions::Dual;
        Eta = iWidth*dx*0.5;
    }
    else if(locResolution == "SINGLE")
    {
        Resolution = Resolutions::Single;
        Eta = iWidth*dx;
    }
    else
    {
        ConsoleOutput::WriteExit("Wrong resolution selected -> " + locResolution, "Settings", "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }

#ifdef MPI_PARALLEL
    MPI_3D_DECOMPOSITION = FileInterface::ReadParameterB(inp, moduleLocation, std::string("MPI3D"), false, false);

    if(MPI_3D_DECOMPOSITION)
    {
        MPI_CART_SIZE[0] = FileInterface::ReadParameterI(inp, moduleLocation, string("Ncx"), false, 1);
        MPI_CART_SIZE[1] = FileInterface::ReadParameterI(inp, moduleLocation, string("Ncy"), false, 1);
        MPI_CART_SIZE[2] = FileInterface::ReadParameterI(inp, moduleLocation, string("Ncz"), false, 1);

        if(MPI_CART_SIZE[0] * MPI_CART_SIZE[1] * MPI_CART_SIZE[2] != MPI_SIZE)
        {
            string message = "The requested number of MPI parallel blocks can not be decomposed by requested number of 3D blocks.";
            ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }
    }
#endif

    SetDimensions(tmpNx,tmpNy,tmpNz);
    assert(tmpNx > 0 || dNx == 0);
    assert(tmpNy > 0 || dNy == 0);
    assert(tmpNz > 0 || dNz == 0);
    
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void GridParameters::ReadJSON(const string InputFileName)
{
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert("Grid Parameters");

    std::ifstream f(InputFileName);
    
    json data = json::parse(f);
    if (data.contains(thisclassname))
    {
        json grid = data[thisclassname];

        int tmpNx = FileInterface::ReadParameter<int>(grid, {"Nx"});
        int tmpNy = FileInterface::ReadParameter<int>(grid, {"Ny"});
        int tmpNz = FileInterface::ReadParameter<int>(grid, {"Nz"});
        dx        = FileInterface::ReadParameter<double>(grid, {"dx"});

        Bcells    = FileInterface::ReadParameter<int>(grid, {"Bcells"},1);
        iWidth    = FileInterface::ReadParameter<int>(grid, {"IWidth"});

        if (tmpNx + tmpNy + tmpNz == 0)
        {
            ConsoleOutput::WriteExit("All dimensions are suppressed: Nx=Ny=Nz=0!\n OpenPhase is not intended for simulations of zero-dimensional systems", thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }

        string locResolution = FileInterface::ReadParameter<std::string>(grid, {"Resolution"},"SINGLE");
        if(locResolution == "DOUBLE" or locResolution == "DUAL")
        {
            Resolution = Resolutions::Dual;
            Eta = iWidth*dx*0.5;
        }
        else if(locResolution == "SINGLE")
        {
            Resolution = Resolutions::Single;
            Eta = iWidth*dx;
        }
        else
        {
            ConsoleOutput::WriteExit("Wrong resolution selected -> " + locResolution, thisclassname, "ReadInput()");
            OP_Exit(EXIT_FAILURE);
        }

    #ifdef MPI_PARALLEL
        MPI_3D_DECOMPOSITION = FileInterface::ReadParameter<bool>(grid, {"MPI3D"},false);

        if(MPI_3D_DECOMPOSITION)
        {
            MPI_CART_SIZE[0] = FileInterface::ReadParameter<int>(grid, {"Ncx"},1);
            MPI_CART_SIZE[1] = FileInterface::ReadParameter<int>(grid, {"Ncy"},1);
            MPI_CART_SIZE[2] = FileInterface::ReadParameter<int>(grid, {"Ncz"},1);

            if(MPI_CART_SIZE[0] * MPI_CART_SIZE[1] * MPI_CART_SIZE[2] != MPI_SIZE)
            {
                string message = "The requested number of MPI parallel blocks can not be decomposed by requested number of 3D blocks.";
                ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
                OP_Exit(EXIT_FAILURE);
            }
        }
    #endif

        SetDimensions(tmpNx,tmpNy,tmpNz);
        assert(tmpNx > 0 || dNx == 0);
        assert(tmpNy > 0 || dNy == 0);
        assert(tmpNz > 0 || dNz == 0);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void GridParameters::Read(const Settings& locSettings, int tStep)
{
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,"GridParameters_", tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        OP_Exit(EXIT_FAILURE);
    };

    GridParameters tmpGrid;

    int tmptStep = 0;
    while(!inp.eof())
    {
        int tmpTotalNx;
        int tmpTotalNy;
        int tmpTotalNz;

        inp >> tmptStep >> tmpTotalNx >> tmpTotalNy >> tmpTotalNz;
        if(tmptStep <= tStep)
        {
            tmpGrid.SetDimensions(tmpTotalNx,tmpTotalNy,tmpTotalNz);
            tmpGrid.Idx = tmptStep;
        }
        else
        {
            break;
        }
    }
    *this = tmpGrid;
}

void GridParameters::Write(const Settings& locSettings, int tStep)
{
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,"GridParameters_", tStep, ".dat");

    ofstream out(FileName.c_str(), ios::app);

    if (!out)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write");
        OP_Exit(EXIT_FAILURE);
    };
    out << Idx << "\t" << TotalNx << "\t" << TotalNy << "\t" << TotalNz << endl;
}

#ifdef MPI_PARALLEL

void GridParameters::SetDimensions(int total_nx, int total_ny, int total_nz)
{
    if(MPI_3D_DECOMPOSITION)
    {
        OP_MPI_Cart_Setup(MPI_CART_SIZE,MPI_RANK,MPI_CART_RANK);

        if(total_nx != 0)
        {
            dNx     = 1;
            Nx      = total_nx / MPI_CART_SIZE[0];
            TotalNx = total_nx;
            OffsetX = MPI_CART_RANK[0] * Nx;

            int blocksize_restX = total_nx % MPI_CART_SIZE[0];

            if (blocksize_restX != 0 and MPI_CART_RANK[0] < MPI_CART_SIZE[0] - 1)
            {
                Nx += 1;
                OffsetX += MPI_CART_RANK[0];
            }
            if(blocksize_restX != 0 and MPI_CART_RANK[0] == MPI_CART_SIZE[0] - 1)
            {
                Nx -= MPI_CART_SIZE[0] - 1 - blocksize_restX;
                OffsetX  += MPI_CART_RANK[0];
            }

            maxNx      = std::max(maxNx, Nx);
            maxTotalNx = std::max(maxTotalNx, TotalNx);

            if(Nx < Bcells)
            {
                ConsoleOutput::WriteExit("Dimension X is too small for " + std::to_string(MPI_CART_SIZE[0]) + " MPI processes", thisclassname, "SetDimensions()");
                OP_Exit(EXIT_FAILURE);
            }
        }
        else if (MPI_CART_SIZE[0] == 1) // reduced X dimension
        {
            dNx        = 0;
            Nx         = 1;
            TotalNx    = 1;
            OffsetX    = 0;
            maxNx      = 1;
            maxTotalNx = 1;
        }
        else
        {
            ConsoleOutput::WriteExit("Dimension X is suppressed and cannot be decomposed into " + std::to_string(MPI_CART_SIZE[0]) + " MPI processes", thisclassname, "SetDimensions()");
            OP_Exit(EXIT_FAILURE);
        }

        if(total_ny != 0)
        {
            dNy     = 1;
            Ny      = total_ny / MPI_CART_SIZE[1];
            TotalNy = total_ny;
            OffsetY = MPI_CART_RANK[1] * Ny;

            int blocksize_restY = total_ny % MPI_CART_SIZE[1];

            if (blocksize_restY != 0 and MPI_CART_RANK[1] < MPI_CART_SIZE[1] - 1)
            {
                Ny += 1;
                OffsetY += MPI_CART_RANK[1];
            }
            if(blocksize_restY != 0 and MPI_CART_RANK[1] == MPI_CART_SIZE[1] - 1)
            {
                Ny -= MPI_CART_SIZE[1] - 1 - blocksize_restY;
                OffsetY  += MPI_CART_RANK[1];
            }

            maxNy      = std::max(maxNy, Ny);
            maxTotalNy = std::max(maxTotalNy, TotalNy);

            if (Ny < Bcells)
            {
                ConsoleOutput::WriteExit("Dimension Y is too small for " + std::to_string(MPI_CART_SIZE[1]) + " MPI processes", thisclassname, "SetDimensionsMPI()");
                OP_Exit(EXIT_FAILURE);
            }
        }
        else if (MPI_CART_SIZE[1] == 1) // reduced Y dimension
        {
            dNy        = 0;
            Ny         = 1;
            TotalNy    = 1;
            OffsetY    = 0;
            maxNy      = 1;
            maxTotalNy = 1;
        }
        else
        {
            ConsoleOutput::WriteExit("Dimension Y is suppressed and cannot be decomposed into " + std::to_string(MPI_CART_SIZE[1]) + " MPI processes", thisclassname, "SetDimensions()");
            OP_Exit(EXIT_FAILURE);
        }

        if(total_nz != 0)
        {
            dNz     = 1;
            Nz      = total_nz / MPI_CART_SIZE[2];
            TotalNz = total_nz;
            OffsetZ = MPI_CART_RANK[2] * Nz;

            int blocksize_restZ = total_nz % MPI_CART_SIZE[2];

            if (blocksize_restZ != 0 and MPI_CART_RANK[2] < MPI_CART_SIZE[2] - 1)
            {
                Nz += 1;
                OffsetZ += MPI_CART_RANK[2];
            }
            if(blocksize_restZ != 0 and MPI_CART_RANK[2] == MPI_CART_SIZE[2] - 1)
            {
                Nz -= MPI_CART_SIZE[2] - 1 - blocksize_restZ;
                OffsetZ  += MPI_CART_RANK[2];
            }

            maxNz      = std::max(maxNz, Nz);
            maxTotalNz = std::max(maxTotalNz, TotalNz);

            if (Nz < Bcells)
            {
                ConsoleOutput::WriteExit("Dimension Z is too small for " + std::to_string(MPI_CART_SIZE[2]) + " MPI processes", thisclassname, "SetDimensions()");
                OP_Exit(EXIT_FAILURE);
            }

            Nz2 = Nz/2 + 1;
        }
        else if (MPI_CART_SIZE[2] == 1) // reduced Z dimension
        {
            dNz        = 0;
            Nz         = 1;
            Nz2        = 1;
            TotalNz    = 1;
            OffsetZ    = 0;
            maxNz      = 1;
            maxTotalNz = 1;
        }
        else
        {
            ConsoleOutput::WriteExit("Dimension Z is suppressed and cannot be decomposed into " + std::to_string(MPI_CART_SIZE[2]) + " MPI processes", thisclassname, "SetDimensions()");
            OP_Exit(EXIT_FAILURE);
        }
        ConsoleOutput::WriteStandard(thisclassname + "(RANK " + std::to_string(MPI_RANK) + ")", "MPI 3D environment is initialized");
    }
    else
    {
        if(total_nx != 0)
        {
            dNx     = 1;
            Nx      = total_nx / MPI_SIZE;
            TotalNx = total_nx;
            OffsetX = Nx * MPI_RANK;

            int blocksize_rest = total_nx % MPI_SIZE;

            if (blocksize_rest != 0)
            {
                if (MPI_RANK < MPI_SIZE - 1)
                {
                    Nx += 1;
                    OffsetX += MPI_RANK;
                }
                else if (MPI_RANK == MPI_SIZE - 1)
                {
                    Nx -= MPI_SIZE - 1 - blocksize_rest;
                    OffsetX += MPI_RANK;
                }
            }
            if (Nx < Bcells)
            {
                ConsoleOutput::WriteExit("Dimension X is too small for the " + std::to_string(MPI_SIZE) + "MPI processes", thisclassname, "SetDimensionsMPI()");
                OP_Exit(EXIT_FAILURE);
            }
        }
        else if(MPI_SIZE == 1) // reduced X dimension
        {
            dNx        = 0;
            Nx         = 1;
            TotalNx    = 1;
            OffsetX    = 0;
            maxNx      = 1;
            maxTotalNx = 1;
        }
        else
        {
            ConsoleOutput::WriteExit("Dimension X is suppressed and cannot be decomposed into " + std::to_string(MPI_SIZE) + " MPI processes", thisclassname, "SetDimensions()");
            OP_Exit(EXIT_FAILURE);
        }

        if(total_ny != 0)
        {
            dNy        = 1;
            Ny         = total_ny;
            TotalNy    = total_ny;
            OffsetY    = 0;
            maxNy      = std::max(maxNy, Ny);
            maxTotalNy = std::max(maxTotalNy, TotalNy);
        }
        else // reduced Y dimension
        {
            dNy        = 0;
            Ny         = 1;
            TotalNy    = 1;
            OffsetY    = 0;
            maxNy      = 1;
            maxTotalNy = 1;
        }

        if(total_nz != 0)
        {
            dNz        = 1;
            Nz         = total_nz;
            Nz2        = Nz/2 + 1;
            TotalNz    = total_nz;
            OffsetZ    = 0;
            maxNz      = std::max(maxNz, Nz);
            maxTotalNz = std::max(maxTotalNz, TotalNz);
        }
        else // reduced Z dimension
        {
            dNz        = 0;
            Nz         = 1;
            Nz2        = 1;
            TotalNz    = 1;
            OffsetZ    = 0;
            maxNz      = 1;
            maxTotalNz = 1;
        }
        ConsoleOutput::WriteStandard(thisclassname + "(RANK " + std::to_string(MPI_RANK) + ")", "MPI 1D environment is initialized");
    }
}

#else

void GridParameters::SetDimensions(int total_nx, int total_ny, int total_nz)
{
    if(total_nx != 0)
    {
        dNx        = 1;
        Nx         = total_nx;
        TotalNx    = total_nx;
        OffsetX    = 0;
        maxNx      = std::max(maxNx, Nx);
        maxTotalNx = std::max(maxTotalNx, TotalNx);

    }
    else // reduced X dimension
    {
        dNx        = 0;
        Nx         = 1;
        TotalNx    = 1;
        OffsetX    = 0;
        maxNx      = 1;
        maxTotalNx = 1;
    }

    if(total_ny != 0)
    {
        dNy        = 1;
        Ny         = total_ny;
        TotalNy    = total_ny;
        OffsetY    = 0;
        maxNy      = std::max(maxNy, Ny);
        maxTotalNy = std::max(maxTotalNy, TotalNy);
    }
    else // reduced Y dimension
    {
        dNy        = 0;
        Ny         = 1;
        TotalNy    = 1;
        OffsetY    = 0;
        maxNy      = 1;
        maxTotalNy = 1;
    }

    if(total_nz != 0)
    {
        dNz        = 1;
        Nz         = total_nz;
        Nz2        = Nz/2 + 1;
        TotalNz    = total_nz;
        OffsetZ    = 0;
        maxNz      = std::max(maxNz, Nz);
        maxTotalNz = std::max(maxTotalNz, TotalNz);
    }
    else // reduced Z dimension
    {
        dNz        = 0;
        Nz         = 1;
        Nz2        = 1;
        TotalNz    = 1;
        OffsetZ    = 0;
        maxNz      = 1;
        maxTotalNz = 1;
    }
}

#endif
} //namespace openphase
