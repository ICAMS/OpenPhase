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
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#include "MassDensity.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Temperature.h"
#include "BoundaryConditions.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

void MassDensity::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Density";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    size_t Bcells = Grid.Bcells;
    Phase.Allocate(Grid, {Nphases}, Bcells);
    Total.Allocate(Grid, Bcells);

    Initial.Allocate({Nphases});

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void MassDensity::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Density input");
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

void MassDensity::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        stringstream converter;
        converter << string("Rho0_") << alpha;

        Initial({alpha}) = FileInterface::ReadParameterD(inp, moduleLocation, converter.str());
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void MassDensity::Set(PhaseField& PF, Composition& Cx, Temperature& Tx)
{
    size_t Ncomp = Cx.Ncomp;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            Total(i,j,k) = 0.0;
        }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
        if(PF.Fields(i,j,k).wide_interface())
        {
            for(size_t comp = 0; comp < Ncomp; comp++)
            for(size_t alpha = 0; alpha != Nphases; alpha++)
            {
                Total(i,j,k) += PF.Fractions(i,j,k,{alpha})*
                                Cx.MoleFractions(i,j,k,{alpha, comp});
            }
        }
        else
        {
            size_t alpha = PF.FieldsProperties[PF.Fields(i,j,k).front().index].Phase;
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                Total(i,j,k) = Cx.MoleFractions(i,j,k,{alpha, comp});
            }
        }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MassDensity::SetInitial(PhaseField& PF)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Total,0,)
    {
        Total(i,j,k) = 0.0;
        Tensor<double,1> locFractions = PF.Fractions(i,j,k);
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            Phase(i,j,k,{alpha}) = Initial({alpha});
            Total(i,j,k) += locFractions({alpha})*Initial({alpha});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void MassDensity::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Total.Remesh(newNx, newNy, newNz);
    Phase.Remesh(newNx, newNy, newNz);

    Grid.SetDimensions(newNx, newNy, newNz);

    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void MassDensity::MoveFrame(const int dx, const int dy, const int dz, const BoundaryConditions& BC)
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

    SetBoundaryConditions(BC);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Phase(i,j,k) = Phase(i + dx, j + dy, k + dz);
        Total(i, j, k) = Total(i + dx, j + dy, k + dz);
    }

    SetBoundaryConditions(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Frame moved.");
}

void MassDensity::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Phase);
    BC.SetY(Phase);
    BC.SetZ(Phase);

    BC.SetX(Total);
    BC.SetY(Total);
    BC.SetZ(Total);
}

bool MassDensity::Write(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+"_"+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+"_", tStep, ".dat");
#endif

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        return false;
    };
    STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells())
        out.write(reinterpret_cast<const char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Total,Total.Bcells())
        out.write(reinterpret_cast<const char*>(&Total(i,j,k)), sizeof(double));
    STORAGE_LOOP_END
    return true;
}

bool MassDensity::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+"_"+std::to_string(MPI_RANK)+"_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname+"_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        return false;
    };

    STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells())
        inp.read(reinterpret_cast<char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Total,Total.Bcells())
        inp.read(reinterpret_cast<char*>(&Total(i,j,k)), sizeof(double));
    STORAGE_LOOP_END
    SetBoundaryConditions(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Binary Input Read Successfully.");
    return true;
}

void MassDensity::WriteVTK(int tStep, const Settings& locSettings, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t){"Total",  [this](int i,int j,int k){return Total(i,j,k);}});
    for(size_t n = 0; n < Nphases; n++)
    {
        ListOfFields.push_back((VTK::Field_t){"Phase_" + std::to_string(n), [n,this](int i,int j,int k){return Phase(i,j,k,{n});}});
    }
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void MassDensity::PrintPointStatistics(int x, int y, int z)
{
    cout << "Density:\t";
    cout << "Total:\t" << Total(x, y, z) << endl;

    for (size_t n = 0; n < Nphases; n++)
    {
        cout << "Phase " << n << ":\t" << Phase(x,y,z,{n}) << endl;
    }
}

MassDensity& MassDensity::operator=(const MassDensity& rhs)
{
    // protect against self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname != "")
    {
        thisclassname = rhs.thisclassname;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;

        if (Phase.IsNotAllocated())
        {
            Phase.Allocate(Grid, {Nphases}, rhs.Phase.Bcells());
            Total.Allocate(Grid, rhs.Total.Bcells());
            Initial.Allocate({Nphases});
        }
        else if (not Phase.IsSize(Grid.Nx, Grid.Ny, Grid.Nz))
        {
            Phase.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Total.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Initial.Reallocate({Nphases});
        }

        Phase = rhs.Phase;
        Total = rhs.Total;

        for (size_t n = 0; n < Nphases; n++)
        {
            Initial({n}) = rhs.Initial({n});
        }
    }
    return *this;
}

}// namespace openphase
