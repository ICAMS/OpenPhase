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
 *   Main contributors :   Efim Borukhovich; Philipp Engels; Oleg Shchyglo
 *
 */

#include "Orientations.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "Tools.h"
#include "Crystallography.h"
#include "Velocities.h"
#include "VTK.h"
#include "AdvectionHR.h"

namespace openphase
{
using namespace std;

Orientations::Orientations(Settings& locSettings)
{
    Initialize(locSettings);
}

void Orientations::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Orientations";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    size_t Bcells = Grid.Bcells;
    Quaternions.Allocate(Grid, Bcells);

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    initialized = true;
    ConsoleOutput::WriteStandard("Orientations", "Initialized");
}

void Orientations::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Quaternions);
    BC.SetY(Quaternions);
    BC.SetZ(Quaternions);
}

void Orientations::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Quaternions.Remesh(newNx, newNy, newNz);

    Grid.SetDimensions(newNx, newNy, newNz);

    SetBoundaryConditions(BC);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Quaternions, Quaternions.Bcells(),)
    {
        Quaternions(i, j, k).setRotationMatrix();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

Orientations& Orientations::operator= (const Orientations& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Orientations")
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;

        Quaternions = rhs.Quaternions;
    }
    return *this;
}

bool Orientations::Write(const Settings& locSettings, const int tStep) const
{
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,"Rotations_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "WriteRotations");
        OP_Exit(EXIT_FAILURE);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        for(int n = 0; n < 4; n++)
        {
            double tmp = Quaternions(i,j,k)[n];
            out.write(reinterpret_cast<char*>(&tmp), sizeof(double));
        }
    }
    return true;
}

bool Orientations::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir,"Rotations_",
                                                  tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be created",
                        thisclassname, "ReadRotations");
        OP_Exit(EXIT_FAILURE);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        double tmp[4];

        for(int n = 0; n < 4; n++)
        {
            inp.read(reinterpret_cast<char*>(&tmp[n]), sizeof(double));
        }
        Quaternions(i,j,k).set(tmp[0],tmp[1],tmp[2],tmp[3]);
    }
    return true;
}

void Orientations::WriteVTK(const Settings& locSettings, const int tStep,
                            const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"RQ_0", [this](int i,int j,int k){return Quaternions(i,j,k)[0];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_1", [this](int i,int j,int k){return Quaternions(i,j,k)[1];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_2", [this](int i,int j,int k){return Quaternions(i,j,k)[2];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_3", [this](int i,int j,int k){return Quaternions(i,j,k)[3];}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Quaternions_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Orientations::WriteTotalVTK(const Settings& locSettings,
                                 const PhaseField& Phase,
                                 const int tStep,
                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"R", [this,&Phase](int i,int j,int k){return getTotalRotation(Phase, i, j, k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "TotalRotations_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Orientations::WriteMisorientationsVTK(const Settings& locSettings,
                                           const int tStep,
                                           const int precision,
                                           const std::string measure) const
{
    double mult = 1.0;
    if (!measure.compare("deg")) mult*=180.0/Pi;

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Omega", [this,mult](int i,int j,int k){return Tools::getMisorientation(dMatrix3x3::UnitTensor(), Quaternions(i,j,k).RotationMatrix)*mult;}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Misorientations_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void Orientations::WriteGrainEBSDDataQuaternions(const Settings& locSettings,
                                                 const PhaseField& Phase,
                                                 const int tStep)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat real\t"
              << "Quat i\t"
              << "Quat j\t"
              << "Quat k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Quaternions,0)
    {
        EulerAngles tempAng;
        dMatrix3x3  tempMatrix;
        Quaternion  tempQuat;

        outbuffer << index << "\t";
        if (Phase.Fields(i,j,k).interface())
        {
            // linear interpolation in interface via quaternions
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                tempQuat += (Quaternions(i,j,k) + Phase.FieldsProperties[alpha->index].Orientation)*alpha->value;
            }
            tempQuat.normalize();
            // Set phase index to 0 in interface
            outbuffer << 0 << "\t";
        }
        else
        {
            tempQuat = Quaternions(i,j,k) + Phase.FieldsProperties[Phase.Fields(i,j,k).front().index].Orientation;
            outbuffer << Phase.FieldsProperties[Phase.Fields(i,j,k).front().index].Phase + 1 << "\t";
        }

        outbuffer << i*Phase.Grid.dx << "\t"
                  << j*Phase.Grid.dx << "\t"
                  << k*Phase.Grid.dx << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();

    ConsoleOutput::WriteStandard("EBSD file", FileName);
}

void Orientations::WriteRotated100VectorVTK(const Settings& locSettings,
                                            const int tStep,
                                            const int precision)
{
    stringstream outbufer;
    dVector3 X;
    X[0] = 1.0;
    X[1] = 0.0;
    X[2] = 0.0;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Rotated100Vector" << "\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET RECTILINEAR_GRID\n";
    outbufer << "DIMENSIONS " << Grid.Nx << " " << Grid.Ny << " " << Grid.Nz << "\n";
    outbufer << "X_COORDINATES " << Grid.Nx << " double\n";
    for (int i = 0; i < Grid.Nx; i++) outbufer << i << " ";
    outbufer << "\n";
    outbufer << "Y_COORDINATES " << Grid.Ny << " double\n";
    for (int j = 0; j < Grid.Ny; j++) outbufer << j << " ";
    outbufer << "\n";
    outbufer << "Z_COORDINATES " << Grid.Nz << " double\n";
    for (int k = 0; k < Grid.Nz; k++) outbufer << k << " ";
    outbufer << "\n";
    outbufer << "POINT_DATA " << Grid.LocalNumberOfCells() << "\n";

    for (int dir = 0; dir < 3; ++dir)
    {
        outbufer << "SCALARS X_" << dir << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Grid.Nz; k++)
        for (int j = 0; j < Grid.Ny; j++)
        for (int i = 0; i < Grid.Nx; i++)
        {
            outbufer << (Quaternions(i,j,k).RotationMatrix*X)[dir] <<"\n";
        }
    }
    outbufer << "VECTORS X double\n";

    for (int k = 0; k < Grid.Nz; k++)
    for (int j = 0; j < Grid.Ny; j++)
    for (int i = 0; i < Grid.Nx; i++)
    {
        outbufer << (Quaternions(i,j,k).RotationMatrix*X)[0] << " "
                 << (Quaternions(i,j,k).RotationMatrix*X)[1] << " "
                 << (Quaternions(i,j,k).RotationMatrix*X)[2] << "\n";
    }

    string FileName = FileInterface::MakeFileName(locSettings.VTKDir,"Rotated100Vector_",
                                                  tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

void Orientations::Advect(AdvectionHR& Adv,
                          const Velocities& Vel,
                          PhaseField& Phi,
                          const BoundaryConditions& BC,
                          const double dt, const double tStep)
{
    if(QuaternionsDot.IsNotAllocated())
    {
        QuaternionsDot.Allocate(Quaternions);
    }
    if(not QuaternionsDot.IsSize(Grid.Nx, Grid.Ny, Grid.Nz))
    {
        QuaternionsDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    }

    Adv.AdvectField(Quaternions, QuaternionsDot, Vel, BC, dt);
    ConsoleOutput::WriteStandard(thisclassname, "Advected");
}

void Orientations::PrintPointStatistics(int x, int y, int z)
{
    cout << "Point:      (" << x << ", " << y << ", " << z << ")" << endl;
    cout << "Quaternions:\n" << Quaternions(x,y,z).print() << endl;

    cout << endl;
}

}// namespace openphase
