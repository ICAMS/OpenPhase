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
 *   Main contributors :   Efim Borukhovich; Oleg Shchyglo
 *
 */

#include "ElectricalPotential.h"
#include "Composition.h"
#include "PhaseField.h"
#include "Settings.h"
#include "BoundaryConditions.h"
#include "VTK.h"

namespace openphase
{

ElectricalPotential::~ElectricalPotential()
{
    if (initialized)
    {
        //fftw_free(ftRHS);
        //fftw_free(ftPotential);
        //fftw_free(RHS);
        //fftw_free(rlPotential);
        delete[] ftRHS;
        delete[] ftPotential;
        delete[] RHS;
        delete[] rlPotential;
        fftw_destroy_plan(ForwardPlan);
        fftw_destroy_plan(BackwardPlan);
        delete[] Q[0];
        delete[] Q[1];
        delete[] Q[2];

#ifdef MPI_PARALLEL
        op_fftw_mpi_cleanup();
#endif

#ifdef _OPENMP
        fftw_cleanup_threads();
#endif

    }
}

using namespace std;
ElectricalPotential::ElectricalPotential(Settings& locSettings,
                                         const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}
void ElectricalPotential::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ElectricalPotential";
    thisobjectname = thisclassname + ObjectNameSuffix;

    EpsInv  = 1.0/8.85418781762e-12;

    CathodeCharge = 0.0;
    AnodeCharge   = 0.0;

    Grid = locSettings.Grid;

    Size    = Grid.LocalNumberOfCells();
    ftSize  = Grid.Nx*Grid.Ny*Grid.Nz2;

    Ncomp   = locSettings.Ncomp;
    Nphases = locSettings.Nphases;

    Q[0] = new double[ftSize];
    Q[1] = new double[ftSize];
    Q[2] = new double[ftSize];
    QXYZ();

    ftRHS       = new complex<double>[ftSize];
    ftPotential = new complex<double>[ftSize];
    RHS         = new double[Size];
    rlPotential = new double[Size];

    ForwardPlan  = fftw_plan_dft_r2c_3d(Grid.Nx, Grid.Ny, Grid.Nz, RHS,reinterpret_cast<fftw_complex*> (ftRHS),FFTW_ESTIMATE);
    BackwardPlan = fftw_plan_dft_c2r_3d(Grid.Nx, Grid.Ny, Grid.Nz, reinterpret_cast<fftw_complex*> (ftPotential), rlPotential, FFTW_ESTIMATE);

//    PotentialDataDir     = "PotentialData/";
//    int ignore = system(string("mkdir " + PotentialDataDir).c_str());
    size_t Bcells = Grid.Bcells;
    Potential.Allocate(Grid, Bcells);
//    Rho.Allocate(Dimensions, Bcells);

//    MolarVolume.Allocate(nPhases, nComp);
//    NuRef.Allocate(nPhases, nComp);
//    RhoRef.Allocate(nPhases, nComp);
//    Lattice.Allocate(nPhases, nComp);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void ElectricalPotential::ReadInput(string InputFileName)
{
    ConsoleOutput::WriteLineInsert("ElectricalPotential input");
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

void ElectricalPotential::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    Nelectrodes = FileInterface::ReadParameterI(inp, moduleLocation, "Nelectrodes");

    electrodes.resize(Nelectrodes);
    for(int i = 0; i < Nelectrodes; i++)
    {
        electrodes[i].index = i+1;
        electrodes[i].charge = 0;
        electrodes[i].interfaceVolume = 0.0;
        for(int dir = 0; dir < 3; dir++)
        {
            electrodes[i].position[dir] = 0;
        }
        stringstream converter;
        converter << "CntElectr" << i+1;
        electrodes[i].counterElectrode = FileInterface::ReadParameterI(inp, moduleLocation, converter.str())-1;
    }

    for (int q = 0 ; q < Ncomp; q++)
    {
        stringstream converter;
        converter << "Z" << q;
        ElementarCharges.push_back(FileInterface::ReadParameterI(inp, moduleLocation, converter.str()));
        MolarCharge.push_back(ElementarCharges[q] * PhysicalConstants::Faraday); // [C/mol]
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

//void ElectricalPotential::CalculateChargeDensity(MolarDensity& Nu, Diffusion& DF)
//{
//    #pragma omp parallel //OMP BEGIN
//    {
//        int nThreads = omp_get_num_threads();
//        int myThread = omp_get_thread_num();
//
//        int BegR=(myThread*(Nx+2))/nThreads;
//        int EndR=((myThread+1)*(Nx+2))/nThreads;
//
//        for (int i = BegR; i < EndR; i++)
//        for (int j =    0; j < Ny+2; j++)
//        for (int k =    0; k < Nz+2; k++)
//        {
//            Rho(i,j,k) = 0.;
//
//            double C[2];
//            C[0]=DF.Cx(i,j,k)/100.;
//            C[1]=1.-C[0];
//
//            for(int q = 0; q < nComp; q++)
//                Rho(i,j,k) += C[q] * Nu.Nu(i,j,k) * MolarCharge[q];
//        }//loop over the domain
//    }//OMP END
//}

//void ElectricalPotential::SetElectrodeCharges(PhaseField& Phase, Composition& Cx)
//{
//    CalculateElectrodeInterfaceVolumes(Phase);
//}

inline double ElectricalPotential::ChargeDensity(int i, int j, int k,
                                                 const PhaseField& Phase,
                                                 const Composition& Cx)
{
//    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
//    {
//        if(i == electrodes[eIdx].position[0] and
//           j == electrodes[eIdx].position[1] and
//           k == electrodes[eIdx].position[2])
//        {
//            return electrodes[eIdx].charge;
//        }
//    }
    if(Phase.Fields(i,j,k).interface())
    {
        int Index = -1;
        for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
        {
            if(Phase.Fields(i,j,k).get_value(eIdx+1) > 0) Index=eIdx;                     ///assumes that there is no interface between electrodes
        }

//        cout<<"charge at "<<i<<" "<<j<<" "<<k<<" is: "
//            <<electrodes[Index].charge/**Phase.Fields(i,j,k)[Index]/
//                       electrodes[Index].interfaceVolume +
//                       Cx.phase(i,j,k,0)*Phase.Fields(i,j,k)[0]*Faraday/(13.*8.472e-5)*/
//           <<endl;
        if(Index<0)
        {
            stringstream message;
            message << "Electrode Index not found!" << endl;
            ConsoleOutput::WriteExit(message.str(),thisclassname, "ChargeDensity()");
            exit(13);
        }
        return electrodes[Index].charge*Phase.Fields(i,j,k).get_value(Index+1)/
               electrodes[Index].interfaceVolume +
               Cx.MoleFractions(i,j,k,{0,0})*Phase.Fields(i,j,k).get_value(0)*PhysicalConstants::Faraday/(13.0*8.472e-5);
    }
    else
    {
        return Cx.MoleFractions(i,j,k,{0,0})*PhysicalConstants::Faraday/(13.0*8.472e-5);                        ///13 atoms in a molecule of PC, 8.472e-5 mol/m^3
    }
}

void ElectricalPotential::CalculateElectrodeInterfaceVolumes(const PhaseField& Phase) //TODO: parallelize
{
    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
    {
        electrodes[eIdx].interfaceVolume = 0;
    }//electrodes loop
    for(int eIdx = 0; eIdx < Nelectrodes; eIdx++)
    {
        for(int i = 0; i < Grid.Nx; i++)
        for(int j = 0; j < Grid.Ny; j++)
        for(int k = 0; k < Grid.Nz; k++)
        if(Phase.Fields(i,j,k).interface())
        {
            electrodes[eIdx].interfaceVolume += Phase.Fields(i,j,k).get_value(eIdx+1);
        }
    }//electrodes loop
}

//inline double ElectricalPotential::ChargeDensity(int i, int j, int k, MolarDensity& Nu, Composition& Cx)
//{
//    double result = 0.;
//
//    double C[2];
//    C[0]=Cx.total(i,j,k)/100.;
//    C[1]=1.-C[0];
//
//    for(int q = 0; q < Ncomp; q++)
//        result += C[q] * Nu.Nu(i,j,k) * MolarCharge[q];
//    return result;
//}

void ElectricalPotential::QXYZ(void)
{
    double DPi_nX = 2.0*Pi/double(Grid.Nx);
    double DPi_nY = 2.0*Pi/double(Grid.Ny);
    double DPi_nZ = 2.0*Pi/double(Grid.Nz);

    ///copied from the spectral solver
    for(int i = 0; i < Grid.Nx ; i++)
    for(int j = 0; j < Grid.Ny ; j++)
    for(int k = 0; k < Grid.Nz2; k++)
    {
         int XYZ = k + Grid.Nz2*(j + Grid.Ny*i);

         Q[0][XYZ] = DPi_nX*(i*(i <= Grid.Nx/2) - (Grid.Nx-i)*(i > Grid.Nx/2))/Grid.dx;
         Q[1][XYZ] = DPi_nY*(j*(j <= Grid.Ny/2) - (Grid.Ny-j)*(j > Grid.Ny/2))/Grid.dx;
         Q[2][XYZ] = DPi_nZ*(k*(k <= Grid.Nz/2) - (Grid.Nz-k)*(k > Grid.Nz/2))/Grid.dx;
    }
}

void ElectricalPotential::Solve(PhaseField& Phase, Composition& Cx, BoundaryConditions& BC)
{
    CalculateElectrodeInterfaceVolumes(Phase);
    for (int i = 0; i < Grid.Nx; i++)
    for (int j = 0; j < Grid.Ny; j++)
    for (int k = 0; k < Grid.Nz; k++)
    {
        RHS[k + Grid.Nz*(j + Grid.Ny*i)] = ChargeDensity(i, j, k, Phase, Cx);//Rho(i+1,j+1,k+1); //sin(i*6.28/Nx);
    }

    fftw_execute(ForwardPlan);

    for(int XYZ = 0; XYZ < ftSize; XYZ++)
    {
        /// ^phi = - ^rho / (i*q*i*q)
        double Q_sqr = Q[0][XYZ]*Q[0][XYZ] + Q[1][XYZ]*Q[1][XYZ] + Q[2][XYZ]*Q[2][XYZ];
        if (Q_sqr == 0)
        {
            ftPotential[XYZ] = 0.0;
        }
        else
        {
            ftPotential[XYZ] = ftRHS[XYZ] / Q_sqr * EpsInv;
        }
    }

    fftw_execute(BackwardPlan);

    for (int i = 0; i < Grid.Nx; i++)
    for (int j = 0; j < Grid.Ny; j++)
    for (int k = 0; k < Grid.Nz; k++)
    {
        Potential(i,j,k) = rlPotential[k + Grid.Nz*(j + Grid.Ny*i)]/Size;
    }

    ///Boundary Conditions (periodic only):
    BC.SetX(Potential);
    BC.SetY(Potential);
    BC.SetZ(Potential);
    // End Boundary Conditions
}

void ElectricalPotential::WriteChargeDensityVTK(const int tStep, const Settings& locSettings, const PhaseField& Phase, const Composition& Cx)
{
    CalculateElectrodeInterfaceVolumes(Phase);

    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ChargeDensity_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"ChargeDensity", [this,&Phase,&Cx](int i,int j,int k){return ChargeDensity(i, j, k, Phase, Cx)/*Rho(i,j,k)*/;}});
    VTK::Write(Filename, locSettings, ListOfFields);

}  //  WriteVTK

void ElectricalPotential::WritePotentialVTK(const int tStep, const Settings& locSettings) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ElectricalPotential_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) {"ElectricalPotential", [this](int i,int j,int k){return Potential(i,j,k);}});
    VTK::Write(Filename, locSettings, ListOfFields);
}  //  WriteVTK

}
