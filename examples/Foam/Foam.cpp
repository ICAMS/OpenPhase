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

 *   File created :   2020
 *   Main contributors :   Samad Vakili; Oleg Shchyglo;
 *
 */


#include "Settings.h"
#include "RunTimeControl.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "AdvectionHR.h"
#include "Velocities.h"
#include "Tools/CSVParser.h"
#include "VTK.h"
#include "Includes.h"
#include "FluidDynamics/D3Q27.h"

using namespace std;
using namespace openphase;

class MPFSolver: public OPObject
{
 public:
    MPFSolver()
    {
    };
    MPFSolver(const Settings& locSettings,
              const string InputFileName = DefaultInputFileName,
              unsigned int boundary = 1)
    {
        Initialize(locSettings, boundary);
        ReadInput(InputFileName);
    };
    void Initialize(const Settings& locSettings, int boundary = 2);
    void ReadInput(string InputFileName);
    int Nx;
    int Ny;
    int Nz;
    int dNx;
    int dNy;
    int dNz;
    int Nphases;
    double dx;
    double dt;

    double mu;
    double specific_v;
    double Cs2_Bub;
    dVector3 Uinit;
    double IntEng;
    double IntEngC;
    double ReactiveMassDensity;                                                 ///< Density of the gas forming agent in liquid [kg/m^3]
    double ActiveLayer;                                                         ///< Thickness of the active (gas producing) layer around the bubble [m]

    dMatrixNxN Wsqr;

    dVectorN   Mass;                                                            ///< Mass of the individual bubbles
    dVectorN   f;                                                               ///< Bulk free energy density
    dVectorN   Rho0;                                                            ///< Initial density of phases

    Storage3D<NodeA<double>,0> RawIntF;
    dMatrix3x3 VelocityGradients(Velocities& Vel, const int i, const int j, const int k) const;  ///< Calculates local velocity gradient tensor
    dMatrix3x3 AverageVelocityGradients;                                        ///< Average velocity gradient in the simulation domain
    void CalculateAverageVelocityGradients(Velocities& Vel);                               ///< Calculates average velocity gradient in the simulation domain
    void MergeBubblesInContact(BoundaryConditions& BC, PhaseField& Phi, int tStep, double scale);
    void DensityUpdate(BoundaryConditions& BC, PhaseField& Phi, const int tStep);
    void PrintLocalData(PhaseField& Phi) const;
    dVector3 AdvectVelocity(Velocities& Vel, const int i, const int j, const int k);             ///< This is changed and it is different than the one in the library
    dVector3 ViscousForce(Velocities& Vel, const int i, const int j, const int k);
    dVector3 InterfacialForce(BoundaryConditions& BC,PhaseField& Phi, const int i, const int j, const int k);
    void SolveVelocity_PhaseField(BoundaryConditions& BC, PhaseField& Phi, Velocities& Vel);         // bulk force (hydrostatic pressure force)
    void InterfaceEnergyCurv(PhaseField& Phi);
    void WriteVTK(Settings& locSettings, PhaseField& Phi, Velocities& Vel,
                  const int tStep, const int precision = 16) const;
    void WriteSupplementary(const int tStep) const;
    void ReadSupplementary(const int tStep);
    string RawDataDir;
};

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
    Settings           OPSettings;
    OPSettings.ReadInput();

    RunTimeControl     RTC(OPSettings);
    PhaseField         Phi(OPSettings);
    Velocities         Vel(OPSettings);
    BoundaryConditions BC(OPSettings);
    MPFSolver          MPFS(OPSettings);

    MPFS.dt = RTC.dt;

    cout << "Initialization stage! Done!" << endl;

    if (RTC.Restart)
    {
        ConsoleOutput::WriteBlankLine();
        cout << "Restart data being read!" << endl;
        Phi.Read(OPSettings, BC, RTC.tStart);
        Vel.Read(OPSettings, BC, RTC.tStart);
        MPFS.ReadSupplementary(RTC.tStart);
        MPFS.DensityUpdate(BC, Phi, RTC.tStart);
        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        // Set initial geometry of the phases
        int LiquidPFIndex = Initializations::Single(Phi, 0, BC);
        int iRadius = 5;
        //random bubble nucleation
        Initializations::FillGrainWithSpheres(Phi,LiquidPFIndex,1,
                                iRadius, 3*iRadius, BC, 500);

        //Single bubble initialization
        //Initializations::Sphere(Phi, 1, iRadius, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);

        //Two bubble initialization
        //Initializations::Sphere(Phi, 1, iRadius, (OPSettings.Nx+1)/4, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);
        //Initializations::Sphere(Phi, 1, iRadius, 3*(OPSettings.Nx+1)/4, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC, OPSettings);

        MPFS.Mass.Allocate(Phi.FieldsProperties.size());
        MPFS.f.Allocate(Phi.FieldsProperties.size());

        for (size_t i = 0; i < Phi.FieldsProperties.size(); i++)
        {
            double Vol   = MPFS.specific_v * Phi.FieldsProperties[i].Volume;
            MPFS.Mass[i] = MPFS.Rho0[Phi.FieldsProperties[i].Phase]*Vol;
        }
        double rho_ref = 0.9*MPFS.Rho0[1];

        for (size_t n = 0; n < Phi.FieldsProperties.size(); n++)
        {
            size_t pIndex = Phi.FieldsProperties[n].Phase;
            if(pIndex == 1)
            {
                Phi.FieldsProperties[n].Density = 0.1*MPFS.Rho0[1]*(rand()%1000)/1000.0 + rho_ref;
            }
            else
            {
                Phi.FieldsProperties[n].Density = MPFS.Rho0[0];
            }
        }

        MPFS.f[0] = -(6.4*MPFS.Rho0[0]/(6000.0-MPFS.Rho0[0]) - 0.00000075*MPFS.Rho0[0]);

        Vel.SetAverage(MPFS.Uinit);
        MPFS.DensityUpdate(BC, Phi, 0);
        Vel.SetBoundaryConditions(BC);
    }

    /**************************************************************start from here*/

    ConsoleOutput::WriteSimple("Starting simulation...");

    ignore_result(system("if [ -d obs ] ; then rm -rf obs; fi"));
    ignore_result(system("mkdir obs"));

    ofstream out_time;
    out_time.open("obs/time_terms.dat");

    //-------------- The Time Loop -------------//
    for (RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        MPFS.InterfaceEnergyCurv(Phi);
        MPFS.SolveVelocity_PhaseField(BC, Phi, Vel);
        MPFS.MergeBubblesInContact(BC, Phi, RTC.tStep, 2.0);
        MPFS.DensityUpdate(BC, Phi, RTC.tStep);

        // Output to files
        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Vel.WriteVTK(OPSettings, RTC.tStep);
            MPFS.WriteVTK(OPSettings, Phi, Vel, RTC.tStep);
        }

        // Write restart output
        if (RTC.WriteRawData())
        {
            if (RTC.tStep == 20000 )  RTC.tFileWrite = 5000;
            if (RTC.tStep == 60000 )  RTC.tFileWrite = 10000;
            if (RTC.tStep == 100000)  RTC.tFileWrite = 20000;
            if (RTC.tStep == 200000)  RTC.tFileWrite = 50000;
            // Output PhaseField in binary format
            Phi.Write(OPSettings, RTC.tStep);
            Vel.Write(OPSettings, RTC.tStep);
            MPFS.WriteSupplementary(RTC.tStep);
        }

        // Output to screens
        if (RTC.WriteToScreen())
        {
            ConsoleOutput::WriteTimeStep(RTC.tStep, RTC.nSteps);
            //  Statistics
            MPFS.PrintLocalData(Phi);
        }

        out_time << RTC.tStep << " ";
        for (size_t n = 0; n < Phi.FieldsProperties.size(); n++)
        {
            out_time << Phi.FieldsProperties[n].Density << " "
                     << Phi.FieldsProperties[n].Volume*MPFS.specific_v << " "
                     << MPFS.f[n] << " ";
        }
        out_time << MPFS.IntEng << " " << MPFS.IntEngC << endl;
    }
    out_time.close();

    return 0;
}

void MPFSolver::Initialize(const Settings& locSettings, int boundary)
{
    thisclassname = "MPFSolver";

    Nx      = locSettings.Grid.Nx;
    Ny      = locSettings.Grid.Ny;
    Nz      = locSettings.Grid.Nz;
    dNx     = locSettings.Grid.dNx;
    dNy     = locSettings.Grid.dNy;
    dNz     = locSettings.Grid.dNz;
    dx      = locSettings.Grid.dx;

    Nphases = locSettings.Nphases;
    RawDataDir = locSettings.RawDataDir;

    specific_v = pow(dx, dNx + dNy + dNz);

    RawIntF.Allocate(Nx, Ny, Nz, dNx, dNy, dNz, 3);

    Rho0.Allocate(Nphases);
    Mass.Allocate(Nphases);
    f.Allocate(Nphases);
    Wsqr.Allocate(Nphases);                                                     //interface energy coefficients

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void MPFSolver::ReadInput(string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Multiphase Flow settings");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());

    fstream inpF(InputFileName.c_str(), ios::in);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    mu      = FileInterface::ReadParameterD(inp, moduleLocation, "Mu");
    Cs2_Bub = FileInterface::ReadParameterD(inp, moduleLocation, "Cs2_Bub");
    ReactiveMassDensity = FileInterface::ReadParameterD(inp, moduleLocation, "ReMass");
    ActiveLayer = FileInterface::ReadParameterD(inp, moduleLocation, "ALayer");

    for (int alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream counter;
        counter << alpha;
        Rho0[alpha] = FileInterface::ReadParameterD(inp, moduleLocation, string("Rho0_") + counter.str(), true, 0.0);
    }

    for (int alpha = 0; alpha < Nphases; alpha++)
    for (int beta = alpha; beta < Nphases; beta++)
    {
        stringstream counter;
        counter << alpha << "_" << beta;
        Wsqr(alpha,beta) = FileInterface::ReadParameterD(inp, moduleLocation, string("Wsqr_") + counter.str(), true, 0.0);
        Wsqr(beta,alpha) = Wsqr(alpha,beta);
    }

    Uinit[0] = FileInterface::ReadParameterD(inp, moduleLocation, string("u"));
    Uinit[1] = FileInterface::ReadParameterD(inp, moduleLocation, string("v"));
    Uinit[2] = FileInterface::ReadParameterD(inp, moduleLocation, string("w"));

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

dMatrix3x3 MPFSolver::VelocityGradients(Velocities& Vel, const int i, const int j, const int k) const
{
    dMatrix3x3 dVel;
    double factor = 0.5/dx;

    dVel(0,0) = (Vel.Average(i+1,j,k)[0] - Vel.Average(i-1,j,k)[0])*factor;
    dVel(0,1) = (Vel.Average(i,j+1,k)[0] - Vel.Average(i,j-1,k)[0])*factor;
    dVel(0,2) = (Vel.Average(i,j,k+1)[0] - Vel.Average(i,j,k-1)[0])*factor;

    dVel(1,0) = (Vel.Average(i+1,j,k)[1] - Vel.Average(i-1,j,k)[1])*factor;
    dVel(1,1) = (Vel.Average(i,j+1,k)[1] - Vel.Average(i,j-1,k)[1])*factor;
    dVel(1,2) = (Vel.Average(i,j,k+1)[1] - Vel.Average(i,j,k-1)[1])*factor;

    dVel(2,0) = (Vel.Average(i+1,j,k)[2] - Vel.Average(i-1,j,k)[2])*factor;
    dVel(2,1) = (Vel.Average(i,j+1,k)[2] - Vel.Average(i,j-1,k)[2])*factor;
    dVel(2,2) = (Vel.Average(i,j,k+1)[2] - Vel.Average(i,j,k-1)[2])*factor;

    return dVel;
}

void MPFSolver::CalculateAverageVelocityGradients(Velocities& Vel)
{
    dMatrix3x3 dVelAVG;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Vel.Average,0,reduction(dMatrix3x3SUM:dVelAVG))
    {
        dVelAVG += VelocityGradients(Vel,i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    AverageVelocityGradients = dVelAVG/(Nx*Ny*Nz);
}

void MPFSolver::MergeBubblesInContact(BoundaryConditions& BC ,PhaseField& Phi, int tStep, double scale)
{
    double p_crit = (dNx+dNy+dNz - 1.0)*Wsqr(0,1)*Pi*Pi/(8.0*Phi.Grid.Eta*Phi.Grid.Eta); //-2.0*f[0] is from average pressure

    int Nthreads = omp_get_max_threads();
    vector<NodeAB<double,double>> Overlap(Nthreads);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    {
        int thread_num = omp_get_thread_num();
        if(Phi.Fields(i,j,k).flag)
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend(); ++alpha)
        for (auto  beta = alpha + 1;
                   beta != Phi.Fields(i,j,k).cend(); beta++)
        if(Phi.FieldsProperties[alpha->index].Phase != 0 and
           Phi.FieldsProperties[ beta->index].Phase != 0)                       // if phases of both phase-fields are right
        {
            Overlap[thread_num].set_sym1(alpha->index, beta->index, 1.0);       // phase-fields overlap is detected

            if (fabs(f[alpha->index] + f[beta->index]) > p_crit*scale)          // if pressure from both bubbles is above critical
            {
                Overlap[thread_num].set_sym2(alpha->index, beta->index, 1.0);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for (int it = 1; it < Nthreads; it++)
    {
        Overlap[0] += Overlap[it];
    }

    for (auto it = Overlap[0].begin(); it < Overlap[0].end(); it++)
    {
        if (it->value1 == 1.0 and it->value2 == 1.0)                            // if both conditions (overlap and pressure) hold -> merge bubbles
        {
            cout << "Selective combine phase-fields: phase-field " << it->indexB
                 << " is merged into phase-field " << it->indexA
                 << " at time step: " << tStep << endl;
            Phi.SelectiveCombinePhaseFields(BC,it->indexA,it->indexB);
            Mass[it->indexA] += Mass[it->indexB];
            Mass[it->indexB] = 0;
        }
    }
}

void MPFSolver::DensityUpdate(BoundaryConditions& BC, PhaseField& Phi, const int tStep)
{
    double ReactiveMassReduction = 0.0;
    double LiquidVolume = 0.0;

    for (size_t Np = 0; Np < Phi.FieldsProperties.size(); Np++)
    {
        if(Phi.FieldsProperties[Np].Phase == 1)
        {
            switch(dNx+dNy+dNz)
            {
                case 2: //2D
                {
                    double locRadius = sqrt(Phi.FieldsProperties[Np].Volume/Pi)*dx;
                    double locSurfaceVolume = 2.0*Pi*locRadius*ActiveLayer;
                    double MassIncrease = locSurfaceVolume*ReactiveMassDensity;
                    ReactiveMassReduction += MassIncrease;
                    Mass[Np] += MassIncrease;
                    //cout << locSurfaceVolume*ReactiveMassDensity << endl;
                    break;
                }
                case 3: //3D
                {
                    double locRadius = pow(0.75*Phi.FieldsProperties[Np].Volume/Pi,1.0/3.0)*dx;
                    double locSurfaceVolume = 4.0*Pi*locRadius*locRadius*ActiveLayer;
                    double MassIncrease = locSurfaceVolume*ReactiveMassDensity;
                    ReactiveMassReduction += MassIncrease;
                    Mass[Np] += MassIncrease;
                    break;
                }
            }
        }
        else
        {
            LiquidVolume += Phi.FieldsProperties[Np].Volume;
        }
    }

    ReactiveMassDensity -= ReactiveMassReduction/(LiquidVolume*dx*dx*dx);

    for (size_t n = 0; n < Phi.FieldsProperties.size(); n++)
    if(Phi.FieldsProperties[n].Exist)
    {
        if(Phi.FieldsProperties[n].Phase == 1)
        {
            double vol = specific_v * Phi.FieldsProperties[n].Volume;
            Phi.FieldsProperties[n].Density = Mass[n] / vol;
            f[n] = -Cs2_Bub*Phi.FieldsProperties[n].Density;
        }
        else
        {
            f[n] = -(6.4*Phi.FieldsProperties[n].Density/(6000.0 - Phi.FieldsProperties[n].Density) - 0.00000075*Phi.FieldsProperties[n].Density);
        }
    }
}

void MPFSolver::PrintLocalData(PhaseField& Phi) const
{
    for(size_t idx = 0; idx < Phi.FieldsProperties.size(); idx++)
    {
        if(Phi.FieldsProperties[idx].Volume > 0.0)
        {
            ConsoleOutput::WriteStandardNarrow( "Phase-Field", idx);
            ConsoleOutput::WriteStandardNarrow( "Density", Phi.FieldsProperties[idx].Density);
            ConsoleOutput::WriteStandardNarrow( "Volume", Phi.FieldsProperties[idx].Volume*specific_v);
            //Info::WriteStandardNarrow( "dVavg", AverageVelocityGradients.print());
        }
    }
}

dVector3 MPFSolver::AdvectVelocity(Velocities& Vel, const int i, const int j, const int k)
{
    const double dxHalf_Inv = 0.5/dx;
    dVector3 locVec;

    for (auto dir = 0; dir < 3; dir++)
    {
        locVec[dir] = dxHalf_Inv*(
                ( Vel.Average(i,j,k)[0] - fabs(Vel.Average(i,j,k)[0]) )*(  Vel.Average(i+1,j,k)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[0] + fabs(Vel.Average(i,j,k)[0]) )*( -Vel.Average(i-1,j,k)[dir] + Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[1] - fabs(Vel.Average(i,j,k)[1]) )*(  Vel.Average(i,j+1,k)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[1] + fabs(Vel.Average(i,j,k)[1]) )*( -Vel.Average(i,j-1,k)[dir] + Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[2] - fabs(Vel.Average(i,j,k)[2]) )*(  Vel.Average(i,j,k+1)[dir] - Vel.Average(i,j,k)[dir] ) +
                ( Vel.Average(i,j,k)[2] + fabs(Vel.Average(i,j,k)[2]) )*( -Vel.Average(i,j,k-1)[dir] + Vel.Average(i,j,k)[dir] ) );
    }
    return locVec;
}

dVector3 MPFSolver::ViscousForce(Velocities& Vel, const int i, const int j, const int k)
{
    dVector3 locVec;
    double dx_sqr = dx*dx;

    for (int ii = -1; ii <= +1; ++ii)
    for (int jj = -1; jj <= +1; ++jj)
    for (int kk = -1; kk <= +1; ++kk)
    {
        locVec += Vel.Average(i+ii,j+jj,k+kk)*LBStencil3D[ii+1][jj+1][kk+1];
    }

    locVec = (locVec - Vel.Average(i,j,k)) * 6.0/dx_sqr;

    locVec[0] += (
                   Vel.Average(i+1,j,k)[0] + Vel.Average(i-1,j,k)[0] - 2.0*Vel.Average(i,j,k)[0] +
            0.25*( Vel.Average(i+1,j+1,k)[1]-Vel.Average(i+1,j-1,k)[1]-Vel.Average(i-1,j+1,k)[1]+Vel.Average(i-1,j-1,k)[1] ) +
            0.25*( Vel.Average(i+1,j,k+1)[2]-Vel.Average(i+1,j,k-1)[2]-Vel.Average(i-1,j,k+1)[2]+Vel.Average(i-1,j,k-1)[2] ) )/dx_sqr;

    locVec[1] += (
            0.25*( Vel.Average(i+1,j+1,k)[0]-Vel.Average(i+1,j-1,k)[0]-Vel.Average(i-1,j+1,k)[0]+Vel.Average(i-1,j-1,k)[0] ) +
                   Vel.Average(i,j+1,k)[1] + Vel.Average(i,j-1,k)[1] - 2.0*Vel.Average(i,j,k)[1] +
            0.25*( Vel.Average(i,j+1,k+1)[2]-Vel.Average(i,j+1,k-1)[2]-Vel.Average(i,j-1,k+1)[2]+Vel.Average(i,j-1,k-1)[2] ) )/dx_sqr;

    locVec[2] += (
            0.25*( Vel.Average(i+1,j,k+1)[0]-Vel.Average(i+1,j,k-1)[0]-Vel.Average(i-1,j,k+1)[0]+Vel.Average(i-1,j,k-1)[0] ) +
            0.25*( Vel.Average(i,j+1,k+1)[1]-Vel.Average(i,j+1,k-1)[1]-Vel.Average(i,j-1,k+1)[1]+Vel.Average(i,j-1,k-1)[1] ) +
                   Vel.Average(i,j,k+1)[2] + Vel.Average(i,j,k-1)[2] - 2.0*Vel.Average(i,j,k)[2] )/dx_sqr;

    return locVec;
}

dVector3 MPFSolver::InterfacialForce(BoundaryConditions& BC,PhaseField& Phi, const int i, const int j, const int k)
{
    dVector3 locVec;
    const double Prefactor = Pi*Pi/(Phi.Grid.Eta*Phi.Grid.Eta);

    double pre_f = 8.0/Pi;

    if (Phi.Fields(i,j,k).flag)
    {
        const NodePF& locPF = Phi.Fields(i,j,k);
        for (auto alpha = locPF.cbegin();
                  alpha != locPF.cend()-1; alpha++)
        for (auto beta = alpha+1;
                  beta != locPF.cend(); beta++)
        {
            int  pIndexA = Phi.FieldsProperties[alpha->index].Phase;
            int  pIndexB = Phi.FieldsProperties[ beta->index].Phase;

            double phi_a = f[beta->index]*pre_f*sqrt(alpha->value*beta->value) -
                           f[alpha->index]*pre_f*sqrt(alpha->value*beta->value) -
                           0.5*Wsqr(pIndexA,pIndexB)*(beta->laplacian - alpha->laplacian) -
                           0.5*Prefactor*Wsqr(pIndexA,pIndexB)*(beta->value - alpha->value);

            if(locPF.size() > 2)
            for(auto ksi  = locPF.cbegin();
                     ksi != locPF.cend(); ++ksi)
            if((ksi != alpha) && (ksi != beta))
            {
                int  pIndexG = Phi.FieldsProperties[ ksi->index].Phase;

                phi_a += - 0.5*(Wsqr(pIndexA,pIndexG) - Wsqr(pIndexB,pIndexG))*ksi->laplacian -
                           0.5*Prefactor*(Wsqr(pIndexA,pIndexG) - Wsqr(pIndexB,pIndexG))*ksi->value;
            }

            phi_a *= 1.0/Phi.Fields(i,j,k).size();

            RawIntF(i,j,k).add_value(alpha->index, phi_a);
            RawIntF(i,j,k).add_value(beta->index, -phi_a);

            locVec += (alpha->gradient - beta->gradient)*phi_a;
        }
    }
    return locVec;
}

void MPFSolver::SolveVelocity_PhaseField(BoundaryConditions& BC, PhaseField& Phi,
                                         Velocities& Vel)
{
    CalculateAverageVelocityGradients(Vel);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    {
        RawIntF(i,j,k).clear();

        dVector3 locVec_Adv = AdvectVelocity(Vel,i,j,k),
                 locVec_Vis = ViscousForce(Vel,i,j,k),
                 locVec_Int = InterfacialForce(BC,Phi,i,j,k);

        double RhoTotal = 0.0;
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend(); alpha++)
        {
            RhoTotal += Phi.FieldsProperties[alpha->index].Density*alpha->value;
        }

        if (RhoTotal < DBL_EPSILON)
        {
            cout << "*********     Error: Density is zero!     **********"<< endl;
            cout << i << " " << j << " " << k << " " << RhoTotal << endl;
            OP_Exit(EXIT_FAILURE);
        }
        double dt_rho = dt/RhoTotal;

        Vel.Average(i,j,k)[0] += -dt*locVec_Adv[0] + dt_rho*( mu*locVec_Vis[0] - locVec_Int[0] );
        Vel.Average(i,j,k)[1] += -dt*locVec_Adv[1] + dt_rho*( mu*locVec_Vis[1] - locVec_Int[1] );
        Vel.Average(i,j,k)[2] += -dt*locVec_Adv[2] + dt_rho*( mu*locVec_Vis[2] - locVec_Int[2] );
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0,)
    if (Phi.Fields(i,j,k).flag)
    {
        const NodePF& locPF = Phi.Fields(i,j,k);

        for (auto alpha  = locPF.cbegin();
                  alpha != locPF.cend(); alpha++)
        {
            double value = dt*( -(Vel.Average(i,j,k)*alpha->gradient)
                           + RawIntF(i,j,k).get_value(alpha->index));
            Phi.Fields(i,j,k).add_value(alpha->index, value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Phi.Finalize(BC);
    Vel.SetBoundaryConditions(BC);
}

void MPFSolver::InterfaceEnergyCurv(PhaseField& Phi)
{
    double sum[2] = {0.0, 0.0};
    int i = Nx/2;
    int j = Ny/2;
    for (int k = Nz/2+1; k < Nz; k++)
    {
        const NodePF& locPF = Phi.Fields(i,j,k);
        sum[0] += Wsqr(0,1)*locPF.get_gradient(1)[2] * locPF.get_gradient(1)[2]*dx;
        sum[1] += Wsqr(0,1)*locPF.get_gradient(1)[2] * locPF.get_gradient(1)[2]/(k-Nz/2);
    }
    IntEng  = sum[0];
    IntEngC = sum[1];
}

void MPFSolver::WriteVTK(Settings& locSettings, PhaseField& Phi, Velocities& Vel, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");

    ListOfFields.push_back((VTK::Field_t) {"RawInterfacial", [this](int i,int j,int k){return RawIntF(i,j,k).get_value(1);}});
    ListOfFields.push_back((VTK::Field_t) {"VelocityGradients_", [&Vel,this](int i,int j,int k){return VelocityGradients(Vel,i, j, k);}});
    ListOfFields.push_back((VTK::Field_t) {"Density",        [this, &Phi](int i,int j,int k)
    {
        double RhoTotal = 0.0;
        for (auto alpha = Phi.Fields(i,j,k).cbegin();
                  alpha != Phi.Fields(i,j,k).cend(); alpha++)
        {
            RhoTotal += Phi.FieldsProperties[alpha->index].Density*alpha->value;
        }
        return RhoTotal;
    }});

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void MPFSolver::WriteSupplementary(const int tStep) const
{
    std::string FileName = FileInterface::MakeFileName(RawDataDir,"SUP_",tStep,".dat");

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::cout << "File \"" << FileName << "\" could not be created! Terminating!!!" << std::endl;
        OP_Exit(EXIT_FAILURE);
    }

    double value = f[0];
    out.write(reinterpret_cast<char*>(&value), sizeof(double));
    for (size_t alpha = 0; alpha < Mass.size(); alpha++)
    {
        double valueR = Mass[alpha];
        out.write(reinterpret_cast<char*>(&valueR), sizeof(double));
    }
    out.close();
}

void MPFSolver::ReadSupplementary(const int tStep)
{
    std::string FileName = FileInterface::MakeFileName(RawDataDir,"SUP_", tStep, ".dat");

    std::fstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        OP_Exit(EXIT_FAILURE);
    }

    inp.read(reinterpret_cast<char*>(&f[0]), sizeof(double));
    for (size_t alpha = 0; alpha < Mass.size(); alpha++)
    {
        inp.read(reinterpret_cast<char*>(&Mass[alpha]), sizeof(double));
    }
    ConsoleOutput::WriteStandardNarrow( "Supplementary data", "Binary Input Read");
}
