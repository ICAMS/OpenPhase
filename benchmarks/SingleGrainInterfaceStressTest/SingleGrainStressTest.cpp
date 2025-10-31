/*
 *   This file is a part of the OpenPhase software project.
 *   For more details visit www.openphase.de
 *
 *   Created:    2016
 *
 *   Authors:    Raphael Schiedung
 *
 *   Copyright (c) 2009-2020 Interdisciplinary Centre for Advanced Materials
 *                 Simulation (ICAMS). Ruhr-Universitaet Bochum. Germany
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "BoundaryConditions.h"
#include "DoubleObstacle.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"
#include "Tools/TimeInfo.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"

namespace op = openphase;

//void WriteSurfaceLinePlot(const op::PhaseField& Phase,
//        const op::ElasticProperties& EP, const int Index,  const int tStep);

void Diagnostics(op::Settings &OPSettings, op::PhaseField &Phase,
        op::InterfaceProperties &IP, op::ElasticProperties &EP,
        op::DoubleObstacle &DO, double &Radius);

#ifdef MPI_PARALLEL
   int MPI_RANK;
   int MPI_SIZE;
#endif

/****************************** <<< The Main >>> ******************************/
int main(int argc, char **argv)
{
    std::string InputFileName = op::DefaultInputFileName;

#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_SINGLE, &provided);
    op_fftw_mpi_init();
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    // Detects floating point errors on Linux machines
    //feenableexcept(FE_INEXACT);  // NOTE: THIS DETECTS ROUNDING ERRORS
    feenableexcept(FE_DIVBYZERO);  // Devision by zero
    feenableexcept(FE_INVALID);    // domain error occurred
    feenableexcept(FE_OVERFLOW);   // the result was too large
    feenableexcept(FE_UNDERFLOW);  // the result was too small

    op::BoundaryConditions       BC;
    op::DoubleObstacle           DO;
    op::ElasticProperties        EP;
    op::ElasticitySolverSpectral ES;
    op::InterfaceProperties      IP;
    op::PhaseField               Phase;
    op::RunTimeControl           RTC;
    op::Settings                 OPSettings;
    op::TimeInfo                 Time;

    OPSettings.Initialize();
    OPSettings.ReadInput(InputFileName);

    BC.Initialize(OPSettings);
    BC.ReadInput(InputFileName);
    //ES.Setup_MPI(OPSettings, BC);

    ES   .Initialize(OPSettings, BC);
    DO   .Initialize(OPSettings);
    EP   .Initialize(OPSettings);
    IP   .Initialize(OPSettings);
    Phase.Initialize(OPSettings);
    RTC  .Initialize(OPSettings);
    Time .Initialize(OPSettings, "Interface stress benchmark");

    ES   .ReadInput(InputFileName);
    DO   .ReadInput(InputFileName);
    EP   .ReadInput(InputFileName);
    IP   .ReadInput(InputFileName);
    Phase.ReadInput(InputFileName);
    RTC  .ReadInput(InputFileName);

    // Start timer
    Time.SetStart();

    /************************** Set initial condition**************************/
    // Set initial conditions

    double Radius = 0.0;
    if (RTC.Restart == true)
    {
        // Read restart input
        Phase.Read(OPSettings, BC, RTC.tStart);
    }
    else
    {
        op::Initializations::Single(Phase, 0, BC);

        // Initialize spherical inclusion
        const int    Nx   = OPSettings.Grid.Nx;
        const int    Ny   = OPSettings.Grid.Ny;
        const int    Nz   = OPSettings.Grid.Nz;
        const double C[3] = {Nx/2.0, Ny/2.0, Nz/2.0};
        const double R    = C[1]/6;

        //const double C[3] = {OPSettings.Nx/1.0,OPSettings.Ny/2.0,OPSettings.Nz/2.0};
        //const double R    = C[1]/6;

        op::ConsoleOutput::Write("Initialize grain with radius",R);
        op::Initializations::Sphere(Phase, 1, R, C[0], C[1], C[2], BC);

        Radius = R * OPSettings.Grid.dx;
    }
    Phase.SetBoundaryConditions(BC);

    // Calculate elastic properties
    IP.Set(Phase, BC);
    EP.SetEffectiveProperties(Phase, IP);
    Time.SetTimeStamp("Set initial condition");

    EP.WriteEigenStrainsVTK(OPSettings, 0);
    EP.WriteEnergyDensityVTK(OPSettings, Phase, 0);
    EP.WriteElasticStrainsVTK(OPSettings, 0);
    EP.WriteStressesVTK(OPSettings, 0);
    Phase.WriteVTK(OPSettings, 0);

    /************************** Analytic calculation **************************/
    /* Check for Eq.(A12) 'Laplace pressure' in
     * Schiedung, R., Steinbach, I., & Varnik, F. (2018).
     * Multi-phase-field method for surface tension induced elasticity.
     * Physical Review B, 97(3), 035410.*/
    // Get Lames's parameter and calculate prefactor of analytic solution
    //const double lambda0 = EP.PhaseElasticConstants[0](0,1);
    const double lambda1 = EP.PhaseElasticConstants[1](0,1);
    const double mu0     = EP.PhaseElasticConstants[0](3,3);
    const double mu1     = EP.PhaseElasticConstants[1](3,3);
          double prefac  = (3.0*lambda1+2.0*mu1)/(4.0*mu0+3.0*lambda1+2.0*mu1);
    // Check if the simulation is 2D
    if (OPSettings.Grid.Nx == 1 or OPSettings.Grid.Ny == 1 or OPSettings.Grid.Nz == 1)
    {
        prefac *= 0.5;
    }
    const double iSigma      = IP.InterfaceEnergy(0,1).Energy;
    const double dP_Analytic = (Radius > DBL_EPSILON) ? 2.0 * prefac * iSigma/Radius : 0;

    /************************** Numeric calculation **************************/
    op::ConsoleOutput::WriteSimple("Solve Problem...");
    size_t Iteration      = 0;
    double dP_Numeric     = dP_Analytic;
    double dP_NumericOld  = dP_Analytic;
    std::stringstream FName;
    FName << OPSettings.TextDir << "DeltaPressure.csv";
    std::fstream log(FName.str(), std::ios::trunc | std::ios::out);
    auto EndSolve = [&Phase, &EP, &Iteration, &dP_Numeric, &dP_NumericOld, &log]
    {
        double PressureInside  = 0.0;
        double PressureOutside = 0.0;
        int    CellsInside     = 0;
        int    CellsOuside     = 0;

        STORAGE_LOOP_BEGIN(i,j,k,EP.Stresses,0)
        if (Phase.Fields(i,j,k).flag == 0)
        {
            if (Phase.Fields(i,j,k).front().index == 1)
            {
                PressureInside += EP.Stresses(i,j,k).Pressure();
                CellsInside++;
            }
            else
            {
                PressureOutside += EP.Stresses(i,j,k).Pressure();
                CellsOuside++;
            }
        }
        STORAGE_LOOP_END

        if (CellsInside > 0) PressureInside  /= CellsInside;
        if (CellsOuside > 0) PressureOutside /= CellsOuside;

        dP_Numeric    = PressureInside - PressureOutside;
        double Delta  = (dP_Numeric - dP_NumericOld)/dP_NumericOld;

        //std::cout << Iteration << "," << dP_Numeric << "," << Delta << std::endl;
        log << Iteration << "," << dP_Numeric << "," << Delta << std::endl;

        dP_NumericOld = dP_Numeric;
        Iteration++;
        if (std::abs(Delta) < 1e-8) return true;
        else return false;
    };
    ES.Solve(EP, BC, RTC.dt, EndSolve);
    Time.SetTimeStamp("Solve for mechanical equilibrium");

    /***************************** Write Output ****************************/
    const double Error      = std::abs(dP_Numeric - dP_Analytic);
    const double ErrorRel   = Error/dP_Analytic;
    const double IntEnergy  = DO.Energy(Phase, IP, EP);
    const double ElaEnergy  = EP.Energy(Phase);
    const double TotEnergy  = IntEnergy + ElaEnergy;

    op::ConsoleOutput::WriteTimeStep(1,RTC.nSteps);
    op::ConsoleOutput::Write("Radius [m]", Radius);
    op::ConsoleOutput::Write("dPressure numeric [Pa]", dP_Numeric);
    op::ConsoleOutput::Write("dPressure analytic [Pa]", dP_Analytic);
    op::ConsoleOutput::Write("Relative error"         , ErrorRel);
    op::ConsoleOutput::WriteLine("-");
    op::ConsoleOutput::Write("Interface Energy [J]", IntEnergy);
    op::ConsoleOutput::Write("Elastic Energy [J]"  , ElaEnergy);
    op::ConsoleOutput::Write("Total Energy [J]"    , TotEnergy);
    op::ConsoleOutput::WriteBlankLine();

    // Write simulation results
    std::ofstream resultsSim ("Results.sim");
    if (resultsSim.is_open())
    {
        resultsSim << 2 << "\n";
        resultsSim << std::scientific << std::setprecision(8);
        resultsSim << "PressureError " << ErrorRel << " " << ErrorRel << "\n";
        resultsSim.close();
    }

    op::ConsoleOutput::WriteSimple("Write Output");
    op::ConsoleOutput::WriteLine();

    EP.WriteEigenStrainsVTK(OPSettings, 1);
    EP.WriteEnergyDensityVTK(OPSettings, Phase, 1);
    EP.WriteElasticStrainsVTK(OPSettings, 1);
    EP.WriteStressesVTK(OPSettings, 1);
    Phase.WriteVTK(OPSettings, 1);

    //WriteSurfaceLinePlot(Phase,EP, 1, 0);
    Time.SetTimeStamp("Write to disk");

    Time.PrintWallClockSummary();
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif
    return 0;
}

//void Diagnostics(op::Settings &OPSettings, op::PhaseField &Phase,
//        op::InterfaceProperties &IP, op::ElasticProperties &EP,
//        op::DoubleObstacle &DO, double &Radius)
//{
//    /***************************** Check Solution ****************************/
//
//    /* Check for Eq.(A12) 'Laplace pressure' in
//     * Schiedung, R., Steinbach, I., & Varnik, F. (2018).
//     * Multi-phase-field method for surface tension induced elasticity.
//     * Physical Review B, 97(3), 035410.*/
//
//    double PressureInside  = 0.0;
//    double PressureOutside = 0.0;
//    int    CellsInside     = 0;
//    int    CellsOuside     = 0;
//
//    STORAGE_LOOP_BEGIN(i,j,k,EP.Stresses,0)
//    if (Phase.Fields(i,j,k).flag == 0)
//    {
//        if (Phase.Fields(i,j,k).front().index == 1)
//        {
//            PressureInside += EP.Stresses(i,j,k).Pressure();
//            CellsInside++;
//        }
//        else
//        {
//            PressureOutside += EP.Stresses(i,j,k).Pressure();
//            CellsOuside++;
//        }
//    }
//    STORAGE_LOOP_END
//
//    if (CellsInside > 0) PressureInside  /= CellsInside;
//    if (CellsOuside > 0) PressureOutside /= CellsOuside;
//
//    // Get Lames's parameter and calculate prefactor of analytic solution
//    //const double lambda0 = EP.PhaseElasticConstants[0](0,1);
//    const double lambda1 = EP.PhaseElasticConstants[1](0,1);
//    const double mu0     = EP.PhaseElasticConstants[0](3,3);
//    const double mu1     = EP.PhaseElasticConstants[1](3,3);
//          double prefac  = (3.0*lambda1+2.0*mu1)/(4.0*mu0+3.0*lambda1+2.0*mu1);
//
//    // Check if the simulation is 2D
//    if (OPSettings.Nx == 1 or OPSettings.Ny == 1 or OPSettings.Nz == 1)
//    {
//        prefac *= 0.5;
//    }
//
//    const double iSigma         = IP.InterfaceEnergy(0,1).Energy;
//    const double dP_Numeric     = PressureInside - PressureOutside;
//    const double dP_Analytic    = (Radius > DBL_EPSILON) ? 2.0 * prefac * iSigma/Radius : 0;
//    const double Error          = std::abs(dP_Numeric - dP_Analytic);
//    const double ErrorRel       = Error/dP_Analytic;
//    const double IntEnergy      = DO.Energy(Phase, IP, EP);
//    const double ElaEnergy      = EP.Energy();
//    const double TotEnergy      = IntEnergy + ElaEnergy;
//
//}

/********************************* Diagnostics ********************************/
/****************************** NOT ESSENTIAL ! *******************************/

//op::vStress BulkStress(const op::ElasticProperties& EP, const int i,
//        const int j, const int k)
//{
//    return EP.EffectiveElasticConstants(i,j,k) * EP.TotalStrains(i,j,k);
//}
//
//op::vStrain IntStress(const op::ElasticProperties& EP, const int i,
//        const int j, const int k)
//{
//    return EP.EffectiveElasticConstants(i,j,k) *
//        EP.EffectiveEigenStrains(i,j,k) * (-1.);
//}
//
//double ForceXInterface(const op::PhaseField& Phase,
//        const op::ElasticProperties& EP, const int i, const int j,
//        const int k)
//{
//    // Calculate divergence of stress tensor in x-direction
//    const op::vStress Stress2PX = IntStress(EP,i+2,j,k);
//    const op::vStress Stress1PX = IntStress(EP,i+1,j,k);
//    const op::vStress Stress1MX = IntStress(EP,i-1,j,k);
//    const op::vStress Stress2MX = IntStress(EP,i-2,j,k);
//
//    const op::vStress Stress2PY = IntStress(EP,i,j+2,k);
//    const op::vStress Stress1PY = IntStress(EP,i,j+1,k);
//    const op::vStress Stress1MY = IntStress(EP,i,j-1,k);
//    const op::vStress Stress2MY = IntStress(EP,i,j-2,k);
//
//    const op::vStress Stress2PZ = IntStress(EP,i,j,k+2);
//    const op::vStress Stress1PZ = IntStress(EP,i,j,k+1);
//    const op::vStress Stress1MZ = IntStress(EP,i,j,k-1);
//    const op::vStress Stress2MZ = IntStress(EP,i,j,k-2);
//
//    const double fx = (- Stress2PX.get_tensor(0,0)*1./12.
//                       + Stress1PX.get_tensor(0,0)*2./3.
//                       - Stress1MX.get_tensor(0,0)*2./3.
//                       + Stress2MX.get_tensor(0,0)*1./12.)/Phase.dx
//                    + (- Stress2PY.get_tensor(1,0)*1./12.
//                       + Stress1PY.get_tensor(1,0)*2./3.
//                       - Stress1MY.get_tensor(1,0)*2./3.
//                       + Stress2MY.get_tensor(1,0)*1./12.)/Phase.dx
//                    + (- Stress2PZ.get_tensor(2,0)*1./12.
//                       + Stress1PZ.get_tensor(2,0)*2./3.
//                       - Stress1MZ.get_tensor(2,0)*2./3.
//                       + Stress2MZ.get_tensor(2,0)*1./12.)/Phase.dx;
//    return fx;
//}
//
//double ForceXBulk(const op::PhaseField& Phase,
//        const op::ElasticProperties& EP, const int i, const int j, const int k)
//{
//    const op::vStress Stress2PX = BulkStress(EP,i+2,j,k);
//    const op::vStress Stress1PX = BulkStress(EP,i+1,j,k);
//    const op::vStress Stress1MX = BulkStress(EP,i-1,j,k);
//    const op::vStress Stress2MX = BulkStress(EP,i-2,j,k);
//
//    const op::vStress Stress2PY = BulkStress(EP,i,j+2,k);
//    const op::vStress Stress1PY = BulkStress(EP,i,j+1,k);
//    const op::vStress Stress1MY = BulkStress(EP,i,j-1,k);
//    const op::vStress Stress2MY = BulkStress(EP,i,j-2,k);
//
//    const op::vStress Stress2PZ = BulkStress(EP,i,j,k+2);
//    const op::vStress Stress1PZ = BulkStress(EP,i,j,k+1);
//    const op::vStress Stress1MZ = BulkStress(EP,i,j,k-1);
//    const op::vStress Stress2MZ = BulkStress(EP,i,j,k-2);
//
//    const double fx = (- Stress2PX.get_tensor(0,0)*1./12.
//                       + Stress1PX.get_tensor(0,0)*2./3.
//                       - Stress1MX.get_tensor(0,0)*2./3.
//                       + Stress2MX.get_tensor(0,0)*1./12.)/Phase.dx
//                    + (- Stress2PY.get_tensor(1,0)*1./12.
//                       + Stress1PY.get_tensor(1,0)*2./3.
//                       - Stress1MY.get_tensor(1,0)*2./3.
//                       + Stress2MY.get_tensor(1,0)*1./12.)/Phase.dx
//                    + (- Stress2PZ.get_tensor(2,0)*1./12.
//                       + Stress1PZ.get_tensor(2,0)*2./3.
//                       - Stress1MZ.get_tensor(2,0)*2./3.
//                       + Stress2MZ.get_tensor(2,0)*1./12.)/Phase.dx;
//
//    return fx;
//}
//
//double CalcCurvature(const op::PhaseField& Phase, const int Index,
//        const int i, const int j, const int k, const double jj)
//{
//    op::Node Kappa10;
//    op::Node Kappa20;
//
//    op::Node Kappa11;
//    op::Node Kappa21;
//
//    std::tie(Kappa10, Kappa20) = Phase.PrincipalCurvatures(i,j,k);
//    std::tie(Kappa11, Kappa21) = Phase.PrincipalCurvatures(i,j-1,k);
//
//    double kappa0 = (Kappa10.get(Index) + Kappa20.get(Index));
//    double kappa1 = (Kappa11.get(Index) + Kappa21.get(Index));
//    double ratio  = (j-jj > 0 and i-jj < 1) ? i-jj : 0;
//
//    double kappa = (1.0-ratio)*kappa0 + ratio*kappa1;
//
//    return kappa;
//}
//
//double CalcContourLine(const op::PhaseField& Phase, const int Index,
//        const int i, const int j, const int k)
//{
//    const double A =
//          - Phase.Fields(i,j-1,k).get(Index)
//          + Phase.Fields(i,j,k).get(Index);
//    const double b = 0.5 - Phase.Fields(i,j,k).get(Index);
//
//    return (std::abs(A) > 1E-9)?(b/A + j):0.0;
//}
//
//void WriteSurfaceLinePlot(const op::PhaseField& Phase,
//        const op::ElasticProperties& EP, const int Index,  const int tStep)
//{
//    std::ofstream log("TextData/tab_data.csv");
//    log << std::scientific
//        << std::setprecision(16)
//        << "i,r,y,kappa,p,sxx,syy,szz,syz,sxz,sxy,gamma,Fx,FxInt,FxBulk"
//        << ",Bulksxx,Bulksyy,Bulkszz,Bulksyz,Bulksxz,Bulksxy"
//        << ",Intsxx,Intsyy,Intszz,Intsyz,Intsxz,Intsxy,I1,I2,I3"<< std::endl;
//
//    long double gamma = 0;
//    const long int i0 = (Phase.Fields.sizeX())/2;
//    const long int j0 = (Phase.Fields.sizeY())/2;
//    const long int k0 = (Phase.Fields.sizeZ())/2;
//    for(long int i = Phase.Fields.sizeX()/2; i <  Phase.Fields.sizeX()-1; i++)
//    {
//        double jj = 0;
//        double H  = 0;
//        //double phi = Phase.Fields(i,j0,k).get(Index);
//        for(long int j = j0; j <= Phase.Fields.sizeY(); j++)
//        if (Phase.Interface(i,j,k0))
//        {
//
//            const double bufferA = Phase.Fields(i,j-1,k0).get(Index) - 0.5;
//            const double bufferB = Phase.Fields(i,j  ,k0).get(Index) - 0.5;
//
//            if ((bufferA > 0 and bufferB < 0))
//            {
//                jj = j;
//                jj = CalcContourLine(Phase, Index, i, jj, k0);
//                H  = CalcCurvature(Phase, Index, i, jj, k0, jj);
//            }
//        }
//        // Calculate contour line in between the grid points
//        const double p_tang = -EP.Stresses(i,j0,k0)[1];
//        const double p_perp = -EP.Stresses(i,j0,k0)[0];
//        gamma += (p_perp - p_tang) * Phase.dx;
//
//        // Calculate stress invariants
//        const op::dVector3 I   = EP.Stresses(i,j0,k0).Invariants();
//
//        const double locForceXInterface = ForceXInterface(Phase,EP,i,j0,k0);
//        const double locForceXBulk      = ForceXBulk     (Phase,EP,i,j0,k0);
//        const double locForceX          = locForceXBulk + locForceXInterface;
//        const op::vStrain locBulkStress = BulkStress(EP,i,j0,k0);
//        const op::vStrain locIntStress  = IntStress(EP,i,j0,k0);
//        log << std::scientific;
//        log << i << "," << (i-i0)*Phase.dx << "," << jj << "," << H << ","
//            << EP.Stresses(i,j0,k0).Pressure() << "," //4
//            << EP.Stresses(i,j0,k0)[0] << "," //  5
//            << EP.Stresses(i,j0,k0)[1] << "," //  6
//            << EP.Stresses(i,j0,k0)[2] << "," //  7
//            << EP.Stresses(i,j0,k0)[3] << "," //  8
//            << EP.Stresses(i,j0,k0)[4] << "," //  9
//            << EP.Stresses(i,j0,k0)[5] << "," // 10
//            << gamma                      << "," // 11
//            << locForceX                  << "," // 12
//            << locForceXInterface         << "," // 13
//            << locForceXBulk              << "," // 14
//            << locBulkStress[0]           << "," // 15
//            << locBulkStress[1]           << "," // 16
//            << locBulkStress[2]           << "," // 17
//            << locBulkStress[3]           << "," // 18
//            << locBulkStress[4]           << "," // 19
//            << locBulkStress[5]           << "," // 20
//            << locIntStress[0]            << "," // 21
//            << locIntStress[1]            << "," // 22
//            << locIntStress[2]            << "," // 23
//            << locIntStress[3]            << "," // 24
//            << locIntStress[4]            << "," // 25
//            << locIntStress[5]            << "," // 26
//            << I[0]   << "," << I[1]   << "," << I[2] << std::endl;// 27 28 29
//    }
//    log.close();
//}
