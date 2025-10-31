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

 *   File created :   2022
 *   Main contributors :   Oleg Shchyglo; Hossein Jafarzadeh, Muhammad Adil Ali
 *
 */

#include "Includes.h"
#include "FractureField.h"
#include "ElasticProperties.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "DoubleObstacle.h"

#include "Composition.h"
#include "Temperature.h"

#include "Settings.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include "Velocities.h"
#include "AdvectionHR.h"
#include "RunTimeControl.h"

namespace openphase
{
using namespace std;

void FractureField::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "FractureField";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid         = locSettings.Grid;
    Ncomp        = locSettings.Ncomp;
    ElementNames = locSettings.ElementNames;

    sWidth = 1;

    Mobility      = 1.0;
    Xi            = 0.9;
    delta_g_0     = 30000;
    PFcutOff      = 0.5;
    InitialSurfaceEnergy = 1.0;

    FractureFieldLaplacianStencil = LaplacianStencils::Isotropic;
    FractureFieldGradientStencil  = GradientStencils::Isotropic;
    FractureModel = FractureModels::Obstacle;

    Fields          .Allocate(Grid, Grid.Bcells);
    Fields_dot      .Allocate(Grid, Grid.Bcells);
    Laplacian       .Allocate(Grid, Grid.Bcells);
    Flag            .Allocate(Grid, Grid.Bcells);
    DisplacementsOLD.Allocate(Grid, Grid.Bcells);
    SurfaceEnergy     .Allocate(Grid, Grid.Bcells);
    TestOutput.Allocate(Grid, Grid.Bcells);
    TestOutput2.Allocate(Grid, Grid.Bcells);

    switch(Grid.Active())
    {
        case 1:
        {
            LStencil.Set(LaplacianStencil1D_3, Grid);
            GStencil.Set(GradientStencil1D, Grid);
            break;
        }
        case 2:
        {
            switch(FractureFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil2D_5, Grid);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil2D_9, Grid);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil2D_LB, Grid);
                    break;
                }
            }

            switch(FractureFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, Grid);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil2D, Grid);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil2D_LB, Grid);
                    break;
                }
            }
            break;
        }
        case 3:
        {
            switch(FractureFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil3D_7, Grid);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil3D_27a, Grid);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil3D_LB, Grid);
                    break;
                }
            }

            switch(FractureFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, Grid);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil3D, Grid);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil3D_LB, Grid);
                    break;
                }
            }
            break;
        }
    }
    initialized = true;
    locSettings.AddForRemeshing(*this);
    locSettings.AddForAdvection(*this);
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void FractureField::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Fracture field input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void FractureField::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    sWidth        = FileInterface::ReadParameterD(inp, moduleLocation, string("SWidth"),    true, 1.0);
    Mobility      = FileInterface::ReadParameterD(inp, moduleLocation, string("Mu"),        true, 1.0);
    delta_g_0     = FileInterface::ReadParameterD(inp, moduleLocation, string("delta_g_0"), false, 30000);
    string mode   = FileInterface::ReadParameterK(inp, moduleLocation, string("FractureModel"), false, string("Obstacle"));
    InitialSurfaceEnergy = FileInterface::ReadParameterD(inp, moduleLocation, string("Sigma"),     true, 1.0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        SurfaceEnergy(i, j, k) = InitialSurfaceEnergy;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    ChemoMechanicalCoupling = FileInterface::ReadParameterB(inp, moduleLocation, "ChemoMechanicalCoupling", false, false);
    if(ChemoMechanicalCoupling)
    {
        Xi = FileInterface::ReadParameterD(inp, moduleLocation, string("Xi"));
        std::string RefCompName = FileInterface::ReadParameterK(inp, moduleLocation, "RefElement");
        for(size_t n = 0; n < Ncomp; n++)
            if(ElementNames[n] != RefCompName) Comp = n;
    }

    string tmp2 = FileInterface::ReadParameterK(inp, moduleLocation, string("FractureFieldLaplacianStencil"), false, string("ISOTROPIC"));
    if(tmp2 == "SIMPLE")
    {
        FractureFieldLaplacianStencil = LaplacianStencils::Simple;
    }
    if(tmp2 == "ISOTROPIC")
    {
        FractureFieldLaplacianStencil = LaplacianStencils::Isotropic;
    }
    if(tmp2 == "LB")
    {
        FractureFieldLaplacianStencil = LaplacianStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong Laplacian stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    string tmp3 = FileInterface::ReadParameterK(inp, moduleLocation, string("FractureFieldGradientStencil"), false, string("ISOTROPIC"));
    if(tmp3 == "SIMPLE")
    {
        FractureFieldGradientStencil = GradientStencils::Simple;
    }
    if(tmp3 == "ISOTROPIC")
    {
        FractureFieldGradientStencil = GradientStencils::Isotropic;
    }
    if(tmp3 == "LB")
    {
        FractureFieldGradientStencil = GradientStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong gradient stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    switch(Grid.Active())
    {
        case 1:
        {
            LStencil.Set(LaplacianStencil1D_3, Grid);
            GStencil.Set(GradientStencil1D, Grid);
            break;
        }
        case 2:
        {
            switch(FractureFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil2D_5, Grid);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil2D_9, Grid);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil2D_LB, Grid);
                    break;
                }
            }

            switch(FractureFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, Grid);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil2D, Grid);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil2D_LB, Grid);
                    break;
                }
            }
            break;
        }
        case 3:
        {
            switch(FractureFieldLaplacianStencil)
            {
                case LaplacianStencils::Simple:
                {
                    LStencil.Set(LaplacianStencil3D_7, Grid);
                    break;
                }
                case LaplacianStencils::Isotropic:
                {
                    LStencil.Set(LaplacianStencil3D_27a, Grid);
                    break;
                }
                case LaplacianStencils::LB:
                {
                    LStencil.Set(LaplacianStencil3D_LB, Grid);
                    break;
                }
            }

            switch(FractureFieldGradientStencil)
            {
                case GradientStencils::Simple:
                {
                    GStencil.Set(GradientStencil1D, Grid);
                    break;
                }
                case GradientStencils::Isotropic:
                {
                    GStencil.Set(GradientStencil3D, Grid);
                    break;
                }
                case GradientStencils::LB:
                {
                    GStencil.Set(GradientStencil3D_LB, Grid);
                    break;
                }
            }
            break;
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void FractureField::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    if (Fields_dot.IsNotAllocated()) Fields_dot.Allocate(Fields);
    if(not Fields_dot.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) Fields_dot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    Adv.AdvectField(Fields, Fields_dot, Vel, BC, dt);
}

void FractureField::CreatePF(PhaseField& Phase, size_t PhaseIndex, BoundaryConditions& BC)
{
    size_t index = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells())
    {
        if(Fields(i, j, k) > PFcutOff)
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i, j, k).set_value(index, 1.0);
        }
    }
    STORAGE_LOOP_END
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    if(Phase.Fields(i,j,k).size() == 1)
    {
        size_t alpha_index = Phase.Fields(i,j,k).begin()->index;

        for(int di = -1; di <= 1; di++)
        for(int dj = -1; dj <= 1; dj++)
        for(int dk = -1; dk <= 1; dk++)
        if(i+di >= 0 and i+di < Grid.Nx and
           j+dj >= 0 and j+dj < Grid.Ny and
           k+dk >= 0 and k+dk < Grid.Nz and
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
}

void FractureField::CreateCrack(dVector3& CrackStart, dVector3& CrackEnd, Settings& OP, BoundaryConditions& BC)
{
    double CrackLine[3]  = {-CrackEnd[1] + CrackStart[1], ///< {a, b, c} if crack line equation is ax+by+c=0
                             CrackEnd[0] - CrackStart[0],
                            (CrackEnd[1] - CrackStart[1]) * CrackStart[0] -
                            (CrackEnd[0] - CrackStart[0]) * CrackStart[1]};

    double CrackLength =
        sqrt((CrackEnd[1] - CrackStart[1]) * (CrackEnd[1] - CrackStart[1]) + ///< Crack length in xy plane
             (CrackEnd[0] - CrackStart[0]) * (CrackEnd[0] - CrackStart[0]));

    STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells())
    {
        Fields(i, j, k) = 0;
        Flag(i, j, k)   = 2;
        double Distance = sqrt((CrackLine[0] * i + CrackLine[1] * j + CrackLine[2]) * ///< distance of the point (i, j, k) from line ax+by+c=0
                               (CrackLine[0] * i + CrackLine[1] * j + CrackLine[2]) /
                               (CrackLine[0] * CrackLine[0] + CrackLine[1] * CrackLine[1]));

        if (Distance < (sWidth * 1))
        {
            double DistanceStart = (i - CrackStart[0]) * (i - CrackStart[0]) + (j - CrackStart[1]) * (j - CrackStart[1]);
            double DistanceEnd   = (i - CrackEnd[0]) * (i - CrackEnd[0]) + (j - CrackEnd[1]) * (j - CrackEnd[1]);
            double DistanceSq    = DistanceStart + DistanceEnd;

            if (DistanceSq < 2 * Distance * Distance + CrackLength * CrackLength)
            {
                Fields(i, j, k) = pow((1 - Distance / (sWidth * 1)), 2);
            }
            else if (DistanceStart < DistanceEnd and sqrt(DistanceStart) < (sWidth * 1))
            {
                Fields(i, j, k) = pow((1 - sqrt(DistanceStart) / (sWidth * 1)), 2);
            }
            else if (DistanceStart > DistanceEnd and sqrt(DistanceEnd) < (sWidth * 1))
            {
                Fields(i, j, k) = pow((1 - sqrt(DistanceEnd) / (sWidth * 1)), 2);
            }
        }
    }
    STORAGE_LOOP_END
    Finalize(BC);
    CorrectFractureProfile(BC, 1e-8, 1000);
}

void FractureField::CorrectFractureProfile( BoundaryConditions& BC, double locMobility, size_t nSteps)
{
    auto CalculateFielddot = [&](int i, int j, int k) -> void
    {
        double Eta = sWidth * Grid.dx;
        switch(FractureModel)
        {
            case FractureModels::Well:
            {
                if (Flag(i, j, k) == 2)
                {
                    Fields_dot(i, j, k) =
                        locMobility * ((0.708 * Laplacian(i, j, k) * Eta - 0.5 * (Fields(i, j, k) + 0.5) / Eta));
                }
                break;
            }
            case FractureModels::Obstacle:
            {
                double epsilon = 3 / 16.0 * Eta;
                double K       = 9 / 64.0;
                if (Flag(i, j, k) == 2)
                {
                    Fields_dot(i, j, k)     =  locMobility * ((2 * epsilon * Laplacian(i, j, k) - K / epsilon));
                }
                break;
            }
            default:
            {
                break;
            }
        }
    };
    for (size_t i = 0; i < nSteps; ++i)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
        {
            CalculateFielddot(i, j, k); 
        }
       OMP_PARALLEL_STORAGE_LOOP_END
       MergeIncrements(BC, 0.01);
    }
}

void FractureField::SetSurfaceEnergy(Composition& Cx, Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        double expx            = exp(-delta_g_0 / (PhysicalConstants::R * Tx(i, j, k)));
        double theta           = Cx.MoleFractions(i, j, k, {Comp}) / (expx + Cx.MoleFractions(i, j, k, {Comp}));
        SurfaceEnergy(i, j, k) = InitialSurfaceEnergy * (1.0 - (Xi * theta));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::Solve(PhaseField& Phase,
                          DoubleObstacle& DO,
                          InterfaceProperties& IP,
                          BoundaryConditions& BC,
                          ElasticProperties& EP,
                          double dt)
{
    switch(FractureModel)
    {
        case FractureModels::Obstacle:
        {
            CalculateIncrementsObstacle(Phase, IP, EP, DO);
            break;
        }
        case FractureModels::Well:
        {
            CalculateIncrementsWell(Phase, IP, EP, DO);
            break;
        }
        default:
        {
            break;
        }
    }
    MergeIncrements(BC, dt);
    //CorrectFractureProfile(BC, 1e-8, 1000);
    //SetFracturedElasticConstants(EP);
}

void FractureField::CalculateIncrementsWell(const PhaseField& Phase, InterfaceProperties& IP, ElasticProperties& EP, DoubleObstacle& DO)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if (Flag(i, j, k))
        {
            Fields_dot(i, j, k) =
                Mobility *
                (SurfaceEnergy(i, j, k) * (0.708 * Laplacian(i, j, k) * Eta - 0.5 * (Fields(i, j, k) + 0.5) / Eta) +
                 6.0 * Fields(i, j, k) * (1.0 - Fields(i, j, k)) * DO.PointEnergy(Phase, IP, i, j, k));
        }

        if (Flag(i, j, k) == 2 and Fields(i, j, k) != 1.0)
        {
             Fields_dot(i, j, k) += 2.0 * Mobility * EP.EnergyDensity(Phase,i, j, k) / (1.0 - Fields(i, j, k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::CalculateIncrementsObstacle(const PhaseField& Phase, InterfaceProperties& IP, ElasticProperties& EP, DoubleObstacle& DO)
{
    double Eta     = sWidth * Grid.dx;
    double epsilon = 3 / 16.0 * Eta;
    double K       = 9 / 64.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if (Flag(i, j, k))
        {
            Fields_dot(i, j, k) = Mobility * (SurfaceEnergy(i, j, k) * (2 * epsilon * Laplacian(i, j, k) - K/epsilon));
            TestOutput(i,j,k) = Fields_dot(i,j,k);
        }

        if (Flag(i, j, k) and Fields(i, j, k) != 1.0)
        {
            Fields_dot(i, j, k) += 2.0 * Mobility * EP.EnergyDensity(Phase,i, j, k) / (1.0 - Fields(i, j, k));
            TestOutput2(i,j,k) =   2.0 * Mobility * EP.EnergyDensity(Phase,i, j, k) / (1.0 - Fields(i, j, k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::ApplyAdvection(Settings &locSettings, PhaseField& Phase, ElasticProperties& EP, Orientations& OR,
                                   AdvectionHR &Adv, Velocities &Vel, BoundaryConditions& BC, RunTimeControl &RTC, double dt)
{
    //CalculateCrackTip();
    SetDisplacementsIncrement(Vel, EP, dt);
    locSettings.AdvectAll    (Adv, Vel, Phase, BC, dt, RTC.TimeStep);
    SaveDisplacements        (EP);
    FixAdvectionFractureProfile       (EP);
    Finalize(BC);
}

void FractureField::SetDisplacementsIncrement(Velocities& Vel, ElasticProperties &EP, double dt)
{
    double maxDisplacement = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Displacements, EP.Displacements.Bcells(), reduction(max:maxDisplacement) )
    {

        EP.Displacements(i, j, k) = EP.Displacements(i, j, k) - DisplacementsOLD(i, j, k);
        EP.Displacements(i, j, k) = EP.Displacements(i, j, k) / 10;
        Vel.Average(i, j, k)      = EP.Displacements(i, j, k) * Grid.dx / dt;

        // int NeibhourPoints = 3;
        // if(CrackTip(i, j, k))
        // for (int ii = -NeibhourPoints * Grid.dNx; ii <= NeibhourPoints * Grid.dNx; ++ii)
        // for (int jj = -NeibhourPoints * Grid.dNy; jj <= NeibhourPoints * Grid.dNy; ++jj)
        // for (int kk = -NeibhourPoints * Grid.dNz; kk <= NeibhourPoints * Grid.dNz; ++kk) 
        // {
        //     Vel.Average(i + ii, j + jj, k + kk).set_to_zero();
        // }
        maxDisplacement           = max(EP.Displacements(i, j, k).abs(), maxDisplacement);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &maxDisplacement, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
    if (maxDisplacement > 1.0)
    {
        ConsoleOutput::WriteLine("=");
        string message = "Max displacement : " + to_string(maxDisplacement) +
                         " > 1.0, which is too large for the CFL condition.";
        ConsoleOutput::WriteStandard(thisclassname + "::SetDisplacementsIncrement()", message);
        ConsoleOutput::WriteLine("=");
    }
}

void FractureField::SaveDisplacements(ElasticProperties &EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Displacements, EP.Displacements.Bcells(), )
    {
        DisplacementsOLD(i, j, k) += EP.Displacements(i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::FixAdvectionFractureProfile(ElasticProperties &EP)
{
    auto maxDisplacement = [&](int i, int j, int k)
    {
        double Displacement = 0.0;
        for (int ii = -Grid.dNx; ii <= +Grid.dNx; ii += 2)
        for (int jj = -Grid.dNy; jj <= +Grid.dNy; jj += 2)
        for (int kk = -Grid.dNz; kk <= +Grid.dNz; kk += 2)
        {
            Displacement = max(EP.Displacements(i + ii, j + jj, k + kk).abs(), Displacement);
        }
        return Displacement;
    };
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if (Flag(i, j, k))
        {
            if (Fields(i, j, k) > 0.96)
                Fields(i, j, k) = 1;
            if (Fields(i, j, k) > 0.5)
            {
                Fields(i, j, k) += maxDisplacement(i, j, k);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::SetMobilityPF(PhaseField &Phase, InterfaceProperties &IP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, IP.Properties, 0, )
    {

        if(Phase.Fields(i,j,k).wide_interface())
        {
            for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            if(alpha != beta)
            {
                if(Fields(i,j,k) < PFcutOff)
                {
                    IP.Properties(i, j, k).set_mobility(alpha->index, beta->index, 0.0);
                }
                else
                {
                    double mobility = IP.Properties(i, j, k).get_mobility(alpha->index, beta->index);
                    IP.Properties(i, j, k).set_mobility(alpha->index, beta->index, mobility * Fields(i,j,k));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::SetDrivingForcePF(const PhaseField &Phase, DrivingForce &dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields,0,)
    if(Phase.Fields(i,j,k).wide_interface())
    {
        for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta  = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            if(Fields(i,j,k) < PFcutOff)
                dGab.Force(i,j,k).set_raw(alpha->index, beta->index, 0.0);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void FractureField::SetFracturedElasticConstants(ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(), )
    {
        EP.EffectiveElasticConstants(i, j, k) *= (1.0 - Fields(i, j, k)) * (1.0 - Fields(i, j, k));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::MergeIncrements(const BoundaryConditions& BC,
                                           const double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if (Flag(i, j, k))
        {
           if (Fields_dot(i,j,k) < 0) Fields_dot(i, j, k) = 0;
           if (Fields(i, j, k) != 1.0)
           {
               Fields(i, j, k) += Fields_dot(i, j, k) * dt;
           }
           Fields_dot(i, j, k) = 0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

void FractureField::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        Fields    (i, j, k) = 0;
        Fields_dot(i, j, k) = 0;
        Laplacian (i, j, k) = 0;
        Flag      (i, j, k) = 0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

dVector3 FractureField::Gradient(const int i, const int j, const int k) const
{
    double GWeights[3] = {-0.5 / Grid.dx, 0.0, 0.5 / Grid.dx};

    dVector3 locGradient;
    for (int ii = -1; ii <= +1; ii += 2)
    {
        locGradient[0] += GWeights[ii + 1] * Fields(i + ii, j, k);
        locGradient[1] += GWeights[ii + 1] * Fields(i, j + ii, k);
        locGradient[2] += GWeights[ii + 1] * Fields(i, j, k + ii);
    }
    return locGradient;
}

dVector3 FractureField::Gradient_v2 (const int i, const int j, const int k) const
{
    dVector3 locGradient;
    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        int ii = gs->di;
        int jj = gs->dj;
        int kk = gs->dk;

        double value_x = gs->weightX * Fields(i + ii, j + jj, k + kk);
        double value_y = gs->weightY * Fields(i + ii, j + jj, k + kk);
        double value_z = gs->weightZ * Fields(i + ii, j + jj, k + kk);

        locGradient += (dVector3){value_x, value_y, value_z};
    }
    return locGradient;
}

void FractureField::CalculateCrackTip()
{
//    auto computeCrackTip =
//        [&](int i, int j, int k)
//    {
//        double localCrackTip = 0.0;
//        int value = 0;
//        if (Fields(i, j, k) > 0.99)
//        {
//            for (int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
//            for (int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
//            for (int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
//            {
//                if (ii == 0 && jj == 0 && kk == 0) continue;
//                if (Fields(i + ii, j + jj, k + kk) > 0.99) {
//                    localCrackTip += 1.0;
//                }
//                value++;
//            }
//        }
//            if (Fields(i, j, k) > 0.99 and localCrackTip > 1.0)
//                return 0.0;
//            else if (Fields(i, j, k) > 0.99)
//                return 1.0;
//            return 0.0;
//    };
//
//    const int offset = Fields.Bcells() - 1;
//    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset, )
//    {
//        double locGradient = computeCrackTip(i, j, k);
//        // CrackTip(i, j, k) = pow((locGradient.abs() / maxGradientLocal), 5);
//        //CrackTip(i, j, k) = locGradient;
//    }
//    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::CalculateLaplacians(void)
{
    const int offset = Fields.Bcells() - 1;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset, )
    {
        if (Flag(i, j, k))
        {
            Laplacian(i, j, k) = 0.0;
            for (auto ls = LStencil.begin(); ls != LStencil.end(); ls++)
            {
                int ii = i + ls->di;
                int jj = j + ls->dj;
                int kk = k + ls->dk;

                Laplacian(i, j, k) += ls->weight * Fields(ii, jj, kk);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::SetFlags(void)
{
    const int offset = Fields.Bcells() - 1;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset,)
    {
        if (Flag(i, j, k) == 2)
        {
            for(int ii = -Grid.dNx; ii <= +Grid.dNx; ++ii)
            for(int jj = -Grid.dNy; jj <= +Grid.dNy; ++jj)
            for(int kk = -Grid.dNz; kk <= +Grid.dNz; ++kk)
            if (!(Flag(i + ii, j + jj, k + kk)))
            {
                Flag(i + ii, j + jj, k + kk) = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FractureField::Finalize(const BoundaryConditions& BC)
{
    const int offset = Fields.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset, )
    {
        if (Flag(i, j, k))
        {
            if (Fields(i, j, k) < 1.0e-10)
                Fields(i, j, k) = 0;
            if (Fields(i, j, k) > 1)
                Fields(i, j, k) = 1;

            if (Fields(i, j, k) > 0)
                Flag(i, j, k) = 2;
            else
                Flag(i, j, k) = 0;

            Fields_dot(i, j, k) = 0;
            Laplacian(i, j, k)  = 0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
    SetFlags();
    CalculateLaplacians();
}

void FractureField::PrintVolumeFractions(void)
{
    double VolumeFraction = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, reduction(+:VolumeFraction))
    {
        if (Fields(i, j, k) > 0.0)
        {
            VolumeFraction += Fields(i, j, k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    VolumeFraction /= double(Grid.TotalNumberOfCells());

    ConsoleOutput::WriteStandard("Fracture Volume Fraction [%]", VolumeFraction * 100.0);
}

void FractureField::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Fields);
    BC.SetY(Fields);
    BC.SetZ(Fields);

    BC.SetX(Fields_dot);
    BC.SetY(Fields_dot);
    BC.SetZ(Fields_dot);

    BC.SetX(Laplacian);
    BC.SetY(Laplacian);
    BC.SetZ(Laplacian);

    BC.SetX(Flag);
    BC.SetY(Flag);
    BC.SetZ(Flag);
}

bool FractureField::Write(const Settings& locSettings, const int tStep) const
{
    string FileName =
        FileInterface::MakeFileName(locSettings.RawDataDir,"FractureField_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened",
                thisclassname, "Write()");
        return false;
    };

    int tmp = Grid.Nx;
    out.write(reinterpret_cast<const char*>(&tmp), sizeof(int));
    tmp = Grid.Ny;
    out.write(reinterpret_cast<const char*>(&tmp), sizeof(int));
    tmp = Grid.Nz;
    out.write(reinterpret_cast<const char*>(&tmp), sizeof(int));

    STORAGE_LOOP_BEGIN(i, j, k, Fields, 0)
    {
        double val = Fields(i, j, k);
        out.write(reinterpret_cast<const char*>(&val), sizeof(double));
    }
    STORAGE_LOOP_END

    out.close();
    return true;
}

bool FractureField::Read(const Settings& locSettings,
                         const BoundaryConditions& BC, const int tStep)
{
    string FileName =
        FileInterface::MakeFileName(locSettings.InputRawDataDir,"FractureField_", tStep, ".dat");

    bool read_success = Read(FileName);
    Finalize(BC);
    return read_success;
}

bool FractureField::Read(string FileName)
{
    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteWarning(FileName + " could not be opened",
                thisclassname, "Read()");
        return false;
    };

    int locNx = Grid.Nx;
    int locNy = Grid.Ny;
    int locNz = Grid.Nz;
    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx
                << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Grid.Nx
                << ", " << Grid.Ny << ", " << Grid.Nz << ") grid points.\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }

    STORAGE_LOOP_BEGIN(i, j, k, Fields, 0)
    {
        double val;
        inp.read(reinterpret_cast<char*>(&val), sizeof(double));            // Field(i,j,k)->value
        Fields(i, j, k) = val;
    }
    STORAGE_LOOP_END

    inp.close();

    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

void FractureField::Remesh(int newNx, int newNy, int newNz,
                           const BoundaryConditions& BC)
{
    Grid.Nx = newNx;
    Grid.Ny = newNy;
    Grid.Nz = newNz;

    Fields.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        Flag(i, j, k) = 2 * (Fields(i, j, k) > 0);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Finalize(BC);
        
    RefVolume = Pi;

    if(Grid.Nx < sWidth) RefVolume *= (Grid.Nx+1)/2.0;
    else RefVolume *= 1.1*sWidth;
    if(Grid.Ny < sWidth) RefVolume *= (Grid.Ny+1)/2.0;
    else RefVolume *= 1.1*sWidth;
    if(Grid.Nz < sWidth) RefVolume *= (Grid.Nz+1)/2.0;
    else RefVolume *= 1.1*sWidth;

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

FractureField& FractureField::operator= (const FractureField& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "FractureField")
    {
        thisclassname  = rhs.thisclassname;
        thisobjectname = rhs.thisobjectname;

        Grid = rhs.Grid;

        sWidth    = rhs.sWidth;
        RefVolume = rhs.RefVolume;
        Mobility  = rhs.Mobility;

        if (Fields.IsNotAllocated())
        {
            Fields    .Allocate(Grid, rhs.Fields.Bcells());
            Fields_dot.Allocate(Grid, rhs.Fields.Bcells());
            Laplacian .Allocate(Grid, rhs.Fields.Bcells());
            Flag      .Allocate(Grid, rhs.Fields.Bcells());
            SurfaceEnergy.Allocate(Grid, rhs.Fields.Bcells());
        }
        else if (not Fields.IsSize(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz))
        {
            Fields    .Reallocate(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz);
            Fields_dot.Reallocate(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz);
            Laplacian .Reallocate(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz);
            Flag      .Reallocate(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz);
            SurfaceEnergy.Reallocate(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
        {
            Fields    (i, j, k) = rhs.Fields(i, j, k);
            Fields_dot(i, j, k) = rhs.Fields_dot(i, j, k);
            Laplacian (i, j, k) = rhs.Laplacian(i, j, k);
            Flag      (i, j, k) = rhs.Flag(i, j, k);
            SurfaceEnergy(i, j, k) = rhs.SurfaceEnergy(i, j, k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

double FractureField::PointEnergy(const InterfaceProperties& IP, int i, int j, int k)
{
    const double Eta = sWidth * Grid.dx;
    double energy = 0;
    double Prefactor = Eta*Eta;

    if (Flag(i, j, k))
    {
        dVector3 locGradient = Gradient(i, j, k);
        energy += 0.5 * SurfaceEnergy(i, j, k) / Eta *
                  (Fields(i, j, k) * Fields(i, j, k) + Prefactor * (locGradient * locGradient));
    }
    return energy;
}

double FractureField::Energy(const InterfaceProperties& IP)
{
    double Energy = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, reduction(+ : Energy))
    {
        Energy += PointEnergy(IP, i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Energy*Grid.CellVolume();
}

double FractureField::AverageEnergyDensity(const InterfaceProperties& IP)
{
    double Energy = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, reduction(+ : Energy))
    {
        Energy += PointEnergy(IP, i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Energy/(Grid.LocalNumberOfCells());
}

void FractureField::WriteVTK(const Settings& locSettings, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + "_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t){"Interface",     [this](int i, int j, int k) { return Interface  (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"Flag",          [this](int i, int j, int k) { return Flag       (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"FractureField", [this](int i, int j, int k) { return Fields     (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"Laplacian",     [this](int i, int j, int k) { return Laplacian  (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"SurfaceEnergy", [this](int i, int j, int k) { return SurfaceEnergy(i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"TestOutPut",    [this](int i, int j, int k) { return TestOutput(i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"TestOutPut2",    [this](int i, int j, int k) { return TestOutput2(i, j, k);}});
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void FractureField::WriteDistortedVTK(
        const Settings& locSettings,
        const ElasticProperties& EP,
        const int tStep,
        const int precision) const
{
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + "Distorted_", tStep, ".vts");

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Interfaces",    [this](int i, int j, int k) { return Interface(i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"Flags",         [this](int i, int j, int k) { return Flag     (i, j, k);}});
    ListOfFields.push_back((VTK::Field_t){"FractureField", [this](int i, int j, int k) { return Fields   (i, j, k);}});
    VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields, precision);
}

void FractureField::WriteEnergyVTK(const Settings& locSettings, const int tStep,
                                   const InterfaceProperties& IP, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname + "Energy_", tStep, ".vts");
    ListOfFields.push_back((VTK::Field_t) { "FractureEnergyDensity", [&IP,this](int i, int j, int k)
    { 
        return PointEnergy(IP, i, j, k); 
    }}); 
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

}
