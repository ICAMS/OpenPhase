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
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#include "Includes.h"
#include "FluidDynamics/InteractionSolidFluid.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Velocities.h"
#include "Tools.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace openphase
{

void InteractionSolidFluid::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
     thisclassname  = "InteractionSolidFluid";
     thisobjectname = thisclassname + ObjectNameSuffix;
}

void InteractionSolidFluid::ReadInput(std::string InputFileName)
{
     ConsoleOutput::WriteLineInsert(thisclassname+"input");
     ConsoleOutput::WriteStandard("Source", InputFileName);

     std::fstream inp(InputFileName.c_str(), std::ios::in);

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

void InteractionSolidFluid::ReadInput(std::stringstream& inp)
{
     int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

     DoInertia                         = FileInterface::ReadParameterB(inp, moduleLocation, "Inertia");
     DoRotations                       = FileInterface::ReadParameterB(inp, moduleLocation, "Rotations");
     DoLimitVelocity                   = FileInterface::ReadParameterB(inp, moduleLocation, "LimitVelocity", false, false);
     DoEnforceZeroTotalMomentum        = FileInterface::ReadParameterB(inp, moduleLocation, "EnforceZeroTotalMomentum", false, false);
     DoEnforceZeroTotalAngularMomentum = FileInterface::ReadParameterB(inp, moduleLocation, "DoEnforceZeroTotalAngularMomentum", false, false);
     TauTranslation                    = FileInterface::ReadParameterD(inp, moduleLocation, "TauTrans",!DoInertia,0.0);
     TauRotation                       = FileInterface::ReadParameterD(inp, moduleLocation, "TauRot",!DoInertia,0.0);
     VelocityLimit                     = FileInterface::ReadParameterD(inp, moduleLocation, "VelocityLimit", DoLimitVelocity,0.0);
}

void InteractionSolidFluid::CalculateCenterOfMass(PhaseField& Phase,
        const BoundaryConditions& BC)
{
    int Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif

    // The old center of mass is used to transform the system coordinates
    // into the old center of mass. The new center of mass is calculated and
    // the system coordinates are transformed back. This solves the problem
    // of the periodic boundary conditions.
    bool setInitialRcm = false;
    std::vector<dVector3> oldRcm(Phase.FieldsProperties.size());
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        oldRcm[idx] = Phase.FieldsProperties[idx].Rcm;
        Phase.FieldsProperties[idx].Rcm.set_to_zero();
        if (applicable(Phase.FieldsProperties[idx]))
        if (oldRcm[idx].abs() < std::numeric_limits<double>::epsilon())
        {
            setInitialRcm = true;
        }
    }

    // If no old or initial center of mass has been set, an arbitrary point in
    // the grain is searched and set as the center of old center of mass.
    if (setInitialRcm)
    {
        std::vector<std::vector<dVector3>> initRcm(Nthreads);
        for (size_t t = 0; t < initRcm.size(); t++)
        {
            initRcm[t].resize(Phase.FieldsProperties.size());
            for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            {
                initRcm[t][idx] = oldRcm[idx];
            }
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
            if (initRcm[omp_get_thread_num()][alpha->index].abs() < std::numeric_limits<double>::epsilon())
            {
                initRcm[omp_get_thread_num()][alpha->index] = Tools::Position(dVector3({double(i),double(j),double(k)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // OMP reduction
        for (size_t idx = 0; idx < oldRcm.size(); idx++)
        {
            for (size_t t = 0; t < initRcm.size(); t++)
            if (initRcm[t][idx].abs() > oldRcm[idx].abs())
            {
                oldRcm[idx] = initRcm[t][idx];
                break;
            }
        }

        // MPI reduction
        #ifdef MPI_PARALLEL
        for (size_t idx = 0; idx < oldRcm.size(); idx++)
        {
            struct {double value; int rank;} in, out;
            in.value = oldRcm[idx].abs();
            in.rank = MPI_RANK;
            OP_MPI_Allreduce(&in, &out, 1, OP_MPI_DOUBLE_INT, OP_MPI_MAXLOC, OP_MPI_COMM_WORLD);
            OP_MPI_Bcast(oldRcm[idx].data(), 3, OP_MPI_DOUBLE, out.rank, OP_MPI_COMM_WORLD);
        }
        #endif
    }
    
    // Calculate new center of mass
    std::vector<std::vector<dVector3>> localRcm(Nthreads);
    for (size_t t = 0; t < localRcm.size(); t++)
    {
        localRcm[t].resize(Phase.FieldsProperties.size());
        for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
        {
            localRcm[t][idx].set_to_zero();
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        if (applicable(Phase.FieldsProperties[alpha->index]))
        {
            const dVector3 position = Tools::Position(dVector3({double(i),double(j),double(k)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
            const dVector3 distance = Tools::Distance(position, oldRcm[alpha->index], Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC);
            localRcm[omp_get_thread_num()][alpha->index] += distance*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (applicable(Phase.FieldsProperties[idx]))
    {
        //OMP Reduction
        for (auto it : localRcm) Phase.FieldsProperties[idx].Rcm += it[idx];

        #ifdef MPI_PARALLEL
        auto localMPI_Rcm = Phase.FieldsProperties[idx].Rcm;
        OP_MPI_Allreduce(&localMPI_Rcm[0], &(Phase.FieldsProperties[idx].Rcm[0]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&localMPI_Rcm[1], &(Phase.FieldsProperties[idx].Rcm[1]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&localMPI_Rcm[2], &(Phase.FieldsProperties[idx].Rcm[2]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        Phase.FieldsProperties[idx].Rcm *= 1.0/double(Phase.FieldsProperties[idx].Volume);
        Phase.FieldsProperties[idx].Rcm += oldRcm[idx];
        Tools::MapIntoSimulationDomain(Phase.FieldsProperties[idx].Rcm, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC);
    }
}

void InteractionSolidFluid::CalculateMomentOfInterita(PhaseField& Phase,
        const BoundaryConditions& BC)
{
    for(unsigned int idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        Phase.FieldsProperties[idx].InertiaM.set_to_zero();
    }

    int Nthreads = 1;
    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif

    std::vector<std::vector<dMatrix3x3>> locInertiaMs(Nthreads);
    for (size_t t = 0; t < size_t(Nthreads); t++)
    {
        locInertiaMs[t].resize(Phase.FieldsProperties.size());
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        const double dV = Phase.Grid.CellVolume(true);
        const dVector3 pos = Tools::Position(dVector3({double(i),double(j),double(k)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
        for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        if (applicable(Phase.FieldsProperties[alpha->index]))
        {
            const Grain& grain = Phase.FieldsProperties[alpha->index];
            const dVector3 r_Rcm = Tools::Distance(pos, Phase.FieldsProperties[alpha->index].Rcm, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Phase.Grid.dx;

            locInertiaMs[omp_get_thread_num()][alpha->index](0,0) += alpha->value*grain.Density*(r_Rcm[1]*r_Rcm[1] + r_Rcm[2]*r_Rcm[2])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](1,1) += alpha->value*grain.Density*(r_Rcm[0]*r_Rcm[0] + r_Rcm[2]*r_Rcm[2])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](2,2) += alpha->value*grain.Density*(r_Rcm[0]*r_Rcm[0] + r_Rcm[1]*r_Rcm[1])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](0,1) -= alpha->value*grain.Density*(r_Rcm[0]*r_Rcm[1])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](1,0) -= alpha->value*grain.Density*(r_Rcm[1]*r_Rcm[0])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](0,2) -= alpha->value*grain.Density*(r_Rcm[0]*r_Rcm[2])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](2,0) -= alpha->value*grain.Density*(r_Rcm[2]*r_Rcm[0])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](1,2) -= alpha->value*grain.Density*(r_Rcm[1]*r_Rcm[2])*dV;
            locInertiaMs[omp_get_thread_num()][alpha->index](2,1) -= alpha->value*grain.Density*(r_Rcm[2]*r_Rcm[1])*dV;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if (applicable(Phase.FieldsProperties[idx]))
    {
        for (size_t t = 0; t < size_t(Nthreads); t++)
        {
            Phase.FieldsProperties[idx].InertiaM += locInertiaMs[t][idx];
        }

        #ifdef MPI_PARALLEL
        auto locInertiaM = Phase.FieldsProperties[idx].InertiaM;
        OP_MPI_Allreduce(&locInertiaM(0,0), &(Phase.FieldsProperties[idx].InertiaM(0,0)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(1,1), &(Phase.FieldsProperties[idx].InertiaM(1,1)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(2,2), &(Phase.FieldsProperties[idx].InertiaM(2,2)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(0,1), &(Phase.FieldsProperties[idx].InertiaM(0,1)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(1,0), &(Phase.FieldsProperties[idx].InertiaM(1,0)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(0,2), &(Phase.FieldsProperties[idx].InertiaM(0,2)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(2,0), &(Phase.FieldsProperties[idx].InertiaM(2,0)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(1,2), &(Phase.FieldsProperties[idx].InertiaM(1,2)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&locInertiaM(2,1), &(Phase.FieldsProperties[idx].InertiaM(2,1)), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif
    }
}

void InteractionSolidFluid::CalculateRigidBodyVelocities(PhaseField& Phase,
        const BoundaryConditions& BC, const double dt) const
{
    const double dV = Phase.Grid.CellVolume(true);
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        Grain& grain = Phase.FieldsProperties[idx];
        if (applicable(grain))
        {
            double GrainMass_1 = 1.0/(grain.Density*grain.Volume*dV);
            grain.Acm  += grain.Force * GrainMass_1;
            grain.Vcm  += grain.Acm * dt;

            if (DoRotations)
            {
                grain.aAcc += grain.InertiaM.inverted() * grain.Torque;
                grain.aVel += grain.aAcc * dt;
                EulerAngles locAngles({grain.aVel[0] * dt,
                                       grain.aVel[1] * dt,
                                       grain.aVel[2] * dt}, XYZ);
                grain.Orientation += locAngles.getQuaternion();
            }
        }
        else
        {
            grain.Acm  = {0,0,0};
            grain.Vcm  = {0,0,0};
            grain.aAcc = {0,0,0};
            grain.aVel = {0,0,0};
        }
    }
}

void InteractionSolidFluid::CalculateRigidBodyVelocitiesWithoutInertia(
        PhaseField& Phase, const BoundaryConditions& BC, const double dt) const
{
    const double dV = Phase.Grid.CellVolume(true);

    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        Grain& grain = Phase.FieldsProperties[idx];
        if (applicable(grain))
        {
            double GrainMass = grain.Density*grain.Volume*dV;
            double MobilityTranslation = TauTranslation/GrainMass;
            grain.Vcm = grain.Force * MobilityTranslation;
            if (DoRotations)
            {
                dMatrix3x3 MobilityRotion = grain.InertiaM.inverted()*TauRotation;
                grain.aVel = MobilityRotion * grain.Torque;
                EulerAngles locAngles({grain.aVel[0] * dt,
                                       grain.aVel[1] * dt,
                                       grain.aVel[2] * dt}, XYZ);
                grain.Orientation += locAngles.getQuaternion();
            }
        }
        else
        {
            grain.Vcm  = {0,0,0};
            grain.aVel = {0,0,0};
        }
    }
}

void InteractionSolidFluid::EnforceZeroTotalMomentum(PhaseField& Phase)
{
    dVector3 Momentum = {0.0,0.0,0.0};
    double   Mass     = 0.0;

    const double dx = Phase.Grid.dx;
    const double dV = dx*dx;
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
       Grain& grain = Phase.FieldsProperties[idx];
       if (applicable(grain))
       {
           const double GrainMass = grain.Density*grain.Volume*dV;
           Mass     += GrainMass;
           Momentum += grain.Vcm * GrainMass ;
       }
    }

    if (Mass > std::numeric_limits<double>::epsilon())
    {
        const dVector3 VcmFix = Momentum/Mass;
        for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
        {
            Grain& grain = Phase.FieldsProperties[idx];
            if (applicable(grain))
            {
                grain.Vcm -= VcmFix;
            }
        }
    }
}

void InteractionSolidFluid::EnforceZeroTotalAngularMomentum(PhaseField& Phase)
{
    // TODO consider Avel
    dVector3 AMomentum = {0.0,0.0,0.0};
    double   Mass      = 0.0;
    double   Radius    = 0.0;

    const double dx = Phase.Grid.dx;
    const double dV = dx*dx;
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        Grain& grain = Phase.FieldsProperties[idx];
        if (applicable(grain))
        {
            const double   GrainMass = grain.Density*grain.Volume*dV;
            const dVector3 Momentum  = grain.Vcm * GrainMass;
            Mass      += GrainMass;
            AMomentum += grain.Rcm.cross(Momentum);
            Radius    += grain.Rcm.abs();
        }
    }

    if (Radius*Mass     > std::numeric_limits<double>::epsilon())
    if (AMomentum.abs() > std::numeric_limits<double>::epsilon())
    {
        const double   Fix  = AMomentum.abs()/(Radius*Mass);
        const dVector3 e_AM = AMomentum/AMomentum.abs();

        for (size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
        {
            Grain& grain = Phase.FieldsProperties[idx];
            if (applicable(grain))
            {
                const dVector3 e_R   = grain.Rcm/grain.Rcm.abs();
                const dVector3 e_Fix = e_AM.cross(e_R);

                Phase.FieldsProperties[idx].Vcm -= e_Fix * Fix;
            }
        }
    }
}

void InteractionSolidFluid::LimitVelocity(PhaseField& Phase, double VLimit)
{
    int VelocityLimitApplied = 0;
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    {
        Grain& grain = Phase.FieldsProperties[idx];
        if (applicable(grain))
        if (grain.Vcm.abs() > VLimit)
        {
            VelocityLimitApplied++;
            grain.Vcm /= grain.Vcm.abs();
            grain.Vcm *= VLimit;
        }
    }

    if (VelocityLimitApplied)
    {
        std::stringstream message;
        message << "Solid velocity limit as applied!";
        ConsoleOutput::WriteWarning(message.str(), "InteractionSolidFluid", "LimitVelocity");
    }
}

double InteractionSolidFluid::LocalSolidFraction(long i, long j, long k, const PhaseField& Phase)
{
    double SolidFraction = 0.0;
    for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
    {
        if (applicable(Phase.FieldsProperties[alpha->index]))
        if (alpha->value > std::numeric_limits<double>::epsilon())
        {
            SolidFraction += alpha->value;
        }
    }
    return SolidFraction;
}

dVector3 InteractionSolidFluid::LocalGrainVelocity(long i, long j, long k,
        const Grain& grain, const BoundaryConditions& BC,
        const PhaseField& Phase)
{
    const dVector3 pos = Tools::Position(dVector3({double(i),double(j),double(k)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
    const dVector3 locR = Tools::Distance(pos, grain.Rcm, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC) * Phase.Grid.dx;
    return grain.Vcm + grain.aVel.cross(locR);
}

void InteractionSolidFluid::CalculateSolidPhaseVelocities(Velocities& Vel,
        const PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        double SolidFraction = LocalSolidFraction(i,j,k,Phase);
        for(size_t pIdx = 0; pIdx != Vel.Nphases; pIdx++)
        {
            Vel.Phase(i,j,k,{pIdx}).set_to_zero();
        }
        for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); alpha++)
        {
            const Grain& grain = Phase.FieldsProperties[alpha->index];
            if (applicable(grain))
            if (alpha->value > std::numeric_limits<double>::epsilon())
            {
                assert(SolidFraction != 0.0);
                double weight = alpha->value/SolidFraction;
                Vel.Phase(i,j,k,{grain.Phase}) += LocalGrainVelocity(i,j,k,grain,BC,Phase)*weight;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Vel.SetBoundaryConditions(BC);
    Vel.CalculateAverage(Phase);
}

void InteractionSolidFluid::CalculateSolidVelocities(PhaseField& Phase,
       Velocities& Vel, const BoundaryConditions& BC, const double dt) const
{
    CalculateCenterOfMass(Phase, BC);
    if (DoRotations) CalculateMomentOfInterita(Phase, BC);

    if (DoInertia) CalculateRigidBodyVelocities(Phase, BC, dt);
    else           CalculateRigidBodyVelocitiesWithoutInertia(Phase, BC, dt);

    if (DoLimitVelocity) LimitVelocity(Phase, VelocityLimit);
    if (DoEnforceZeroTotalMomentum) EnforceZeroTotalMomentum(Phase);
    if (DoEnforceZeroTotalAngularMomentum) EnforceZeroTotalAngularMomentum(Phase);

    CalculateSolidPhaseVelocities(Vel,Phase,BC);
}

} //namespace openphase
