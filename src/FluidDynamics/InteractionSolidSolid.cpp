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

 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev;
 *                         Marvin Tegeler; Raphael Schiedung
 *
 */

#include "FluidDynamics/InteractionSolidSolid.h"
#include "Tools.h"

namespace openphase
{

void InteractionSolidSolid::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
     thisclassname  = "InteractionSolidSolid";
     thisobjectname = thisclassname + ObjectNameSuffix;

     Grid = locSettings.Grid;

     order    = 0;
     cutoff   = 0;
     strength = 0.0;
     strength_wang = 0.0;

     Model = SolidSolidInteractionModel::Standard;
     Nphases = locSettings.Nphases;
     PhaseAggregateStates = locSettings.PhaseAggregateStates;

     elastic.resize(Nphases);
     for (size_t n = 0; n < Nphases; ++n)
     {
          elastic[n].resize(Nphases);
          std::fill(elastic[n].begin(),elastic[n].end(),0.0);
     }
}

void InteractionSolidSolid::ReadInput(std::string InputFileName)
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

void InteractionSolidSolid::ReadInput(std::stringstream& inp)
{
     int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

     std::string tmp  = FileInterface::ReadParameterK(inp, moduleLocation, "Model");

     if      (tmp == "STANDARD") Model = SolidSolidInteractionModel::Standard;
     else if (tmp == "WANG")     Model = SolidSolidInteractionModel::Wang;
     else
     {
         ConsoleOutput::WriteExit("Nonexistent solid solid interaction model, "+tmp+", selected", thisclassname, "ReadInput()");
         OP_Exit(EXIT_FAILURE);
     }

    switch (Model)
    {
        case SolidSolidInteractionModel::Wang:
            strength_wang = FileInterface::ReadParameterD(inp, moduleLocation, "SWang");
            break;

        default:
            order    = FileInterface::ReadParameterI(inp, moduleLocation, "Order",  false, 4);
            cutoff   = FileInterface::ReadParameterI(inp, moduleLocation, "Cutoff", false, 3);
            strength = FileInterface::ReadParameterD(inp, moduleLocation, "Strength");

            for (size_t n = 0; n < Nphases; ++n)
            if (PhaseAggregateStates[n] == AggregateStates::Solid)
            for (size_t m = n; m < Nphases; ++m)
            if (PhaseAggregateStates[m] == AggregateStates::Solid)
            {
                std::stringstream TEMPCstr;
                TEMPCstr << "Elastic[" << n << "][" << m << "]";
                elastic[n][m] = FileInterface::ReadParameterD(inp, moduleLocation, TEMPCstr.str());
            }
    }
}

void InteractionSolidSolid::CalculateLocal(
        const int i, const int j, const int k,
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const double dt) const
{
    const double dV = Grid.CellVolume(true);
    //const NodeAB<dVector3,dVector3> locNormal = Phase.Normals(i,j,k);

    for (auto alpha  = Phase.Fields(i,j,k).cbegin();
              alpha != Phase.Fields(i,j,k).cend(); ++alpha)
    if (applicable(Phase.FieldsProperties[alpha->index]))
    {
        const size_t pidxA = Phase.FieldsProperties[alpha->index].Phase;
        const dVector3& pos_cm_alpha = Phase.FieldsProperties[alpha->index].Rcm;

        for (int ii = -cutoff*Grid.dNx; ii <= cutoff*Grid.dNx; ++ii)
        for (int jj = -cutoff*Grid.dNy; jj <= cutoff*Grid.dNy; ++jj)
        for (int kk = -cutoff*Grid.dNz; kk <= cutoff*Grid.dNz; ++kk)
        if (ii*ii+jj*jj+kk*kk > 0)
        if (ii*ii+jj*jj+kk*kk < cutoff*cutoff)
        {
            // Calculate position of collision event
            const dVector3 loc_pos_c = {i+0.5*ii,j+0.5*jj,k+0.5*kk};
            const dVector3 pos_c = Tools::Position(loc_pos_c, Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);

            for(auto beta  = Phase.Fields(i+ii,j+jj,k+kk).cbegin();
                     beta != Phase.Fields(i+ii,j+jj,k+kk).cend(); ++beta)
            if (applicable(Phase.FieldsProperties[beta->index]))
            if(alpha->index != beta->index)
            {
                const double r = sqrt(ii*ii+jj*jj+kk*kk); // Distance between (i,j,k) and (i+ii,j+jj,k+kk).
                const size_t pidxB = Phase.FieldsProperties[beta->index].Phase;
                const dVector3& pos_cm_beta = Phase.FieldsProperties[beta->index].Rcm;

                // Calculate collision velocity v_c projected on connection between (i,j,k) and (i+ii,j+jj,k+kk).
                const dVector3 r_c_alpha = Tools::Distance(pos_c, pos_cm_alpha, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
                const dVector3 r_c_beta  = Tools::Distance(pos_c, pos_cm_beta,  Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
                const dVector3 v_alpha = Phase.FieldsProperties[alpha->index].Vcm + Phase.FieldsProperties[alpha->index].aVel.cross(r_c_alpha);
                const dVector3 v_beta  = Phase.FieldsProperties[beta ->index].Vcm + Phase.FieldsProperties[beta ->index].aVel.cross(r_c_beta );
                const dVector3 normal = {ii/r,jj/r,kk/r};
                const double v_c = normal*(v_alpha-v_beta);

                // Calculate force density
                const double Vol1 = 4.0/3.0*Pi*alpha->value*alpha->value*alpha->value;
                const double Vol2 = 4.0/3.0*Pi*beta ->value*beta ->value*beta ->value;
                const double v0   = Grid.dx/dt; // Reference velocity to get physical dimensions right
                double dEnergyDensity_dr = 0.5*strength*Vol1*Vol2*order*std::pow((r-cutoff)/cutoff,order-1)/cutoff;
                dEnergyDensity_dr = elastic[pidxA][pidxB]*dEnergyDensity_dr + (1.0-elastic[pidxA][pidxB])*v_c/v0*dEnergyDensity_dr;
                const dVector3 force_density = {dEnergyDensity_dr*ii/r,dEnergyDensity_dr*jj/r,dEnergyDensity_dr*kk/r};

                // Apply force density
                const dVector3 pos_alpha = Tools::Position(dVector3({double(i),    double(j),    double(k)}),    Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
                const dVector3 pos_beta  = Tools::Position(dVector3({double(i+ii), double(j+jj), double(k+kk)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);
                const dVector3 r_alpha = Tools::Distance(pos_alpha, pos_cm_alpha, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
                const dVector3 r_beta  = Tools::Distance(pos_beta,  pos_cm_beta,  Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    Phase.FieldsProperties[alpha->index].Force  += force_density*dV;
                    Phase.FieldsProperties[beta ->index].Force  -= force_density*dV;
                    Phase.FieldsProperties[alpha->index].Torque += r_alpha.cross(force_density)*dV;
                    Phase.FieldsProperties[beta ->index].Torque -= r_beta .cross(force_density)*dV;
                }
            }
        }
    }
}

void InteractionSolidSolid::CalculateLocalWang(
        const int i, const int j, const int k,
        PhaseField& Phase,
        const BoundaryConditions& BC,
        const std::function<double(int,int,int)>& MassDensity) const
{
    const double dV = Grid.CellVolume(true);
    //The calculation of the force density is based on, Wang, Y. U. (2006).
    //Acta materialia, 54(4), 953-961.

    const dVector3 pos = Tools::Position(dVector3({double(i), double(j), double(k)}), Phase.Grid.OffsetX, Phase.Grid.OffsetY, Phase.Grid.OffsetZ);

    double locMassDensity = MassDensity(i,j,k); // actual local mass
    double locMassDensity_Uncompressed = 0.0;   // local equilibrium mass density without pressure
    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
    {
        locMassDensity_Uncompressed += Phase.FieldsProperties[alpha->index].Density*alpha->value;
    }

    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha != Phase.Fields(i,j,k).cend(); ++alpha)
    if (applicable(Phase.FieldsProperties[alpha->index]))
    {
        const dVector3 gradient_alpha = Phase.Fields(i,j,k).get_gradient(alpha->index);
        const dVector3 pos_cm_alpha   = Phase.FieldsProperties[alpha->index].Rcm;
        const dVector3 r_alpha        = Tools::Distance(pos, pos_cm_alpha, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
        for (auto beta = alpha + 1; beta != Phase.Fields(i,j,k).cend(); ++beta)
        if (applicable(Phase.FieldsProperties[beta->index]))
        {
            const dVector3 gradient_beta  = Phase.Fields(i,j,k).get_gradient(beta->index);
            const dVector3 ForceDensity   = (gradient_alpha-gradient_beta)*(locMassDensity-locMassDensity_Uncompressed)*strength_wang;
            const dVector3 pos_cm_beta    = Phase.FieldsProperties[beta->index].Rcm;
            const dVector3 r_beta         = Tools::Distance(pos, pos_cm_beta, Phase.Grid.TotalNx, Phase.Grid.TotalNy, Phase.Grid.TotalNz, BC)*Grid.dx;
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                Phase.FieldsProperties[alpha->index].Force  += ForceDensity*dV;
                Phase.FieldsProperties[beta ->index].Force  -= ForceDensity*dV;
                Phase.FieldsProperties[alpha->index].Torque += r_alpha.cross(ForceDensity)*dV;
                Phase.FieldsProperties[beta ->index].Torque -= r_beta .cross(ForceDensity)*dV;
            }
        }
    }
}

void InteractionSolidSolid::Calculate(PhaseField& Phase,
        const BoundaryConditions& BC,
        const std::function<double(int,int,int)>& MassDensity,
        const double dt) const
{
    switch (Model)
    {
        case SolidSolidInteractionModel::Wang:
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
            {
                CalculateLocalWang(i,j,k,Phase,BC,MassDensity);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;

        default:
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
            {
                CalculateLocal(i,j,k,Phase,BC,dt);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
    }
}
} //namespace openphase
