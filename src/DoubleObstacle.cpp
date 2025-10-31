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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Johannes Goerler; Raphael Schiedung; Stephan Hubig
 *
 */

#include "DoubleObstacle.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "InterfaceRegularization.h"
#include "PhaseField.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

void DoubleObstacle::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "DoubleObstacle";
    thisobjectname = thisclassname + ObjectNameSuffix;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSharp(PhaseField& Phase, InterfaceProperties& IP)
{
    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSharpSR(Phase, IP);
            break;
        }
        case Resolutions::Dual:
        {
            CalculatePhaseFieldIncrementsSharpDR(Phase, IP);
            break;
        }
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSharpSR(PhaseField& Phase,
                                                          InterfaceProperties& IP)
{
    dVector3 n_vector = {1.0, 0.0, 0.0};
    n_vector.normalize();

    double dx = Phase.Grid.dx;
    vector<PotentialCorrections> StencilDirections;

    for(auto ds = Phase.LStencil.cbegin(); ds != Phase.LStencil.cend(); ds++)
    {
        double d_x = ds->di;
        double d_y = ds->dj;
        double d_z = ds->dk;

        PotentialCorrections dirNew;

        dirNew.d_x = d_x;
        dirNew.d_y = d_y;
        dirNew.d_z = d_z;

        dirNew.stencil_weight = ds->weight;

        dirNew.scal_prod = (d_x*n_vector[0] + d_y*n_vector[1] + d_z*n_vector[2])*dx;

        dirNew.alpha     = cos(Pi*(dirNew.scal_prod)/Phase.Grid.Eta);
        dirNew.beta      = sin(Pi*(dirNew.scal_prod)/Phase.Grid.Eta);
        dirNew.phi_tl    = 0.5*cos(Pi*(1.0 - (dirNew.scal_prod)/Phase.Grid.Eta)) + 0.5;
        dirNew.phi_tu    = 0.5*cos(Pi*(dirNew.scal_prod)/Phase.Grid.Eta) + 0.5;

        StencilDirections.push_back(dirNew);
    }

    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.Fields(i,j,k);

            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double dPhi_dt = 0.0;

                double pot_term_a = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                {
                    double fik_term = 0;
                    if(StDir->scal_prod < 0 && alpha->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && alpha->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(alpha->value - 0.5) - StDir->beta*sqrt(alpha->value*(1.0 - alpha->value));
                    }
                    pot_term_a += -StDir->stencil_weight*(fik_term - alpha->value);
                }

                double pot_term_b = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend() ; StDir++)
                {
                    double fik_term = 0.0;
                    if(StDir->scal_prod < 0 && beta->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && beta->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(beta->value - 0.5) - StDir->beta*sqrt(beta->value*(1.0 - beta->value));
                    }
                    pot_term_b += -StDir->stencil_weight*(fik_term - beta->value);
                }

                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio*
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                dPhi_dt = IP.Properties(i,j,k).get_energy(alpha->index, beta->index)
                         *((alpha->laplacian - beta->laplacian) + 0.5*(1.0 + scale)*(pot_term_a - pot_term_b));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    double pot_term_c = 0.0;
                    for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                    {
                        double fik_term = 0.0;
                        if(StDir->scal_prod < 0 && gamma->value > StDir->phi_tu)
                        {
                            fik_term = 1.0;
                        }
                        else if(StDir->scal_prod > 0 && gamma->value < StDir->phi_tl)
                        {
                            fik_term = 0.0;
                        }
                        else
                        {
                            fik_term = 0.5 + StDir->alpha*(gamma->value - 0.5) - StDir->beta*sqrt(gamma->value*(1.0 - gamma->value));
                        }
                        pot_term_c += -StDir->stencil_weight*(fik_term - gamma->value);
                    }

                    dPhi_dt += (IP.Properties(i,j,k).get_energy( beta->index, gamma->index) -
                                IP.Properties(i,j,k).get_energy(alpha->index, gamma->index))
                              *(gamma->laplacian + pot_term_c + Prefactor*0.5);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            dPhi_dt += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                dPhi_dt *= IP.Properties(i,j,k).get_mobility(alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void DoubleObstacle::CalculatePhaseFieldIncrementsSharpDR(PhaseField& Phase,
                                                          InterfaceProperties& IP)
{
    dVector3 n_vector = {1.0, 0.0, 0.0};
    n_vector.normalize();

    double dx = Phase.Grid.dx;
    vector<PotentialCorrections> StencilDirections;

    for(auto ds = Phase.LStencil.cbegin(); ds != Phase.LStencil.cend(); ds++)
    {
        double d_x = ds->di;
        double d_y = ds->dj;
        double d_z = ds->dk;

        PotentialCorrections dirNew;

        dirNew.d_x = d_x;
        dirNew.d_y = d_y;
        dirNew.d_z = d_z;

        dirNew.stencil_weight = ds->weight;

        dirNew.scal_prod = (d_x*n_vector[0] + d_y*n_vector[1] + d_z*n_vector[2])*dx;

        dirNew.alpha     = cos(Pi*(dirNew.scal_prod)/Phase.Grid.Eta);
        dirNew.beta      = sin(Pi*(dirNew.scal_prod)/Phase.Grid.Eta);
        dirNew.phi_tl    = 0.5*cos(Pi*(1.0 - (dirNew.scal_prod)/Phase.Grid.Eta)) + 0.5;
        dirNew.phi_tu    = 0.5*cos(Pi*(dirNew.scal_prod)/Phase.Grid.Eta) + 0.5;

        StencilDirections.push_back(dirNew);
    }

    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,0,)
    {
        if(Phase.FieldsDR(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.FieldsDR(i,j,k);

            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double dPhi_dt = 0.0;

                double pot_term_a = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                {
                    double fik_term = 0;
                    if(StDir->scal_prod < 0 && alpha->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && alpha->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(alpha->value - 0.5) - StDir->beta*sqrt(alpha->value*(1.0 - alpha->value));
                    }
                    pot_term_a += -StDir->stencil_weight*(fik_term - alpha->value);
                }

                double pot_term_b = 0.0;
                for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend() ; StDir++)
                {
                    double fik_term = 0.0;
                    if(StDir->scal_prod < 0 && beta->value > StDir->phi_tu)
                    {
                        fik_term = 1.0;
                    }
                    else if(StDir->scal_prod > 0 && beta->value < StDir->phi_tl)
                    {
                        fik_term = 0.0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(beta->value - 0.5) - StDir->beta*sqrt(beta->value*(1.0 - beta->value));
                    }
                    pot_term_b += -StDir->stencil_weight*(fik_term - beta->value);
                }

                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio*
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                dPhi_dt = 0.5*(1.0 + scale)*IP.PropertiesDR(i,j,k).get_energy(alpha->index, beta->index)
                         *((alpha->laplacian - beta->laplacian) + (pot_term_a - pot_term_b));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    double pot_term_c = 0.0;
                    for(auto StDir = StencilDirections.cbegin(); StDir != StencilDirections.cend(); StDir++)
                    {
                        double fik_term = 0.0;
                        if(StDir->scal_prod < 0 && gamma->value > StDir->phi_tu)
                        {
                            fik_term = 1.0;
                        }
                        else if(StDir->scal_prod > 0 && gamma->value < StDir->phi_tl)
                        {
                            fik_term = 0.0;
                        }
                        else
                        {
                            fik_term = 0.5 + StDir->alpha*(gamma->value - 0.5) - StDir->beta*sqrt(gamma->value*(1.0 - gamma->value));
                        }
                        pot_term_c += -StDir->stencil_weight*(fik_term - gamma->value);
                    }

                    dPhi_dt += (IP.PropertiesDR(i,j,k).get_energy( beta->index, gamma->index) -
                                IP.PropertiesDR(i,j,k).get_energy(alpha->index, gamma->index))
                              *(gamma->laplacian + pot_term_c + Prefactor*0.5);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            dPhi_dt += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                dPhi_dt *= IP.PropertiesDR(i,j,k).get_mobility(alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculateCurvatureDrivingForce(PhaseField& Phase,
                                                    InterfaceProperties& IP,
                                                    DrivingForce& dG)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    const double Prefactor2 = Phase.Grid.Eta/(2.0*Pi);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        {
            const NodePF& locPF = Phase.Fields(i,j,k);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta and alpha->value*beta->value > DBL_EPSILON)
            {
                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio *
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                double loc_dG = IP.Properties(i,j,k).get_energy(alpha->index,beta->index)
                                 *((alpha->laplacian + 0.5*(1.0 + scale)*Prefactor*alpha->value) -
                                   ( beta->laplacian + 0.5*(1.0 + scale)*Prefactor* beta->value));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    loc_dG += (IP.Properties(i,j,k).get_energy( beta->index, gamma->index) -
                               IP.Properties(i,j,k).get_energy(alpha->index, gamma->index))
                             *(gamma->laplacian + 0.5*(1.0 + scale)*Prefactor*gamma->value);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                loc_dG += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            loc_dG += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                loc_dG *= Prefactor2 * scale / sqrt(alpha->value*beta->value);
                loc_dG *= IP.CurvatureFactor;
                dG.Force(i,j,k).add_raw(alpha->index, beta->index, loc_dG);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase,
                                                   InterfaceProperties& IP)
{
    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSR(Phase, IP);
            break;
        }
        case Resolutions::Dual:
        {
            CalculatePhaseFieldIncrementsDR(Phase, IP);
            break;
        }
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase,
                                                   InterfaceProperties& IP,
                                                   DrivingForce& dG)
{
    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSR(Phase, IP);
            break;
        }
        case Resolutions::Dual:
        {
            CalculatePhaseFieldIncrementsDR(Phase, IP);
            break;
        }
    }
    if (IP.FullAnisotropy)
    {
        CalculateFullAnisotropy(Phase,IP,dG);
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.Fields(i,j,k);
            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio *
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                double dPhi_dt = IP.Properties(i,j,k).get_energy(alpha->index,beta->index)
                                 *((alpha->laplacian + 0.5*(1.0 + scale)*Prefactor*alpha->value) -
                                   ( beta->laplacian + 0.5*(1.0 + scale)*Prefactor* beta->value));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    dPhi_dt += (IP.Properties(i,j,k).get_energy( beta->index, gamma->index) -
                                IP.Properties(i,j,k).get_energy(alpha->index, gamma->index))
                              *(gamma->laplacian + 0.5*(1.0 + scale)*Prefactor*gamma->value);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            dPhi_dt += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                dPhi_dt *= IP.Properties(i,j,k).get_mobility(alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrementsDR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,0,)
    {
        if(Phase.FieldsDR(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.FieldsDR(i,j,k);
            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio *
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                double dPhi_dt = 0.5*(1.0 + scale)*IP.PropertiesDR(i,j,k).get_energy(alpha->index,beta->index)
                                 *((alpha->laplacian + Prefactor*alpha->value) -
                                   ( beta->laplacian + Prefactor* beta->value));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    dPhi_dt += 0.5*(1.0 + scale)*(IP.PropertiesDR(i,j,k).get_energy( beta->index, gamma->index) -
                                IP.PropertiesDR(i,j,k).get_energy(alpha->index, gamma->index))
                              *(gamma->laplacian + Prefactor*gamma->value);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                dPhi_dt += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            dPhi_dt += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                dPhi_dt *= IP.PropertiesDR(i,j,k).get_mobility(alpha->index, beta->index)*norm_1*scale;
                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase,
                                                   InterfaceProperties& IP,
                                                   InterfaceRegularization& Kappa)
{
    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSR(Phase, IP, Kappa);
            break;
        }
        case Resolutions::Dual:
        {
            CalculatePhaseFieldIncrementsDR(Phase, IP, Kappa);
            break;
        }
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase,
                                                   InterfaceProperties& IP,
                                                   DrivingForce& dG,
                                                   InterfaceRegularization& Kappa)
{
    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            CalculatePhaseFieldIncrementsSR(Phase, IP, Kappa);
            break;
        }
        case Resolutions::Dual:
        {
            CalculatePhaseFieldIncrementsDR(Phase, IP, Kappa);
            break;
        }
    }
    if (IP.FullAnisotropy)
    {
        CalculateFullAnisotropy(Phase,IP,dG);
    }
}

void DoubleObstacle::CalculatePhaseFieldIncrementsSR(PhaseField& Phase,
                                                     InterfaceProperties& IP,
                                                     InterfaceRegularization& Kappa)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    const double Prefactor2 = Phase.Grid.Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.Fields(i,j,k);
            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio *
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                double loc_increment = IP.Properties(i,j,k).get_energy(alpha->index,beta->index)
                                 *((alpha->laplacian + 0.5*(1.0 + scale)*Prefactor*alpha->value) -
                                   ( beta->laplacian + 0.5*(1.0 + scale)*Prefactor* beta->value));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) and (gamma != beta))
                {
                    loc_increment += (IP.Properties(i,j,k).get_energy( beta->index, gamma->index) -
                                      IP.Properties(i,j,k).get_energy(alpha->index, gamma->index))*
                                     (gamma->laplacian  + 0.5*(1.0 + scale)*Prefactor*gamma->value);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                loc_increment += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            loc_increment += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                double dPhi_dt = loc_increment * IP.Properties(i,j,k).get_mobility(alpha->index, beta->index) * norm_1 * scale;
                if(IP.RegularizationFactor != 1.0 or IP.CurvatureFactor != 1.0)
                {
                    dPhi_dt *= IP.RegularizationFactor;

                    if(alpha->value != 0.0 and beta->value != 0.0)
                    {
                        double loc_curvature = loc_increment * Prefactor2 * scale * norm_1 / sqrt(alpha->value*beta->value);
                        loc_curvature *= IP.CurvatureFactor - IP.RegularizationFactor;
                        Kappa.Curvature(i,j,k).add_raw(alpha->index, beta->index, loc_curvature);
                    }
                }
                Phase.FieldsDot(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::CalculatePhaseFieldIncrementsDR(PhaseField& Phase,
                                                     InterfaceProperties& IP,
                                                     InterfaceRegularization& Kappa)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    const double Prefactor2 = Phase.Grid.Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.FieldsDR,0,)
    {
        if(Phase.FieldsDR(i,j,k).wide_interface())
        {
            const NodePF& locPF = Phase.FieldsDR(i,j,k);
            double norm_1 = 1.0/Phase.LocalN(locPF);

            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != locPF.cend(); ++beta)
            if(alpha != beta)
            {
                double scale = sqrt(Phase.FieldsProperties[alpha->index].VolumeRatio *
                                    Phase.FieldsProperties[ beta->index].VolumeRatio);

                double loc_increment = 0.5*(1.0 + scale)*IP.PropertiesDR(i,j,k).get_energy(alpha->index,beta->index)
                                 *((alpha->laplacian + Prefactor*alpha->value) -
                                   ( beta->laplacian + Prefactor* beta->value));

                if(locPF.size() > 2)
                for(auto gamma  = locPF.cbegin();
                         gamma != locPF.cend(); ++gamma)
                if((gamma != alpha) && (gamma != beta))
                {
                    loc_increment += 0.5*(1.0 + scale)*(IP.PropertiesDR(i,j,k).get_energy( beta->index, gamma->index) -
                                      IP.PropertiesDR(i,j,k).get_energy(alpha->index, gamma->index))
                                    *(gamma->laplacian  + Prefactor*gamma->value);

                    switch(IP.TripleJunctionModel)
                    {
                        case TripleJunctionModels::Model:
                        {
                            if(IP.TripleJunctionFactor != 0.0 or scale < 1.0)
                            {
                                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;
                                size_t pIndexG = Phase.FieldsProperties[gamma->index].Phase;

                                loc_increment += max(IP.TripleJunctionFactor, (1.0 - scale))*Prefactor
                                               *(IP.InterfaceEnergy(pIndexA, pIndexB).Energy +
                                                 IP.InterfaceEnergy(pIndexB, pIndexG).Energy +
                                                 IP.InterfaceEnergy(pIndexA, pIndexG).Energy)
                                               *(alpha->value*gamma->value - beta->value*gamma->value);
                            }
                            break;
                        }
                        case TripleJunctionModels::Value:
                        {
                            loc_increment += IP.TripleJunctionEnergy*Prefactor
                                           *(alpha->value*gamma->value - beta->value*gamma->value);
                            break;
                        }
                        case TripleJunctionModels::None:
                        {
                            break;
                        }
                    }
                }
                double dPhi_dt = loc_increment * IP.PropertiesDR(i,j,k).get_mobility(alpha->index, beta->index) * norm_1 * scale;
                if(IP.RegularizationFactor != 1.0 or IP.CurvatureFactor != 1.0)
                {
                    dPhi_dt *= IP.RegularizationFactor;

                    if(alpha->value*beta->value > DBL_EPSILON)
                    {
                        double loc_curvature = loc_increment * Prefactor2 * scale * norm_1 / sqrt(alpha->value*beta->value);
                        loc_curvature *= IP.CurvatureFactor - IP.RegularizationFactor;
                        Kappa.Curvature(i,j,k).add_raw(alpha->index, beta->index, loc_curvature);
                    }
                }
                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DoubleObstacle::StabilizeThinChannels(PhaseField& Phase, InterfaceProperties& IP, double penalty)
{
    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);
    //const double Prefactor2 = Phase.Eta/(2.0*Pi);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Fields(i,j,k).interface_halo())
        {
            std::vector<NodeAB<dVector3,dVector3>> allNormals;
            std::vector<iVector3> directions;
            // Collect all normals and directions on both sides of the thin channel
            for(auto it = Phase.LStencil.begin(); it != Phase.LStencil.end(); it++)
            if(Phase.Fields(i + it->di,j + it->dj,k + it->dk).interface())
            {
                allNormals.push_back(Phase.Normals(i + it->di,j + it->dj,k + it->dk));
                directions.push_back((iVector3){it->di,it->dj,it->dk});
            }

            std::vector<NodeAB<double,double>> allProducts;

            // Find opposite directions and normals products
            for(size_t n =   0; n < directions.size() - 1; n++)
            for(size_t m = n+1; m < directions.size(); m++)
            if(directions[n] == directions[m]*(-1))
            {
                NodeAB<double,double> locProduct;
                for(auto nm = allNormals[n].begin(); nm != allNormals[n].end(); nm++)
                if(allNormals[m].present(nm->indexA, nm->indexB))
                {
                    locProduct.set_asym1(nm->indexA, nm->indexB, nm->value1*allNormals[m].get_asym1(nm->indexA, nm->indexB));
                }
                allProducts.push_back(locProduct);
            }

            // Analyze normal products and detect opposite normal directions
            for(size_t p = 0; p < allProducts.size(); p++)
            for(auto it = allProducts[p].begin(); it != allProducts[p].end(); it++)
            if(it->value1 < 0.0)
            {
                Phase.FieldsDot(i,j,k).add_asym1(it->indexA, it->indexB, it->value1*penalty*Prefactor);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double DoubleObstacle::PointEnergy(const PhaseField& Phase,
                                   const InterfaceProperties& IP,
                                   int i, int j, int k) const
{
    double energy = 0.0;
    double Prefactor = Phase.Grid.Eta*Phase.Grid.Eta/(Pi*Pi);

    if(Phase.Fields(i,j,k).wide_interface())
    {
        for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        for(auto  beta  = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        if(alpha != beta)
        {
            energy += 4.0*IP.Properties(i,j,k).get_energy(alpha->index, beta->index)*
                          (alpha->value*beta->value - Prefactor*(alpha->gradient*beta->gradient))/Phase.Grid.Eta;
        }
    }
    return energy;
}

double DoubleObstacle::PointEnergy(const PhaseField& Phase,
                                   const InterfaceProperties& IP,
                                   const ElasticProperties& EP,
                                   const int i, const int j, const int k) const
{
    double energy = 0.0;
    double Prefactor = Phase.Grid.Eta*Phase.Grid.Eta/(Pi*Pi);

    if(Phase.Fields(i,j,k).wide_interface())
    {
        for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        for(auto  beta  = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        if(alpha != beta)
        {
            energy += 4.0*IP.Properties(i,j,k).get_energy(alpha->index, beta->index)*
                          EP.DeformationGradientsElastic(i,j,k).determinant()*
                          (alpha->value*beta->value - Prefactor*(alpha->gradient*beta->gradient))/Phase.Grid.Eta;
        }
    }
    return energy;
}

double DoubleObstacle::Energy(const PhaseField& Phase,
                              const InterfaceProperties& IP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    OP_MPI_Allreduce(&loc_Energy, &Energy, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    return Energy*Phase.Grid.CellVolume();
}

double DoubleObstacle::Energy(const PhaseField& Phase,
                              const InterfaceProperties& IP,
                              const ElasticProperties& EP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, EP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    OP_MPI_Allreduce(&loc_Energy, &Energy, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    return Energy*Phase.Grid.CellVolume();
}

double DoubleObstacle::AverageEnergyDensity(const PhaseField& Phase,
                                            const InterfaceProperties& IP) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, IP, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    double loc_Energy = Energy;
    OP_MPI_Allreduce(&loc_Energy, &Energy, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif

    return Energy/double(Phase.Grid.TotalNumberOfCells());
}

void DoubleObstacle::WriteEnergyVTK(const int tStep,
                                    const Settings& locSettings,
                                    const PhaseField& Phase,
                                    const InterfaceProperties& IP) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Interface Energy Density [J/m^3]", [this, &Phase, &IP](int i,int j,int k){return double(PointEnergy(Phase, IP, i,j,k));}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "InterfaceEnergyDensity_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void DoubleObstacle::WriteEnergyVTK(const int tStep,
                                    const Settings& locSettings,
                                    const PhaseField& Phase,
                                    const InterfaceProperties& IP,
                                    const ElasticProperties& EP) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t){"Interface Energy Density [J/m^3]", [this, &Phase, &IP, &EP](int i,int j,int k){return double(PointEnergy(Phase, IP, EP, i,j,k));}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "InterfaceEnergyDensity_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void DoubleObstacle::CalculateFullAnisotropy(const PhaseField& Phase,
                                             InterfaceProperties& IP,
                                             DrivingForce& dG)
{
    IP.InterfaceStiffnessTMP.Clear();

    const double Prefactor1 = 4.0/Phase.Grid.Eta;
    const double Prefactor2 = Phase.Grid.Eta*Phase.Grid.Eta/Pi/Pi;
    const double Prefactor3 = - 4.0*Phase.Grid.Eta*Phase.Grid.Eta/Pi/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    if(Phase.Fields(i,j,k).interface())
    {
        for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        for(auto  beta  = alpha + 1;
                  beta != Phase.Fields(i,j,k).cend(); ++beta)
        if(alpha != beta)
        {
            const size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
            const size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

            if(IP.InterfaceEnergy(pIndexA, pIndexB).Model != InterfaceEnergyModels::Iso)
            if(IP.InterfaceEnergy(pIndexA, pIndexB).Type  == InterfaceEnergyModelType::Energy)
            {
                const double I_AB = Prefactor1*(alpha->value*beta->value - Prefactor2*(alpha->gradient*beta->gradient));
                const dVector3 value = (IP.dEnergy_dGradientAlpha(Phase, alpha, beta) - IP.dEnergy_dGradientAlpha(Phase, beta, alpha))*(-I_AB);

                if (value.abs() > 0.0)
                {
                    IP.InterfaceStiffnessTMP(i,j,k).add_asym1(alpha->index, beta->index, value);
                }

                for(auto gs = Phase.GStencil.cbegin(); gs != Phase.GStencil.cend(); ++gs)
                {
                    const double sigma     = IP.Properties(i+gs->di, j+gs->dj, k+gs->dk).get_energy(alpha->index, beta->index);
                    const double dsigma_dx = gs->weightX * sigma;
                    const double dsigma_dy = gs->weightY * sigma;
                    const double dsigma_dz = gs->weightZ * sigma;
                    const dVector3 sigma_gradient ={dsigma_dx,dsigma_dy,dsigma_dz};
                    const double value = Prefactor3*(sigma_gradient*(alpha->gradient - beta->gradient));
                    dG.Force(i,j,k).add_raw(alpha->index, beta->index, -value);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    const double DWeights[3] = {-0.5/Phase.Grid.dx, 0.0, 0.5/Phase.Grid.dx};
    if (Phase.Grid.dNx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP.InterfaceStiffnessTMP,0,)
        {
            if(Phase.Fields(i,j,k).interface())
            for (int ii = -1; ii <= +1; ii+=2)
            for (auto it = IP.InterfaceStiffnessTMP(i+ii,j,k).cbegin();
                      it < IP.InterfaceStiffnessTMP(i+ii,j,k).cend(); ++it)
            {
                dG.Force(i,j,k).add_raw(it->indexA, it->indexB,-DWeights[ii+1]*it->value1[0]);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    if (Phase.Grid.dNy)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP.InterfaceStiffnessTMP,0,)
        {
            if(Phase.Fields(i,j,k).interface())
            for (int ii = -1; ii <= +1; ii+=2)
            for (auto it = IP.InterfaceStiffnessTMP(i,j+ii,k).cbegin();
                      it < IP.InterfaceStiffnessTMP(i,j+ii,k).cend(); ++it)
            {
                dG.Force(i,j,k).add_raw(it->indexA, it->indexB,-DWeights[ii+1]*it->value1[1]);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    if (Phase.Grid.dNz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP.InterfaceStiffnessTMP,0,)
        {
            if(Phase.Fields(i,j,k).interface())
            for (int ii = -1; ii <= +1; ii+=2)
            for (auto it = IP.InterfaceStiffnessTMP(i,j,k+ii).cbegin();
                      it < IP.InterfaceStiffnessTMP(i,j,k+ii).cend(); ++it)
            {
                dG.Force(i,j,k).add_raw(it->indexA, it->indexB,-DWeights[ii+1]*it->value1[2]);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

}// namespace openphase
