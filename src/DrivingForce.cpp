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

#include "DrivingForce.h"
#include "Settings.h"
#include "VTK.h"
#include "PhaseField.h"
#include "InterfaceProperties.h"
#include "BoundaryConditions.h"
#include "Temperature.h"

namespace openphase
{

using namespace std;
DrivingForce::DrivingForce(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void DrivingForce::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "DrivingForce";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    Limiting = true;
    Limit.Allocate(Nphases,Nphases);
    Averaging = false;

    WeightsMode = AveragingWeightsModes::PhaseFields;

    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            Range = Grid.iWidth;

            if(Grid.iWidth < 5.0)
            {
                PhiThreshold = Grid.iWidth/15.0;
            }
            else
            {
                PhiThreshold = 1.0/3.0;
            }
            break;
        }
        case Resolutions::Dual:
        {
            Range = (Grid.iWidth)/2;

            if(Grid.iWidth < 5.0)
            {
                PhiThreshold = Grid.iWidth/30.0;
            }
            else
            {
                PhiThreshold = 1.0/6.0;
            }
            break;
        }
    }

    Force.Allocate(Grid, Range);

    OvershootCounter = 0;
    MAXOvershootPOS.Allocate(Nphases,Nphases);
    MAXOvershootNEG.Allocate(Nphases,Nphases);
    MAXDrivingForcePOS.Allocate(Nphases,Nphases);
    MAXDrivingForceNEG.Allocate(Nphases,Nphases);

    locSettings.AddForRemeshing(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void DrivingForce::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("DrivingForce input");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());

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

void DrivingForce::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    Averaging    = FileInterface::ReadParameterB(inp, moduleLocation, string("Average"), false, true);
    Unifying     = FileInterface::ReadParameterB(inp, moduleLocation, string("bUnify"), false, false);
    Range        = FileInterface::ReadParameterI(inp, moduleLocation, string("Range"), false, Range);
    PhiThreshold = FileInterface::ReadParameterD(inp, moduleLocation, string("Threshold"), false, PhiThreshold);
    Limiting     = FileInterface::ReadParameterB(inp, moduleLocation, vector<string>{"Limiting","dGcut"}, false, true);
    string tmp2  = FileInterface::ReadParameterK(inp, moduleLocation, "WeightsMode", false, "PHASEFIELDS");

    if(tmp2 == string("RANGE"))
    {
        WeightsMode = AveragingWeightsModes::Range;
    }
    else if(tmp2 == string("PHASEFIELDS"))
    {
        WeightsMode = AveragingWeightsModes::PhaseFields;
    }
    else if(tmp2 == string("COUNTER"))
    {
        WeightsMode = AveragingWeightsModes::Counter;
    }
    else
    {
        ConsoleOutput::WriteWarning("Wrong driving force weights mode specified!\nThe default \"PHASEFIELDS\" mode is used!", thisclassname, "ReadInput()");
    }

    for(size_t m = 0; m < Nphases; ++m)
    for(size_t n = m; n < Nphases; ++n)
    {
        stringstream converter;
        converter << m << "_" << n;
        string counter = converter.str();

        Limit(m,n) = FileInterface::ReadParameterD(inp, moduleLocation, vector<string>{"Limit_" + counter, "CutOff_" + counter}, false, 0.95);
        if(Limit(m,n) == 0.0)
        {
            stringstream message;
            message << "Limit(" << m << ", " << n << ") = 0.0, please correct the driving force input!";
            ConsoleOutput::WriteExit(message.str(), thisclassname, "ReadInput()");
            exit(13);
        }
        Limit(n,m) = Limit(m,n);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void DrivingForce::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,Force.Bcells(),)
    {
        Force(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Force);
    BC.SetY(Force);
    BC.SetZ(Force);
}

NodeDF DrivingForce::CalcUnified(const PhaseField& Phase, const BoundaryConditions& BC)
{
    NodeDF Unified;
    //OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0, reduction(+:UnifiedDrivingForce),reduction(+:Counter))
    STORAGE_LOOP_BEGIN(i,j,k,Force,0)
    {
        if (Phase.Fields(i,j,k).interface())
        for(const auto& it : Force(i,j,k))
        {
            const double PhiAlpha = Phase.Fields(i,j,k).get_value(it.indexA);
            const double PhiBeta  = Phase.Fields(i,j,k).get_value(it.indexB);
            const double weight   = std::sqrt(PhiAlpha*PhiBeta);

            Unified.add_raw(it.indexA, it.indexB, it.raw*weight);
            Unified.add_average(it.indexA, it.indexB, it.average*weight);
            Unified.add_weight(it.indexA, it.indexB, weight);
        }
    }
    STORAGE_LOOP_END

    for(auto& it : Unified)
    {
        const double norm = Unified.get_weight(it.indexA,it.indexB);
        if (norm > DBL_EPSILON)
        {
            it.raw     /= norm;
            it.average /= norm;
        }
    }
    return Unified;
}

void DrivingForce::Unify(const PhaseField& Phase, const BoundaryConditions& BC)
{
    if (Unifying)
    {
        NodeDF Unified = CalcUnified(Phase, BC);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0,)
        {
            if (Phase.Fields(i,j,k).interface())
            {
                Force(i,j,k) = Unified;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        SetBoundaryConditions(BC);
    }
}

void DrivingForce::Average(const PhaseField& Phase, const BoundaryConditions& BC)
{
    if(Averaging)
    {
        SetWeights(Phase);
        SetBoundaryConditions(BC);
        CollectAverage(Phase);
        SetBoundaryConditions(BC);
        DistributeAverage(Phase);
    }
    else
    {
        SkipAverage(Phase);
    }
}

void DrivingForce::SkipAverage(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Force(i,j,k).begin();
                 it != Force(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get_value(it->indexB);

            it->weight  = sqrt(PhiAlpha*PhiBeta);
            it->average = it->raw;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::SetWeights(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Force(i,j,k).begin();
                 it != Force(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get_value(it->indexB);

            it->weight = sqrt(PhiAlpha*PhiBeta);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::CollectAverage(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Force(i,j,k).begin();
                 it != Force(i,j,k).end(); ++it)
        {
            double scale = Phase.FieldsProperties[it->indexA].VolumeRatio *
                           Phase.FieldsProperties[it->indexB].VolumeRatio;

            if (it->weight > weight_threshold*scale)
            {
                double value       = 0.0;
                double sum_weights = 0.0;
                int    counter     = 0;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double weight_dist = Range - sqrt(ii*ii + jj*jj + kk*kk);
                    double weight_phi  = Force(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB);

                    if(weight_dist > 0.0 and weight_phi > 0.0)
                    {
                        double weight = 0.0;
                        switch(WeightsMode)
                        {
                            case AveragingWeightsModes::Range:
                            {
                                weight = weight_dist;
                                sum_weights += weight;
                                break;
                            }
                            case AveragingWeightsModes::PhaseFields:
                            {
                                weight = weight_phi;
                                sum_weights += weight;
                                break;
                            }
                            case AveragingWeightsModes::Counter:
                            {
                                counter++;
                                weight = 1.0;
                                break;
                            }
                        }
                        value += weight*Force(i+ii, j+jj, k+kk).get_raw(it->indexA, it->indexB);
                    }
                }
                switch(WeightsMode)
                {
                    case AveragingWeightsModes::Range:
                    case AveragingWeightsModes::PhaseFields:
                    {
                        if(sum_weights > 0.0)
                        {
                            it->tmp = value/sum_weights;
                        }
                        else
                        {
                            it->tmp = it->raw;
                        }
                        break;
                    }
                    case AveragingWeightsModes::Counter:
                    {
                        if(counter != 0)
                        {
                            it->tmp = value/counter;
                        }
                        else
                        {
                            it->tmp = it->raw;
                        }
                        break;
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::DistributeAverage(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Force,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Force(i,j,k).begin();
                 it != Force(i,j,k).end(); ++it)
        {
            if (it->weight > 0.0)
            {
                double value       = 0.0;
                double sum_weights = 0.0;
                int    counter     = 0;

                double scale = Phase.FieldsProperties[it->indexA].VolumeRatio *
                               Phase.FieldsProperties[it->indexB].VolumeRatio;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double weight_dist = Range - sqrt(ii*ii + jj*jj + kk*kk);
                    double weight_phi  = Force(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB);
                    if(weight_dist > 0.0 and weight_phi > weight_threshold*scale)
                    {
                        double weight = 0.0;
                        switch(WeightsMode)
                        {
                            case AveragingWeightsModes::Range:
                            {
                                weight = weight_dist;
                                sum_weights += weight;
                                break;
                            }
                            case AveragingWeightsModes::PhaseFields:
                            {
                                weight = weight_phi;
                                sum_weights += weight;
                                break;
                            }
                            case AveragingWeightsModes::Counter:
                            {
                                counter++;
                                weight = 1.0;
                                break;
                            }
                        }
                        value += weight*Force(i+ii, j+jj, k+kk).get_tmp(it->indexA, it->indexB);
                    }
                }

                switch(WeightsMode)
                {
                    case AveragingWeightsModes::Range:
                    case AveragingWeightsModes::PhaseFields:
                    {
                        if(sum_weights > 0.0)
                        {
                            it->average = value/sum_weights;
                        }
                        else
                        {
                            it->average = it->tmp;
                        }
                        break;
                    }
                    case AveragingWeightsModes::Counter:
                    {
                        if(counter != 0)
                        {
                            it->average = value/counter;
                        }
                        else
                        {
                            it->average = it->tmp;
                        }
                        break;
                    }
                }
            }
            else if(Phase.FieldsProperties[it->indexA].Stage == GrainStages::Seed or
                    Phase.FieldsProperties[it->indexB].Stage == GrainStages::Seed)
            {
                it->average = it->raw;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            MergePhaseFieldIncrementsSR(Phase, IP);
            break;
        }
        case Resolutions::Dual:
        {
            MergePhaseFieldIncrementsDR(Phase, IP);
            break;
        }
    }
}

void DrivingForce::MergePhaseFieldIncrementsSR(PhaseField& Phase,
                                               InterfaceProperties& IP)
{
    double locMaxPsi = 0.0;
    long int locOvershootCounter = 0;
    Matrix<double> locMAXOvershootPOS = MAXOvershootPOS;
    Matrix<double> locMAXOvershootNEG = MAXOvershootNEG;
    Matrix<double> locMAXDrivingForcePOS = MAXDrivingForcePOS;
    Matrix<double> locMAXDrivingForceNEG = MAXDrivingForceNEG;

    const double Prefactor = Pi/Phase.Grid.Eta;
    //const double Prefactor = 2.0*Pi/(Phase.Eta*Phase.LocalN(Phase.Fields(i,j,k)));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0, \
            reduction(+:locOvershootCounter) \
            reduction(MatrixDMAX:locMAXOvershootPOS, locMAXDrivingForcePOS) \
            reduction(MatrixDMIN:locMAXOvershootNEG, locMAXDrivingForceNEG) \
            reduction(max: locMaxPsi))
    {
        if (Phase.Fields(i,j,k).interface())
        {
            for(auto it  = Force(i,j,k).begin();
                     it != Force(i,j,k).end(); ++it)
            {
                double locdG = 0.0;
                if(Averaging)
                {
                    locdG = it->average;
                }
                else
                {
                    locdG = it->raw;
                }

                double norm = it->weight;

                size_t indexA = it->indexA;
                size_t indexB = it->indexB;
                size_t pIndexA = Phase.FieldsProperties[indexA].Phase;
                size_t pIndexB = Phase.FieldsProperties[indexB].Phase;

                if(Phase.FieldsProperties[indexA].Stage == GrainStages::Seed or
                   Phase.FieldsProperties[indexB].Stage == GrainStages::Seed)
                {
                    /* In the case of a freshly nucleated grain, a small value
                    is assigned to the normalization coefficient to allow the
                    initial growth of the nucleus. */

                    norm += 1.0e-6;
                }

                if(Limiting)
                {
                    /* If the local driving force exceeds the interface stabilizing
                    force, the interface profile can distort leading to simulation
                    artifacts. Therefore the driving force is limited to a safe
                    value controlled by the user specified Limit(pIndexA,pIndexB)
                    parameter.*/

                    double absDG = fabs(locdG);
                    double allowedDG = Limit(pIndexA,pIndexB)*Prefactor
                                      *IP.InterfaceEnergy(pIndexA,pIndexB).MaxEnergy;

                    allowedDG *= IP.RegularizationFactor;

                    /* Collecting statistics for later output */
                    if(absDG > 0.4*allowedDG)
                    {
                        locOvershootCounter++;
                    }

                    double tmpMAXOvershoot  = locdG/allowedDG;

                    if(tmpMAXOvershoot > locMAXOvershootPOS(pIndexA,pIndexB))
                    {
                        locMAXOvershootPOS(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    if(tmpMAXOvershoot < locMAXOvershootNEG(pIndexA,pIndexB))
                    {
                        locMAXOvershootNEG(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    locMAXDrivingForcePOS(pIndexA,pIndexB) = max(locdG, locMAXDrivingForcePOS(pIndexA,pIndexB));
                    locMAXDrivingForceNEG(pIndexA,pIndexB) = min(locdG, locMAXDrivingForceNEG(pIndexA,pIndexB));

                    /* End collecting statistics for later output */
                    locdG = allowedDG*tanh(locdG/allowedDG);
                }

                double dPsi_dt = locdG*IP.Properties(i,j,k).get_mobility(indexA,indexB)*norm*Prefactor;

                locMaxPsi = max(fabs(dPsi_dt), locMaxPsi);

                Phase.FieldsDot(i,j,k).add_asym1(indexA,indexB,dPsi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    maxPsi = locMaxPsi;

    OvershootCounter  += locOvershootCounter;
    MAXOvershootPOS    = locMAXOvershootPOS;
    MAXOvershootNEG    = locMAXOvershootNEG;
    MAXDrivingForcePOS = locMAXDrivingForcePOS;
    MAXDrivingForceNEG = locMAXDrivingForceNEG;
}

void DrivingForce::MergePhaseFieldIncrementsDR(PhaseField& Phase,
                                               InterfaceProperties& IP)
{
    double locMaxPsi = 0.0;
    long int locOvershootCounter = 0;
    Matrix<double> locMAXOvershootPOS = MAXOvershootPOS;
    Matrix<double> locMAXOvershootNEG = MAXOvershootNEG;
    Matrix<double> locMAXDrivingForcePOS = MAXDrivingForcePOS;
    Matrix<double> locMAXDrivingForceNEG = MAXDrivingForceNEG;

    const double Prefactor = Pi/Phase.Grid.Eta;
    //const double Prefactor = 2.0*Pi/(Phase.Eta*Phase.LocalN(Phase.FieldsDR(i,j,k)));

    const double di = Grid.dNx*0.25;
    const double dj = Grid.dNy*0.25;
    const double dk = Grid.dNz*0.25;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.FieldsDR, 0, \
            reduction(+:locOvershootCounter) \
            reduction(MatrixDMAX:locMAXOvershootPOS, locMAXDrivingForcePOS) \
            reduction(MatrixDMIN:locMAXOvershootNEG, locMAXDrivingForceNEG) \
            reduction(max: locMaxPsi))
    {
        if (Phase.FieldsDR(i,j,k).interface())
        {
            NodeDF locForce = Force_at(i/2 - di, j/2 - dj, k/2 - dk);

            for(auto alpha  = Phase.FieldsDR(i,j,k).begin();
                     alpha != Phase.FieldsDR(i,j,k).end(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.FieldsDR(i,j,k).end(); ++beta)
            {
                /* For each pair of grains the driving force is treated
                individually. */

                double locdG = 0.0;
                if(Averaging)
                {
                    locdG = locForce.get_average(alpha->index, beta->index);
                }
                else
                {
                    locdG = locForce.get_raw(alpha->index, beta->index);
                }

                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                double norm = sqrt(alpha->value * beta->value);

                if(Phase.FieldsProperties[alpha->index].Stage == GrainStages::Seed or
                   Phase.FieldsProperties[ beta->index].Stage == GrainStages::Seed)
                {
                    /* In the case of a freshly nucleated grain, a small value
                    is assigned to the normalization coefficient to allow the
                    initial growth of the nucleus. */

                    norm += 1.0e-6;
                }

                if(Limiting)
                {
                    /* If the local driving force exceeds the interface stabilizing
                    force, the interface profile can distort leading to simulation
                    artifacts. Therefore the driving force is limited to a safe
                    value controlled by the user specified Limit(pIndexA,pIndexB)
                    parameter.*/

                    double absDG = fabs(locdG);
                    double allowedDG = Limit(pIndexA,pIndexB)*Prefactor
                                      *IP.InterfaceEnergy(pIndexA,pIndexB).MaxEnergy;

                    allowedDG *= IP.RegularizationFactor;

                    /* Collecting statistics for later output */
                    if(absDG > 0.4*allowedDG)
                    {
                        locOvershootCounter++;
                    }

                    double tmpMAXOvershoot  = locdG/allowedDG;

                    if(tmpMAXOvershoot > locMAXOvershootPOS(pIndexA,pIndexB))
                    {
                        locMAXOvershootPOS(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }
                    if(tmpMAXOvershoot < locMAXOvershootNEG(pIndexA,pIndexB))
                    {
                        locMAXOvershootNEG(pIndexA,pIndexB) = tmpMAXOvershoot;
                    }

                    locMAXDrivingForcePOS(pIndexA,pIndexB) = max(locdG, locMAXDrivingForcePOS(pIndexA,pIndexB));
                    locMAXDrivingForceNEG(pIndexA,pIndexB) = min(locdG, locMAXDrivingForceNEG(pIndexA,pIndexB));

                    /* End collecting statistics for later output */

                    locdG = allowedDG*tanh(locdG/allowedDG);
                }

                double dPsi_dt = locdG*IP.PropertiesDR(i,j,k).get_mobility(alpha->index, beta->index)*norm*Prefactor;

                locMaxPsi = max(fabs(dPsi_dt), locMaxPsi);

                Phase.FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index, dPsi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    maxPsi = locMaxPsi;

    OvershootCounter  += locOvershootCounter;
    MAXOvershootPOS    = locMAXOvershootPOS;
    MAXOvershootNEG    = locMAXOvershootNEG;
    MAXDrivingForcePOS = locMAXDrivingForcePOS;
    MAXDrivingForceNEG = locMAXDrivingForceNEG;
}

double DrivingForce::MaxTimeStep(PhaseField& Phase, InterfaceProperties& IP,
                     Settings& OPSettings, double TheorLimit, double NumerLimit)
{
    /** This function calculates the maximum time step for the phase field
    equation time integration.

    This function uses two different stability criteria:

    1) Neumann stability criterion for phase-field equation:
       dt < 0.5*dx^2/(InterfaceEnergy*InterfaceMobility)
       (TheorLimit is set to 1.0 by default);
    2) Numerical criterion, where the maximum phase-field increment
       in the whole simulation domain has to be lower than 1E-3
       (NumerLimit is set to 1E-3 by default).

    If phase field gets unstable, decrease TheorLimit, if diffusion field
    becomes unstable due to too high phase-field increments, decrease
    NumerLimit.*/

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &maxPsi, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif

    // Calculate theoretical maximum allowed time step
    double maxTheorTimeStep = TheorLimit*IP.ReportMaximumTimeStep();

    // Calculate maximum numerical allowed time step
    double maxNumerTimeStep = 0.0;
    if(maxPsi > DBL_EPSILON)
    {
        maxNumerTimeStep = NumerLimit/maxPsi;
    }

    // Return maximum allowed time step
    return min(maxTheorTimeStep, maxNumerTimeStep);
}

void DrivingForce::PrintDiagnostics()
{
#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &OvershootCounter, 1, OP_MPI_LONG, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, MAXOvershootPOS.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, MAXOvershootNEG.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MIN, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, MAXDrivingForcePOS.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, MAXDrivingForceNEG.data(), Nphases*Nphases, OP_MPI_DOUBLE, OP_MPI_MIN, OP_MPI_COMM_WORLD);
#endif

    if (OvershootCounter)
    {
        std::string message = "To prevent interface distortion the driving force has been limited " +
                               std::to_string(OvershootCounter) + " times!\n";

        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = n; m < Nphases; m++)
        if(MAXOvershootPOS(n, m) > DBL_EPSILON or MAXOvershootNEG(n, m) < -DBL_EPSILON or
           MAXOvershootPOS(m, n) > DBL_EPSILON or MAXOvershootNEG(m, n) < -DBL_EPSILON)
        {
            message += "   Phase pair (" + std::to_string(n) + ", " + std::to_string(m) + "): \n";
            if(n != m)
            {
                message += "       Max (positive) driving force value (overshoot ratio): "
                        + std::to_string(max(MAXDrivingForcePOS(n, m), -MAXDrivingForceNEG(m, n)))
                        + " ("
                        + std::to_string(max(MAXOvershootPOS(n, m), -MAXOvershootNEG(m, n)))
                        + ")\n";
                message += "       Min (negative) driving force value (overshoot ratio): "
                        + std::to_string(min(MAXDrivingForceNEG(n, m), -MAXDrivingForcePOS(m, n)))
                        + " ("
                        + std::to_string(min(MAXOvershootNEG(n, m), -MAXOvershootPOS(m, n)))
                        + ")\n";
            }
            else
            {
                message += "      Max (positive) driving force value (overshoot ratio): "
                        + std::to_string(MAXDrivingForcePOS(n, m))
                        + " ("
                        + std::to_string(MAXOvershootPOS(n, m))
                        + ")\n";
                message += "      Min (negative) driving force value (overshoot ratio): "
                        + std::to_string(MAXDrivingForceNEG(n, m))
                        + " ("
                        + std::to_string(MAXOvershootNEG(n, m))
                        + ")\n";
            }
        }
        ConsoleOutput::WriteStandard(thisclassname + "::PrintDiagnostics()", message);
        OvershootCounter = 0;
        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = 0; m < Nphases; m++)
        {
            MAXOvershootPOS(n,m) = 0.0;
            MAXOvershootNEG(n,m) = 0.0;
            MAXDrivingForcePOS(n,m) = 0.0;
            MAXDrivingForceNEG(n,m) = 0.0;
        }
    }
}

void DrivingForce::PrintPointStatistics(const int x, const int y, const int z) const
{
    std::stringstream pointstat;
    pointstat << "DrivingForce Indices:\t";

    for (auto alpha  = Force(x,y,z).cbegin();
              alpha != Force(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->indexA << ", " << alpha->indexB << "\t\t";
    }
    pointstat << endl;
    pointstat << "DrivingForce Values:\t";

    for (auto alpha  = Force(x,y,z).cbegin();
              alpha != Force(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->average << "\t\t";
    }
    pointstat << endl;

    ConsoleOutput::WriteSimple(pointstat.str());
}

void DrivingForce::WriteVTK(const Settings& locSettings, const int tStep,
                            const size_t indexA, const size_t indexB,
                            const int precision) const
{
    stringstream converter;
    converter << indexA << "," << indexB;
    string phases = converter.str();

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"dGraw(" + phases + ")", [indexA,indexB,this](int i,int j,int k){if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN(); return Force(i,j,k).get_raw(indexA, indexB);}});
    ListOfFields.push_back((VTK::Field_t) {"dGtmp(" + phases + ")", [indexA,indexB,this](int i,int j,int k){if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN(); return Force(i,j,k).get_tmp(indexA, indexB);}});
    ListOfFields.push_back((VTK::Field_t) {"dGavg(" + phases + ")", [indexA,indexB,this](int i,int j,int k){if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN(); return Force(i,j,k).get_average(indexA, indexB);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "DrivingForce_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void DrivingForce::WriteVTK(const Settings& locSettings,
                            const PhaseField& Phase,
                            const int tStep,
                            const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;

    for(size_t alpha =     0; alpha < Phase.Nphases; alpha++)
    for(size_t beta  = alpha;  beta < Phase.Nphases; beta++)
    {
        stringstream converter;
        converter << alpha << "," << beta;
        string phases = converter.str();

        ListOfFields.push_back((VTK::Field_t) {"dGraw(" + phases + ")", [alpha,beta,&Phase,this](int i,int j,int k){
            if(Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();
            double tempdG = 0.0;
            int counter = 0;
            for(auto it1  = Phase.Fields(i,j,k).cbegin();
                     it1 != Phase.Fields(i,j,k).cend(); ++it1)
            for(auto it2  = Phase.Fields(i,j,k).cbegin();
                     it2 != Phase.Fields(i,j,k).cend(); ++it2)
            if(Phase.FieldsProperties[it1->index].Phase == alpha and
               Phase.FieldsProperties[it2->index].Phase == beta)
            {
                tempdG += Force(i,j,k).get_raw(it1->index, it2->index);
                counter ++;
            }
            if(counter > 1.0)
            {
                tempdG /= counter;
            }
            return tempdG;
            }});
        ListOfFields.push_back((VTK::Field_t) {"dGtmp(" + phases + ")", [alpha,beta,&Phase,this](int i,int j,int k){
            if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();
            double tempdG = 0.0;
            int counter = 0;
            for(auto it1  = Phase.Fields(i,j,k).cbegin();
                     it1 != Phase.Fields(i,j,k).cend(); ++it1)
            for(auto it2  = Phase.Fields(i,j,k).cbegin();
                     it2 != Phase.Fields(i,j,k).cend(); ++it2)
            if(Phase.FieldsProperties[it1->index].Phase == alpha and
               Phase.FieldsProperties[it2->index].Phase == beta)
            {
                tempdG += Force(i,j,k).get_tmp(it1->index, it2->index);
                counter ++;
            }
            if(counter > 1.0)
            {
                tempdG /= counter;
            }
            return tempdG;
            }});
        ListOfFields.push_back((VTK::Field_t) {"dGavg(" + phases + ")", [alpha,beta,&Phase,this](int i,int j,int k){
            if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();
            double tempdG = 0.0;
            int counter = 0;
            for(auto it1  = Phase.Fields(i,j,k).cbegin();
                     it1 != Phase.Fields(i,j,k).cend(); ++it1)
            for(auto it2  = Phase.Fields(i,j,k).cbegin();
                     it2 != Phase.Fields(i,j,k).cend(); ++it2)
            if(Phase.FieldsProperties[it1->index].Phase == alpha and
               Phase.FieldsProperties[it2->index].Phase == beta)
            {
                tempdG += Force(i,j,k).get_average(it1->index, it2->index);
                counter ++;
            }
            if(counter > 1.0)
            {
                tempdG /= counter;
            }
            return tempdG;
            }});
        ListOfFields.push_back((VTK::Field_t) {"dGwgt(" + phases + ")", [alpha,beta,&Phase,this](int i,int j,int k){
            if (Force(i,j,k).size() == 0) return std::numeric_limits<double>::quiet_NaN();
            double tempdG = 0.0;
            int counter = 0;
            for(auto it1  = Phase.Fields(i,j,k).cbegin();
                     it1 != Phase.Fields(i,j,k).cend(); ++it1)
            for(auto it2  = Phase.Fields(i,j,k).cbegin();
                     it2 != Phase.Fields(i,j,k).cend(); ++it2)
            if(Phase.FieldsProperties[it1->index].Phase == alpha and
               Phase.FieldsProperties[it2->index].Phase == beta)
            {
                tempdG += Force(i,j,k).get_weight(it1->index, it2->index);
                counter ++;
            }
            if(counter > 1.0)
            {
                tempdG /= counter;
            }
            return tempdG;
            }});

    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "DrivingForcePhases_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void DrivingForce::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    Force.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

DrivingForce& DrivingForce::operator= (const DrivingForce& rhs)
{
    // protect against self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "DrivingForce")
    {
        thisclassname = rhs.thisclassname;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;
        Range   = rhs.Range;
        PhiThreshold = rhs.PhiThreshold;

        OvershootCounter   = rhs.OvershootCounter;
        MAXOvershootPOS    = rhs.MAXOvershootPOS;
        MAXOvershootNEG    = rhs.MAXOvershootNEG;
        MAXDrivingForcePOS = rhs.MAXDrivingForcePOS;
        MAXDrivingForceNEG = rhs.MAXDrivingForceNEG;

        Limit      = rhs.Limit;
        Limiting   = rhs.Limiting;
        Averaging  = rhs.Averaging;

        Force = rhs.Force;
    }
    return *this;
}

NodeDF DrivingForce::Force_at(const double x, const double y, const double z) const
{
#ifdef DEBUG
    if(x > Grid.Nx + Force.Bcells()*Grid.dNx - 1 or
       y > Grid.Ny + Force.Bcells()*Grid.dNy - 1 or
       z > Grid.Nz + Force.Bcells()*Grid.dNz - 1 or
       x < -Force.Bcells()*Grid.dNx or
       y < -Force.Bcells()*Grid.dNy or
       z < -Force.Bcells()*Grid.dNz)
    {
        std::stringstream message;
        message << "ERROR: DrivingForce::Force_at()\n"
                << "Access beyond storage range -> ("
                << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                << -Force.Bcells()*Grid.dNx << ", " << -Force.Bcells()*Grid.dNy << ", " << -Force.Bcells()*Grid.dNz << "] and ("
                << Grid.Nx + Force.Bcells()*Grid.dNx << ", " << Grid.Ny + Force.Bcells()*Grid.dNy << ", " << Grid.Nz + Force.Bcells()*Grid.dNz << ")\n"
                << "Terminating!!!\n";
        throw std::logic_error(message.str());
    }
#endif

    long int x0 = floor(x)*Grid.dNx;
    long int y0 = floor(y)*Grid.dNy;
    long int z0 = floor(z)*Grid.dNz;
    double dx = fabs(x - x0)*Grid.dNx;
    double dy = fabs(y - y0)*Grid.dNy;
    double dz = fabs(z - z0)*Grid.dNz;

    NodeDF loc_dG;

    double weight_x0y0z0 = ((1.0 - dx)*(1.0 - dy)*(1.0 - dz));
    for(auto it  = Force(x0,y0,z0).cbegin();
             it != Force(x0,y0,z0).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x0y0z0,
                                               it->tmp*weight_x0y0z0,
                                               it->average*weight_x0y0z0,
                                               weight_x0y0z0);
    }

    double weight_x1y0z0 = ((dx)*(1.0 - dy)*(1.0 - dz));
    for(auto it  = Force(x0+1,y0,z0).cbegin();
             it != Force(x0+1,y0,z0).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x1y0z0,
                                               it->tmp*weight_x1y0z0,
                                               it->average*weight_x1y0z0,
                                               weight_x1y0z0);
    }

    double weight_x0y1z0 = ((1.0 - dx)*(dy)*(1.0 - dz));
    for(auto it  = Force(x0,y0+1,z0).cbegin();
             it != Force(x0,y0+1,z0).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x0y1z0,
                                               it->tmp*weight_x0y1z0,
                                               it->average*weight_x0y1z0,
                                               weight_x0y1z0);
    }

    double weight_x0y0z1 = ((1.0 - dx)*(1.0 - dy)*(dz));
    for(auto it  = Force(x0,y0,z0+1).cbegin();
             it != Force(x0,y0,z0+1).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x0y0z1,
                                               it->tmp*weight_x0y0z1,
                                               it->average*weight_x0y0z1,
                                               weight_x0y0z1);
    }

    double weight_x1y1z0 = ((dx)*(dy)*(1.0 - dz));
    for(auto it  = Force(x0+1,y0+1,z0).cbegin();
             it != Force(x0+1,y0+1,z0).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x1y1z0,
                                               it->tmp*weight_x1y1z0,
                                               it->average*weight_x1y1z0,
                                               weight_x1y1z0);
    }

    double weight_x1y0z1 = ((dx)*(1.0 - dy)*(dz));
    for(auto it  = Force(x0+1,y0,z0+1).cbegin();
             it != Force(x0+1,y0,z0+1).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x1y0z1,
                                               it->tmp*weight_x1y0z1,
                                               it->average*weight_x1y0z1,
                                               weight_x1y0z1);
    }

    double weight_x0y1z1 = ((1.0 - dx)*(dy)*(dz));
    for(auto it  = Force(x0,y0+1,z0+1).cbegin();
             it != Force(x0,y0+1,z0+1).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x0y1z1,
                                               it->tmp*weight_x0y1z1,
                                               it->average*weight_x0y1z1,
                                               weight_x0y1z1);
    }

    double weight_x1y1z1 = ((dx)*(dy)*(dz));
    for(auto it  = Force(x0+1,y0+1,z0+1).cbegin();
             it != Force(x0+1,y0+1,z0+1).cend(); ++it)
    {
        loc_dG.add_all(it->indexA, it->indexB, it->raw*weight_x1y1z1,
                                               it->tmp*weight_x1y1z1,
                                               it->average*weight_x1y1z1,
                                               weight_x1y1z1);
    }
    // normalize entries
    for(auto it  = loc_dG.begin();
             it != loc_dG.end(); ++it)
    {
        if(it->weight > DBL_EPSILON)
        {
            double weight_1 = 1.0/it->weight;

            it->raw     *= weight_1;
            it->tmp     *= weight_1;
            it->average *= weight_1;
        }
    }

    return loc_dG;
}

double DrivingForce::GetDrivingForce(PhaseField& Phi,
                                     const int i, const int j, const int k,
                                     const size_t alpha, const size_t beta) const
{
    double tempdG = 0.0;
    for(auto it1  = Phi.Fields(i,j,k).cbegin();
             it1 != Phi.Fields(i,j,k).cend(); ++it1)
    for(auto it2  = Phi.Fields(i,j,k).cbegin();
             it2 != Phi.Fields(i,j,k).cend(); ++it2)
    if((Phi.FieldsProperties[it1->index].Phase == alpha)
    and(Phi.FieldsProperties[it2->index].Phase == beta))
    {
        tempdG += Force(i,j,k).get_average(it1->index, it2->index);
    }
    return tempdG;
}

void DrivingForce::AverageGlobal(PhaseField& Phase, double time)
{
    AverageDG2.clear();
    AverageDG2 = AverageDG1;
    AverageDG1.clear();
    AverageDG1 = AverageDG;
    AverageDG.clear();
    STORAGE_LOOP_BEGIN(i,j,k,Force,0)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Force(i,j,k).begin();
                 it != Force(i,j,k).end(); ++it)
        {
            size_t pIndexA = Phase.FieldsProperties[it->indexA].Phase;
            size_t pIndexB = Phase.FieldsProperties[it->indexB].Phase;

            double weight = sqrt(Phase.Fractions(i,j,k, {pIndexA})*Phase.Fractions(i,j,k, {pIndexB}));

            AverageDG.add_raw(it->indexA, it->indexB, weight*Force(i, j, k).get_raw(it->indexA, it->indexB));
            AverageDG.add_weight(it->indexA,it->indexB, weight);
        }
    }
    STORAGE_LOOP_END
    std::ofstream outfile("AvDg.dat", std::ofstream::app);
    outfile << time << " ";
    for(auto it  = AverageDG.begin();
             it != AverageDG.end(); ++it)
    {
        it->raw /= it->weight;
    }
    for (size_t i = 0; i < Nphases; ++i)
    for (size_t j = 0; j < Nphases; ++j)
    {
        outfile << AverageDG.get_raw(i,j) << " ";
    }
    outfile << std::endl;
    outfile.close();
}

}// namespace openphase

