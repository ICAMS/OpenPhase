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
 *   Main contributors :   Raphael Schiedung
 *
 */

#include "Includes.h"
#include "InterfaceProperties.h"
#include "Settings.h"
#include "InterfaceDiffusion.h"
#include "PhaseField.h"
#include "VTK.h"

namespace openphase
{
using namespace std;

InterfaceDiffusion::InterfaceDiffusion(Settings& locSettings,
                                       const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void InterfaceDiffusion::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "InterfaceDiffusion";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    Grid    = locSettings.Grid;

    Coefficients.Allocate(Nphases, Nphases);
    DiffusionPotential.Allocate(Grid, Grid.Bcells);

    switch(locSettings.Grid.Active())
    {
        case 1:
        {
            LStencil.Set(LaplacianStencil1D_3, Grid);
            break;
        }
        case 2:
        {
            LStencil.Set(LaplacianStencil2D_9, Grid);
            break;
        }
        case 3:
        {
            LStencil.Set(LaplacianStencil3D_27a, Grid);
            break;
        }
    }

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void InterfaceDiffusion::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("InterfaceDiffusion input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);
    inp.close();
}

void InterfaceDiffusion::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    DoubleObstacleSmoothnessRange =
        FileInterface::ReadParameterD(inp, moduleLocation, string("dSR"));

    for(unsigned int alpha = 0;     alpha < Nphases; alpha++)
    for(unsigned int beta  = alpha; beta  < Nphases;  beta++)
    {
        double value = 0.0;
        stringstream converter;
        converter << alpha << "_" << beta;
        string counter = converter.str();
        value = FileInterface::ReadParameterD(inp, moduleLocation, string("dIDC_") + counter);
        Coefficients(alpha, beta) = value;
        Coefficients(beta, alpha) = value;
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void InterfaceDiffusion::CalculatePhaseFieldIncrements(PhaseField& Phase,
        const InterfaceProperties& IP)
{
    const unsigned int PhaseFieldBCells = Phase.Fields.Bcells();
    if (PhaseFieldBCells <= 2)
    {
        ConsoleOutput::WriteExit("Not enough PhaseField boundary cells (minimum 3 cells!)",
                thisclassname, "Calculate");
        OP_Exit(EXIT_FAILURE);
    }

    const unsigned int InterfaceEnergyBCells = IP.Properties.Bcells();
    if (InterfaceEnergyBCells <= 1)
    {
        ConsoleOutput::WriteExit("Not enough InterfaceEnergy boundary cells (minimum 2 cells!)",
                thisclassname, "Calculate");
        OP_Exit(EXIT_FAILURE);
    }

    // Override flag functionality
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 3,)
        Phase.Fields(i,j,k).flag = 2;
    OMP_PARALLEL_STORAGE_LOOP_END

    CalculateDiffusionPotential          (Phase, IP);
    CalculateDiffusionPotentialLaplacian (Phase);
}

void InterfaceDiffusion::CalculateDiffusionPotentialLaplacian(
        PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 0,)
    {
        for (auto ls = LStencil.begin(); ls != LStencil.end(); ls++)
        {
            const int ii = ls->di;
            const int jj = ls->dj;
            const int kk = ls->dk;

            for (auto it  = DiffusionPotential(i+ii,j+jj,k+kk).cbegin();
                      it != DiffusionPotential(i+ii,j+jj,k+kk).cend(); ++it)
            {
                const int pIndexA = Phase.FieldsProperties[it->indexA].Phase;
                const int pIndexB = Phase.FieldsProperties[it->indexB].Phase;

                double value = - Coefficients(pIndexA,pIndexB)*ls->weight*it->value1;
                Phase.FieldsDot(i,j,k).add_asym1(it->indexA, it->indexB, value);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 1,)
        DiffusionPotential(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}

double InterfaceDiffusion::PotentialDerivative(const double alpha,
        const double beta) const
{
    double value = 0.0;
    // This method calulates the
    if (DoubleObstacleSmoothnessRange == 0.0)
    {
        // Normal double-obstacle potential
        double Pot = alpha * beta;
        if (Pot < 0.0)
            value = - beta;
        else if (Pot > 0.0)
            value = beta;
    }
    else if (DoubleObstacleSmoothnessRange == 1.0)
    {
        // Normal double-well potential
        value =  2.0 * alpha * pow(beta,2);
    }
    else
    {
        // Bent-cable model for the potential
        // Calculate bent-cable function of phi-beta
        double BentCableBeta = 0.0;
        if (abs(beta) < DoubleObstacleSmoothnessRange)
        {
            BentCableBeta -= 1.0 * pow(beta,4)
                / (16.0 * pow(DoubleObstacleSmoothnessRange,3));
            BentCableBeta += 3.0 * pow(beta,2)
                / ( 8.0 * pow(DoubleObstacleSmoothnessRange,1));
            BentCableBeta += 1.0 * pow(beta,1) /   2.0;
            BentCableBeta += 3.0 * DoubleObstacleSmoothnessRange / 16.0;
        }
        else if (beta >= DoubleObstacleSmoothnessRange) BentCableBeta = beta;

        // Calculate smoothed absolute value of phi-alpha
        double SAbsBeta  = - beta  + 2.0 * BentCableBeta;

        // Calculate bent-cable function of phi-alpha
        double dBentCableAlpha_dAlpha = 0.0;
        if (abs(alpha) < DoubleObstacleSmoothnessRange)
        {
            dBentCableAlpha_dAlpha -= 1.0 * pow(alpha,3)
                / (4.0 * pow(DoubleObstacleSmoothnessRange,3));
            dBentCableAlpha_dAlpha += 6.0 * pow(alpha,1)
                / (8.0 * pow(DoubleObstacleSmoothnessRange,1));
            dBentCableAlpha_dAlpha += 0.5;
        }
        else if (alpha >= DoubleObstacleSmoothnessRange)
            dBentCableAlpha_dAlpha = 1.0;

        // Calculate smoothed absolute value of phi-alpha
        double dSAbsAlpha_dAlpha = - 1.0 + 2.0 * dBentCableAlpha_dAlpha;

        return dSAbsAlpha_dAlpha * SAbsBeta;
    }
    return value;
}

void InterfaceDiffusion::CalculateDiffusionPotential(PhaseField& Phase,
        const InterfaceProperties& IP)
{
    const double Prefactor2 = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 2,)
    {
        double norm_1 = 1.0/double(Phase.Fields(i,j,k).size());

        const NodePF& locPF = Phase.Fields(i,j,k);

        for(auto alpha = locPF.cbegin();
                 alpha < locPF.cend(); ++alpha)
        for(auto  beta = alpha + 1;
                  beta < locPF.cend(); ++beta)
        {
            double dDiffusionPotential_dt =
                IP.Properties(i,j,k).get_energy(alpha->index, beta->index) *
                ((alpha->laplacian + Prefactor2 *
                  PotentialDerivative(beta->value,  alpha->value)) -
                 (beta->laplacian + Prefactor2 *
                  PotentialDerivative(alpha->value, beta->value)));

            if(locPF.size() > 2)
            for(auto gamma = locPF.cbegin();
                     gamma < locPF.cend(); ++gamma)
            if((gamma != alpha) && (gamma != beta))
            {
                dDiffusionPotential_dt +=
                     IP.Properties(i,j,k).get_energy( beta->index, gamma->index) *
                    (gamma->laplacian + Prefactor2 *
                     PotentialDerivative(beta->value,  gamma->value)) -
                     IP.Properties(i,j,k).get_energy(alpha->index, gamma->index) *
                    (gamma->laplacian + Prefactor2 *
                     PotentialDerivative(alpha->value, gamma->value));
            }

            dDiffusionPotential_dt *= norm_1;

            DiffusionPotential(i,j,k).add_asym1(alpha->index, beta->index,
                    dDiffusionPotential_dt);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

InterfaceDiffusionAnisotropic::InterfaceDiffusionAnisotropic(Settings& locSettings,
                                                    const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void InterfaceDiffusionAnisotropic::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    InterfaceDiffusion::Initialize(locSettings, ObjectNameSuffix);

    thisclassname = "InterfaceDiffusionAnisotropic";
    thisobjectname = thisclassname + ObjectNameSuffix;

    DiffusionFlux.Allocate (Grid, Grid.Bcells);
}

void InterfaceDiffusionAnisotropic::ReadInput(const std::string InputFileName)
{
    InterfaceDiffusion::ReadInput(InputFileName);
    if (DoubleObstacleSmoothnessRange < 0.04)
    {
        ConsoleOutput::WriteWarning("DoubleObstacleSmoothnessRange might be too small to the algorithm to work!",
                thisclassname, "ReadInput");
    }
}

void InterfaceDiffusionAnisotropic::CalculatePhaseFieldIncrements(PhaseField& Phase,
        const InterfaceProperties& IP)
{
    const int PhaseFieldBCells = Phase.Fields.Bcells();
    if (PhaseFieldBCells <= 2)
    {
        ConsoleOutput::WriteExit("Not enough PhaseField boundary cell (min 3 cells!)",
                thisclassname, "CalculateAnsiotropic");
        exit(EXIT_FAILURE);
    }

    const int InterfaceEnergyBCells = IP.Properties.Bcells();
    if (InterfaceEnergyBCells <= 1)
    {
        ConsoleOutput::WriteExit("Not enough InterfaceEnergy boundary cell (min 2 cells!)",
                thisclassname, "CalculateAnsiotropic");
        exit(EXIT_FAILURE);
    }

    // Override flag functionality
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 3,)
        Phase.Fields(i,j,k).flag = 2;
    OMP_PARALLEL_STORAGE_LOOP_END

    // This method solves the anisotropic diffusion-equation
    // PhiDot = Nabla ( Diffusion-Tensor * Nabla( Diffusion-potential))
    CalculateDiffusionPotential          (Phase, IP);
    CalculateDiffusionPotentialGradients ();
    CalculateDiffusionFlux               (Phase);
    CalculateDiffusionFluxDivergence     (Phase);
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionFlux(PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, 1,)
    {
        NodeAB<dVector3,dVector3> locNormals       = Phase.Normals(i,j,k);
        NodeA<dVector3>           locGradients     = Phase.Fields(i,j,k).get_gradients();
        NodeAB<dVector3,dVector3> locDiffusionFlux = DiffusionFlux(i,j,k);

        for (auto it = locDiffusionFlux.cbegin();
                  it < locDiffusionFlux.cend(); ++ it)
        {
            dMatrix3x3 locDT;
            locDT.set_to_unity();

            // Calculation of the anisotropy:
            // The ratio of tangential to normal Interface Diffusion will be
            // calculated
            //TODO Note: The bent-cable potential is not used !!
            double Ratio     = 0.0; // Ration of normal to tangential diffusion
            double Potential = 0.0; // Stores value of double obstacle potential
            Potential  = abs(Phase.Fields(i,j,k).get_value(it->indexA));
            Potential *= abs(Phase.Fields(i,j,k).get_value(it->indexB));
            if (Potential > 1.0e-200)
            {
               Ratio  = abs(locGradients.get_value(it->indexA)
                          * locGradients.get_value(it->indexB));
               Ratio /= Potential;
               Ratio *= (Phase.Grid.Eta * Phase.Grid.Eta)/(Pi * Pi);
            }

            // Limit Ratio to meaningful values
            if (Ratio > 1.0) Ratio = 1.0; // Extinguishes normal diffusion completely
            if (Ratio < 0.0) Ratio = 0.0; // Results in the scalar model

            // Calculate local interface diffusion mobility tensor
            dVector3 locNormal = locNormals.get_asym1(it->indexA, it->indexB);
            locDT(0,0) -= Ratio * locNormal.getX()*locNormal.getX();
            locDT(0,1) -= Ratio * locNormal.getX()*locNormal.getY();
            locDT(0,2) -= Ratio * locNormal.getX()*locNormal.getZ();
            locDT(1,0) -= Ratio * locNormal.getY()*locNormal.getX();
            locDT(1,1) -= Ratio * locNormal.getY()*locNormal.getY();
            locDT(1,2) -= Ratio * locNormal.getY()*locNormal.getZ();
            locDT(2,0) -= Ratio * locNormal.getZ()*locNormal.getX();
            locDT(2,1) -= Ratio * locNormal.getZ()*locNormal.getY();
            locDT(2,2) -= Ratio * locNormal.getZ()*locNormal.getZ();

            // Apply diffuse interface interpolation
            int pIndexA = Phase.FieldsProperties[it->indexA].Phase;
            int pIndexB = Phase.FieldsProperties[it->indexB].Phase;
            locDT *= - Coefficients(pIndexA,pIndexB);

            // Apply Projection Operator
            dVector3 value = locDT * locDiffusionFlux.get_asym1(it->indexA,it->indexB);
            DiffusionFlux(i,j,k).set_asym1(it->indexA, it->indexB, value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionFluxDivergence(PhaseField& Phase)
{
    const double DWeights[3] = {-0.5/Grid.dx, 0.0, 0.5/Grid.dx};

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, 0,)
    for (int ii = -1; ii <= +1; ii+=2)
    {
        double DWeight = DWeights[ii+1];
        if (Grid.dNx)
        for (auto it = DiffusionFlux(i+ii,j,k).cbegin();
                  it < DiffusionFlux(i+ii,j,k).cend(); ++it)
        {
            double valueX = DWeight*(it->value1[0]);
            Phase.FieldsDot(i,j,k).add_asym1(it->indexA, it->indexB, valueX);
        }
        if (Grid.dNy)
        for (auto it = DiffusionFlux(i,j+ii,k).cbegin();
                  it < DiffusionFlux(i,j+ii,k).cend(); ++it)
        {
            double valueY = DWeight*(it->value1[1]);
            Phase.FieldsDot(i,j,k).add_asym1(it->indexA, it->indexB, valueY);
        }
        if (Grid.dNz)
        for (auto it = DiffusionFlux(i,j,k+ii).cbegin();
                  it < DiffusionFlux(i,j,k+ii).cend(); ++it)
        {
            double valueZ = DWeight*(it->value1[2]);
            Phase.FieldsDot(i,j,k).add_asym1(it->indexA, it->indexB, valueZ);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    const int BCells = DiffusionFlux.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, BCells,)
        DiffusionFlux(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionPotentialGradients()
{
    const double GWeights[3] = {-0.5/Grid.dx, 0.0, 0.5/Grid.dx};

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 1,)
        DiffusionFlux(i,j,k).clear();
        for (int ii = -1; ii <= +1; ii+=2)
        {
            double GWeight = GWeights[ii+1];
            if (Grid.dNx)
            for (auto it = DiffusionPotential(i+ii,j,k).cbegin();
                      it < DiffusionPotential(i+ii,j,k).cend(); ++it)
            {
                double value = GWeight*(it->value1);
                DiffusionFlux(i,j,k).add_asym1(it->indexA, it->indexB, (dVector3){value, 0.0, 0.0});
            }
            if (Grid.dNy)
            for (auto it = DiffusionPotential(i,j+ii,k).cbegin();
                      it < DiffusionPotential(i,j+ii,k).cend(); ++it)
            {
                double value = GWeight*(it->value1);
                DiffusionFlux(i,j,k).add_asym1(it->indexA, it->indexB, (dVector3){0.0, value, 0.0});
            }
            if (Grid.dNz)
            for (auto it = DiffusionPotential(i,j,k+ii).cbegin();
                      it < DiffusionPotential(i,j,k+ii).cend(); ++it)
            {
                double value = GWeight*(it->value1);
                DiffusionFlux(i,j,k).add_asym1(it->indexA, it->indexB, (dVector3){0.0, 0.0, value});
            }
        }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 2,)
        DiffusionPotential(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}
}
