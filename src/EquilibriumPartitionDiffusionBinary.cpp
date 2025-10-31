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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#include "EquilibriumPartitionDiffusionBinary.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Settings.h"
#include "Temperature.h"
#include "VTK.h"

namespace openphase
{

using namespace std;

inline pair<int, int> three_state_bounds_selector(int selector,
                                                  int lower_bound,
                                                  int upper_bound)
{
    pair<int, int> result;
    switch (selector)
    {
        case -1: // lower bound only
        {
            result = make_pair(lower_bound, lower_bound + 1);
            break;
        }
        case 0: // both bounds
        {
            result = make_pair(lower_bound, upper_bound);
            break;
        }
        case 1: // upper bound only
        {
            result = make_pair(upper_bound - 1, upper_bound);
            break;
        }
    }
    return result;
}

/*************************************************************************/
void EquilibriumPartitionDiffusionBinary::Initialize(Settings& locSettings,
                                                     std::string ObjectNameSuffix)
{
    thisclassname = "EquilibriumPartitionDiffusion";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;

    Grid = locSettings.Grid;

    ElementNames = locSettings.ElementNames;

    Comp = 0;
    R = PhysicalConstants::R;

    EnableAntiTrapping = true;

    DrivingForceModel = BinaryDrivingForceModels::Standard;

    size_t Bcells = Grid.Bcells;
    dMu.Allocate(Grid, {Nphases, Ncomp}, Bcells);

    Tcross.Allocate({Nphases, Nphases});
    Ccross.Allocate({Nphases, Nphases});
    Slope.Allocate({Nphases, Nphases});

    IDC.resize(Nphases,1.0);
    DC.resize(Nphases,0.0);
    DC0.resize(Nphases,0.0);
    AE.resize(Nphases,0.0);
    Stoichiometric.resize(Nphases,0);
    Entropy.resize(Nphases,0);
    DiffusionStencil = LaplacianStencils::Isotropic;

    switch(locSettings.Grid.Active())
    {
        case 1:
        {
            DStencil.SetNoCenter(LaplacianStencil1D_3, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            break;
        }
        case 2:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil2D_5, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil2D_9, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            break;
        }
        case 3:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil3D_7, Grid.dx);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil3D_27a, Grid.dx);
            }
            break;
        }
        default:
        {

        }
    }
    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void EquilibriumPartitionDiffusionBinary::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("EquilibriumPartitionDiffusionBinary input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };

    ConsoleOutput::WriteLineInsert("EquilibriumPartitionDiffusion properties");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void EquilibriumPartitionDiffusionBinary::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    string RefCompName = FileInterface::ReadParameterK(inp, moduleLocation, string("RefElement"), true, "NN");
    EnableAntiTrapping = FileInterface::ReadParameterB(inp, moduleLocation, string("EnableAntiTrapping"), false, true);
    ConsiderChemicalPotential = FileInterface::ReadParameterB(inp, moduleLocation, string("ConsiderChemicalPotential"), false, false);

    std::string dGmodel = FileInterface::ReadParameterK(inp, moduleLocation,"DrivingForceModel", false, "STANDARD");
    if (dGmodel == "STANDARD")
    {
        DrivingForceModel = BinaryDrivingForceModels::Standard;
    }
    else if (dGmodel == "LOWERSLOPE")
    {
        DrivingForceModel = BinaryDrivingForceModels::LowerSlope;
    }
    else if (dGmodel == "WEIGHTED")
    {
        DrivingForceModel = BinaryDrivingForceModels::Weighted;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong driving force model specified, the standard model used instead!", thisclassname, "ReadInput()");
        DrivingForceModel = BinaryDrivingForceModels::Standard;
    }

    bool comp_set = false;
    for(size_t n = 0; n < Ncomp; n++)
    if(ElementNames[n] == RefCompName)
    {
        RefComp = n;
        comp_set = true;
    }
    else
    {
        Comp = n;
    }
    if(!comp_set)
    {
        ConsoleOutput::WriteExit("Nonexistent reference element is selected", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }
    for(size_t n = 0; n < Nphases; ++n)
    {
        stringstream converter;
        converter << n;
        string counter = converter.str();
        Stoichiometric[n] = FileInterface::ReadParameterB(inp, moduleLocation, string("Flag_") + counter);
        IDC[n]            = FileInterface::ReadParameterD(inp, moduleLocation, string("IDC_") + counter, false, 1.0);
        DC0[n]            = FileInterface::ReadParameterD(inp, moduleLocation, string("DC_") + counter);
        AE[n]             = FileInterface::ReadParameterD(inp, moduleLocation, string("AE_") + counter);
        Entropy[n]        = FileInterface::ReadParameterD(inp, moduleLocation, string("EF_") + counter);
    }

    for(size_t n =   0; n < Nphases-1; ++n)
    for(size_t m = n+1; m < Nphases  ; ++m)
    {
        stringstream converter;
        converter << n << "_" << m;
        string counter = converter.str();
        Ccross({n, m}) = Ccross({m, n}) = FileInterface::ReadParameterD(inp, moduleLocation, string("Cs_") + counter);
        Tcross({n, m}) = Tcross({m, n}) = FileInterface::ReadParameterD(inp, moduleLocation, string("Ts_") + counter);
        Slope({n, m})  = FileInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter);
        stringstream converter2;
        converter2 << m << "_" << n;
        string counter2 = converter2.str();
        Slope({m, n})  = FileInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter2);
    }

    for(size_t n =   0; n < Nphases; ++n)
    {
        Slope({n, n}) = 1.0;
    }

    string tmp1 = FileInterface::ReadParameterK(inp, moduleLocation, string("DiffusionStencil"), false, string("ISOTROPIC"));
    if(tmp1 == "SIMPLE")
    {
        DiffusionStencil = LaplacianStencils::Simple;
    }
    else if(tmp1 == "ISOTROPIC")
    {
        DiffusionStencil = LaplacianStencils::Isotropic;
    }
    else if(tmp1 == "LB")
    {
        DiffusionStencil = LaplacianStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong diffusion stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    switch(Grid.Active())
    {
        case 1:
        {
            DStencil.SetNoCenter(LaplacianStencil1D_3, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            break;
        }
        case 2:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil2D_5, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil2D_9, Grid.dx,Grid.dNx,Grid.dNy,Grid.dNz);
            }
            break;
        }
        case 3:
        {
            if(DiffusionStencil == LaplacianStencils::Simple)
            {
                DStencil.SetNoCenter(LaplacianStencil3D_7, Grid.dx);
            }
            else
            {
                DStencil.SetNoCenter(LaplacianStencil3D_27a, Grid.dx);
            }
            break;
        }
        default:
        {

        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void EquilibriumPartitionDiffusionBinary::CalculateDiffusionCoefficients(Temperature& Tx)
{
    double invRT  = 1.0/(R*Tx.Tavg);
    maxDC = 0.0;

    for (size_t n = 0; n < Nphases; ++n)
    {
        DC[n] = DC0[n] * exp(-AE[n] * invRT);
        maxDC = max(DC[n], maxDC);
    }
#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &maxDC, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
}

double EquilibriumPartitionDiffusionBinary::ReportMaximumTimeStep(Temperature& Tx)
{
    CalculateDiffusionCoefficients(Tx);

    double max_dt = DBL_MAX;
    if(maxDC != 0.0)
    {
        switch(Grid.Active())
        {
            case 1:
            {
                max_dt = 0.50*Grid.dx*Grid.dx/maxDC;
                break;
            }
            case 2:
            {
                max_dt = 0.25*Grid.dx*Grid.dx/maxDC;
                break;
            }
            case 3:
            {
                max_dt = 0.16*Grid.dx*Grid.dx/maxDC;
                break;
            }
        }
    }
    return max_dt;
}

double EquilibriumPartitionDiffusionBinary::EquilibriumComposition(size_t n,
                                                                   size_t m,
                                                                   double Temp)
{
    double eqC = Ccross({n, m});
    if(Slope({n, m}) != 0.0)
    {
        eqC += (Temp - Tcross({n, m}))/Slope({n, m});
    }
    eqC = std::clamp(eqC, 0.0, 1.0);
    return eqC;
}

double EquilibriumPartitionDiffusionBinary::PartitioningCoefficient(size_t n,
                                                                    size_t m)
{
    double pt_coefficient = 0.0;
    if(Slope({m, n}) != 0.0)
    {
        pt_coefficient = fabs(Slope({n, m}) / Slope({m, n}));
    }
    return pt_coefficient;
}

void EquilibriumPartitionDiffusionBinary::CalculateLocalPhaseConcentrations(
                                                        PhaseField& Phase,
                                                        Composition& Cx,
                                                        Temperature& Tx,
                                                        int i, int j, int k)
{
    if(Phase.Fields(i,j,k).bulk())
    {
        Cx.MoleFractions(i,j,k).set_to_zero();
        size_t pIndex = Phase.FieldsProperties[Phase.Fields(i,j,k).front().index].Phase;

        Cx.MoleFractions(i,j,k,{pIndex,   Comp}) = Cx.MoleFractionsTotal(i,j,k,{   Comp});
        Cx.MoleFractions(i,j,k,{pIndex,RefComp}) = Cx.MoleFractionsTotal(i,j,k,{RefComp});
    }
    else
    {
        vector<double> partitioning_reduction(Nphases, 0.0);
        bool iterate = false;

        size_t MaxIterations = 10;
        size_t iteration = 0;

        do
        {
            iterate = false;

            Cx.MoleFractions(i,j,k).set_to_zero();

            // Calculate equilibrium concentrations
            for(size_t n = 0; n < Nphases; ++n)
            {
                double eqCx    = 0.0;
                double divisor = 0.0;
                for(size_t m = 0; m < Nphases; ++m)
                if(n != m and Phase.Fractions(i,j,k,{m}) != 0.0)
                {
                    eqCx += Phase.Fractions(i,j,k,{m})*EquilibriumComposition(n,m,Tx(i,j,k))*(1.0 - partitioning_reduction[n]);
                    if(partitioning_reduction[n] > DBL_EPSILON)
                    {
                        eqCx += Phase.Fractions(i,j,k,{m})*Cx.MoleFractionsTotal(i,j,k,{Comp})*partitioning_reduction[n];
                    }
                    divisor += Phase.Fractions(i,j,k,{m});
                }

                if(divisor != 0.0)
                {
                    Cx.MoleFractions(i,j,k,{n, Comp}) = eqCx / divisor;
                }
                else
                {
                    Cx.MoleFractions(i,j,k,{n, Comp}) = Cx.MoleFractionsTotal(i,j,k,{Comp});
                }
            }

            // Calculate concentrations delta
            double dCx = Cx.MoleFractionsTotal(i,j,k,{Comp});
            for(size_t n = 0; n < Nphases; ++n)
            {
                dCx -= Cx.MoleFractions(i,j,k,{n, Comp})*Phase.Fractions(i,j,k,{n});
            }

            // Calculate phase concentrations
            for(size_t n = 0; n < Nphases; ++n)
            if (!Stoichiometric[n])
            {
                double numerator   = Phase.Fractions(i,j,k,{n});
                double denominator = Phase.Fractions(i,j,k,{n});

                for(size_t m = 0; m < Nphases; ++m)
                if(n != m and !Stoichiometric[m])
                {
                    numerator   += Phase.Fractions(i,j,k,{m});
                    denominator += Phase.Fractions(i,j,k,{m})*
                                      (partitioning_reduction[n]*1.0 +
                                      (1.0 - partitioning_reduction[n])*PartitioningCoefficient(n, m));
                }
                if(denominator != 0.0)
                {
                    numerator /= denominator;
                }
                else
                {
                    numerator = 0.0;
                }
                Cx.MoleFractions(i,j,k,{n, Comp}) += dCx*numerator;
            }

//            // Calculate concentration deviation
//            double deviation = Cx.MoleFractionsTotal(i,j,k,{Comp});
//            for(size_t n = 0; n < Cx.Nphases; n++)
//            {
//                deviation -= Cx.MoleFractions(i,j,k,{n, Comp})*Phase.Fractions(i,j,k,{n});
//            }
//
//            // Check and adjust mass conservation
//            if (fabs(deviation) > 10.0*DBL_EPSILON)
//            {
//                //cout << deviation << endl;
//                double deltaCstoich = 0.0;
//                for(size_t n = 0; n < Cx.Nphases; n++)
//                if(Stoichiometric[n])
//                {
//                    deltaCstoich += Cx.MoleFractions(i,j,k,{n,Comp})*Phase.Fractions(i,j,k,{n});
//                }
//                double c_total_sans_stoich = Cx.MoleFractionsTotal(i,j,k,{Comp}) - deltaCstoich;
//                double correction_factor = c_total_sans_stoich - deviation;
//                for(size_t n = 0; n < Cx.Nphases; n++)
//                if(Phase.Fractions(i,j,k,{n}) != 0.0 and
//                   !Stoichiometric[n] and
//                   fabs(correction_factor) > 10.0*DBL_EPSILON)
//                {
//                    Cx.MoleFractions(i,j,k,{n, Comp}) *= c_total_sans_stoich/correction_factor;
//                }
//            }

            // Reduce partitioning coefficients if phase concentrations go out of range
            for(size_t n = 0; n < Nphases; n++)
            if(Phase.Fractions(i,j,k,{n}) != 0.0 and
               (Cx.MoleFractions(i,j,k,{n, Comp}) < Cx.MinMoleFraction(n, Comp) or
                Cx.MoleFractions(i,j,k,{n, Comp}) > Cx.MaxMoleFraction(n, Comp)))
            {
                partitioning_reduction[n] += 0.1;
                iterate = true;
            }
            iteration++;
        }
        while(iterate and iteration < MaxIterations);

        // Set reference element concentrations
        for(size_t n = 0; n < Nphases; ++n)
        {
            Cx.MoleFractions(i,j,k,{n,RefComp}) = 1.0 - Cx.MoleFractions(i,j,k,{n,Comp});
        }

#ifdef DEBUG
        if (iterate)
        {
            ConsoleOutput::WriteWarning(thisclassname,"CalculateLocalPhaseConcentrations","Calculation of mole fractions did not converge");
            for(size_t n = 0; n < Nphases; n++)
            {
                std::stringstream ss;
                ss << "Cx.MoleFractions(" << i << "," << j << "," << k << ",{" << n << "," << ElementNames[Comp] << "})";
                ConsoleOutput::Write(ss.str(),Cx.MoleFractions(i,j,k,{n, Comp}));
            }
            ConsoleOutput::Write("");
        }
#endif
    }
}

void EquilibriumPartitionDiffusionBinary::CalculatePhaseConcentrations(PhaseField& Phase,
                                                                       Composition& Cx,
                                                                       Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,Cx.MoleFractionsTotal.Bcells(),)
    {
        CalculateLocalPhaseConcentrations(Phase, Cx, Tx, i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::CalculateDrivingForce(PhaseField& Phase,
                                                                Composition& Cx,
                                                                Temperature& Tx,
                                                                DrivingForce& dG)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            CalculateLocalPhaseConcentrations(Phase, Cx, Tx, i, j, k);

            for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                if(pIndexA != pIndexB)
                {
                    double counter = 0.0;

                    double locdG_AB  = 0.0;
                    if(not Stoichiometric[pIndexA])
                    {
                        locdG_AB = (Tcross({pIndexA, pIndexB}) + Slope({pIndexA, pIndexB})*
                                   (Cx.MoleFractions(i,j,k,{pIndexA, Comp}) - Ccross({pIndexA, pIndexB})) - Tx(i,j,k))*(Entropy[pIndexB] - Entropy[pIndexA]);
                        counter++;
                    }
                    double locdG_BA  = 0.0;
                    if(not Stoichiometric[pIndexB])
                    {
                        locdG_BA = (Tcross({pIndexB, pIndexA}) + Slope({pIndexB, pIndexA})*
                                   (Cx.MoleFractions(i,j,k,{pIndexB, Comp}) - Ccross({pIndexB, pIndexA})) - Tx(i,j,k))*(Entropy[pIndexB] - Entropy[pIndexA]);
                        counter++;
                    }

                    double dG_AB = 0.0;

                    switch(DrivingForceModel)
                    {
                        case BinaryDrivingForceModels::Standard:
                        {
                            if(counter != 0.0)
                            {
                                dG_AB = (locdG_AB + locdG_BA)/counter;
                            }
                            break;
                        }
                        case BinaryDrivingForceModels::LowerSlope:
                        {
                            if (std::fabs(Slope({pIndexA, pIndexB})) < std::fabs(Slope({pIndexB, pIndexA})))
                            {
                                dG_AB = locdG_AB;
                            }
                            else
                            {
                                dG_AB = locdG_BA;
                            }
                            break;
                        }
                        case BinaryDrivingForceModels::Weighted:
                        {
                            double SlopeAB_1 = (Slope({pIndexA, pIndexB}) != 0.0) ? 1.0/Slope({pIndexA, pIndexB}) : 0.0;
                            double SlopeBA_1 = (Slope({pIndexB, pIndexA}) != 0.0) ? 1.0/Slope({pIndexB, pIndexA}) : 0.0;

                            double total = std::fabs(SlopeAB_1) +
                                           std::fabs(SlopeBA_1);

                            if(total > DBL_EPSILON)
                            {
                                dG_AB = (std::fabs(SlopeAB_1)*locdG_AB +
                                         std::fabs(SlopeBA_1)*locdG_BA)/total;
                            }
                            break;
                        }
                    }
                    dG.Force(i,j,k).add_raw(alpha->index, beta->index, dG_AB);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::Remesh(int newNx,
                                                 int newNy,
                                                 int newNz,
                                                 const BoundaryConditions& BC)
{
    dMu.Reallocate(newNx, newNy, newNz);

    Grid.SetDimensions(newNx, newNy, newNz);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed!");
}

void EquilibriumPartitionDiffusionBinary::CalculateInterfaceMobility(PhaseField& Phase,
                                                                     Temperature& Tx,
                                                                     BoundaryConditions& BC,
                                                                     InterfaceProperties& IP)
{
    CalculateDiffusionCoefficients(Tx);

    Matrix<double> IntMob(Nphases,Nphases);

    for(size_t pIndexA = 0; pIndexA < Nphases; pIndexA++)
    for(size_t pIndexB = 0; pIndexB < Nphases; pIndexB++)
    {
        if(IP.InterfaceMobility(pIndexA,pIndexB).Model == InterfaceMobilityModels::Ext)
        {
            double mob = 0.0;
            if(!Stoichiometric[pIndexA] or !Stoichiometric[pIndexB])
            {
                double deltaS = (Entropy[pIndexA] - Entropy[pIndexB]);
                double deltaC = 0.0;

                if(Slope({pIndexA, pIndexB}) != 0.0)
                {
                    deltaC += (Tx.Tavg - Tcross({pIndexA, pIndexB}))/Slope({pIndexA, pIndexB});
                }

                if(Slope({pIndexB, pIndexA}) != 0.0)
                {
                    deltaC -= (Tx.Tavg - Tcross({pIndexB, pIndexA}))/Slope({pIndexB, pIndexA});
                }

                double Slope_tmp = 0.0;
                if(!Stoichiometric[pIndexA])
                {
                    Slope_tmp = min(fabs(Slope({pIndexA, pIndexB})), fabs(Slope({pIndexB, pIndexA})));
                }
                else
                {
                    Slope_tmp = Slope({pIndexB, pIndexA});
                }

                double denominator = Slope_tmp*Phase.Grid.Eta*deltaS*deltaC;

                if(denominator != 0.0)
                {
                    mob = fabs(8.0*(DC[pIndexA] + DC[pIndexB])/denominator);
                }
            }
            IntMob(pIndexA, pIndexB) = mob;
            IntMob(pIndexB, pIndexA) = mob;
        }
    }

    switch(Phase.Grid.Resolution)
    {
        case Resolutions::Single:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP.Properties,0,)
            {
                if(Phase.Fields(i,j,k).wide_interface())
                {
                    for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                             alpha != Phase.Fields(i,j,k).cend(); ++alpha)
                    for(auto  beta  = alpha + 1;
                              beta != Phase.Fields(i,j,k).cend(); ++beta)
                    {
                        size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                        size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                        if(IP.InterfaceMobility(pIndexA,pIndexB).Model == InterfaceMobilityModels::Ext)
                        {
                            double mob = IntMob(pIndexA, pIndexB);
                            IP.Properties(i,j,k).set_mobility(alpha->index, beta->index, mob);
                        }
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        case Resolutions::Dual:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,IP.PropertiesDR,0,)
            {
                if(Phase.FieldsDR(i,j,k).wide_interface())
                {
                    for(auto alpha  = Phase.FieldsDR(i,j,k).cbegin();
                             alpha != Phase.FieldsDR(i,j,k).cend(); ++alpha)
                    for(auto  beta  = alpha + 1;
                              beta != Phase.FieldsDR(i,j,k).cend(); ++beta)
                    {
                        size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;
                        size_t pIndexB = Phase.FieldsProperties[ beta->index].Phase;

                        if(IP.InterfaceMobility(pIndexA,pIndexB).Model == InterfaceMobilityModels::Ext)
                        {
                            double mob = IntMob(pIndexA, pIndexB);
                            IP.PropertiesDR(i,j,k).set_mobility(alpha->index, beta->index, mob);
                        }
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            IP.SetBoundaryConditionsDR(BC);
            IP.Coarsen(Phase);
            break;
        }
    }
    IP.SetBoundaryConditions(BC);
}

void EquilibriumPartitionDiffusionBinary::CalculateAntitrappingIncrements(PhaseField& Phase,
                                                                          Composition& Cx)
{
    const double coef = Phase.Grid.Eta/Pi;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;
                if(Phase.Fields(i+di, j+dj, k+dk).interface())
                {
                    NodePF locPF2 = Phase.Fields(i,j,k);
                    locPF2.add_existing_all(Phase.Fields(i+di, j+dj, k+dk));

                    if(locPF2.size() > 1)
                    {
                        locPF2 *= 0.5;
                        dVector3 StencilDir{double(di),double(dj),double(dk)};
                        StencilDir.normalize();

                        for(auto alpha  = locPF2.cbegin();
                                 alpha != locPF2.cend(); ++alpha)
                        {
                            size_t pIndexA = Phase.FieldsProperties[alpha->index].Phase;

                            double locDelta = 0.0;
                            double locMoleFractionA = 0.5*(Cx.MoleFractions(i   , j   , k   ,{pIndexA, Comp}) +
                                                           Cx.MoleFractions(i+di, j+dj, k+dk,{pIndexA, Comp}));

                            for(auto beta  = locPF2.cbegin();
                                     beta != locPF2.cend(); ++beta)
                            if(alpha != beta)
                            {
                                size_t pIndexB = Phase.FieldsProperties[beta->index].Phase;
                                if (pIndexA != pIndexB and fabs(DC0[pIndexA] - DC0[pIndexB]) > DBL_EPSILON)
                                {
                                    double locMoleFractionB = 0.5*(Cx.MoleFractions(i   , j   , k   ,{pIndexB, Comp}) +
                                                                   Cx.MoleFractions(i+di, j+dj, k+dk,{pIndexB, Comp}));

                                    double scale = 1.0;

                                    if(Phase.FieldsProperties[alpha->index].Stage != GrainStages::Stable)
                                    {
                                        scale *= Phase.FieldsProperties[alpha->index].VolumeRatio;
                                    }
                                    if(Phase.FieldsProperties[ beta->index].Stage != GrainStages::Stable)
                                    {
                                        scale *= Phase.FieldsProperties[ beta->index].VolumeRatio;
                                    }

                                    dVector3 locNormal = Phase.Normal(alpha, beta);

                                    double Projection = locNormal*StencilDir;

                                    double locPsi  = 0.5*(Phase.FieldsDot(i   , j   , k   ).get_asym1(alpha->index, beta->index) +
                                                          Phase.FieldsDot(i   , j   , k   ).get_asym2(alpha->index, beta->index) +
                                                          Phase.FieldsDot(i+di, j+dj, k+dk).get_asym1(alpha->index, beta->index) +
                                                          Phase.FieldsDot(i+di, j+dj, k+dk).get_asym2(alpha->index, beta->index));

                                    locDelta += Grid.dx * ds->weight * locPsi * Projection * coef * scale *
                                                (alpha->value + beta->value) * 0.5 *
                                                (locMoleFractionB - locMoleFractionA) *

                                                (DC[pIndexA] - DC[pIndexB])/(DC[pIndexA] + DC[pIndexB]);
                                }
                            }
                            Cx.MoleFractionsTotalDot(i,j,k,{Comp}) += locDelta;
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::CalculateDiffusionIncrements(
                                                        PhaseField& Phase,
                                                        Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        for(size_t n = 0; n < Nphases; n++)
        if(Phase.Fractions(i,j,k,{n}) != 0.0)
        {
            double locDelta = 0.0;
            for(auto ds = DStencil.cbegin(); ds != DStencil.cend(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                if(Phase.Fractions(i+di, j+dj, k+dk, {n}) != 0.0)
                {
                    locDelta += ds->weight*DC[n]*IDC[n]*
                              0.5*(Phase.Fractions(i   , j   , k   , {n}) +
                                   Phase.Fractions(i+di, j+dj, k+dk, {n}))*
                                  (Cx.MoleFractions(i+di, j+dj, k+dk, {n, Comp}) -
                                   Cx.MoleFractions(i   , j   , k   , {n, Comp}));
                }
            }
            Cx.MoleFractionsTotalDot(i,j,k,{Comp}) += locDelta;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::CalculateChemicalPotentialContribution(
                                                            PhaseField& Phase,
                                                            Composition& Cx,
                                                            Temperature& Tx)
{
    vector<double> Mob(Nphases,0.0);
    for(size_t n = 0; n < Nphases; n++)
    {
        Mob[n] = DC[n]*Cx.TotalMolarVolume/(R*Tx.Tavg);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        for(size_t n = 0; n < Nphases; n++)
        if(Phase.Fractions(i, j, k, {n}) != 0.0 and
           !Stoichiometric[n])
        {
            double locDelta = 0.0;
            for(auto ds = DStencil.cbegin(); ds != DStencil.cend(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                if(Phase.Fractions(i+di, j+dj, k+dk, {n}) != 0.0)
                {
                    locDelta += ds->weight*Mob[n]*IDC[n]*
                                0.5*(Phase.Fractions(i+di, j+dj, k+dk, {n}) +
                                     Phase.Fractions(i   , j   , k   , {n}))*
                                    (dMu(i+di, j+dj, k+dk, {n, Comp}) -
                                     dMu(i   , j   , k   , {n, Comp}));
                }
            }
            Cx.MoleFractionsTotalDot(i, j, k, {Comp}) += locDelta;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::ClearChemicalPotential(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,dMu,dMu.Bcells(),)
    {
        for(size_t alpha = 0; alpha != Nphases; ++alpha)
        {
            dMu(i,j,k,{alpha, Comp}) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

bool EquilibriumPartitionDiffusionBinary::CalculateLimits(PhaseField& Phase,
                                                          Composition& Cx,
                                                          double dt)
{
    bool LimitingNeeded = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.NormTotal,Cx.NormTotal.Bcells(),)
    {
        Cx.NormTotal(i,j,k,{Comp}) = 1.0;

        Cx.Limiting(i,j,k) = false;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotalDot,0,reduction(||:LimitingNeeded))
    {
        double dTotal = Cx.MoleFractionsTotalDot(i,j,k,{Comp})*dt;

        if(fabs(dTotal) > DBL_EPSILON)
        {
            double locMIN = 0.0;
            double locMAX = 0.0;

            if(Phase.Fields(i,j,k).wide_interface())
            {
                Tensor<double,1> locPhaseFractions = Phase.Fractions(i,j,k);

                for(size_t alpha = 0; alpha < Nphases; alpha++)
                {
                    locMIN += Cx.MinMoleFraction(alpha, Comp)*locPhaseFractions({alpha});
                    locMAX += Cx.MaxMoleFraction(alpha, Comp)*locPhaseFractions({alpha});
                }
            }
            else
            {
                size_t alpha = Phase.FieldsProperties[Phase.Fields(i,j,k).front().index].Phase;

                locMIN = Cx.MinMoleFraction(alpha, Comp);
                locMAX = Cx.MaxMoleFraction(alpha, Comp);
            }

            double TotalOld = Cx.MoleFractionsTotal(i,j,k,{Comp});
            double TotalNew = TotalOld + dTotal;
            if(TotalNew > locMAX)
            {
                Cx.NormTotal(i,j,k,{Comp}) = (locMAX - TotalOld)/dTotal;
                LimitingNeeded = true;

                for (int x = -Grid.dNx; x <= Grid.dNx; ++x)
                for (int y = -Grid.dNy; y <= Grid.dNy; ++y)
                for (int z = -Grid.dNz; z <= Grid.dNz; ++z)
                {
                    Cx.Limiting(i+x, j+y, k+z) = true;
                }
            }
            if(TotalNew < locMIN)
            {
                Cx.NormTotal(i,j,k,{Comp}) = (locMIN - TotalOld)/dTotal;
                LimitingNeeded = true;

                for (int x = -Grid.dNx; x <= Grid.dNx; ++x)
                for (int y = -Grid.dNy; y <= Grid.dNy; ++y)
                for (int z = -Grid.dNz; z <= Grid.dNz; ++z)
                {
                    Cx.Limiting(i+x, j+y, k+z) = true;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &LimitingNeeded, 1, OP_MPI_CXX_BOOL, OP_MPI_LOR, OP_MPI_COMM_WORLD);
#endif
    return LimitingNeeded;
}

void EquilibriumPartitionDiffusionBinary::LimitDiffusionIncrements(PhaseField& Phase,
                                                                   Composition& Cx,
                                                                   Temperature& Tx)
{
    vector<double> Mob(Nphases,0.0);
    for(size_t n = 0; n < Nphases; n++)
    {
        Mob[n] = DC[n]*Cx.TotalMolarVolume/(R*Tx.Tavg);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotalDot,0,)
    {
        if(Cx.Limiting(i, j, k))
        {
            for(size_t alpha = 0; alpha < Nphases; ++alpha)
            if(Phase.Fractions(i, j, k, {alpha}) != 0.0)
            for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
            {
                int di = ds->di;
                int dj = ds->dj;
                int dk = ds->dk;

                if (Phase.Fractions(i+di, j+dj, k+dk, {alpha}) != 0.0)
                {
                    double localDelta = ds->weight*Mob[alpha]*IDC[alpha]*
                            0.5*(Phase.Fractions(i, j, k, {alpha}) +
                                 Phase.Fractions(i+di, j+dj, k+dk,{alpha}))*

                                (dMu(i+di, j+dj, k+dk)({alpha}) -
                                 dMu(i   , j   , k   )({alpha}));

                    if (Cx.NormTotal(i   , j   , k   )({Comp}) <=
                        Cx.NormTotal(i+di, j+dj, k+dk)({Comp}))
                    {
                        double localNorm = Cx.NormTotal(i, j, k)({Comp});
                        double localIncrement = localDelta*(1.0 - localNorm);
                        if(localIncrement > DBL_EPSILON)
                        {
                            Cx.MoleFractionsTotalDot(i, j, k)({Comp}) -=
                                                              localIncrement;
                        }
                    }
                    if (Cx.NormTotal(i   , j   , k   )({Comp}) >
                        Cx.NormTotal(i+di, j+dj, k+dk)({Comp}))
                    {
                        double localNorm = Cx.NormTotal(i+di, j+dj, k+dk)({Comp});
                        double localIncrement = localDelta*(1.0 - localNorm);
                        if(localIncrement > DBL_EPSILON)
                        {
                            Cx.MoleFractionsTotalDot(i, j, k)({Comp}) -=
                                                              localIncrement;
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::ApplyIncrements(PhaseField& Phase,
                                                          Composition& Cx,
                                                          double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        Cx.MoleFractionsTotal(i,j,k,{Comp})   += Cx.MoleFractionsTotalDot(i,j,k,{Comp})*dt;
        Cx.MoleFractionsTotal(i,j,k,{RefComp}) = 1.0 - Cx.MoleFractionsTotal(i,j,k,{Comp});
        Cx.MoleFractionsTotalDot(i,j,k,{Comp}) = 0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitionDiffusionBinary::SolveDiffusion(PhaseField& Phase,
                                                         Composition& Cx,
                                                         Temperature& Tx,
                                                         BoundaryConditions& BC,
                                                         double dt)
{
    CalculateDiffusionCoefficients(Tx);

    CalculatePhaseConcentrations(Phase, Cx, Tx);

    CalculateDiffusionIncrements(Phase, Cx);

    if(EnableAntiTrapping)
    {
        CalculateAntitrappingIncrements(Phase, Cx);
    }

    ApplyIncrements(Phase, Cx, dt);
    Cx.SetBoundaryConditions(BC);

    if(ConsiderChemicalPotential)
    {
        CalculateChemicalPotentialContribution(Phase, Cx, Tx);
        bool LimitingNeeded = CalculateLimits(Phase, Cx, dt);

        if(LimitingNeeded)
        {
            Cx.SetTotalLimitsBoundaryConditions(BC);
            LimitDiffusionIncrements(Phase, Cx, Tx);
        }
        ApplyIncrements(Phase, Cx, dt);
        ClearChemicalPotential();
        Cx.SetBoundaryConditions(BC);
    }

    // Solve diffusion in 1D extensions
    if(Cx.ExtensionsActive)
    {
#ifdef MPI_PARALLEL
        if(MPI_RANK == 0)
#endif
        if(Cx.ExtensionX0.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionX0, BC.BC0X, dt);}
#ifdef MPI_PARALLEL
        if(MPI_RANK == MPI_SIZE - 1)
#endif
        if(Cx.ExtensionXN.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionXN, BC.BCNX, dt);}

        if(Cx.ExtensionY0.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionY0, BC.BC0Y, dt);}
        if(Cx.ExtensionYN.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionYN, BC.BCNY, dt);}
        if(Cx.ExtensionZ0.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionZ0, BC.BC0Z, dt);}
        if(Cx.ExtensionZN.isActive()) {Solve1Dextension(Cx, Tx, Cx.ExtensionZN, BC.BCNZ, dt);}

        Cx.SetBoundaryConditions(BC);
    }
}

void EquilibriumPartitionDiffusionBinary::Solve1Dextension(Composition& Cx,
                                              Temperature& Tx,
                                              Composition1Dextension& CxExt,
                                              BoundaryConditionTypes& extBC,
                                              double dt)
{
    std::pair<int, int> Xbounds = three_state_bounds_selector(CxExt.Direction[0], 0, Cx.Grid.Nx);
    std::pair<int, int> Ybounds = three_state_bounds_selector(CxExt.Direction[1], 0, Cx.Grid.Ny);
    std::pair<int, int> Zbounds = three_state_bounds_selector(CxExt.Direction[2], 0, Cx.Grid.Nz);

    double boundaryDC = 0.0;
    Tensor<double,1> boundaryValue({Cx.Ncomp});
    double area = 0.0;

    // calculate boundary values of composition and diffusion coefficients
    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    {
        for (size_t comp = 0; comp < Cx.Ncomp; comp++)
        {
            boundaryValue({comp}) += Cx.MoleFractionsTotal(i,j,k,{comp});
        }
        boundaryDC += DC[CxExt.PhaseIndex];
        area++;
    }

#ifdef MPI_PARALLEL
    if(CxExt.Direction[0] == 0)
    {
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, boundaryValue.data(), Cx.Ncomp, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &boundaryDC, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &area, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
#endif

    // set boundary value of composition and diffusion coefficients in the extension
    CxExt.Data(-1) = boundaryValue/area;
    double DC      = boundaryDC/area;
    double dx_2    = 1.0/pow(Cx.Grid.dx,2);

    // calculate diffusion increments
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        CxExt.DataDot(x,{Comp}) += dx_2*DC*(CxExt.Data(x-1,{Comp}) - 2.0*CxExt.Data(x,{Comp}) + CxExt.Data(x+1,{Comp}));
    }

    // apply diffusion increments
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        CxExt.Data(x,{Comp}) += CxExt.DataDot(x,{Comp})*dt;
        CxExt.DataDot(x).set_to_zero();
    }

    // calculate mole fraction of reference element
    for(long int x = 0; x < (long int)CxExt.size(); x++)
    {
        CxExt.Data(x,{RefComp}) = 1.0 - CxExt.Data(x,{Comp});
    }

    // set extension's boundary conditions
    CxExt.setBC(extBC);

    // set boundary values in the simulation domain
    for (int i = Xbounds.first; i < Xbounds.second; ++i)
    for (int j = Ybounds.first; j < Ybounds.second; ++j)
    for (int k = Zbounds.first; k < Zbounds.second; ++k)
    for (size_t comp = 0; comp < Cx.Ncomp; comp++)
    {
        Cx.MoleFractionsTotal(i+CxExt.Direction[0],j+CxExt.Direction[1],k+CxExt.Direction[2],{comp}) = CxExt.Data(0)({comp});
    }
}

}// namespace openphase
