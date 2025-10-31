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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Johannes Goerler
 *
 */

#include "ElasticProperties.h"
#include "Containers/dMatrix3x3.h"
#include "Settings.h"
#include "SymmetryVariants.h"
#include "PhaseField.h"
#include "Composition.h"
#include "DrivingForce.h"
#include "BoundaryConditions.h"
#include "Temperature.h"
#include "Velocities.h"
#include "VTK.h"
#include "InterfaceProperties.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "Tools.h"
#include "AdvectionHR.h"

namespace openphase
{
using namespace std;

ElasticProperties::ElasticProperties(Settings& locSettings,
                                     const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void ElasticProperties::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "ElasticProperties";
    thisobjectname = thisclassname + ObjectNameSuffix;

    LargeDeformations = false;
    AnyPlasticity = true;
    KeepAspectRatio = false;
    PreventShear = false;

    Grid = locSettings.Grid;

    ElasticityModel = ElasticityModels::Khachaturyan;
    StrainModel  = StrainModels::Small;

    Nphases = locSettings.Nphases;
    Ncomp = locSettings.Ncomp;
    ElementNames = locSettings.ElementNames;

    AppliedStrainMask.set_to_zero();
    AppliedStrainRateMask.set_to_zero();
    AppliedStressMask.set_to_zero();

    AppliedStrain.set_to_zero();
    AppliedStress.set_to_zero();
    AppliedStrainRate.set_to_zero();

    AverageStrain.set_to_zero();
    RemeshedStrain.set_to_zero();
    StrainToRemesh.set_to_zero();
    //AverageDeformationGradient.set_to_unity();

    size_t Bcells = Grid.Bcells;

    Stresses.Allocate(Grid, Bcells);

    DeformationJumps.Allocate(Grid, Bcells);
    NonLocalDeformationJumpTerm.Allocate(Grid, Bcells);

    DeformationGradientsTotal.Allocate(Grid, Bcells);
    DeformationGradientsEigen.Allocate(Grid, Bcells);
    EffectiveElasticConstants.Allocate(Grid, Bcells);
    Displacements.Allocate(Grid, Bcells);
    ElasticConstants.Allocate(Grid, Bcells);
    TransformationStretches.Allocate(Grid, Bcells);

    // TODO: Plasticity related storages should be allocated only if plasticity is active.
    DeformationGradientsPlastic.Allocate(Grid, Bcells);
    VelocityGradientsPlastic.Allocate(Grid, Bcells);

    PhaseElasticConstants.Allocate(Nphases);
    PhaseTransformationStretches.Allocate(Nphases);

    PhaseAlpha.Allocate(Nphases);
    PhaseGamma.Allocate(Nphases);
    Tref.Allocate(Nphases);
    PoissonRatio.Allocate(Nphases);

    GrainElasticConstants.Allocate(Nphases);
    GrainTransformationStretches.Allocate(Nphases);

    GrainAlpha.Allocate(Nphases);
    GrainGamma.Allocate(Nphases);

    if(Ncomp > 0)
    {
        PhaseKappa.Allocate({Nphases, Ncomp});
        PhaseLambda.Allocate({Nphases, Ncomp});
        Cref.Allocate({Nphases, Ncomp});

        GrainKappa.Allocate({Nphases, Ncomp});
        GrainLambda.Allocate({Nphases, Ncomp});
    }

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseAlpha[alpha].set_to_zero();
        PhaseGamma[alpha].set_to_zero();

        Tref[alpha] = 0.0;

        PhaseTransformationStretches[alpha].set_to_unity();
        PhaseElasticConstants[alpha].set_to_zero();

        for(size_t comp = 0; comp != Ncomp; comp++)
        {
            Cref({alpha, comp}) = 0.0;
            PhaseLambda({alpha, comp}).set_to_zero();
            PhaseKappa({alpha, comp}).set_to_zero();
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Stresses,Bcells,)
    {
        DeformationGradientsTotal(i,j,k).set_to_unity();                        //NOTE: zero strain
        DeformationGradientsEigen(i,j,k).set_to_unity();                        //NOTE: zero strain
        Stresses(i,j,k).set_to_zero();
        // TODO: Plasticity related storages should be allocated only if plasticity is active.
        DeformationGradientsPlastic(i,j,k).set_to_unity();                      //NOTE: zero strain
        VelocityGradientsPlastic(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Variants.Initialize(locSettings);

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    ConsoleOutput::Write(thisclassname, "Initialized");
}

void ElasticProperties::InitializeLD(Settings& locSettings)
{
    if (!Stresses.IsAllocated())
    {
        ConsoleOutput::WriteWarning("ElasticProperties should be initialized before LargeDeformations", thisclassname, "InitializedLD");
        OP_Exit(EXIT_FAILURE);
    }

    LargeDeformations = true;

    size_t Bcells = Grid.Bcells;

    StressIncrements.Allocate(Grid, Bcells);
    VelocityGradientsTotal.Allocate(Grid, Bcells);
    LocalRotations.Allocate(Grid, Bcells);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Stresses, Bcells, )
    {
        StressIncrements(i,j,k).set_to_zero();
        VelocityGradientsTotal(i,j,k).set_to_zero();
        LocalRotations(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    ConsoleOutput::Write("ElasticProperties", "Initialized for Large Deformations");
}

void ElasticProperties::ReadInput(const string InputFileName)
{
    ConsoleOutput::Write("Source", InputFileName);
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File " + InputFileName + " could not be opened",thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ReadInput(data);
}

void ElasticProperties::ReadInput(stringstream& inp)
{
    ConsoleOutput::WriteLineInsert("ElasticProperties input");
    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseElasticConstants[alpha].set_to_zero();
    }

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    string tmp1 = FileInterface::ReadParameterK(inp, moduleLocation, "EModel", false, "KHACHATURYAN");

    if (tmp1 == "KHACHATURYAN")
    {
        ElasticityModel = ElasticityModels::Khachaturyan;
    }
    else if (tmp1 == "STEINBACH")
    {
        ElasticityModel = ElasticityModels::Steinbach;
    }
    else if (tmp1 == "VOIGT")
    {
        ElasticityModel = ElasticityModels::Voigt;
    }
    else if (tmp1 == "REUSS")
    {
        ElasticityModel = ElasticityModels::Reuss;
    }
    else if (tmp1 == "RANK1")
    {
        ElasticityModel = ElasticityModels::Rank1;
    }
    else if (tmp1 == "RANK1NL")
    {
        ElasticityModel = ElasticityModels::Rank1NL;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong elasticity model specified!\nThe default \"Khachaturyan\" model is used!", thisclassname, "ReadInput()");
    }

    /// Mechanical boundary conditions

    std::vector<std::string> BCnames {"BCX", "BCY","BCZ","BCYZ","BCXZ","BCXY"};
    std::vector<std::string> BCvalueNames {"BCValueX", "BCValueY","BCValueZ","BCValueYZ","BCValueXZ","BCValueXY"};

    for(int n = 0; n < 6; n++)
    {
        string tmp2 = FileInterface::ReadParameterK(inp, moduleLocation, BCnames[n], 0);

        if (tmp2 == "FREEBOUNDARIES")
        {
            AppliedStressMask[n] = 1.0;
            AppliedStress[n] = 0.0;
        }
        if (tmp2 == "APPLIEDSTRAIN")
        {
            AppliedStrain[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
            AppliedStrainMask[n] = 1.0;
        }
        if (tmp2 == "APPLIEDSTRAINRATE")
        {
            AppliedStrainRate[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
            AppliedStrainRateMask[n] = 1.0;
        }
        if (tmp2 == "APPLIEDSTRESS")
        {
            AppliedStressMask[n] = 1.0;
            AppliedStress[n] = FileInterface::ReadParameterD(inp, moduleLocation, BCvalueNames[n]);
        }
    }

    string tmp3 = FileInterface::ReadParameterK(inp, moduleLocation, "Restrict", 0);

    if (tmp3 == "ASPECTRATIO")
    {
        KeepAspectRatio = true;
    }

    if (tmp3 == "SHEAR")
    {
        PreventShear = true;
    }

    for(size_t alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseElasticConstants[alpha].set_to_zero();
    }
    
    string tmp4 = FileInterface::ReadParameterK(inp, moduleLocation, "ElasticConstantMode", false, "ELASTICMODULI");
    
    if (tmp4 != "ELASTICMODULI" and tmp4 != "TENSOR")
    {
        std::cerr << "WARNING: ElasticConstantMode should be ELASTICMODULI or TENSOR, trying ELASTICMODULI." << std::endl;
        tmp4 = "ELASTICMODULI";        
    }

    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        // Read elastic moduli if provided
        std::stringstream KConv;
        std::stringstream EConv;
        std::stringstream LConv;
        std::stringstream GConv;
        std::stringstream NConv;
        std::stringstream MConv;

        KConv << "K_" << pIndex;
        EConv << "E_" << pIndex;
        LConv << "L_" << pIndex;
        GConv << "G_" << pIndex;
        NConv << "Nu_"<< pIndex;
        MConv << "M_" << pIndex;

        double K      = FileInterface::ReadParameterD(inp, moduleLocation, KConv.str(), false, 0); ///< Bulk modulus
        double E      = FileInterface::ReadParameterD(inp, moduleLocation, EConv.str(), false, 0); ///< Young's modulus
        double lambda = FileInterface::ReadParameterD(inp, moduleLocation, LConv.str(), false, 0); ///< First Lame
        double G      = FileInterface::ReadParameterD(inp, moduleLocation, GConv.str(), false, 0); ///< Shear modulus
        double nu     = FileInterface::ReadParameterD(inp, moduleLocation, NConv.str(), false, 0); ///< Poisson's ratio
        double M      = FileInterface::ReadParameterD(inp, moduleLocation, MConv.str(), false, 0); ///< P-wave modulus

        if      (K      != 0.0 and E      != 0.0) {lambda = 3.0*K*(3.0*K-E)/(9.0*K-E);    G = 3.0*K*E/(9.0*K-E);}
        else if (K      != 0.0 and lambda != 0.0) {                                       G = 3*(K-lambda)/2.0;}
        else if (K      != 0.0 and G      != 0.0) {lambda = K-2.0*G/3.0;}
        else if (K      != 0.0 and nu     != 0.0) {lambda = 3.0*K*nu/(1+nu);              G = 3*K*(1-2.0*nu)/(2.0*(1.0+nu));}
        else if (K      != 0.0 and M      != 0.0) {lambda = (3.0*K-M)/2;                  G = 3.0*(M-K)/4.0;}
        else if (E      != 0.0 and lambda != 0.0) {                                       G = (E-3.0*lambda+std::sqrt(E*E+9.0*lambda*lambda+2.0*E*lambda))/4.0;}
        else if (E      != 0.0 and G      != 0.0) {lambda = G*(E-2.0*G)/(3.0*G-E);}
        else if (E      != 0.0 and nu     != 0.0) {lambda = E*nu/((1.0+nu)*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu));}
        else if (E      != 0.0 and M      != 0.0) {ConsoleOutput::WriteExit("Choice of elastic moduli is not unique", thisclassname, "ReadInput"); std::exit(EXIT_FAILURE);}
        else if (lambda != 0.0 and nu     != 0.0) {                                       G = lambda*(1.0-2.0*nu)/(2.0*nu);}
        else if (lambda != 0.0 and M      != 0.0) {                                       G = (M-lambda)/2.0;}
        else if (G      != 0.0 and nu     != 0.0) {lambda = 2.0*G*nu/(1.0-2.0*nu);}
        else if (G      != 0.0 and M      != 0.0) {lambda = M-2.0*G;}
        else if (nu     != 0.0 and M      != 0.0) {lambda = M*nu/(1.0-nu);                G = M*(1.0-2.0*nu)/(2.0*(1.0-nu));}

        if (std::fabs(lambda) > DBL_EPSILON and std::fabs(G) > DBL_EPSILON and tmp4 == "ELASTICMODULI")
        {
            PhaseElasticConstants[pIndex](0,0) = 2.0*G + lambda;
            PhaseElasticConstants[pIndex](1,1) = 2.0*G + lambda;
            PhaseElasticConstants[pIndex](2,2) = 2.0*G + lambda;

            PhaseElasticConstants[pIndex](0,1) = lambda;
            PhaseElasticConstants[pIndex](0,2) = lambda;
            PhaseElasticConstants[pIndex](1,2) = lambda;

            PhaseElasticConstants[pIndex](1,0) = lambda;
            PhaseElasticConstants[pIndex](2,0) = lambda;
            PhaseElasticConstants[pIndex](2,1) = lambda;

            PhaseElasticConstants[pIndex](3,3) = G;
            PhaseElasticConstants[pIndex](4,4) = G;
            PhaseElasticConstants[pIndex](5,5) = G;

        }
        else
        {
            stringstream converter;
            converter << "C" << "_" << pIndex;
            PhaseElasticConstants[pIndex] = FileInterface::ReadParameterM6x6(inp, moduleLocation, converter.str(),false,dMatrix6x6::ZeroTensor());
            if(PhaseElasticConstants[pIndex] == dMatrix6x6::ZeroTensor())
            {
                for(int ii =  1; ii <= 6; ii++)
                for(int jj = ii; jj <= 6; jj++)
                {
                    stringstream Cij;
                    Cij << "C" << ii << jj << "_" << pIndex;

                    PhaseElasticConstants[pIndex](ii-1,jj-1) = FileInterface::ReadParameterD(inp, moduleLocation, Cij.str(),false,0);
                    if(ii != jj)
                    {
                        PhaseElasticConstants[pIndex](jj-1,ii-1) = PhaseElasticConstants[pIndex](ii-1,jj-1);
                    }
                }
            }
        }
        PoissonRatio[pIndex] = PhaseElasticConstants[pIndex](0,1) / (PhaseElasticConstants[pIndex](0,1) + PhaseElasticConstants[pIndex](0,0));

        // Check if stiffness constants are read.
        if (PhaseElasticConstants[pIndex].norm() <= DBL_EPSILON)
        {
            std::string message = "Elastic constants for phase " + to_string(pIndex) + " not set.";
            ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
            exit(3);
        }
    }

    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "U_" << pIndex;
        PhaseTransformationStretches[pIndex] = FileInterface::ReadParameterM3x3(inp, moduleLocation, converter.str(),true,dMatrix3x3::UnitTensor());
    }
    // Considering external forces
    ConsiderExternalForces = FileInterface::ReadParameterB(inp, moduleLocation, "ConsiderExternalForces", false, false);
    if(ConsiderExternalForces)
    {
        ForceDensity.Allocate(Grid, 1);
    }

    // Reading chemo-mechanical coupling parameters if considered
    ChemoMechanicalCoupling = FileInterface::ReadParameterB(inp, moduleLocation, "ChemoMechanicalCoupling",false,false);
    if(ChemoMechanicalCoupling)
    {
        int counter = 0;
        for(size_t comp = 0; comp < Ncomp; comp++)
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            stringstream converter;
            converter << ElementNames[comp] << "_" << alpha;
            string counter = converter.str();
            Cref({alpha, comp}) = FileInterface::ReadParameterD(inp, moduleLocation, string("Cref_") + counter, false, 0.0);
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            for (size_t comp = 0; comp < Ncomp; comp++)
            {
                stringstream converter;
                converter << "Kappa_" << pIndex << "_" << ElementNames[comp];
                PhaseKappa({ pIndex, comp }) = FileInterface::ReadParameterM6x6(inp, moduleLocation, converter.str(),false,dMatrix6x6::ZeroTensor());

                if (PhaseKappa({pIndex, comp}).norm() != 0.0)
                {
                    counter++;
                }
            }
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        for (size_t comp = 0; comp < Ncomp; comp++)
        {
            stringstream converter;
            converter << "Lambda_" << pIndex << "_" << ElementNames[comp];
            PhaseLambda({ pIndex, comp }) = FileInterface::ReadParameterM3x3(inp, moduleLocation, converter.str(),false,dMatrix3x3::ZeroTensor());

            if (PhaseLambda({pIndex, comp}).norm() != 0.0)
            {
                counter++;
            }
        }

        if(counter == 0)
        {
            ConsoleOutput::WriteWarning("ChemoMechanicalCoupling is ON but no coupling parameters specified", thisclassname, "ReadInput()");
        }
    }
    // Reading thermo-mechanical coupling parameters if considered
    ThermoMechanicalCoupling = FileInterface::ReadParameterB(inp, moduleLocation, "ThermoMechanicalCoupling",false,false);
    if(ThermoMechanicalCoupling)
    {
        int counter = 0;
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            stringstream converter;
            converter << "Tref_" << alpha;
            Tref[alpha] = FileInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            stringstream converter;
            converter << "Gamma_" << pIndex;

            double ScalarGamma = FileInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);

            if (ScalarGamma != 0.0)
            {
                for (int ii = 0; ii < 6; ii++)
                for (int jj = 0; jj < 6; jj++)
                {
                    PhaseGamma[pIndex](ii,jj) = ScalarGamma;
                }
            }
            else
            {
                for (int ii = 1; ii <= 6; ii++)
                for (int jj = ii; jj <= 6; jj++)
                {
                    stringstream GammaString;
                    GammaString << "Gamma" << ii << jj << "_" << pIndex;
                    PhaseGamma[pIndex](ii - 1, jj - 1) =
                        FileInterface::ReadParameterD(inp, moduleLocation, GammaString.str(), false, 0.0);
                    if (ii != jj)
                    {
                        PhaseGamma[pIndex](jj - 1, ii - 1) = PhaseGamma[pIndex](ii - 1, jj - 1);
                    }
                }
            }

            if (PhaseGamma[pIndex].norm() != 0.0)
            {
                counter++;

                stringstream GammaNX;
                GammaNX << "Gamma_" << pIndex;
                ConsoleOutput::Write(GammaNX.str(),PhaseGamma[pIndex],6);
            }
        }

        for (size_t pIndex = 0; pIndex < Nphases; pIndex++)
        {
            for (int ii = 1; ii <= 3; ii++)
            for (int jj = ii; jj <= 3; jj++)
            {
                stringstream converter;
                converter << "Alpha" << ii << jj << "_" << pIndex;

                PhaseAlpha[pIndex](ii - 1, jj - 1) =
                    FileInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
                if (ii != jj)
                {
                    PhaseAlpha[pIndex](jj - 1, ii - 1) = PhaseAlpha[pIndex](ii - 1, jj - 1);
                }
            }
            if(PhaseAlpha[pIndex].norm() != 0.0)
            {
                counter++;
                stringstream AlphaNX;
                AlphaNX << "Alpha_" << pIndex;
                ConsoleOutput::Write(AlphaNX.str(),PhaseAlpha[pIndex],6);
            }
        }

        if(counter == 0)
        {
            ConsoleOutput::WriteWarning("ThermoMechanicalCoupling is ON but no coupling parameters specified", thisclassname, "ReadInput()");
        }
    }
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();

    Variants.ReadInput(inp);
}

void ElasticProperties::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);
    ///Changing box size
    dMatrix3x3 stretch;
    stretch.set_to_zero();
    stretch(0,0) = double(newNx)/double(Grid.Nx);
    stretch(1,1) = double(newNy)/double(Grid.Ny);
    stretch(2,2) = double(newNz)/double(Grid.Nz);

    RemeshedStrain[0] += stretch(0,0) - 1.0;
    RemeshedStrain[1] += stretch(1,1) - 1.0;
    RemeshedStrain[2] += stretch(2,2) - 1.0;

    StrainToRemesh[0] -= stretch(0,0) - 1.0;
    StrainToRemesh[1] -= stretch(1,1) - 1.0;
    StrainToRemesh[2] -= stretch(2,2) - 1.0;

    Grid.SetDimensions(newNx, newNy, newNz);

    Stresses.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    StressIncrements.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    DeformationGradientsTotal.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    DeformationGradientsEigen.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    DeformationGradientsPlastic.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    VelocityGradientsPlastic.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    Displacements.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    EffectiveElasticConstants.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    ElasticConstants.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    TransformationStretches.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    if(VelocityGradientsTotal.IsAllocated()) VelocityGradientsTotal.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    if(LocalRotations.IsAllocated()) LocalRotations.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    DeformationJumps.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
    NonLocalDeformationJumpTerm.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);

    SetBoundaryConditions(BC);

    //WriteRemeshingData(tStep, sim_time);
    ConsoleOutput::Write(thisclassname, "Remeshed");
}

void ElasticProperties::SetGrainsProperties(const PhaseField& Phase)
{
    size_t size = Phase.FieldsProperties.size();
    if(GrainTransformationStretches.size() != size)
    {
        GrainTransformationStretches.Reallocate(size);
        GrainElasticConstants.Reallocate(size);
        if(Ncomp)
        {
            GrainLambda.Reallocate({size, Ncomp});
            GrainKappa.Reallocate({size, Ncomp});
        }
        GrainAlpha.Reallocate(size);
        GrainGamma.Reallocate(size);
    }
    for(size_t alpha = 0; alpha != size; alpha++)
    if(Phase.FieldsProperties[alpha].Exist)
    {
        size_t pIndex = Phase.FieldsProperties[alpha].Phase;
        size_t vIndex = Phase.FieldsProperties[alpha].Variant;

        GrainTransformationStretches[alpha] = PhaseTransformationStretches[pIndex];
        GrainElasticConstants[alpha] = PhaseElasticConstants[pIndex];
        GrainAlpha[alpha] = PhaseAlpha[pIndex];
        GrainGamma[alpha] = PhaseGamma[pIndex];

        if(Variants.set)
        {
            GrainTransformationStretches[alpha].rotate(Variants(pIndex, vIndex));
            GrainElasticConstants[alpha].rotate(Variants(pIndex, vIndex));
            GrainAlpha[alpha].rotate(Variants(pIndex, vIndex));
            GrainGamma[alpha].rotate(Variants(pIndex, vIndex));
        }

        GrainTransformationStretches[alpha].rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);
        GrainElasticConstants[alpha].rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);
        GrainAlpha[alpha].rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);
        GrainGamma[alpha].rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);

        for(size_t comp = 0; comp != Ncomp; comp++)
        {
            GrainLambda({alpha, comp}) = PhaseLambda({pIndex, comp});
            GrainKappa({alpha, comp})  = PhaseKappa({pIndex, comp});

            if(Variants.set)
            {
                GrainLambda({alpha, comp}).rotate(Variants(pIndex, vIndex));
                GrainKappa({alpha, comp}).rotate(Variants(pIndex, vIndex));
            }

            GrainLambda({alpha, comp}).rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);
            GrainKappa({alpha, comp}).rotate(Phase.FieldsProperties[alpha].Orientation.RotationMatrix);
        }
    }
}

void ElasticProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    // Only periodic BC are correct. For non periodic boundary conditions gradients should be treated differently.

    BC.SetX(DeformationGradientsTotal);
    BC.SetY(DeformationGradientsTotal);
    BC.SetZ(DeformationGradientsTotal);

    BC.SetX(DeformationGradientsEigen);
    BC.SetY(DeformationGradientsEigen);
    BC.SetZ(DeformationGradientsEigen);

    BC.SetX(Displacements);
    BC.SetY(Displacements);
    BC.SetZ(Displacements);

    if(AnyPlasticity)
    {
        SetBoundaryConditionsPlastic(BC);
    }

    BC.SetX(Stresses);
    BC.SetY(Stresses);
    BC.SetZ(Stresses);

    BC.SetX(EffectiveElasticConstants);
    BC.SetY(EffectiveElasticConstants);
    BC.SetZ(EffectiveElasticConstants);

    if(LargeDeformations)
    {
        BC.SetX(StressIncrements);
        BC.SetY(StressIncrements);
        BC.SetZ(StressIncrements);

        BC.SetX(VelocityGradientsTotal);
        BC.SetY(VelocityGradientsTotal);
        BC.SetZ(VelocityGradientsTotal);
    }
}

void ElasticProperties::SetBoundaryConditionsPlastic(const BoundaryConditions& BC)
{
    BC.SetX(DeformationGradientsPlastic);
    BC.SetY(DeformationGradientsPlastic);
    BC.SetZ(DeformationGradientsPlastic);

    BC.SetX(VelocityGradientsPlastic);
    BC.SetY(VelocityGradientsPlastic);
    BC.SetZ(VelocityGradientsPlastic);
}

ElasticProperties& ElasticProperties::operator= (const ElasticProperties& rhs)
{
    // protect against self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "ElasticProperties")
    {
        thisclassname = rhs.thisclassname;

        LargeDeformations = rhs.LargeDeformations;
        AnyPlasticity = rhs.AnyPlasticity;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;
        Ncomp   = rhs.Ncomp;
        ElementNames = rhs.ElementNames;
        ElasticityModel = rhs.ElasticityModel;
        StrainModel = rhs.StrainModel;

        if (DeformationGradientsTotal.IsNotAllocated())
        {
            DeformationGradientsTotal.Allocate(Grid, rhs.DeformationGradientsTotal.Bcells());
            DeformationGradientsEigen.Allocate(Grid, rhs.DeformationGradientsEigen.Bcells());
            Stresses.Allocate(Grid, rhs.Stresses.Bcells());
            EffectiveElasticConstants.Allocate(Grid, rhs.EffectiveElasticConstants.Bcells());
            Displacements.Allocate(Grid, 0);
            ElasticConstants.Allocate(Grid, rhs.ElasticConstants.Bcells());
            TransformationStretches.Allocate(Grid, rhs.TransformationStretches.Bcells());

            if (AnyPlasticity)
            {
                VelocityGradientsPlastic.Allocate(Grid, rhs.VelocityGradientsTotal.Bcells());
                DeformationGradientsPlastic.Allocate(Grid, rhs.DeformationGradientsEigen.Bcells());
            }

            if (LargeDeformations)
            {
                VelocityGradientsTotal.Allocate(Grid, rhs.VelocityGradientsTotal.Bcells());
                StressIncrements.Allocate(Grid, rhs.StressIncrements.Bcells());
            }
        }
        else if (not DeformationGradientsTotal.IsSize(rhs.Grid.Nx, rhs.Grid.Ny, rhs.Grid.Nz))
        {
            DeformationGradientsTotal.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            DeformationGradientsEigen.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Stresses.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            EffectiveElasticConstants.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Displacements.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            ElasticConstants.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            TransformationStretches.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

            if (LargeDeformations)
            {
                StressIncrements.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
                VelocityGradientsTotal.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            }
        }

        if (PhaseElasticConstants.IsNotAllocated())
        {
            PhaseElasticConstants.Allocate(Nphases);
            PhaseTransformationStretches.Allocate(Nphases);

            PhaseKappa.Allocate({Nphases, Ncomp});
            PhaseLambda.Allocate({Nphases, Ncomp});
            PhaseGamma.Allocate(Nphases);
            PhaseAlpha.Allocate(Nphases);
            Tref.Allocate(Nphases);
            Cref.Allocate({Nphases, Ncomp});
        }
        else if (PhaseElasticConstants.size() != rhs.Nphases)
        {
            PhaseElasticConstants.Reallocate(Nphases);
            PhaseTransformationStretches.Reallocate(Nphases);
            PhaseKappa.Reallocate({Nphases, Ncomp});
            PhaseLambda.Reallocate({Nphases, Ncomp});
            PhaseGamma.Reallocate(Nphases);
            PhaseAlpha.Reallocate(Nphases);
            Tref.Reallocate(Nphases);
            Cref.Reallocate({Nphases, Ncomp});
        }
        for(size_t alpha = 0; alpha != Nphases; alpha++)
        {
            PhaseAlpha[alpha] = rhs.PhaseAlpha[alpha];
            PhaseGamma[alpha] = rhs.PhaseGamma[alpha];

            Tref[alpha] = rhs.Tref[alpha];

            PhaseTransformationStretches[alpha] = rhs.PhaseTransformationStretches[alpha];
            PhaseElasticConstants[alpha] = rhs.PhaseElasticConstants[alpha];

            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                PhaseLambda({alpha, comp}) = rhs.PhaseLambda({alpha, comp});
                PhaseKappa({alpha, comp}) = rhs.PhaseKappa({alpha, comp});
                Cref({alpha, comp}) = rhs.Cref({alpha, comp});
            }
        }
        GrainElasticConstants = rhs.GrainElasticConstants;
        GrainTransformationStretches = rhs.GrainTransformationStretches;

        GrainAlpha = rhs.GrainAlpha;
        GrainGamma = rhs.GrainGamma;
        GrainLambda = rhs.GrainLambda;
        GrainKappa = rhs.GrainKappa;

        RemeshedStrain = rhs.RemeshedStrain;
        StrainToRemesh = rhs.StrainToRemesh;
        AverageStrain = rhs.AverageStrain;
        AppliedStress = rhs.AppliedStress;
        AppliedStrain = rhs.AppliedStrain;
        AppliedStrainRate = rhs.AppliedStrainRate;

        AppliedStressMask = rhs.AppliedStressMask;
        AppliedStrainMask = rhs.AppliedStrainMask;
        AppliedStrainRateMask = rhs.AppliedStrainRateMask;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells(),)
        {
            DeformationGradientsTotal(i,j,k) = rhs.DeformationGradientsTotal(i,j,k);
            DeformationGradientsEigen(i,j,k) = rhs.DeformationGradientsEigen(i,j,k);
            Stresses(i,j,k) = rhs.Stresses(i,j,k);

            EffectiveElasticConstants(i,j,k) = rhs.EffectiveElasticConstants(i,j,k);
            Displacements(i,j,k) = rhs.Displacements(i,j,k);

            if (AnyPlasticity)
            {
                VelocityGradientsPlastic(i,j,k) = rhs.VelocityGradientsPlastic(i,j,k);
                DeformationGradientsPlastic(i,j,k) = rhs.DeformationGradientsPlastic(i,j,k);
            }
            if (LargeDeformations)
            {
                StressIncrements(i,j,k) = rhs.StressIncrements(i,j,k);
                VelocityGradientsTotal(i,j,k) = rhs.VelocityGradientsTotal(i,j,k);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void ElasticProperties::SetBaseTransformationStretches(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, TransformationStretches, TransformationStretches.Bcells(),)
    {
        TransformationStretches(i,j,k).clear();
        for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            TransformationStretches(i,j,k).set_value(alpha->index, GrainTransformationStretches[alpha->index]);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ElasticProperties::ApplyLocalRotationsToDeformationGradientsEigen(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DeformationGradientsEigen, DeformationGradientsEigen.Bcells(),)
    {
        DeformationGradientsEigen(i,j,k).rotate(LocalRotations(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetDeformationGradientsEigen(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EffectiveElasticConstants, EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            DeformationGradientsEigen(i,j,k).set_to_zero();
            for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                DeformationGradientsEigen(i, j, k) += TransformationStretches(i, j, k).get_value(alpha->index) * alpha->value;
            }
        }
        else
        {
            DeformationGradientsEigen(i,j,k) = TransformationStretches(i, j, k).front().value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::AddVegardsExpansion(const PhaseField& Phase, const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, TransformationStretches, TransformationStretches.Bcells(),)
    {
        for(auto alpha  = TransformationStretches(i,j,k).begin();
                 alpha != TransformationStretches(i,j,k).end(); ++alpha)
        {
            size_t pIndex = Phase.FieldsProperties[alpha->index].Phase;
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.MoleFractions(i,j,k,{pIndex, comp}) - Cref({pIndex, comp}));
                alpha->value += alpha->value.Hadamard_product(GrainLambda({alpha->index, comp}))*delta;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::AddThermalExpansion(const PhaseField& Phase, const Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,TransformationStretches,TransformationStretches.Bcells(),)
    {
        for(auto alpha  = TransformationStretches(i,j,k).begin();
                 alpha != TransformationStretches(i,j,k).end(); ++alpha)
        {
            int pIndex = Phase.FieldsProperties[alpha->index].Phase;
            double delta = (Tx(i,j,k) - Tref[pIndex]);
            alpha->value += alpha->value.Hadamard_product(GrainAlpha[alpha->index])*delta;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::AddInterfaceStressContribution(const PhaseField& Phase, const InterfaceProperties& IP)
{
    //Raphael Schiedung, Ingo Steinbach, and Fathollah Varnik.
    //"Multi-phase-field method for surface tension induced elasticity."
    //Physical Review B 97.3 (2018): 035410.
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsEigen,0,)
    if (Phase.Fields(i,j,k).interface())
    {
        const double pre1 = 4.0/Phase.Grid.Eta;
        const double pre2 = Phase.Grid.Eta * Phase.Grid.Eta/Pi;

        vStress locInterfaceStress;
        // Calculate phase gradients and normals
        const NodeAB<dVector3,dVector3>& locNormals   = Phase.Normals(i,j,k);
        const NodeA<dVector3>& locGradients = Phase.Fields(i,j,k).get_gradients();

        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            // Calculate phase field part of the interface stress
            const double PhiA = abs(Phase.Fields(i,j,k).get_value(it->indexA));
            const double PhiB = abs(Phase.Fields(i,j,k).get_value(it->indexB));

            if (PhiA * PhiB < 0.02) continue;

            const double interpol = pre1 * (pre2 *
                    (locGradients.get_value(it->indexA) *
                     locGradients.get_value(it->indexB)) + PhiA * PhiB);

            // Calculate local projection matrix
            dMatrix3x3 Projection;
            Projection.set_to_unity();
            Projection(0,0) -= it->value1[0]*it->value1[0];
            Projection(1,0) -= it->value1[1]*it->value1[0];
            Projection(2,0) -= it->value1[2]*it->value1[0];
            Projection(0,1) -= it->value1[0]*it->value1[1];
            Projection(1,1) -= it->value1[1]*it->value1[1];
            Projection(2,1) -= it->value1[2]*it->value1[1];
            Projection(0,2) -= it->value1[0]*it->value1[2];
            Projection(1,2) -= it->value1[1]*it->value1[2];
            Projection(2,2) -= it->value1[2]*it->value1[2];

            // Calculate local interface stress
            locInterfaceStress -= VoigtStress(Projection) * interpol *
                    IP.Properties(i,j,k).get_energy(it->indexA,it->indexB);
        }

        //TODO account for derivatives Interface energy with respect to the strain
        vStrain locEigenStrainInterface =
            EffectiveElasticConstants(i,j,k).inverted() * locInterfaceStress *(-1.);

        if(StrainModel == StrainModels::Small)
        {
            DeformationGradientsEigen(i,j,k) += locEigenStrainInterface.tensor();
        }
//        else
//        {
//            vStrain locEigenStrain = StrainGreenLagrange(DeformationGradientsEigen(i,j,k));
//            locEigenStrain += locEigenStrainInterface;
//
//            DeformationGradientsEigen(i,j,k) =
//            sqrtM3x3(locEigenStrain.tensor()*2.0 + dMatrix3x3::UnitTensor());
//        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetBaseElasticConstants(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ElasticConstants, ElasticConstants.Bcells(),)
    {
        ElasticConstants(i,j,k).clear();
        for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            ElasticConstants(i,j,k).set_value(alpha->index, GrainElasticConstants[alpha->index]);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ElasticProperties::AddThermalSoftening(const PhaseField& Phase, const Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ElasticConstants, ElasticConstants.Bcells(),)
    {
        for(auto alpha  = ElasticConstants(i,j,k).begin();
                 alpha != ElasticConstants(i,j,k).end(); ++alpha)
        {
            size_t pIndex = Phase.FieldsProperties[alpha->index].Phase;
            double delta = (Tx(i,j,k) - Tref[pIndex]);
            alpha->value += alpha->value.Hadamard_product(GrainGamma[alpha->index])*delta;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::AddChemicalSoftening(const PhaseField& Phase, const Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, ElasticConstants, ElasticConstants.Bcells(),)
    {
        for(auto alpha  = ElasticConstants(i,j,k).begin();
                 alpha != ElasticConstants(i,j,k).end(); ++alpha)
        {
            size_t pIndex = Phase.FieldsProperties[alpha->index].Phase;
            for(size_t comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.MoleFractions(i,j,k,{pIndex, comp}) - Cref({pIndex, comp}));
                alpha->value += alpha->value.Hadamard_product(GrainKappa({alpha->index, comp}))*delta;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ElasticProperties::SetEffectiveProperties(const PhaseField& Phase)
{
    SetGrainsProperties(Phase);
    SetBaseElasticConstants(Phase);
    SetBaseTransformationStretches(Phase);

    SetEffectiveElasticConstants(Phase);
    SetDeformationGradientsEigen(Phase);
    if(LargeDeformations)
    {
        ApplyLocalRotationsToElasticConstants();
        ApplyLocalRotationsToDeformationGradientsEigen();
    }
}

void ElasticProperties::SetEffectiveProperties(const PhaseField& Phase, const Composition& Cx)
{
    SetGrainsProperties(Phase);
    SetBaseElasticConstants(Phase);
    SetBaseTransformationStretches(Phase);
    if(ChemoMechanicalCoupling)
    {
        AddChemicalSoftening(Phase, Cx);
        AddVegardsExpansion(Phase, Cx);
    }

    SetEffectiveElasticConstants(Phase);
    SetDeformationGradientsEigen(Phase);
    if(LargeDeformations)
    {
        ApplyLocalRotationsToElasticConstants();
        ApplyLocalRotationsToDeformationGradientsEigen();
    }
}

void ElasticProperties::SetEffectiveProperties(const PhaseField& Phase, const Temperature& Tx)
{
    SetGrainsProperties(Phase);
    SetBaseElasticConstants(Phase);
    SetBaseTransformationStretches(Phase);
    if(ThermoMechanicalCoupling)
    {
        AddThermalSoftening(Phase, Tx);
        AddThermalExpansion(Phase, Tx);
    }

    SetEffectiveElasticConstants(Phase);
    SetDeformationGradientsEigen(Phase);
    if(LargeDeformations)
    {
        ApplyLocalRotationsToElasticConstants();
        ApplyLocalRotationsToDeformationGradientsEigen();
    }
}

void ElasticProperties::SetEffectiveProperties(const PhaseField& Phase, const Composition& Cx, const Temperature& Tx)
{
    SetGrainsProperties(Phase);
    SetBaseElasticConstants(Phase);
    SetBaseTransformationStretches(Phase);
    if(ChemoMechanicalCoupling)
    {
        AddChemicalSoftening(Phase, Cx);
        AddVegardsExpansion(Phase, Cx);
    }
    if(ThermoMechanicalCoupling)
    {
        AddThermalSoftening(Phase, Tx);
        AddThermalExpansion(Phase, Tx);
    }
    SetEffectiveElasticConstants(Phase);
    SetDeformationGradientsEigen(Phase);
    if(LargeDeformations)
    {
        ApplyLocalRotationsToElasticConstants();
        ApplyLocalRotationsToDeformationGradientsEigen();
    }
}
void ElasticProperties::SetEffectiveProperties(const PhaseField& Phase, const InterfaceProperties& IP)
{
    SetGrainsProperties(Phase);
    SetBaseElasticConstants(Phase);
    SetBaseTransformationStretches(Phase);

    SetEffectiveElasticConstants(Phase);
    SetDeformationGradientsEigen(Phase);
    if(LargeDeformations)
    {
        ApplyLocalRotationsToElasticConstants();
        ApplyLocalRotationsToDeformationGradientsEigen();
    }
}

dMatrix6x6 ElasticProperties::MAXElasticConstants(void) const
{
    dMatrix6x6 Cij_max;
    for(size_t alpha = 0; alpha != GrainElasticConstants.size(); alpha++)
    for(int n = 0; n < 6; n++)
    for(int m = 0; m < 6; m++)
    {
        Cij_max(n,m) = std::max(fabs(GrainElasticConstants[alpha](n,m)), Cij_max(n,m));
    }
    return Cij_max;
}

dMatrix6x6 ElasticProperties::AverageElasticConstants(void) const
{
    dMatrix6x6 locAverageElasticConstants;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EffectiveElasticConstants, 0, reduction(dMatrix6x6SUM: locAverageElasticConstants))
    {
        locAverageElasticConstants += EffectiveElasticConstants(i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE,locAverageElasticConstants.data(),36,OP_MPI_DOUBLE,OP_MPI_SUM,OP_MPI_COMM_WORLD);
#endif

    return locAverageElasticConstants /= double(Grid.TotalNumberOfCells());
}

void ElasticProperties::ApplyLocalRotationsToElasticConstants(void)
{
    if(LargeDeformations)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EffectiveElasticConstants, EffectiveElasticConstants.Bcells(),)
    {
        EffectiveElasticConstants(i,j,k).rotate(LocalRotations(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetEffectiveElasticConstants(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EffectiveElasticConstants, EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            switch(ElasticityModel)
            {
                case ElasticityModels::Khachaturyan:
                case ElasticityModels::Voigt:
                case ElasticityModels::Rank1:
                case ElasticityModels::Rank1NL:
                {
                    EffectiveElasticConstants(i,j,k).set_to_zero();
                    for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                             alpha != Phase.Fields(i,j,k).cend(); ++alpha)
                    {
                        EffectiveElasticConstants(i, j, k) += ElasticConstants(i, j, k).get_value(alpha->index) * alpha->value;
                    }
                    break;
                }
                case ElasticityModels::Steinbach:
                case ElasticityModels::Reuss:
                {
                    EffectiveElasticConstants(i,j,k).set_to_zero();
                    for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                             alpha != Phase.Fields(i,j,k).cend(); ++alpha)
                    {
                        EffectiveElasticConstants(i, j, k) += ElasticConstants(i, j, k).get_value(alpha->index).inverted() * alpha->value;
                    }
                    EffectiveElasticConstants(i, j, k).invert();
                    break;
                }
            }
        }
        else
        {
            EffectiveElasticConstants(i,j,k) = ElasticConstants(i, j, k).front().value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateDrivingForce(const PhaseField& Phase, DrivingForce& dGab)
{
    if(ElasticityModel == ElasticityModels::Rank1)
    {
        CalculateDeformationJumps(Phase);
    }
    if(ElasticityModel == ElasticityModels::Rank1NL)
    {
        CalculateDeformationJumps(Phase);
        CalculateNonlocalJumpContribution(Phase);
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal,0,)
    if(Phase.Fields(i,j,k).interface())
    {
        const vStrain ElasticStrains = EffectiveElasticConstants(i,j,k).inverted()*Stresses(i, j, k);
        NodeAB<dVector3,dVector3> locNormals = Phase.Normals(i,j,k);

        for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta  = alpha + 1;
                  beta != Phase.Fields(i, j, k).cend();  ++beta)
        {
            dMatrix3x3 locStretchesAlpha = TransformationStretches(i,j,k).get_value(alpha->index);
            dMatrix3x3 locStretchesBeta  = TransformationStretches(i,j,k).get_value( beta->index);

            vStrain locEigenStrainDifference =
                EigenStrainDifference(locStretchesAlpha, locStretchesBeta, i,j,k);

            double dG_AB = 0.0;
            switch(ElasticityModel)
            {
                case ElasticityModels::Khachaturyan:
                {
                    dG_AB += (ElasticStrains*
                             ((ElasticConstants(i,j,k).get_value( beta->index) -
                               ElasticConstants(i,j,k).get_value(alpha->index))*
                              ElasticStrains))*0.5;

                    dG_AB -= Stresses(i, j, k)*locEigenStrainDifference;
                    break;
                }
                case ElasticityModels::Steinbach:
                case ElasticityModels::Reuss:
                {
                    dG_AB += (Stresses(i, j, k)*
                             ((ElasticConstants(i,j,k).get_value(alpha->index).inverted() -
                               ElasticConstants(i,j,k).get_value( beta->index).inverted())*
                              Stresses(i, j, k)))*0.5;

                    dG_AB -= Stresses(i, j, k)*locEigenStrainDifference;
                    break;
                }
                case ElasticityModels::Voigt:                                   // Defined in small strain setting
                {
                    vStrain locStrain = StrainSmall(DeformationGradientsTotal(i,j,k))
                                      - StrainSmall(DeformationGradientsPlastic(i,j,k));

                    vStrain locEigenStrainAlpha = StrainSmall(TransformationStretches(i,j,k).get_value(alpha->index));
                    vStrain locEigenStrainBeta  = StrainSmall(TransformationStretches(i,j,k).get_value( beta->index));

                    dG_AB += ((locStrain - locEigenStrainBeta)*
                              (ElasticConstants(i,j,k).get_value(beta->index)*
                              (locStrain - locEigenStrainBeta)) -

                              (locStrain - locEigenStrainAlpha)*
                              (ElasticConstants(i,j,k).get_value(alpha->index)*
                              (locStrain - locEigenStrainAlpha)))*0.5;
                    break;
                }
                case ElasticityModels::Rank1:                                   // Defined in small strain setting
                case ElasticityModels::Rank1NL:                                 // Defined in small strain setting
                {
                    if(alpha->value > DBL_EPSILON and beta->value > DBL_EPSILON)
                    {
                        vStrain locStrain = StrainSmall(DeformationGradientsTotal(i,j,k))
                                          - StrainSmall(DeformationGradientsPlastic(i,j,k));

                        double scale = 1.0;
                        if(Phase.FieldsProperties[alpha->index].Stage != GrainStages::Stable)
                        {
                            scale *= Phase.FieldsProperties[alpha->index].VolumeRatio;
                        }
                        if(Phase.FieldsProperties[beta->index].Stage != GrainStages::Stable)
                        {
                            scale *= Phase.FieldsProperties[beta->index].VolumeRatio;
                        }

                        dMatrix3x3 locFjumpAlpha = DeformationJumps(i,j,k).get_sym2(alpha->index,  beta->index);
                        dMatrix3x3 locFjumpBeta  = DeformationJumps(i,j,k).get_sym2( beta->index, alpha->index);

//                        for(auto gamma  = beta + 1;
//                                 gamma != Phase.Fields(i, j, k).cend(); ++gamma)
//                        if(gamma->value > DBL_EPSILON)
//                        {
//                            locFjumpAlpha += DeformationJumps(i,j,k).get_asym2(alpha->index, gamma->index);
//                            locFjumpBeta  += DeformationJumps(i,j,k).get_asym2( beta->index, gamma->index);
//                        }

                        double sign = 1.0;
                        if(alpha->index > beta->index)
                        {
                            sign = -1.0;
                        }

                        vStrain locStrainJumpA = StrainSmall(locFjumpAlpha + dMatrix3x3::UnitTensor())*scale*sign;
                        vStrain locStrainJumpB = StrainSmall(locFjumpBeta  + dMatrix3x3::UnitTensor())*scale*(-sign);

                        vStrain StrainAlpha = locStrain - locStrainJumpA* beta->value/(alpha->value + beta->value);// TODO: check effect of normalization by (alpha->value + beta->value)
                        vStrain StrainBeta  = locStrain - locStrainJumpB*alpha->value/(alpha->value + beta->value);// TODO: check effect of normalization by (alpha->value + beta->value)

                        vStrain locEigenStrainAlpha = StrainSmall(locStretchesAlpha);
                        vStrain locEigenStrainBeta  = StrainSmall(locStretchesBeta);

                        dG_AB += ((StrainBeta - locEigenStrainBeta)*
                                  (ElasticConstants(i,j,k).get_value(beta->index)*
                                  (StrainBeta - locEigenStrainBeta)) -

                                  (StrainAlpha - locEigenStrainAlpha)*
                                  (ElasticConstants(i,j,k).get_value(alpha->index)*
                                  (StrainAlpha - locEigenStrainAlpha)))*0.5*scale;

                        dG_AB -= Stresses(i,j,k)*(locStrainJumpA - locStrainJumpB)*0.5*scale;

                        if(ElasticityModel == ElasticityModels::Rank1NL)
                        {
                            if(Grid.dNx)
                            {
                                dG_AB -= (NonLocalDeformationJumpTerm(i+1,j,k).get_asym1(alpha->index, beta->index)[0] - NonLocalDeformationJumpTerm(i-1,j,k).get_asym1(alpha->index, beta->index)[0])*scale/(2.0*Grid.dx);
                            }
                            if(Grid.dNy)
                            {
                                dG_AB -= (NonLocalDeformationJumpTerm(i,j+1,k).get_asym1(alpha->index, beta->index)[1] - NonLocalDeformationJumpTerm(i,j-1,k).get_asym1(alpha->index, beta->index)[1])*scale/(2.0*Grid.dx);
                            }
                            if(Grid.dNz)
                            {
                                dG_AB -= (NonLocalDeformationJumpTerm(i,j,k+1).get_asym1(alpha->index, beta->index)[2] - NonLocalDeformationJumpTerm(i,j,k-1).get_asym1(alpha->index, beta->index)[2])*scale/(2.0*Grid.dx);
                            }
                        }

                        //Use Khachaturyan's driving force for small nuclei
                        if((1.0 - scale) > FLT_EPSILON)
                        {
                            dG_AB += (ElasticStrains*
                                     ((ElasticConstants(i,j,k).get_value( beta->index) -
                                       ElasticConstants(i,j,k).get_value(alpha->index))*
                                      ElasticStrains))*0.5*(1.0 - scale);

                            dG_AB -= Stresses(i, j, k)*locEigenStrainDifference*(1.0 - scale);
                        }
                    }
                    break;
                }
            }
            dGab.Force(i,j,k).add_raw(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateDeformationJumps(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationJumps,0,)
    {
        DeformationJumps(i,j,k).clear();

        if(Phase.Fields(i,j,k).interface())
        {
            vStrain locStrain = StrainSmall(DeformationGradientsTotal(i,j,k))
                              - StrainSmall(DeformationGradientsPlastic(i,j,k));
            NodeAB<dVector3,dVector3> locNormals = Phase.Normals(i,j,k);
            double residual = 0.0;
            do
            {
                NodeAB<dVector3,dMatrix3x3> locJumps;
                residual = 0.0;
                for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                         alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
                for(auto  beta  = alpha + 1;
                          beta != Phase.Fields(i, j, k).cend();  ++beta)
                if(alpha->value > DBL_EPSILON and beta->value > DBL_EPSILON)
                {
                    dVector3 locNormalAB = locNormals.get_asym1(alpha->index, beta->index);

                    vStrain dStrainA = locStrain - StrainSmall(TransformationStretches(i,j,k).get_value(alpha->index));
                    vStrain dStrainB = locStrain - StrainSmall(TransformationStretches(i,j,k).get_value( beta->index));

                    /*for(auto gamma  = beta + 1;
                             gamma != Phase.Fields(i, j, k).cend(); ++gamma)
                    if(gamma->value > DBL_EPSILON)
                    {
                        dVector3 locNormalAG = locNormals.get_asym(alpha->index, gamma->index);
                        dVector3 locNormalBG = locNormals.get_asym( beta->index, gamma->index);

                        dVector3 locJumpA = DeformationJumps(i,j,k).get_asym1(alpha->index, gamma->index);
                        dVector3 locJumpB = DeformationJumps(i,j,k).get_asym1( beta->index, gamma->index);

                        dMatrix3x3 locFjumpA = locJumpA.dyadic(locNormalAG);
                        dMatrix3x3 locFjumpB = locJumpB.dyadic(locNormalBG);

                        dStrainA -= VoigtStrain(locFjumpA + locFjumpA.transposed())*gamma->value*0.5;
                        dStrainB -= VoigtStrain(locFjumpB + locFjumpB.transposed())*gamma->value*0.5;
                    }*/

                    dMatrix6x6 Cij_alpha = ElasticConstants(i,j,k).get_value(alpha->index);
                    dMatrix6x6 Cij_beta  = ElasticConstants(i,j,k).get_value( beta->index);
                    dMatrix6x6 Cij = Cij_alpha*(1.0 - alpha->value) +
                                     Cij_beta *(1.0 -  beta->value);

                    if(locNormalAB.abs() > DBL_EPSILON)
                    {
                        dMatrix3x3 projCij = Cij.project(locNormalAB);
                        dVector3 locJump = projCij.inverted()*(Cij_alpha*dStrainA - Cij_beta*dStrainB).tensor()*locNormalAB;
                        dMatrix3x3 locDefJump = locJump.dyadic(locNormalAB);
                        locJumps.set_asym1(alpha->index, beta->index, locJump);
                        locJumps.set_asym2(alpha->index, beta->index, locDefJump);

                        dVector3 locJumpOLD = DeformationJumps(i,j,k).get_asym1(alpha->index, beta->index);
                        residual = max((locJumpOLD - locJump).abs()*beta->value, residual);
                    }
                }
                DeformationJumps(i,j,k) = locJumps;
            }
            while(residual > 1.0e-6);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::CalculateNonlocalJumpContribution(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationJumps,0,)
    {
        NonLocalDeformationJumpTerm(i,j,k).clear();

        if(Phase.Fields(i,j,k).interface())
        {
            vStrain locStrain = StrainSmall(DeformationGradientsTotal(i,j,k))
                              - StrainSmall(DeformationGradientsPlastic(i,j,k));
            NodeAB<dVector3,dVector3> locNormals = Phase.Normals(i,j,k);

            for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i, j, k).cend();  ++beta)
            if(alpha->value > DBL_EPSILON and beta->value > DBL_EPSILON)
            {
                dVector3 locNormalAB = locNormals.get_asym1(alpha->index, beta->index);

                dVector3   defJumpV = DeformationJumps(i,j,k).get_asym1(alpha->index, beta->index);
                dMatrix3x3 defJumpT = DeformationJumps(i,j,k).get_asym2(alpha->index, beta->index);
                dMatrix3x3 strainJumpAB = (defJumpT.transposed() + defJumpT)*0.5;

                vStrain dStrainA = (locStrain - StrainSmall(TransformationStretches(i,j,k).get_value(alpha->index)) - VoigtStrain(strainJumpAB)*beta->value);
                vStrain dStrainB = (locStrain - StrainSmall(TransformationStretches(i,j,k).get_value( beta->index)) + VoigtStrain(strainJumpAB)*alpha->value);

                vStress stressAlpha = ElasticConstants(i,j,k).get_value(alpha->index) * dStrainA;
                vStress stressBeta  = ElasticConstants(i,j,k).get_value( beta->index) * dStrainB;

                dMatrix3x3 locUnity = dMatrix3x3::UnitTensor() - locNormalAB.dyadic(locNormalAB);

                dVector3 nonLocalTerm = (((stressAlpha - stressBeta).tensor() * locUnity) * defJumpV) * (alpha->value*beta->value)/alpha->gradient.length();

                NonLocalDeformationJumpTerm(i,j,k).set_asym1(alpha->index, beta->index, nonLocalTerm);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void ElasticProperties::CalculateChemicalPotentialContribution(
                                        const PhaseField& Phase,
                                        EquilibriumPartitionDiffusionBinary& DF)
{
    /** This function calculates the partial derivative of the mechanical
    energy with respect to the composition. This is later used in the EQPD-model
    to calculate an additional diffusion flux. */

    const size_t comp = DF.Comp;
    if(ChemoMechanicalCoupling)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DF.dMu, DF.dMu.Bcells(),)
    {
        vStrain locStrain = ElasticStrains(i,j,k);

        for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        {
            double locdMuEl = 0.0;
            vStrain locLambda = VoigtStrain(GrainLambda({alpha->index, comp}));

            for(int ii = 0; ii < 6; ii++)
            {
                locdMuEl -= locLambda[ii]*Stresses(i, j, k)[ii];
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * GrainKappa({alpha->index, comp})(ii,jj)*
                                      locStrain[ii]*
                                      GrainElasticConstants[alpha->index](ii,jj)*
                                      locStrain[jj];
                }
            }
            size_t pIndex = Phase.FieldsProperties[alpha->index].Phase;
            DF.dMu(i,j,k,{pIndex,comp}) += locdMuEl*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi, const BoundaryConditions& BC, const double dt, const double tStep)
{
    if(LargeDeformations)
    {
    //    if (DeformationGradientsTotalAdvectionDot.IsNotAllocated()) DeformationGradientsTotalAdvectionDot.Allocate(DeformationGradientsTotal);
    //    if (DeformationGradientsEigenAdvectionDot.IsNotAllocated()) DeformationGradientsEigenAdvectionDot.Allocate(DeformationGradientsEigen);
    //    if (DeformationGradientsPlasticAdvectionDot.IsNotAllocated()) DeformationGradientsPlasticAdvectionDot.Allocate(DeformationGradientsPlastic);
        if (StressesAdvectionDot.IsNotAllocated()) StressesAdvectionDot.Allocate(Grid, Grid.Bcells);
        if (LocalRotationsAdvectionDot.IsNotAllocated()) LocalRotationsAdvectionDot.Allocate(Grid, Grid.Bcells);

    //    if(not DeformationGradientsTotalAdvectionDot.IsSize(Nx, Ny, Nz)) DeformationGradientsTotalAdvectionDot.Reallocate(Nx, Ny, Nz);
    //    if(not DeformationGradientsEigenAdvectionDot.IsSize(Nx, Ny, Nz)) DeformationGradientsEigenAdvectionDot.Reallocate(Nx, Ny, Nz);
    //    if(not DeformationGradientsPlasticAdvectionDot.IsSize(Nx, Ny, Nz)) DeformationGradientsPlasticAdvectionDot.Reallocate(Nx, Ny, Nz);
        if(not StressesAdvectionDot.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) StressesAdvectionDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
        if(not LocalRotationsAdvectionDot.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) LocalRotationsAdvectionDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    //    Adv.AdvectField(DeformationGradientsTotal, DeformationGradientsTotalAdvectionDot, Vel, BC, dx, dt);
    //    Adv.AdvectField(DeformationGradientsEigen, DeformationGradientsEigenAdvectionDot, Vel, BC, dx, dt);
    //    Adv.AdvectField(DeformationGradientsPlastic, DeformationGradientsPlasticAdvectionDot, Vel, BC, dx, dt);

        Adv.AdvectField(Stresses, StressesAdvectionDot, Vel, BC, dt);
        Adv.AdvectField(LocalRotations, LocalRotationsAdvectionDot, Vel, BC, dt);
    
        ConsoleOutput::WriteStandard(thisclassname, "Advected");
    }
}

void ElasticProperties::WriteTotalRotationsVTK(Settings& locSettings, const int tStep, const int precision) const
{
    auto CalculateRotations = [&](int i, int j, int k)
    {
        double locAngle = 0.0;
        dVector3 locAxis;
        Tools::getAxisAngle(DeformationGradientsTotal(i,j,k),locAxis,locAngle);
        locAngle *= 180.0 / Pi;
        return std::make_pair(locAngle, locAxis);
    };

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Angle", [=](int i,int j,int k){return CalculateRotations(i,j,k).first;}});
    ListOfFields.push_back((VTK::Field_t) {"Axis",  [=](int i,int j,int k){return CalculateRotations(i,j,k).second;}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Rotations_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with rotations data", "application/xml");
    #endif
}

vStrain ElasticProperties::AveragePlasticStrain(void) const
{
    vStrain locAveragePlasticStrain;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsPlastic,0, reduction(vStrainSUM: locAveragePlasticStrain))
    {
        locAveragePlasticStrain += PlasticStrains(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    vStrain tmpAveragePlasticStrain = locAveragePlasticStrain;
    OP_MPI_Allreduce(tmpAveragePlasticStrain.data(), locAveragePlasticStrain.data(), 6, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
#endif
    return locAveragePlasticStrain/(Grid.TotalNumberOfCells());
}

bool ElasticProperties::Write(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_", tStep, ".dat");
#endif
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteExit("File " + FileName + " could not be created",thisclassname);
        OP_Exit(EXIT_FAILURE);
    };

    out.write(reinterpret_cast<const char*>(&(AverageStrain)), sizeof(vStrain));

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,0)
    {
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
        out.write(reinterpret_cast<const char*>(&(DeformationGradientsPlastic(i,j,k))), sizeof(dMatrix3x3));
        out.write(reinterpret_cast<const char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END

    out.close();
    return true;
}

bool ElasticProperties::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
#ifdef MPI_PARALLEL
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_" + std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, thisclassname + "_", tStep, ".dat");
#endif

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        ConsoleOutput::WriteWarning("File " + FileName + " could not be opened",thisclassname);
        return false;
    };

    inp.read(reinterpret_cast<char*>(&(AverageStrain)), sizeof(vStrain));

    STORAGE_LOOP_BEGIN(i,j,k,DeformationGradientsTotal,DeformationGradientsTotal.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsTotal(i,j,k))), sizeof(dMatrix3x3));
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsEigen(i,j,k))), sizeof(dMatrix3x3));
        inp.read(reinterpret_cast<char*>(&(DeformationGradientsPlastic(i,j,k))), sizeof(dMatrix3x3));
        inp.read(reinterpret_cast<char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END
    SetBoundaryConditions(BC);
    ConsoleOutput::Write(thisclassname, "Binary input loaded");
    return true;
}

void ElasticProperties::WriteTotalStrainsVTK(Settings& locSettings,
                                             const int tStep,
                                             const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"TotalStrain", [this](int i,int j,int k){return TotalStrains(i, j, k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "TotalStrain_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with total strain data", "application/xml");
    #endif
}

void ElasticProperties::WritePlasticStrainsVTK(Settings& locSettings,
                                               const int tStep,
                                               const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"PlasticStrain", [this](int i,int j,int k){return PlasticStrains(i, j, k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "PlasticStrain_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with plastic strain data", "application/xml");
    #endif
}

void ElasticProperties::WriteStressesVTK(Settings& locSettings,
                                         const int tStep,
                                         const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Stresses", [this](int i,int j,int k){return Stresses(i, j, k);}});
    ListOfFields.push_back((VTK::Field_t) {"von Mises", [this](int i,int j,int k){return Stresses(i, j, k).Mises();}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return Stresses(i, j, k).Pressure();}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Stresses_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with stress data", "application/xml");
    #endif
}

void ElasticProperties::WriteElasticStrainsVTK(Settings& locSettings,
                                               const int tStep,
                                               const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"ElasticStrains", [this](int i,int j,int k){return ElasticStrains(i,j,k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ElasticStrains_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with elastic strain data", "application/xml");
    #endif
}

void ElasticProperties::WriteEigenStrainsVTK(Settings& locSettings,
                                             const int tStep,
                                             const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"EigenStrains", [this](int i,int j,int k){return EigenStrains(i,j,k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "EigenStrains_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with Eigen strain data", "application/xml");
    #endif
}

void ElasticProperties::WriteDeformationGradientsTotalVTK(Settings& locSettings,
                                                          const int tStep,
                                                          const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"DeformationGradientsTotal", [this](int i,int j,int k){return DeformationGradientsTotal(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "DeformationGradientsTotal_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with deformation gradient data", "application/xml");
    #endif
}

void ElasticProperties::WriteEffectiveElasticConstantsVTK(Settings& locSettings,
                                                          const int tStep,
                                                          const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"EffectiveElasticConstants", [this](int i,int j,int k){return EffectiveElasticConstants(i, j, k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "EffectiveElasticConstants_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
    #ifndef WIN32       
    locSettings.Meta.AddPart(Filename, "File", "VTK file with elastic constants data", "application/xml");
    #endif
}

void ElasticProperties::WriteForceDensityVTK(Settings& locSettings,
                                             const int tStep,
                                             const int precision) const
{
    if(ConsiderExternalForces)
    {
        std::vector<VTK::Field_t> ListOfFields;
        ListOfFields.push_back((VTK::Field_t) {"ForceDensity", [this](int i,int j,int k){return ForceDensity(i, j, k);}});

        std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ForceDensity_", tStep, ".vts");
        VTK::Write(Filename, locSettings, ListOfFields, precision);
        #ifndef WIN32       
        locSettings.Meta.AddPart(Filename, "File", "VTK file with force density data", "application/xml");
        #endif
    }
    else
    {
        ConsoleOutput::WriteExit("ForceDensity Storage not allocated!", "ElasticProperties", "WriteForceDensityVTK");
        exit(EXIT_FAILURE);
    }
}

void ElasticProperties::WriteForceDensityDistortedVTK(Settings& locSettings,
                                                      const int tStep,
                                                      const int precision) const
{
    if(ConsiderExternalForces)
    {
        std::vector<VTK::Field_t> ListOfFields;
        ListOfFields.push_back((VTK::Field_t) {"ForceDensity", [this](int i,int j,int k){return ForceDensity(i, j, k);}});

        std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ForceDensity_", tStep, ".vts");
        VTK::WriteDistorted(Filename, locSettings, *this, ListOfFields, precision);
        #ifndef WIN32       
        locSettings.Meta.AddPart(Filename, "File", "VTK file with distorted force density data", "application/xml");
        #endif
    }
    else
    {
        ConsoleOutput::WriteExit("ForceDensity Storage not allocated!", "ElasticProperties", "WriteForceDensityVTK");
        exit(EXIT_FAILURE);
    }
}
void ElasticProperties::WriteStressIncrementsVTK(Settings& locSettings,
                                                 const int tStep,
                                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"StressIncrements", [this](int i,int j,int k){return StressIncrements(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t) {"Pressure", [this](int i,int j,int k){return StressIncrements(i,j,k).Pressure();}});
    ListOfFields.push_back((VTK::Field_t) {"vonMises", [this](int i,int j,int k){return StressIncrements(i,j,k).Mises();}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "StressIncrements_", tStep, ".vts");
    VTK::WriteDistorted(Filename, locSettings, *this, ListOfFields, precision);
}

void ElasticProperties::WriteVelocityGradientVTK(Settings& locSettings,
                                                 const int tStep,
                                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"VelocityGradients_", [this](int i,int j,int k){return VelocityGradientsTotal(i, j, k);}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "VelocityGradients_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void ElasticProperties::WriteDisplacementsVTK(Settings& OPSettings,
                                              const int tStep,
                                              const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"U", [this](int i, int j, int k){return Displacements(i,j,k)*Grid.dx;}});

    std::string Filename = FileInterface::MakeFileName(OPSettings.VTKDir, "Displacements_", tStep, ".vts");
    VTK::Write(Filename, OPSettings, ListOfFields, precision);
}

void ElasticProperties::WriteDeformationJumpsVTK(Settings& OPSettings,
                                                 const int tStep,
                                                 const int idxA,
                                                 const int idxB,
                                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"DeformationJumps", [this, idxA, idxB](int i, int j, int k){return DeformationJumps(i,j,k).get_asym1(idxA,idxB);}});

    std::string Filename = FileInterface::MakeFileName(OPSettings.VTKDir, "DeformationJumps_", tStep, ".vts");
    VTK::Write(Filename, OPSettings, ListOfFields, precision);
}

void ElasticProperties::WriteLocalRotationsVTK(Settings& locSettings, const int tStep,
                            const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"RQ_0", [this](int i,int j,int k){return LocalRotations(i,j,k)[0];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_1", [this](int i,int j,int k){return LocalRotations(i,j,k)[1];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_2", [this](int i,int j,int k){return LocalRotations(i,j,k)[2];}});
    ListOfFields.push_back((VTK::Field_t) {"RQ_3", [this](int i,int j,int k){return LocalRotations(i,j,k)[3];}});

    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "LocalRotations_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields);
}

void ElasticProperties::PrintPointStatistics(const int x, const int y, const int z) const
{
    ConsoleOutput::WriteCoordinate(x, y, z, Grid.dx);
    stringstream outstream;
    outstream << "Stresses    : " <<     Stresses(x, y, z).print() << endl;
    outstream << "TotalStrains: " << TotalStrains(x, y, z).print() << endl;
    outstream << "EigenStrains: " << EigenStrains(x, y, z).print() << endl;
    if(AnyPlasticity)
    {
        outstream << "PlasticStrains: " << PlasticStrains(x, y, z).print() << endl;
    }
    outstream << "ElasticConstants: " << endl
              << EffectiveElasticConstants(x, y, z).print() << endl;

    ConsoleOutput::WriteSimple(outstream.str());
}

void ElasticProperties::WriteStressStrainData(int time_step, double time, string filename) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    string separator = " ";

    if(!std::filesystem::exists(filename) or time_step <= 0)
    {
        ofstream file(filename, ios::out);
        file << "TimeStep"   << separator
             << "Time"       << separator
             << "Epsilon_xx" << separator
             << "Epsilon_yy" << separator
             << "Epsilon_zz" << separator
             << "Epsilon_yz" << separator
             << "Epsilon_xz" << separator
             << "Epsilon_xy" << separator
             << "Sigma_xx"   << separator
             << "Sigma_yy"   << separator
             << "Sigma_zz"   << separator
             << "Sigma_yz"   << separator
             << "Sigma_xz"   << separator
             << "Sigma_xy"   << separator
             << "Mises"      << separator
             << "Pressure"   << separator
             << "Norm"       << endl;
        file.close();
    }

    ofstream outputFile(filename, ios::app);
    outputFile << time_step        << separator
               << time             << separator

               << AverageStrain[0] << separator
               << AverageStrain[1] << separator
               << AverageStrain[2] << separator
               << AverageStrain[3] << separator
               << AverageStrain[4] << separator
               << AverageStrain[5] << separator

               << AverageStress[0] << separator
               << AverageStress[1] << separator
               << AverageStress[2] << separator
               << AverageStress[3] << separator
               << AverageStress[4] << separator
               << AverageStress[5] << separator

               << AverageStress.Mises()    << separator
               << AverageStress.Pressure() << separator
               << AverageStress.norm()     << endl;

    outputFile.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void ElasticProperties::WriteDilatometerCurve(std::string filename, double temperature, double real_size) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    ofstream outputFile;
    outputFile.open(filename, ios::app|ios::out);
    outputFile << temperature  << ", " << (pow(AverageDeformationGradient().determinant(), 1.0/(Grid.Active()))-1.0)*real_size << endl;
    outputFile.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void ElasticProperties::SetVelocityGradients(double dt)
{
    if(VelocityGradientsTotal.IsAllocated())
    {
        double velnorm = 1.0/dt;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityGradientsTotal,0,)
        {
            VelocityGradientsTotal(i, j, k) = (DeformationGradientsTotal(i, j, k) - dMatrix3x3::UnitTensor()) * velnorm;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        ConsoleOutput::WriteExit("VelocityGradient is not allocated!", "ElasticProperties",
                                                       "SetVelocityGradient()");
        exit(13);
    }
}

void ElasticProperties::SetLocalRotations(double dt)
{
    if(VelocityGradientsTotal.IsAllocated())
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityGradientsTotal,0,)
        {
            // Update rotation via Euler-Rodrigues formula
            dMatrix3x3 W = (VelocityGradientsTotal(i,j,k) - VelocityGradientsTotal(i,j,k).transposed())*0.5;

            dMatrix3x3 RotMat;
            RotMat.set_to_unity();

            double wnorm = 0.0;
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            {
                wnorm += W(ii,jj)*W(ii,jj);
            }
            wnorm = sqrt(wnorm*0.5);

            if (wnorm > 0.0)
            {
                RotMat += W*(sin(wnorm*dt)/wnorm)+W*W*((1.0-cos(wnorm*dt))/(wnorm*wnorm)); //Euler-Rodrigues formula
            }

            LocalRotations(i,j,k).set(RotMat*LocalRotations(i,j,k).RotationMatrix);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        ConsoleOutput::WriteExit("VelocityGradient or OR.Rotations is not allocated!",
                        "ElasticProperties", "SetRotations()");
        exit(13);
    }
}

void ElasticProperties::SetStressesRX(PhaseField& Phase, BoundaryConditions& BC, std::vector<int> targetPhases)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Stresses,0,)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsProperties[phaseIndex].Phase;
                for (auto jt = targetPhases.begin(); jt < targetPhases.end(); jt++)
                {
                    if (*jt == phaseIndex)
                    {
                        vStress temp = Stresses(i,j,k) * alpha->value;
                        Stresses(i,j,k) *= (1 - alpha->value);
                        vStrain tempStrain = EffectiveElasticConstants(i,j,k).inverted()*temp;
                        temp = PhaseElasticConstants[thPhaseIndex] * tempStrain;
//                        double hydroStress = 0.0;
                        for(int kt = 0; kt < 3; kt++)
                        {
//                            hydroStress += temp[kt];
                            temp[kt+3] = 0.0;
                        }
//                        hydroStress /= 3;
//                        for(int kt = 0; kt < 3; kt++)
//                        {
//                            temp[kt] = hydroStress;
//                        }
//                        vStress temp = EffectiveElasticConstants(i,j,k) * (EigenStrains[thPhaseIndex]) * alpha->value;
                        Stresses(i,j,k) += temp * 0.1;
                    }
                }
            } // end alpha
        }
        else
        {
            int locIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsProperties[locIndex].Phase;
            for (auto jt = targetPhases.begin(); jt < targetPhases.end(); jt++)
            {
                if (*jt == locIndex)
                {
                    vStress temp = Stresses(i,j,k);
//                    Stresses(i,j,k) *= (1 - alpha->value);
                    vStrain tempStrain = EffectiveElasticConstants(i,j,k).inverted()*temp;
                    temp = PhaseElasticConstants[thPhaseIndex] * tempStrain;
//                    double hydroStress = 0.0;
                    for(int kt = 0; kt < 3; kt++)
                    {
//                        hydroStress += temp[kt];
                        temp[kt+3] = 0.0;
                    }
//                    hydroStress /= 3;
//                    for(int kt = 0; kt < 3; kt++)
//                    {
//                        temp[kt] *= 0.5;
////                        temp[kt] = hydroStress;
//                    }
//                    vStress temp = EffectiveElasticConstants(i,j,k) * (EigenStrains(i,j,k));
                    Stresses(i,j,k) = temp * 0.1;
                }
            }
        } // end interface
    } // end loop space
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}

double ElasticProperties::Energy(const PhaseField& Phase) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0, reduction(+:Energy))
    {
        Energy += EnergyDensity(Phase,i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE,&Energy,1,OP_MPI_DOUBLE,OP_MPI_SUM,OP_MPI_COMM_WORLD);
#endif

    return Energy*Grid.CellVolume();
}

double ElasticProperties::AverageEnergyDensity(const PhaseField& Phase) const
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0, reduction(+:Energy))
    {
        Energy += EnergyDensity(Phase,i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

#ifdef MPI_PARALLEL
    OP_MPI_Allreduce(OP_MPI_IN_PLACE,&Energy,1,OP_MPI_DOUBLE,OP_MPI_SUM,OP_MPI_COMM_WORLD);
#endif

    return Energy/double(Grid.TotalNumberOfCells());
}

void ElasticProperties::CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceProperties& IP) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DeformationGradientsTotal, 0,)
    {
        if(Phase.Fields(i,j,k).wide_interface())
        {
            for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != Phase.Fields(i, j, k).cend(); ++beta)
            {
                IP.Properties(i,j,k).add_energy(alpha->index, beta->index, EnergyDensity(Phase,i,j,k)*Phase.Grid.Eta);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double ElasticProperties::EnergyDensity(const PhaseField& Phase, int i, int j, int k) const
{
    double energy = 0.0;

    switch(ElasticityModel)
    {
        case ElasticityModels::Khachaturyan:
        case ElasticityModels::Steinbach:
        case ElasticityModels::Voigt:
        case ElasticityModels::Reuss:
        {
            vStrain locElasticStrains = ElasticStrains(i,j,k);
            energy = 0.5*(locElasticStrains*(EffectiveElasticConstants(i,j,k)*locElasticStrains));

            break;
        }
        case ElasticityModels::Rank1:
        case ElasticityModels::Rank1NL:
        {
            if(Phase.Fields(i,j,k).interface())
            {
                vStrain locStrain = TotalStrains(i,j,k);

                for(auto alpha  = Phase.Fields(i, j, k).cbegin();
                         alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
                {
                    vStrain StrainAlpha = locStrain;

                    for(auto  beta  = alpha + 1;
                              beta != Phase.Fields(i, j, k).cend();  ++beta)
                    if(alpha->value > DBL_EPSILON and beta->value > DBL_EPSILON)
                    {
                        double scale = 1.0;
                        if(Phase.FieldsProperties[alpha->index].Stage != GrainStages::Stable)
                        {
                            scale = Phase.FieldsProperties[alpha->index].VolumeRatio;
                        }
                        if(Phase.FieldsProperties[beta->index].Stage != GrainStages::Stable)
                        {
                            scale = Phase.FieldsProperties[beta->index].VolumeRatio;
                        }

                        dMatrix3x3 locFjumpAlpha = DeformationJumps(i,j,k).get_asym2(alpha->index, beta->index);

                        dMatrix3x3 locStrainJumpA = (locFjumpAlpha + locFjumpAlpha.transposed())*0.5;

                        int sign = 1.0;
                        if(alpha->index > beta->index)
                        {
                            sign = -1.0;
                        }
                        StrainAlpha -= VoigtStrain(locStrainJumpA)*scale*sign*beta->value/(alpha->value + beta->value);
                    }
                    vStrain locEigenStrainAlpha = StrainSmall(GrainTransformationStretches[alpha->index]);

                    energy += ((StrainAlpha - locEigenStrainAlpha)*
                               (GrainElasticConstants[alpha->index]*
                               (StrainAlpha - locEigenStrainAlpha)))*0.5*alpha->value;
                }
            }
            else
            {
                vStrain locElasticStrains = ElasticStrains(i,j,k);
                energy = 0.5*(locElasticStrains*(EffectiveElasticConstants(i,j,k)*locElasticStrains));
            }
            break;
        }
    }
    return energy;
}

void ElasticProperties::WriteEnergyDensityVTK(Settings& locSettings, const PhaseField& Phase, const int tStep, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Elastic Energy Density", [this, &Phase] (int i, int j, int k){return EnergyDensity(Phase,i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "ElasticEnergyDensity_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

vStress ElasticProperties::ShearStresses(int i, int j, int k) const
{
    vStress tmpStress = Stresses(i,j,k);

    double locTrace = tmpStress.trace()/(Grid.Active());
    if(Grid.dNx) tmpStress[0] -= locTrace;
    if(Grid.dNy) tmpStress[1] -= locTrace;
    if(Grid.dNz) tmpStress[2] -= locTrace;

    return tmpStress;
}

dMatrix3x3 ElasticProperties::AverageDeformationGradient(double accuracy) const
{
    return AverageStrain.tensor() + dMatrix3x3::UnitTensor();
}

vStrain ElasticProperties::EigenStrainDifference(
        const dMatrix3x3 locStretchesAlpha,
        const dMatrix3x3 locStretchesBeta,
        const int i, const int j, const int k) const
{
    dMatrix3x3 dFeigenAlphaBeta = locStretchesBeta - locStretchesAlpha;
    return VoigtStrain(dFeigenAlphaBeta.transposed() + dFeigenAlphaBeta)*0.5;
}

vStrain ElasticProperties::TotalStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsTotal(i,j,k));
}

vStrain ElasticProperties::EigenStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsEigen(i,j,k));
}

vStrain ElasticProperties::PlasticStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsPlastic(i,j,k));
}

vStrain ElasticProperties::StressFreeStrains(int i, int j, int k) const
{
    vStrain locStrain = StrainSmall(DeformationGradientsEigen(i,j,k));
    if(AnyPlasticity)
    {
        locStrain += StrainSmall(DeformationGradientsPlastic(i,j,k));
    }
    return locStrain;
}

vStrain ElasticProperties::ElasticStrains(int i, int j, int k) const
{
    return StrainSmall(DeformationGradientsElastic(i,j,k));
}

vStress ElasticProperties::Stress2PK(int i, int j, int k, const dMatrix3x3& locDefGrad) const
{
    vStress locStress = EffectiveElasticConstants(i, j, k)*
                        (StrainSmall(locDefGrad) - StressFreeStrains(i, j, k));
    return locStress;
}

dMatrix3x3 ElasticProperties::Stress2PKpullBack(int i, int j, int k) const
{
    dMatrix3x3 Elastic2PKstress = (EffectiveElasticConstants(i,j,k)*ElasticStrains(i,j,k)).tensor();

    return Elastic2PKstress;
}

dMatrix3x3 ElasticProperties::Stress1PK(int i, int j, int k) const
{
    return Stress2PKpullBack(i,j,k);
}

dMatrix3x3 ElasticProperties::DeformationGradientsElastic(int i, int j, int k) const
{
    dMatrix3x3 locDeformationGradientsElastic = DeformationGradientsTotal(i,j,k);
    locDeformationGradientsElastic -= DeformationGradientsEigen(i,j,k) - dMatrix3x3::UnitTensor();
    if(AnyPlasticity)
    {
        locDeformationGradientsElastic -= DeformationGradientsPlastic(i,j,k) - dMatrix3x3::UnitTensor();
    }
    return locDeformationGradientsElastic;
}

}// namespace openphase
