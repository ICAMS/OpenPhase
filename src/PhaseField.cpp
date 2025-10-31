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

 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev;
 *                         Reza Darvishi Kamachali; Raphael Schiedung;
 *                         Johannes Goerler; Marvin Tegeler
 *
 */

#include "BoundaryConditions.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "ElasticProperties.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"
#include "Velocities.h"
#include "AdvectionHR.h"
#include "H5Interface.h"

namespace openphase
{
using namespace std;
using json = nlohmann::json;

void PhaseField::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "PhaseField";
    thisobjectname = thisclassname + ObjectNameSuffix;

    NucleationPresent = false;
    InterfaceNormalModel = InterfaceNormalModels::AverageGradient;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;
    PhaseNames = locSettings.PhaseNames;
    PhaseAggregateStates = locSettings.PhaseAggregateStates;
    PhaseEquilibriumDensities = locSettings.PhaseEquilibriumDensities;

    NucleusVolumeFactor = 1.0;

    ConsiderNucleusVolume = true;
    Combine.resize(Nphases, false);

    PhaseFieldLaplacianStencil = LaplacianStencils::Isotropic;
    PhaseFieldGradientStencil = GradientStencils::Isotropic;

    locSettings.AddForAdvection(*this);
    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void PhaseField::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("PhaseField input");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());
    std::string filetype = FileInterface::getFileExtension(InputFileName);
    if (filetype == "opi")
    {
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
    else
    if (filetype == "json")
    {
        ReadJSON(InputFileName);
    }
    else
    {
        std::cerr << "Filetype " << filetype << " not recognized. Filetype must be opi or json." << std::endl;
        OP_Exit(1);
    }
}

void PhaseField::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    ConsiderNucleusVolume = FileInterface::ReadParameterB(inp, moduleLocation, string("ConsiderNucleusVolume"), false, true);
    NucleusVolumeFactor   = FileInterface::ReadParameterD(inp, moduleLocation, string("NucleusVolumeFactor"), false, 1.0);

    // Reading combine phase fields conditions for all phases
    for(size_t pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "CombinePhaseFields_" << pIndex;
        Combine[pIndex] = FileInterface::ReadParameterB(inp, moduleLocation, converter.str(),false,false);
    }

    string tmp1 = FileInterface::ReadParameterK(inp, moduleLocation, string("InterfaceNormalModel"), false, string("AVERAGEGRADIENT"));
    if(tmp1 == "AVERAGEGRADIENT")
    {
        InterfaceNormalModel = InterfaceNormalModels::AverageGradient;
    }
    else if(tmp1 == "WEIGHTEDGRADIENT")
    {
        InterfaceNormalModel = InterfaceNormalModels::WeightedGradient;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong interface normal model specified!\nThe default \"AVERAGEGRADIENT\" model is used!", thisclassname, "ReadInput()");
    }

    string tmp2 = FileInterface::ReadParameterK(inp, moduleLocation, string("PhaseFieldLaplacianStencil"), false, string("ISOTROPIC"));
    if(tmp2 == "SIMPLE")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::Simple;
    }
    else if(tmp2 == "ISOTROPIC")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::Isotropic;
    }
    else if(tmp2 == "LB")
    {
        PhaseFieldLaplacianStencil = LaplacianStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong Laplaciant stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    string tmp3 = FileInterface::ReadParameterK(inp, moduleLocation, string("PhaseFieldGradientStencil"), false, string("ISOTROPIC"));
    if(tmp3 == "SIMPLE")
    {
        PhaseFieldGradientStencil = GradientStencils::Simple;
    }
    else if(tmp3 == "ISOTROPIC")
    {
        PhaseFieldGradientStencil = GradientStencils::Isotropic;
    }
    else if(tmp3 == "LB")
    {
        PhaseFieldGradientStencil = GradientStencils::LB;
    }
    else
    {
        ConsoleOutput::WriteWarning("No or wrong gradient stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
    }

    SetStencils(Grid);
    AllocateStorages(Grid);

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void PhaseField::ReadJSON(const string InputFileName)
{
    std::ifstream f(InputFileName);
    json data = json::parse(f);
    if (data.contains(thisclassname))
    {
        json phasefield = data[thisclassname];
        ConsiderNucleusVolume = FileInterface::ReadParameter<bool>(phasefield, {"ConsiderNucleusVolume"}, false);
        NucleusVolumeFactor   = FileInterface::ReadParameter<double>(phasefield, {"NucleusVolumeFactor"}, 1.0);

        string tmp1 = FileInterface::ReadParameter<std::string>(phasefield, {"InterfaceNormalModel"}, "AVERAGEGRADIENT");
        if(tmp1 == "AVERAGEGRADIENT")
        {
            InterfaceNormalModel = InterfaceNormalModels::AverageGradient;
        }
        else if(tmp1 == "WEIGHTEDGRADIENT")
        {
            InterfaceNormalModel = InterfaceNormalModels::WeightedGradient;
        }
        else
        {
            ConsoleOutput::WriteWarning("No or wrong interface normal model specified!\nThe default \"AVERAGEGRADIENT\" model is used!", thisclassname, "ReadJSON()");
        }

        string tmp2 = FileInterface::ReadParameter<std::string>(phasefield, {"PhaseFieldLaplacianStencil"}, "ISOTROPIC");
        if(tmp2 == "SIMPLE")
        {
            PhaseFieldLaplacianStencil = LaplacianStencils::Simple;
        }
        else if(tmp2 == "ISOTROPIC")
        {
            PhaseFieldLaplacianStencil = LaplacianStencils::Isotropic;
        }
        else if(tmp2 == "LB")
        {
            PhaseFieldLaplacianStencil = LaplacianStencils::LB;
        }
        else
        {
            ConsoleOutput::WriteWarning("No or wrong Laplaciant stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadJSON()");
        }

        string tmp3 = FileInterface::ReadParameter<std::string>(phasefield, {"PhaseFieldGradientStencil"}, "ISOTROPIC");
        if(tmp3 == "SIMPLE")
        {
            PhaseFieldGradientStencil = GradientStencils::Simple;
        }
        else if(tmp3 == "ISOTROPIC")
        {
            PhaseFieldGradientStencil = GradientStencils::Isotropic;
        }
        else if(tmp3 == "LB")
        {
            PhaseFieldGradientStencil = GradientStencils::LB;
        }
        else
        {
            ConsoleOutput::WriteWarning("No or wrong gradient stencil specified!\nThe default \"ISOTROPIC\" model is used!", thisclassname, "ReadInput()");
        }
    }
    SetStencils(Grid);
    AllocateStorages(Grid);
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void PhaseField::SetStencils(GridParameters& Grid)
{
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
            switch(PhaseFieldLaplacianStencil)
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

            switch(PhaseFieldGradientStencil)
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
            switch(PhaseFieldLaplacianStencil)
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

            switch(PhaseFieldGradientStencil)
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
}

void PhaseField::AllocateStorages(GridParameters& Grid)
{
    // extra cells needed for driving force averaging
    size_t Bcells = max(Grid.Bcells, int(Grid.iWidth)-1);

    if(Grid.Resolution == Resolutions::Dual)
    {
        Bcells = max(Grid.Bcells, int(Grid.iWidth/2));
    }

    Fields   .Allocate(Grid, Bcells);
    FieldsDot.Allocate(Grid, Bcells);

    Fractions.Allocate(Grid,{Nphases}, Bcells);

    FractionsTotal.resize(Nphases, 0.0);

    if(Grid.Resolution == Resolutions::Dual)
    {
        GridParameters DoubleDimensions = Grid.DoubleResolution();

        FieldsDR   .Allocate(DoubleDimensions, Bcells*2);
        FieldsDotDR.Allocate(DoubleDimensions, Bcells*2);
    }
}

double PhaseField::CalculateReferenceVolume(double radius)
{
    double ref_volume = 0.0;

    switch(Grid.Active())
    {
        case 1:
        {
            ref_volume = radius;
            break;
        }
        case 2:
        {
            ref_volume = Pi*radius*radius;

            break;
        }
        case 3:
        {
            ref_volume = (4.0/3.0)*Pi*radius*radius*radius;
            break;
        }
    }
    return ref_volume*NucleusVolumeFactor;
}

void PhaseField::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        Fields(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if(Grid.Resolution == Resolutions::Dual)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
        {
            FieldsDR(i,j,k).clear();
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void PhaseField::CalculateFractions(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        Fractions(i,j,k).set_to_zero();
        for (auto it  = Fields(i,j,k).cbegin();
                  it != Fields(i,j,k).cend(); ++it)
        {
            Fractions(i,j,k,{FieldsProperties[it->index].Phase}) += it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

//Tensor<double,1> PhaseField::Fractions(const int x, const int y, const int z) const
//{
//    Tensor<double,1> locFractions({Nphases});
//    for (auto it  = Fields(x,y,z).cbegin();
//              it != Fields(x,y,z).cend(); ++it)
//    {
//        locFractions({FieldsProperties[it->index].Phase}) += it->value;
//    }
//    return locFractions;
//}

//double PhaseField::Fraction(const int x, const int y, const int z, const size_t idx) const
//{
//    double locFraction = 0.0;
//    for (auto it  = Fields(x,y,z).cbegin();
//              it != Fields(x,y,z).cend(); ++it)
//    {
//        if(FieldsProperties[it->index].Phase == idx)
//        {
//            locFraction += it->value;
//        }
//    }
//    return locFraction;
//}
/*
NodePF PhaseField::FieldsDerivatives(const int i, const int j, const int k) const
{
    NodePF locPF = Fields(i,j,k);

    for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
    {
        int ii = ls->di;
        int jj = ls->dj;
        int kk = ls->dk;

        for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                  it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double laplacian = ls->weight * it->value;
            locPF.add_laplacian(it->index, laplacian);
        }
    }
    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        int ii = gs->di;
        int jj = gs->dj;
        int kk = gs->dk;

        for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                  it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double value_x = gs->weightX * it->value;
            double value_y = gs->weightY * it->value;
            double value_z = gs->weightZ * it->value;
            locPF.add_gradient(it->index, (dVector3){value_x,value_y,value_z});
        }
    }
    return locPF;
}*/
/*
NodePF PhaseField::FieldsDerivativesDR(const int i, const int j, const int k) const
{
    NodePF locPF = FieldsDR(i,j,k);

    for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
    {
        int ii = ls->di;
        int jj = ls->dj;
        int kk = ls->dk;

        for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                  it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double laplacian = 4.0*ls->weight * it->value;
            locPF.add_laplacian(it->index, laplacian);
        }
    }
    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        int ii = gs->di;
        int jj = gs->dj;
        int kk = gs->dk;

        for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                  it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double value_x = 2.0*gs->weightX * it->value;
            double value_y = 2.0*gs->weightY * it->value;
            double value_z = 2.0*gs->weightZ * it->value;
            locPF.add_gradient(it->index, (dVector3){value_x,value_y,value_z});
        }
    }
    return locPF;
}*/
/*
NodePF PhaseField::FieldsLaplacians(const int i, const int j, const int k) const
{
    NodePF locPF = Fields(i,j,k);

    for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
    {
        int ii = ls->di;
        int jj = ls->dj;
        int kk = ls->dk;

        for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                  it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double laplacian = ls->weight * it->value;
            locPF.add_laplacian(it->index, laplacian);
        }
    }
    return locPF;
}*/
/*
NodePF PhaseField::FieldsLaplaciansDR(const int i, const int j, const int k) const
{
    NodePF locPF = FieldsDR(i,j,k);

    for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
    {
        int ii = ls->di;
        int jj = ls->dj;
        int kk = ls->dk;

        for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                  it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double laplacian = 4.0*ls->weight * it->value;
            locPF.add_laplacian(it->index, laplacian);
        }
    }
    return locPF;
}*/

void PhaseField::CalculateDerivativesSR(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells()-1,)
    {
        if(Fields(i,j,k).wide_interface())
        {
            Fields(i,j,k).set_temporary();

            for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
            {
                int ii = ls->di;
                int jj = ls->dj;
                int kk = ls->dk;

                for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                          it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double laplacian = ls->weight * it->value;
                    Fields(i,j,k).add_laplacian_tmp(it->index, laplacian);
                }
            }
            for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
            {
                int ii = gs->di;
                int jj = gs->dj;
                int kk = gs->dk;

                for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                          it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double value_x = gs->weightX * it->value;
                    double value_y = gs->weightY * it->value;
                    double value_z = gs->weightZ * it->value;
                    Fields(i,j,k).add_gradient_tmp(it->index, (dVector3){value_x,value_y,value_z});
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells()-1,)
    {
        if(Fields(i,j,k).wide_interface())
        {
            Fields(i,j,k).copy_from_temporary();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CalculateDerivativesDR(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, FieldsDR.Bcells()-1,)
    {
        if(FieldsDR(i,j,k).wide_interface())
        {
            FieldsDR(i,j,k).set_temporary();

            for (auto ls = LStencil.cbegin(); ls != LStencil.cend(); ls++)
            {
                int ii = ls->di;
                int jj = ls->dj;
                int kk = ls->dk;

                for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double laplacian = 4.0*ls->weight * it->value;
                    FieldsDR(i,j,k).add_laplacian_tmp(it->index, laplacian);
                }
            }
            for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
            {
                int ii = gs->di;
                int jj = gs->dj;
                int kk = gs->dk;

                for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                          it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
                if (it->value != 0.0)
                {
                    double value_x = 2.0*gs->weightX * it->value;
                    double value_y = 2.0*gs->weightY * it->value;
                    double value_z = 2.0*gs->weightZ * it->value;
                    FieldsDR(i,j,k).add_gradient_tmp(it->index, (dVector3){value_x,value_y,value_z});
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, FieldsDR.Bcells()-1,)
    {
        if(FieldsDR(i,j,k).wide_interface())
        {
            FieldsDR(i,j,k).copy_from_temporary();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
/*
NodePF PhaseField::FieldsGradients(const int i, const int j, const int k) const
{
    NodePF locPF = Fields(i,j,k);

    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        int ii = gs->di;
        int jj = gs->dj;
        int kk = gs->dk;

        for (auto it  = Fields(i + ii, j + jj, k + kk).cbegin();
                  it != Fields(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double value_x = gs->weightX * it->value;
            double value_y = gs->weightY * it->value;
            double value_z = gs->weightZ * it->value;
            locPF.add_gradient(it->index, (dVector3){value_x,value_y,value_z});
        }
    }
    return locPF;
}*/
/*
NodePF PhaseField::FieldsGradientsDR(const int i, const int j, const int k) const
{
    NodePF locPF = FieldsDR(i,j,k);

    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        int ii = gs->di;
        int jj = gs->dj;
        int kk = gs->dk;

        for (auto it  = FieldsDR(i + ii, j + jj, k + kk).cbegin();
                  it != FieldsDR(i + ii, j + jj, k + kk).cend(); ++it)
        if (it->value != 0.0)
        {
            double value_x = 2.0*gs->weightX * it->value;
            double value_y = 2.0*gs->weightY * it->value;
            double value_z = 2.0*gs->weightZ * it->value;
            locPF.add_gradient(it->index, (dVector3){value_x,value_y,value_z});
        }
    }
    return locPF;
}*/

dVector3 PhaseField::Normal(NodePF::citerator alpha, NodePF::citerator beta) const
{
    dVector3 locNormal;

    switch(InterfaceNormalModel)
    {
        case InterfaceNormalModels::AverageGradient:
        {
            locNormal = beta->gradient - alpha->gradient;
            break;
        }
        case InterfaceNormalModels::WeightedGradient:
        {
            locNormal = beta->gradient*alpha->value - alpha->gradient*beta->value;
            break;
        }
    }
    locNormal.normalize();

    return locNormal;
}

NodeAB<dVector3,dVector3> PhaseField::Normals(const int i, const int j, const int k) const
{
    NodeAB<dVector3,dVector3> locNormals;

    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    for (auto  beta  = alpha + 1;
               beta != Fields(i,j,k).cend(); ++beta)
    {
        dVector3 value = Normal(alpha,beta);
        locNormals.set_asym1(alpha->index, beta->index, value);
    }
    return locNormals;
}

NodeAB<dVector3,dVector3> PhaseField::NormalsDR(const int i, const int j, const int k) const
{
    NodeAB<dVector3,dVector3> locNormals;

    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    for (auto  beta  = alpha + 1;
               beta != FieldsDR(i,j,k).cend(); ++beta)
    {
        dVector3 value = Normal(alpha,beta);
        locNormals.set_asym1(alpha->index, beta->index, value);
    }
    return locNormals;
}

NodeA<dVector3> PhaseField::NormalsPhase(const int i, const int j, const int k) const
{
    NodeA<dVector3> locNormals = Fields(i,j,k).get_gradients();

    for (auto alpha  = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        alpha->value.normalize();
    }
    return locNormals;
}

NodeA<dVector3> PhaseField::NormalsPhaseDR(const int i, const int j, const int k) const
{
    NodeA<dVector3> locNormals = FieldsDR(i,j,k).get_gradients();

    for (auto alpha  = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        alpha->value.normalize();
    }
    return locNormals;
}

dMatrix3x3 PhaseField::EffectiveOrientation(const int i, const int j, const int k) const
{
    Quaternion locQuaternion;

    if(Fields(i,j,k).interface())
    {
        int index = Fields(i,j,k).majority_index();
        locQuaternion = FieldsProperties[index].Orientation;
    }
    else
    {
        locQuaternion.set(0.0,0.0,0.0,0.0);

        for(auto alpha  = Fields(i,j,k).cbegin();
                 alpha != Fields(i,j,k).cend(); alpha++)
        {
            locQuaternion += FieldsProperties[alpha->index].Orientation*alpha->value;
        }
    }
    return locQuaternion.RotationMatrix;
}

void PhaseField::CalculateGrainsVolume(void)
{
    int Nthreads = 1;

    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif
    const size_t size = FieldsProperties.size();
    vector<vector<double>> Volume(Nthreads);
    for(int t = 0; t < Nthreads; t++)
    {
        Volume[t].resize(size, 0.0);
    }
    // Calculate volumes in each OpenMP chunk
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        int thread = 0;

        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it  = Fields(i,j,k).cbegin();
                 it != Fields(i,j,k).cend(); ++it)
        {
            Volume[thread][it->index] += it->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    // Add volumes from different OpenMP chunks
    for(size_t idx = 0; idx < size; idx++)
    {
        FieldsProperties[idx].Volume = 0.0;

        for(int t = 0; t < Nthreads; t++)
        {
            FieldsProperties[idx].Volume += Volume[t][idx];
        }
    }

    // Update FieldsProperties across MPI domains
#ifdef MPI_PARALLEL
    unsigned long loc_size = FieldsProperties.size();
    unsigned long max_size = loc_size;

    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &max_size, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    if (max_size > loc_size)
    {
        FieldsProperties.Resize(max_size);
    }

    for(size_t idx = 0; idx < FieldsProperties.size(); idx++)
    {
        double loc_volume = FieldsProperties[idx].Volume;
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_volume, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].Volume = loc_volume;

        double loc_maxvolume = FieldsProperties[idx].MAXVolume;
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_maxvolume, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].MAXVolume = loc_maxvolume;

        double loc_refvolume = FieldsProperties[idx].RefVolume;
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_refvolume, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].RefVolume = loc_refvolume;

        int loc_stage = static_cast<int>(FieldsProperties[idx].Stage);
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_stage, 1, OP_MPI_INT, OP_MPI_MAX, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].Stage = static_cast<openphase::GrainStages>(loc_stage);

        unsigned long loc_variant = FieldsProperties[idx].Variant;
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_variant, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].Variant = loc_variant;

        unsigned long loc_phase = FieldsProperties[idx].Phase;
        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &loc_phase, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
        FieldsProperties[idx].Phase = loc_phase;

        //TODO: add other missing reductions
    }
#endif

    // Calculate MAXVolume and VolumeRatio for all phase fields
    for(size_t idx = 0; idx < FieldsProperties.size(); idx++)
    {
        FieldsProperties[idx].MAXVolume = max(FieldsProperties[idx].Volume,
                                              FieldsProperties[idx].MAXVolume);

        FieldsProperties[idx].VolumeRatio = FieldsProperties[idx].MAXVolume/
                                            FieldsProperties[idx].RefVolume;
    }

    // Update statistics based on actual grains volume
    int NumberOfNuclei = 0;

    for(size_t idx = 0; idx < size; idx++)
    if(FieldsProperties[idx].Exist and FieldsProperties[idx].Volume <= 0.0)
    {
        FieldsProperties[idx].Exist  = false;
        FieldsProperties[idx].Stage  = GrainStages::Stable;
        FieldsProperties[idx].Volume = 0.0;
        FieldsProperties[idx].MAXVolume = 0.0;
        FieldsProperties[idx].VolumeRatio = 1.0;
    }
    else if (FieldsProperties[idx].Volume > 0.0)
    {
        FieldsProperties[idx].Exist = true;
        if(FieldsProperties[idx].Stage == GrainStages::Seed)
        {
            FieldsProperties[idx].Stage = GrainStages::Nucleus;
        }

        if(FieldsProperties[idx].VolumeRatio > 1.0)
        {
            FieldsProperties[idx].Stage = GrainStages::Stable;
            FieldsProperties[idx].VolumeRatio = 1.0;
        }
        if(FieldsProperties[idx].Stage != GrainStages::Stable)
        {
            NumberOfNuclei ++;
        }
    }
    NucleationPresent = (NumberOfNuclei > 0);

    for(size_t n = 0; n < Nphases; n++)
    {
        FractionsTotal[n] = 0.0;
    }

    for(size_t idx = 0; idx < FieldsProperties.size(); idx++)
    if(FieldsProperties[idx].Exist)
    {
        FractionsTotal[FieldsProperties[idx].Phase] += FieldsProperties[idx].Volume;
    }

    for(size_t n = 0; n < Nphases; n++)
    {
        FractionsTotal[n] /= double(Grid.TotalNumberOfCells());
    }
}

void PhaseField::SetFlagsSR(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells()-1,)
    {
        if(Fields(i,j,k).flag == 2)
        {
            for(int ii = -Grid.dNx; ii <= +Grid.dNx; ++ii)
            for(int jj = -Grid.dNy; jj <= +Grid.dNy; ++jj)
            for(int kk = -Grid.dNz; kk <= +Grid.dNz; ++kk)
            if(!(Fields(i+ii, j+jj, k+kk).flag))
            {
                Fields(i+ii, j+jj, k+kk).flag = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SetFlagsDR(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, FieldsDR.Bcells()-1,)
    {
        if(FieldsDR(i,j,k).flag == 2)
        {
            for(int ii = -Grid.dNx; ii <= +Grid.dNx; ++ii)
            for(int jj = -Grid.dNy; jj <= +Grid.dNy; ++jj)
            for(int kk = -Grid.dNz; kk <= +Grid.dNz; ++kk)
            if(!(FieldsDR(i+ii, j+jj, k+kk).flag))
            {
                FieldsDR(i+ii, j+jj, k+kk).flag = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::Finalize(const BoundaryConditions& BC, bool finalize)
{
    CombinePhaseFields();

    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            FinalizeSR(BC, finalize);
            break;
        }
        case Resolutions::Dual:
        {
            FinalizeDR(BC, finalize);
            break;
        }
    }
}

void PhaseField::FinalizeSR(const BoundaryConditions& BC, bool finalize)
{
    if(finalize)
    {
        #ifdef MPI_PARALLEL
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        #else
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        #endif
        {
            if(Fields(i,j,k).wide_interface())
            {
                Fields(i,j,k).finalize();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    SetBoundaryConditionsSR(BC);
    SetFlagsSR();
    SetBoundaryConditionsSR(BC);
    CalculateDerivativesSR();
    SetBoundaryConditionsSR(BC);
    CalculateFractions();
    CalculateGrainsVolume();
}

void PhaseField::FinalizeDR(const BoundaryConditions& BC, bool finalize)
{
    if(finalize)
    {
        #ifdef MPI_PARALLEL
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
        #else
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
        #endif
        {
            if(FieldsDR(i,j,k).wide_interface())
            {
                FieldsDR(i,j,k).finalize();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    SetBoundaryConditionsDR(BC);
    SetFlagsDR();
    CalculateDerivativesDR();
    SetBoundaryConditionsDR(BC);

    Coarsen();

    SetBoundaryConditionsSR(BC);
    SetFlagsSR();
    SetBoundaryConditionsSR(BC);
    CalculateDerivativesSR();
    CalculateFractions();
    CalculateGrainsVolume();
}

void PhaseField::FinalizeInitialization(const BoundaryConditions& BC)
{
    FinalizeSR(BC);

    if(Grid.Resolution == Resolutions::Dual)
    {
        Refine();
        FinalizeDR(BC);
    }
}

void PhaseField::FixSpreading(const BoundaryConditions& BC, double cutoff)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        for (auto alpha  = Fields(i,j,k).begin();
                  alpha != Fields(i,j,k).end(); alpha++)
        {
            if(alpha->value < cutoff and FieldsProperties[alpha->index].Stage == GrainStages::Stable)
            {
                alpha->value = 0.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

void PhaseField::Coarsen(void)
{
    double norm = 1.0/pow(2.0,Grid.Active());
    long int fx = 1 + Grid.dNx;
    long int fy = 1 + Grid.dNy;
    long int fz = 1 + Grid.dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if(Fields(i,j,k).wide_interface())
        {
            Fields(i,j,k).clear();

            for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
            for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
            for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
            {
                Fields(i,j,k) += FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2);
            }
            Fields(i,j,k) *= norm;
            Fields(i,j,k).finalize();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CoarsenDot(void)
{
    double norm = 1.0/pow(2.0,Grid.Active());
    long int fx = 1 + Grid.dNx;
    long int fy = 1 + Grid.dNy;
    long int fz = 1 + Grid.dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
    {
        if(Fields(i,j,k).wide_interface())
        {
            FieldsDot(i,j,k).clear();

            for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
            for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
            for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
            {
                FieldsDot(i,j,k).add_asym_pairs(FieldsDotDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2));
            }
            FieldsDot(i,j,k) *= norm;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::Refine(void)
{
    long int fx = 1 + Grid.dNx;
    long int fy = 1 + Grid.dNy;
    long int fz = 1 + Grid.dNz;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        //if(Fields(i,j,k).flag)
        for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
        for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
        for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
        {
            FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2) = Fields.at(i+di*0.25,j+dj*0.25,k+dk*0.25);
            FieldsDR(fx*i+(di+1)/2,fy*j+(dj+1)/2,fz*k+(dk+1)/2).finalize();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

NodeA<double> PhaseField::Dot(const int i, const int j, const int k, const double dt) const
{
    NodeA<double> value;
    if(Fields(i,j,k).wide_interface())
    {
        NodePF OldFields = Fields(i,j,k);
        NodePF NewFields = Fields(i,j,k);
        NodePF tmp;

        for(auto psi = FieldsDot(i,j,k).cbegin(); psi != FieldsDot(i,j,k).cend(); ++psi)
        {
            if(psi->value1 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value1*dt);
                NewFields.add_value(psi->indexB, -psi->value1*dt);
            }
            if(psi->value2 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value2*dt);
                NewFields.add_value(psi->indexB, -psi->value2*dt);
            }
        }
        NodePF NewFieldsFinalized = NewFields;
        NewFieldsFinalized.finalize();

        for (auto alpha = NewFields.cbegin(); alpha != NewFields.end(); ++alpha)
        {
            value.add_value(alpha->index, (NewFieldsFinalized.get_value(alpha->index)-OldFields.get_value(alpha->index))/dt);
        }
    }
    return value;
}

NodeA<double> PhaseField::Dot1(const int i, const int j, const int k, const double dt) const
{
    NodeA<double> value;
    if(Fields(i,j,k).wide_interface())
    {
        NodePF OldFields = Fields(i,j,k);
        NodePF NewFields = Fields(i,j,k);
        //NodePF tmp;

        for(auto psi = FieldsDot(i,j,k).cbegin(); psi != FieldsDot(i,j,k).cend(); ++psi)
        {
            if(psi->value1 != 0.0)
            {
                NewFields.add_value(psi->indexA,  psi->value1*dt);
                NewFields.add_value(psi->indexB, -psi->value1*dt);
                //tmp.add_value(psi->indexA, 1);
                //tmp.add_value(psi->indexB, 1);
            }
        }
        //NewFields.finalize();

        //for (auto alpha = tmp.cbegin(); alpha != tmp.end(); ++alpha)
        for(auto alpha = NewFields.cbegin(); alpha != NewFields.cend(); ++alpha)
        {
            value.add_value(alpha->index, (NewFields.get_value(alpha->index)-OldFields.get_value(alpha->index))/dt);
        }
    }
    return value;
}

NodeA<double> PhaseField::Dot2(const int i, const int j, const int k, const double dt) const
{
    NodeA<double> value;
    if(Fields(i,j,k).wide_interface())
    {
        NodeA<double> locDot1 = Dot1(i,j,k,dt);
        NodeA<double> locDot  = Dot (i,j,k,dt);

        NodeA<double> tmp;
        for (auto it : locDot1) if (it.value != 0.0) tmp.add_value(it.index, 1);
        for (auto it : locDot ) if (it.value != 0.0) tmp.add_value(it.index, 1);

        for (auto alpha = tmp.cbegin(); alpha != tmp.end(); ++alpha)
        {
            double locDot2Alpha = locDot.get_value(alpha->index)-locDot1.get_value(alpha->index);
            value.add_value(alpha->index,locDot2Alpha);
        }
    }
    return value;
}

void PhaseField::MergeIncrements(const BoundaryConditions& BC,
                                 const double dt,
                                 const bool finalize,
                                 const bool clear)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            MergeIncrementsSR(BC, dt, finalize, clear);
            break;
        }
        case Resolutions::Dual:
        {
            MergeIncrementsDR(BC, dt, finalize, clear);
            break;
        }
    }
}

void PhaseField::MergeIncrementsSR(const BoundaryConditions& BC,
                                   const double dt,
                                   const bool finalize,
                                   const bool clear)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if(Fields(i,j,k).wide_interface())
        {
            for(auto psi  = FieldsDot(i,j,k).cbegin();
                     psi != FieldsDot(i,j,k).cend(); ++psi)
            {
                double factor = 1.0;
                if(FieldsProperties[psi->indexA].GrowthConstraintsViolation != GrowthConstraintsViolations::None and
                   FieldsProperties[psi->indexB].GrowthConstraintsViolation != GrowthConstraintsViolations::None)
                {
                    factor *= PairwiseGrowthFactors.get_sym1(psi->indexA,psi->indexB);
                }
            
                double value = factor*(psi->value1 + psi->value2)*dt;
                if(fabs(value) >= DBL_EPSILON)
                {
                    Fields(i,j,k).add_value(psi->indexA,  value);
                    Fields(i,j,k).add_value(psi->indexB, -value);
                }
            }
            if(clear) FieldsDot(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC, finalize);
}

void PhaseField::MergeIncrementsDR(const BoundaryConditions& BC,
                                   const double dt,
                                   const bool finalize,
                                   const bool clear)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
    {
        if(FieldsDR(i,j,k).wide_interface())
        {
            for(auto psi  = FieldsDotDR(i,j,k).cbegin();
                     psi != FieldsDotDR(i,j,k).cend(); ++psi)
            {
                double factor = 1.0;
                if(FieldsProperties[psi->indexA].GrowthConstraintsViolation != GrowthConstraintsViolations::None and
                   FieldsProperties[psi->indexB].GrowthConstraintsViolation != GrowthConstraintsViolations::None)
                {
                    factor *= PairwiseGrowthFactors.get_sym1(psi->indexA,psi->indexB);
                }

                double value = factor*(psi->value1 + psi->value2)*dt;
                if(fabs(value) >= DBL_EPSILON)
                {
                    FieldsDR(i,j,k).add_value(psi->indexA,  value);
                    FieldsDR(i,j,k).add_value(psi->indexB, -value);
                }
            }
            if(clear) FieldsDotDR(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC, finalize);
}

void PhaseField::NormalizeIncrements(const BoundaryConditions& BC, const double dt)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            NormalizeIncrementsSR(BC, dt);
            break;
        }
        case Resolutions::Dual:
        {
            NormalizeIncrementsDR(BC, dt);
            break;
        }
    }
}

void PhaseField::NormalizeIncrementsSR(const BoundaryConditions& BC, const double dt)
{
    /** This function limits phase-field increments for all present phase-field
    pairs, so that the actual phase-field values are within their natural
    limits of 0.0 and 1.0.*/

    double precision = FLT_EPSILON;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if (Fields(i,j,k).wide_interface())
    {
        switch(FieldsDot(i,j,k).size())
        {
            case 0:
            {
                break;
            }
            case 1:
            {
                size_t indexA = FieldsDot(i,j,k).front().indexA;
                double valueA = (FieldsDot(i,j,k).front().value1
                              +  FieldsDot(i,j,k).front().value2)*dt;
                double old_value = Fields(i,j,k).get_value(indexA);

                if((old_value == 0.0 and valueA < 0.0) or
                   (old_value == 1.0 and valueA > 0.0))
                {
                    FieldsDot(i,j,k).clear();
                }
                else
                {
                    double new_value = old_value + valueA;

                    double norm = 1.0;
                    if(new_value < 0.0)
                    {
                        norm *= -old_value/valueA;
                    }
                    else if(new_value > 1.0)
                    {
                        norm *= (1.0 - old_value)/valueA;
                    }

                    if(norm > DBL_EPSILON)
                    {
                        FieldsDot(i,j,k) *= norm;
                    }
                    else
                    {
                        FieldsDot(i,j,k).clear();
                    }
                }
                break;
            }
            default:
            {
                for(auto alpha  = Fields(i,j,k).cbegin();
                         alpha != Fields(i,j,k).cend(); ++alpha)
                {
                    double dPsiAlpha = 0.0;
                    for(auto beta  = Fields(i,j,k).cbegin();
                             beta != Fields(i,j,k).cend(); ++beta)
                    if(alpha != beta)
                    {
                        dPsiAlpha += FieldsDot(i,j,k).get_asym1(alpha->index,
                                                                 beta->index);
                        dPsiAlpha += FieldsDot(i,j,k).get_asym2(alpha->index,
                                                                 beta->index);
                    }
                    /* Set out of bounds sets to zero */
                    if((alpha->value == 0.0 and dPsiAlpha < 0.0)
                    or (alpha->value == 1.0 and dPsiAlpha > 0.0))
                    {
                        for(auto beta  = Fields(i,j,k).cbegin();
                                 beta != Fields(i,j,k).cend(); ++beta)
                        if(alpha != beta)
                        {
                            FieldsDot(i,j,k).set_sym_pair(alpha->index, beta->index, 0.0, 0.0);
                        }
                    }
                }

                for(auto alpha  = FieldsDot(i,j,k).begin();
                         alpha != FieldsDot(i,j,k).end();)
                {
                    /* Remove zero-sets from the storage */
                    if(alpha->value1 == 0.0 and alpha->value2 == 0.0)
                    {
                        alpha = FieldsDot(i,j,k).erase(alpha);
                    }
                    else
                    {
                        ++alpha;
                    }
                }

                if(FieldsDot(i,j,k).size())
                {
                    /* Limit increments! This is done in a while loop to account for all
                    existing pair-contributions.*/
                    int number_of_iterations = 0;
                    bool LimitingNeeded = true;
                    while (LimitingNeeded)
                    {
                        number_of_iterations++;
                        LimitingNeeded = false;
                        NodeAB<double,double> locIncrements;
                        for(auto it  = FieldsDot(i,j,k).cbegin();
                                 it != FieldsDot(i,j,k).cend(); ++it)
                        {
                            /* Collect increments */
                            if(it->value1 < 0.0)
                            {
                                locIncrements.add_sym2(it->indexA,0, it->value1);
                                locIncrements.add_sym1(it->indexB,0, -it->value1);
                            }
                            if(it->value1 > 0.0)
                            {
                                locIncrements.add_sym1(it->indexA,0, it->value1);
                                locIncrements.add_sym2(it->indexB,0, -it->value1);
                            }
                        }
                        NodeAB<double,double> locLimits;
                        for(auto alpha  = Fields(i,j,k).cbegin();
                                 alpha != Fields(i,j,k).cend(); ++alpha)
                        {
                            /* Calculate limits */
                            double posIncrement = locIncrements.get_sym1(alpha->index,0);
                            double negIncrement = locIncrements.get_sym2(alpha->index,0);
                            double newPFvalue = alpha->value+(posIncrement+negIncrement)*dt;
                            locLimits.set_sym2(alpha->index,0, 1.0);
                            locLimits.set_sym1(alpha->index,0, 1.0);
                            if(newPFvalue < 0.0)
                            {
                                double tmpLim = locLimits.get_sym2(alpha->index,0);
                                double tmpLim2 = min(tmpLim,-(alpha->value+posIncrement*dt)
                                                              /(negIncrement*dt));
                                locLimits.set_sym2(alpha->index,0, tmpLim2);
                            }
                            if(newPFvalue > 1.0)
                            {
                                double tmpLim = locLimits.get_sym1(alpha->index,0);
                                double tmpLim2 = min(tmpLim,(1.0-(alpha->value+negIncrement
                                                             *dt))/(posIncrement*dt));
                                locLimits.set_sym1(alpha->index,0, tmpLim2);
                            }
                        }
                        for(auto it  = FieldsDot(i,j,k).begin();
                                 it != FieldsDot(i,j,k).end(); ++it)
                        {
                            /* Limit increments */
                            if(it->value1 < 0.0)
                            {
                                double tmpLim = min(locLimits.get_sym2(it->indexA,0),
                                                    locLimits.get_sym1(it->indexB,0));
                                it->value1 *= tmpLim;
                                if (tmpLim < 1.0) LimitingNeeded = true;
                            }
                            if(it->value1 > 0.0)
                            {
                                double tmpLim = min(locLimits.get_sym1(it->indexA,0),
                                                    locLimits.get_sym2(it->indexB,0));
                                it->value1 *= tmpLim;
                                if (tmpLim < 1.0) LimitingNeeded = true;
                            }
                        }
                        /* Exit limiting loop, if no convergence after 24 iterations*/
                        if (number_of_iterations > 24) LimitingNeeded = false;
                    }

                    for(auto alpha  = FieldsDot(i,j,k).begin();
                             alpha != FieldsDot(i,j,k).end(); )
                    {
                        /* Clean up values which are too small */
                        if(fabs(alpha->value1*dt) < DBL_EPSILON)
                        {
                            alpha = FieldsDot(i,j,k).erase(alpha);
                        }
                        else
                        {
                            ++alpha;
                        }
                    }

                    /* End plausibility check */
                    NodePF tmpPF = Fields(i,j,k);
                    for(auto psi  = FieldsDot(i,j,k).cbegin();
                             psi != FieldsDot(i,j,k).cend(); ++psi)
                    if(psi->value1 != 0.0)
                    {
                        tmpPF.add_value(psi->indexA,  psi->value1 * dt);
                        tmpPF.add_value(psi->indexB, -psi->value1 * dt);
                    }
                    for(auto it  = tmpPF.cbegin();
                             it != tmpPF.cend(); ++it)
                    if(it->value < -precision or it->value > 1.0 + precision)
                    {
                        double oldFields = Fields(i,j,k).get_value(it->index);
                        string msg = "Normalizing of phase field increments failed in point ("
                                   + to_string(i) + "," + to_string(j) + "," + to_string(k)
                                   + "). " + to_string(Fields(i,j,k).size())
                                   + " fields present. Grain "
                                   + to_string(it->index) + " with a fields-value of "
                                   + to_string(oldFields)
                                   + " is incremented by "
                                   + to_string(it->value-oldFields)
                                   + ", which results in a fields-value of "
                                   + to_string(it->value)
                                   + ". This will result in undefined behavior!";
                        ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                    }
                }
                break;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    SetIncrementsBoundaryConditionsSR(BC);

#ifdef DEBUG

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if(Fields(i,j,k).flag)
    {
        double zero = 1.0;

        for(size_t n = 0; n < Nphases; n++)
        {
            double fraction = Fractions(i,j,k)({n});
            double newfraction = fraction;

            zero -= fraction;

            for(auto it = FieldsDot(i,j,k).begin();
                     it != FieldsDot(i,j,k).end(); ++it)
            {
                if((FieldsProperties[it->indexA].Phase == n)
                and(FieldsProperties[it->indexB].Phase != n))
                {
                    newfraction += it->value1*dt;
                }
                else if((FieldsProperties[it->indexA].Phase != n)
                     and(FieldsProperties[it->indexB].Phase == n))
                {
                    newfraction -= it->value1*dt;
                }
            }

            if((newfraction < -precision)
            or (newfraction >  precision + 1.0)
            or (fraction < -precision)
            or (fraction >  precision + 1.0))
            {
                string msg = "Redistribution of composition not possible. "
                             "Phase-field broke in point ("
                           + to_string(i) + "," + to_string(j) + ","
                           + to_string(k) + ") for phase " + to_string(n)
                           + ". Old fraction "
                           + to_string(fraction)
                           + " New fraction " + to_string(newfraction)
                           + ". This will break mass conservation!";

                ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                //OP_Exit(EXIT_FAILURE);
            }
        }

        if(zero < -precision or zero > precision)
        {
            string msg = "Redistribution of composition not possible. "
                         "Phase-field broke in point ("
                       + to_string(i) + "," + to_string(j) + ","
                       + to_string(k) + "). Phase fractions don't add up to "
                         "unity, difference is " + to_string(zero);

           ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
           OP_Exit(EXIT_FAILURE);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}

void PhaseField::NormalizeIncrementsDR(const BoundaryConditions& BC, const double dt)
{
    /** This function limits phase-field increments for all present phase-field
    pairs, so that the actual phase-field values are within their natural
    limits of 0.0 and 1.0.*/

    double precision = FLT_EPSILON;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDotDR,0,)
    if (FieldsDR(i,j,k).wide_interface())
    {
        switch(FieldsDotDR(i,j,k).size())
        {
            case 0:
            {
                break;
            }
            case 1:
            {
                size_t indexA = FieldsDotDR(i,j,k).front().indexA;
                double valueA = (FieldsDotDR(i,j,k).front().value1
                              +  FieldsDotDR(i,j,k).front().value2)*dt;
                double old_value = FieldsDR(i,j,k).get_value(indexA);

                if((old_value == 0.0 and valueA < 0.0) or
                   (old_value == 1.0 and valueA > 0.0))
                {
                    FieldsDotDR(i,j,k).clear();
                }
                else
                {
                    double new_value = old_value + valueA;

                    double norm = 1.0;
                    if(new_value < 0.0)
                    {
                        norm *= -old_value/valueA;
                    }
                    else if(new_value > 1.0)
                    {
                        norm *= (1.0 - old_value)/valueA;
                    }

                    if(norm > DBL_EPSILON)
                    {
                        FieldsDotDR(i,j,k) *= norm;
                    }
                    else
                    {
                        FieldsDotDR(i,j,k).clear();
                    }
                }
                break;
            }
            default:
            {

                for(auto alpha  = FieldsDR(i,j,k).cbegin();
                         alpha != FieldsDR(i,j,k).cend(); ++alpha)
                {
                    double dPsiAlpha = 0.0;
                    for(auto beta  = FieldsDR(i,j,k).cbegin();
                             beta != FieldsDR(i,j,k).cend(); ++beta)
                    if(alpha != beta)
                    {
                        dPsiAlpha += FieldsDotDR(i,j,k).get_asym1(alpha->index,
                                                                   beta->index);
                        dPsiAlpha += FieldsDotDR(i,j,k).get_asym2(alpha->index,
                                                                   beta->index);
                    }
                    /* Set out of bounds sets to zero */
                    if((alpha->value == 0.0 and dPsiAlpha < 0.0)
                    or (alpha->value == 1.0 and dPsiAlpha > 0.0))
                    {
                        for(auto beta  = FieldsDR(i,j,k).cbegin();
                                 beta != FieldsDR(i,j,k).cend(); ++beta)
                        if(alpha != beta)
                        {
                            FieldsDotDR(i,j,k).set_sym_pair(alpha->index, beta->index, 0.0, 0.0);
                        }
                    }
                }
                for(auto alpha  = FieldsDotDR(i,j,k).begin();
                         alpha != FieldsDotDR(i,j,k).end();)
                {
                    /* Remove zero-sets from the storage */
                    if(alpha->value1 == 0.0 and alpha->value2 == 0.0)
                    {
                        alpha = FieldsDotDR(i,j,k).erase(alpha);
                    }
                    else
                    {
                        ++alpha;
                    }
                }
                if(FieldsDotDR(i,j,k).size())
                {
                    //cout << "while loop!" << endl;
                    /* Limit increments! This is done in a while loop, to acknowledge all
                    existing pair-contributions.*/
                    int number_of_iterations = 0;
                    bool LimitingNeeded = true;
                    while (LimitingNeeded)
                    {
                        number_of_iterations++;
                        LimitingNeeded = false;
                        NodeAB<double,double> locIncrements;
                        for(auto it  = FieldsDotDR(i,j,k).cbegin();
                                 it != FieldsDotDR(i,j,k).cend(); ++it)
                        {
                            /* Collect increments */
                            if(it->value1 < 0.0)
                            {
                                locIncrements.add_sym2(it->indexA,0,  it->value1);
                                locIncrements.add_sym1(it->indexB,0, -it->value1);
                            }
                            if(it->value1 > 0.0)
                            {
                                locIncrements.add_sym1(it->indexA,0,  it->value1);
                                locIncrements.add_sym2(it->indexB,0, -it->value1);
                            }
                        }
                        NodeAB<double,double> locLimits;
                        for(auto alpha  = FieldsDR(i,j,k).cbegin();
                                 alpha != FieldsDR(i,j,k).cend(); ++alpha)
                        {
                            /* Calculate limits */
                            double posIncrement = locIncrements.get_sym1(alpha->index,0);
                            double negIncrement = locIncrements.get_sym2(alpha->index,0);
                            double newPFvalue = alpha->value+(posIncrement+negIncrement)*dt;
                            locLimits.set_sym2(alpha->index,0,1.0);
                            locLimits.set_sym1(alpha->index,0,1.0);
                            if(newPFvalue < 0.0)
                            {
                                double tmpLim = locLimits.get_sym2(alpha->index,0);
                                double tmpLim2 = min(tmpLim,-(alpha->value+posIncrement*dt)
                                                              /(negIncrement*dt));
                                locLimits.set_sym2(alpha->index,0, tmpLim2);
                            }
                            if(newPFvalue > 1.0)
                            {
                                double tmpLim = locLimits.get_sym1(alpha->index,0);
                                double tmpLim2 = min(tmpLim,(1.0-(alpha->value+negIncrement*dt))
                                                                /(posIncrement*dt));
                                locLimits.set_sym1(alpha->index,0,tmpLim2);
                            }
                        }
                        for(auto it  = FieldsDotDR(i,j,k).begin();
                                 it != FieldsDotDR(i,j,k).end(); ++it)
                        {
                            /* Limit increments */
                            if(it->value1 < 0.0)
                            {
                                double tmpLim = min(locLimits.get_sym2(it->indexA,0),
                                                    locLimits.get_sym1(it->indexB,0));
                                it->value1 *= tmpLim;
                                if (tmpLim < 1.0)
                                LimitingNeeded = true;
                            }
                            if(it->value1 > 0.0)
                            {
                                double tmpLim = min(locLimits.get_sym1(it->indexA,0),
                                                    locLimits.get_sym2(it->indexB,0));
                                it->value1 *= tmpLim;
                                if (tmpLim < 1.0)
                                LimitingNeeded = true;
                            }
                        }
                        /* Exit limiting loop, if no convergence after 24 iterations*/
                        if (number_of_iterations > 24)
                        LimitingNeeded = false;
                    }
                    for(auto alpha  = FieldsDotDR(i,j,k).begin();
                             alpha != FieldsDotDR(i,j,k).end(); )
                    {
                        /* Clean up values which are too small */
                        if(fabs(alpha->value1) < DBL_EPSILON)
                        {
                            alpha = FieldsDotDR(i,j,k).erase(alpha);
                        }
                        else
                        {
                            ++alpha;
                        }
                    }

                    /* End plausibility check */
                    NodePF tmpPF = FieldsDR(i,j,k);
                    for(auto psi  = FieldsDotDR(i,j,k).cbegin();
                             psi != FieldsDotDR(i,j,k).cend(); ++psi)
                    if(psi->value1 != 0.0)
                    {
                        tmpPF.add_value(psi->indexA,  psi->value1 * dt);
                        tmpPF.add_value(psi->indexB, -psi->value1 * dt);
                    }
                    for(auto it  = tmpPF.cbegin();
                             it != tmpPF.cend(); ++it)
                    if(it->value < -precision or it->value > 1.0 + precision)
                    {
                        double oldFields = FieldsDR(i,j,k).get_value(it->index);
                        string msg = "Normalizing of phase field increments failed in point ("
                                   + to_string(i) + "," + to_string(j) + "," + to_string(k)
                                   + "). " + to_string(FieldsDR(i,j,k).size())
                                   + " fields present. Grain "
                                   + to_string(it->index) + " with a fields-value of "
                                   + to_string(oldFields)
                                   + " is incremented by "
                                   + to_string(it->value-oldFields)
                                   + ", which results in a fields-value of "
                                   + to_string(it->value)
                                   + ". This will result in undefined behavior!";
                        ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                    }
                }
                break;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    //SetIncrementsBoundaryConditionsDR(BC);
    CoarsenDot();
    SetIncrementsBoundaryConditionsSR(BC);

#ifdef DEBUG

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if(Fields(i,j,k).flag)
    {
        double zero = 1.0;

        for(size_t n = 0; n < Nphases; n++)
        {
            double fraction = Fractions(i,j,k,{n});
            double newfraction = fraction;

            zero -= fraction;

            for(auto it = FieldsDot(i,j,k).begin();
                     it != FieldsDot(i,j,k).end(); ++it)
            {
                if((FieldsProperties[it->indexA].Phase == n)
                and(FieldsProperties[it->indexB].Phase != n))
                {
                    newfraction += it->value1*dt;
                }
                else if((FieldsProperties[it->indexA].Phase != n)
                     and(FieldsProperties[it->indexB].Phase == n))
                {
                    newfraction -= it->value1*dt;
                }
            }

            if((newfraction < -precision)
            or (newfraction >  precision+1.0)
            or (fraction < -precision)
            or (fraction >  precision+1.0))
            {
                string msg = "Redistribution of composition not possible. "
                             "Phase-field broke in point ("
                           + to_string(i) + "," + to_string(j) + ","
                           + to_string(k) + ") for phase " + to_string(n)
                           + ". Old fraction "
                           + to_string(fraction)
                           + " New fraction " + to_string(newfraction)
                           + ". This will break mass conservation!";

                ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
                //OP_Exit(EXIT_FAILURE);
            }
        }

        if(zero < -precision or zero > precision)
        {
            string msg = "Redistribution of composition not possible. "
                         "Phase-field broke in point ("
                       + to_string(i) + "," + to_string(j) + ","
                       + to_string(k) + "). Phase fractions don't add up to "
                         "unity, difference is " + to_string(zero);

           ConsoleOutput::WriteWarning(msg,thisclassname,"NormalizeIncrements");
           OP_Exit(EXIT_FAILURE);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}

void PhaseField::KeepPhaseFieldsVolume(void)
{
    Tensor<bool,2> AllowedTransitions({Nphases,Nphases});
    // Disallow all transitions
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        AllowedTransitions({n,m}) = false;
    }

    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            KeepPhaseVolumeSR(AllowedTransitions);
            break;
        }
        case Resolutions::Dual:
        {
            KeepPhaseVolumeDR(AllowedTransitions);
            break;
        }
    }
}

void PhaseField::KeepPhaseVolumeSR(Tensor<bool,2> AllowedTransitions)
{
    vector<NodeAB<double,double>> PairwiseVolumeChange;

    size_t Nthreads = 1;
#ifdef _OPENMP
    Nthreads = omp_get_max_threads();
#endif
    PairwiseVolumeChange.resize(Nthreads);

    // Collect pairwise increments and corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        int thNum = 0;
#ifdef _OPENMP
        thNum = omp_get_thread_num();
#endif
        if(Fields(i,j,k).wide_interface())
        for(auto alpha = Fields(i,j,k).cbegin();
                 alpha != Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != Fields(i,j,k).cend(); ++beta)
        {
            size_t pIndexA = FieldsProperties[alpha->index].Phase;
            size_t pIndexB = FieldsProperties[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                double locPhiDot = FieldsDot(i,j,k).get_asym1(alpha->index, beta->index);
                PairwiseVolumeChange[thNum].add_asym1(alpha->index, beta->index, locPhiDot);
                PairwiseVolumeChange[thNum].add_sym2(alpha->index, beta->index, sqrt(alpha->value*beta->value));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    for(size_t th = 1; th < Nthreads; th++)
    {
        PairwiseVolumeChange[0].add_asym1(PairwiseVolumeChange[th]);
        PairwiseVolumeChange[0].add_sym2(PairwiseVolumeChange[th]);
    }
    // Distribute pairwise increments using corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, 0, )
    {
        if(Fields(i,j,k).interface())
        for(auto alpha = Fields(i,j,k).cbegin();  alpha != Fields(i,j,k).cend(); ++alpha)
        for(auto  beta = alpha + 1; beta != Fields(i,j,k).cend(); ++beta)
        if(alpha->value*beta->value != 0.0)
        {
            size_t pIndexA = FieldsProperties[alpha->index].Phase;
            size_t pIndexB = FieldsProperties[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                if(PairwiseVolumeChange[0].get_sym2(alpha->index,beta->index) != 0.0)
                {
                    FieldsDot(i,j,k).add_asym1(alpha->index, beta->index,
                        -PairwiseVolumeChange[0].get_asym1(alpha->index,beta->index)*sqrt(alpha->value*beta->value)/
                         PairwiseVolumeChange[0].get_sym2(alpha->index,beta->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::KeepPhaseVolumeDR(Tensor<bool,2> AllowedTransitions)
{
    vector<NodeAB<double,double>> PairwiseVolumeChange;

    size_t Nthreads = 1;
#ifdef _OPENMP
    Nthreads = omp_get_max_threads();
#endif
    PairwiseVolumeChange.resize(Nthreads);

    // Collect pairwise increments and corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, 0, )
    {
        int thNum = 0;
#ifdef _OPENMP
        thNum = omp_get_thread_num();
#endif
        if(FieldsDR(i,j,k).wide_interface())
        for(auto alpha = FieldsDR(i,j,k).cbegin();
                 alpha != FieldsDR(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != FieldsDR(i,j,k).cend(); ++beta)
        {
            size_t pIndexA = FieldsProperties[alpha->index].Phase;
            size_t pIndexB = FieldsProperties[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                double locPhiDot = FieldsDotDR(i,j,k).get_asym1(alpha->index, beta->index);
                PairwiseVolumeChange[thNum].add_asym1(alpha->index, beta->index, locPhiDot);
                PairwiseVolumeChange[thNum].add_sym2(alpha->index, beta->index, sqrt(alpha->value*beta->value));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    for(size_t th = 1; th < Nthreads; th++)
    {
        PairwiseVolumeChange[0].add_asym1(PairwiseVolumeChange[th]);
        PairwiseVolumeChange[0].add_sym2(PairwiseVolumeChange[th]);
    }
    // Distribute pairwise increments using corresponding weights
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, FieldsDR, 0, )
    {
        if(FieldsDR(i,j,k).interface())
        for(auto alpha = FieldsDR(i,j,k).cbegin();
                 alpha != FieldsDR(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta != FieldsDR(i,j,k).cend(); ++beta)
        if(alpha->value*beta->value != 0.0)
        {
            size_t pIndexA = FieldsProperties[alpha->index].Phase;
            size_t pIndexB = FieldsProperties[ beta->index].Phase;
            if(not AllowedTransitions({pIndexA, pIndexB}))
            {
                if(PairwiseVolumeChange[0].get_sym2(alpha->index,beta->index) != 0.0)
                {
                    FieldsDotDR(i,j,k).add_asym1(alpha->index, beta->index,
                        -PairwiseVolumeChange[0].get_asym1(alpha->index,beta->index)*sqrt(alpha->value*beta->value)/
                         PairwiseVolumeChange[0].get_sym2(alpha->index,beta->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(Grid.Resolution == Resolutions::Dual)
    {
        SetBoundaryConditionsDR(BC);
    }
    SetBoundaryConditionsSR(BC);
}

void PhaseField::SetBoundaryConditionsSR(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(Fields);
    if(Grid.dNy) BC.SetY(Fields);
    if(Grid.dNz) BC.SetZ(Fields);
}

void PhaseField::SetBoundaryConditionsDR(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(FieldsDR);
    if(Grid.dNy) BC.SetY(FieldsDR);
    if(Grid.dNz) BC.SetZ(FieldsDR);
}

void PhaseField::SetIncrementsBoundaryConditionsSR(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(FieldsDot);
    if(Grid.dNy) BC.SetY(FieldsDot);
    if(Grid.dNz) BC.SetZ(FieldsDot);
}

void PhaseField::SetIncrementsBoundaryConditionsDR(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(FieldsDotDR);
    if(Grid.dNy) BC.SetY(FieldsDotDR);
    if(Grid.dNz) BC.SetZ(FieldsDotDR);
}

void PhaseField::PrintPointStatistics(const int x, const int y, const int z) const
{
    ConsoleOutput::WriteStandard("Point", iVector3{x,y,z});
    for (auto alpha  = Fields(x,y,z).cbegin();
              alpha != Fields(x,y,z).cend(); ++alpha)
    {
        std::string vname = "Phase Field " + std::to_string(alpha->index);
        ConsoleOutput::WriteStandard(vname, alpha->value);
    }
    for (size_t PhaseIdx = 0; PhaseIdx != Nphases; PhaseIdx++)
    {
        std::string vname = "Phase Fraction " + PhaseNames[PhaseIdx];
        ConsoleOutput::WriteStandard(vname, Fractions(x,y,z,{PhaseIdx}));
    }
    ConsoleOutput::WriteBlankLine();
}

void PhaseField::WriteVTK(Settings& locSettings, const int tStep, const bool CurvatureOutput, const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"_", tStep, ".vts");
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return Interfaces(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return Fields(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return Fields(i,j,k).majority_index();}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){return Fractions(i,j,k,{n});}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions",   [this](int i,int j,int k){return Fields(i,j,k).size();}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",    [this](int i,int j,int k){return Variants(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"ParentGrain", [this](int i,int j,int k){return ParentGrain(i,j,k);}});

            if(CurvatureOutput)
            for (size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"Curvature_" + std::to_string(n),            [n,this](int i,int j,int k){return CurvaturePhase(i,j,k,n);}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_1_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[0];}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_2_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[1];}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision);
#ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with phase field data", "application/xml");
#endif
            break;
        }
        case Resolutions::Dual:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return InterfacesDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return FieldsDR(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return FieldsDR(i,j,k).majority_index();}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){
                double FractionsValue = 0.0;
                for (auto it = FieldsDR(i,j,k).cbegin();
                          it != FieldsDR(i,j,k).cend(); ++it)
                {
                    size_t pIndex = FieldsProperties[it->index].Phase;
                    if(pIndex == n)
                    {
                        FractionsValue += it->value;
                    }
                }
                return FractionsValue;}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions",   [this](int i,int j,int k){return FieldsDR(i,j,k).size();}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",    [this](int i,int j,int k){return VariantsDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"ParentGrain", [this](int i,int j,int k){return ParentGrainDR(i,j,k);}});
            if(CurvatureOutput)
            for (size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"Curvature_" + std::to_string(n),[n,this](int i,int j,int k){return CurvaturePhaseDR(i,j,k,n);}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
#ifndef WIN32
            locSettings.Meta.AddPart(Filename, "File", "VTK file with phase field data", "application/xml");
#endif
            break;
        }
    }
}

void PhaseField::WriteDistortedVTK(
        Settings& locSettings,
        const ElasticProperties& EP,
        const int tStep,
        const bool CurvatureOutput,
        const int precision) const
{
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, thisclassname+"Distorted_", tStep, ".vts");
    std::vector<VTK::Field_t> ListOfFields;
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return Interfaces(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return Fields(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return Fields(i,j,k).majority_index();}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){return Fractions(i,j,k,{n});}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions",   [this](int i,int j,int k){return Fields(i,j,k).size();}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",    [this](int i,int j,int k){return Variants(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"ParentGrain", [this](int i,int j,int k){return ParentGrain(i,j,k);}});
            if(CurvatureOutput)
            for (size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"Curvature_" + std::to_string(n),            [n,this](int i,int j,int k){return CurvaturePhase(i,j,k,n);}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_1_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[0];}});
                ListOfFields.push_back((VTK::Field_t) {"PrincipalCurvature_2_" + std::to_string(n), [n,this](int i,int j,int k){return PrincipalCurvatures(i,j,k,n)[1];}});
            }
            VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields, precision);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with distorted phase field data", "application/xml");
            #endif
            break;
        }
        case Resolutions::Dual:
        {
            ListOfFields.push_back((VTK::Field_t) {"Interfaces",  [this](int i,int j,int k){return InterfacesDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"Flags",       [this](int i,int j,int k){return FieldsDR(i,j,k).flag;}});
            ListOfFields.push_back((VTK::Field_t) {"PhaseFields", [this](int i,int j,int k){return FieldsDR(i,j,k).majority_index();}});
            for(size_t n = 0; n < Nphases; n++)
            {
                ListOfFields.push_back((VTK::Field_t) {"PhaseFraction_" + std::to_string(n), [n,this](int i,int j,int k){
                double FractionsValue = 0.0;
                for (auto it = FieldsDR(i,j,k).cbegin();
                          it != FieldsDR(i,j,k).cend(); ++it)
                {
                    size_t pIndex = FieldsProperties[it->index].Phase;
                    if(pIndex == n)
                    {
                        FractionsValue += it->value;
                    }
                }
                return FractionsValue;}});
            }
            ListOfFields.push_back((VTK::Field_t) {"Junctions",   [this](int i,int j,int k){return FieldsDR(i,j,k).size();}});
            ListOfFields.push_back((VTK::Field_t) {"Variants",    [this](int i,int j,int k){return VariantsDR(i,j,k);}});
            ListOfFields.push_back((VTK::Field_t) {"ParentGrain", [this](int i,int j,int k){return ParentGrainDR(i,j,k);}});

            VTK::WriteDistorted(Filename, locSettings, EP, ListOfFields, precision, 2);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with distorted phase field data", "application/xml");
            #endif
            break;
        }
    }
}

void PhaseField::WriteLaplacianVTK(Settings& locSettings,
                                   const int tStep,
                                   size_t PhiIndex,
                                   const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Laplacian_" + to_string(PhiIndex) + "_", tStep, ".vts");
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            ListOfFields.push_back((VTK::Field_t) {"Laplacian_" + std::to_string(PhiIndex), [PhiIndex,this](int i,int j,int k){return Fields(i,j,k).get_laplacian(PhiIndex);}});
            VTK::Write(Filename, locSettings, ListOfFields, precision);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with phase field laplacian data", "application/xml");
            #endif
            break;
        }
        case Resolutions::Dual:
        {
            ListOfFields.push_back((VTK::Field_t) {"Laplacian_" + std::to_string(PhiIndex), [PhiIndex,this](int i,int j,int k){return FieldsDR(i,j,k).get_laplacian(PhiIndex);}});
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with phase field laplacian data", "application/xml");
            #endif
            break;
        }
    }
}

void PhaseField::WriteIndividualPhaseFieldValuesVTK(Settings& locSettings,
                                                    const int tStep,
                                                    const std::initializer_list<size_t> FieldIndices,
                                                    const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "PhaseFieldValues_", tStep, ".vts");
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            for (auto IteratorIndex = FieldIndices.begin();
                      IteratorIndex != FieldIndices.end(); IteratorIndex++)
            {
                ListOfFields.push_back((VTK::Field_t) {"FieldValue_" + std::to_string(*IteratorIndex), [IteratorIndex,this](int i,int j,int k){return Fields(i,j,k).get_value(*IteratorIndex);}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with individual phase field data", "application/xml");
            #endif
            break;
        }
        case Resolutions::Dual:
        {
            for (auto IteratorIndex = FieldIndices.begin();
                      IteratorIndex != FieldIndices.end(); IteratorIndex++)
            {
                ListOfFields.push_back((VTK::Field_t) {"FieldValue_" + std::to_string(*IteratorIndex), [IteratorIndex,this](int i,int j,int k){return FieldsDR(i,j,k).get_value(*IteratorIndex);}});
            }
            VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
            #ifndef WIN32       
            locSettings.Meta.AddPart(Filename, "File", "VTK file with individual phase field data", "application/xml");
            #endif
            break;
        }
    }
}

bool PhaseField::WriteMPI(const std::string& FileName)
{
	Fields.WriteToFile(FileName+"/Fields.data");
    return true;
}

bool PhaseField::Write(const std::string& FileName) const
{
    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        ConsoleOutput::WriteWarning("File \"" + FileName + "\" could not be opened",
                thisclassname, "Write()");
        return false;
    };

    int Nx = Grid.Nx;
    int Ny = Grid.Ny;
    int Nz = Grid.Nz;

    out.write(reinterpret_cast<const char*>(&Nx), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Ny), sizeof(int));
    out.write(reinterpret_cast<const char*>(&Nz), sizeof(int));
	
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).write(out);
            }
            STORAGE_LOOP_END
            break;
        }
        case Resolutions::Dual:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).write(out);
            }
            STORAGE_LOOP_END

            STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0)
            {
                FieldsDR(i,j,k).write(out);
            }
            STORAGE_LOOP_END
            break;
        }
    }
    out.close();
    return true;
}

bool PhaseField::Write(const Settings& locSettings, const int tStep) const
{
    #ifdef MPI_PARALLEL
    string FileName =
        FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
    #else
    string FileName =
        FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_", tStep, ".dat");
    #endif
    bool write_success = Write(FileName);
    write_success = write_success && FieldsProperties.Write(locSettings, tStep);
    return write_success;
}

void PhaseField::WriteH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    dbuffer.push_back(Grid.Nx);
    dbuffer.push_back(Grid.Ny);
    dbuffer.push_back(Grid.Nz);
    H5.WriteCheckPoint(tStep, "PFDomain", dbuffer);
    dbuffer.clear();
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            dbuffer = Fields.pack();
            H5.WriteCheckPoint(tStep, "Fields", dbuffer);
            break;
        }
        case Resolutions::Dual:
        {
            dbuffer = FieldsDR.pack();
            H5.WriteCheckPoint(tStep, "FieldsDR", dbuffer);
            break;
        }
    }
    FieldsProperties.WriteH5(tStep,H5);
    #else
    ConsoleOutput::WriteWarning("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteH5()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

bool PhaseField::ReadH5(const BoundaryConditions& BC, H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "PFDomain", dbuffer);
    int locNx = dbuffer[0];
    int locNy = dbuffer[1];
    int locNz = dbuffer[2];
    dbuffer.clear();
    if(locNx != Grid.Nx or locNy != Grid.Ny or locNz != Grid.Nz)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                << "Required data dimensions: (" << Grid.Nx << ", " << Grid.Ny << ", " << Grid.Nz << ") grid points.\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }
    switch (Resolution)
    {
        case Resolutions::Single:
        {
            H5.ReadCheckPoint(tStep, "Fields", dbuffer);
            Fields.unpack(dbuffer);
            break;
        }
        case Resolutions::Dual:
        {
            H5.ReadCheckPoint(tStep, "FieldsDR", dbuffer);
            FieldsDR.unpack(dbuffer);
            break;
        }
    }
    FieldsProperties.ReadH5(tStep,H5);
    Finalize(BC);
    #else
    ConsoleOutput::WriteWarning("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "ReadH5()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
    return true;
}

bool PhaseField::Read(string FileName)
{
    ifstream inp(FileName.c_str(), ios::in | ios::binary);

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
                << "Input data dimensions: ("
                << locNx << ", "
                << locNy << ", "
                << locNz << ") grid points.\n"
                << "Required data dimensions: ("
                << Grid.Nx << ", "
                << Grid.Ny << ", "
                << Grid.Nz << ") grid points.\n";
        ConsoleOutput::WriteWarning(message.str(), thisclassname, "Read()");
        return false;
    }
    switch (Grid.Resolution)
    {
        case Resolutions::Single:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).read(inp);
            }
            STORAGE_LOOP_END
            break;
        }
        case Resolutions::Dual:
        {
            STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
            {
                Fields(i,j,k).read(inp);
            }
            STORAGE_LOOP_END

            STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0)
            {
                FieldsDR(i,j,k).read(inp);
            }
            STORAGE_LOOP_END
            break;
        }
    }
    inp.close();
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

bool PhaseField::Read(const Settings& locSettings, const BoundaryConditions& BC, int tStep)
{
#ifdef MPI_PARALLEL
    string FileName =
        FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_"+ std::to_string(MPI_RANK) + "_", tStep, ".dat");
#else
    string FileName =
        FileInterface::MakeFileName(locSettings.InputRawDataDir,thisclassname+"_", tStep, ".dat");
#endif

    bool read_status = Read(FileName);
    read_status = read_status && FieldsProperties.Read(locSettings, tStep);
    Finalize(BC);
    return read_status;
}
void PhaseField::WriteAverageVolume(const int tStep, const size_t PhaseIndex) const
{
    stringstream converter;
    converter << PhaseIndex;

    string FileName = string("AverageVolumeOfPhase_")
        + converter.str() + string(".txt");

    if(!tStep)
    {
        fstream tout(FileName.c_str(), ios::out);
        tout << "Time\tAvgVolume" << endl;
        tout.close();
    }

    double avgVol = 0;
    size_t count = 0;
    for(size_t n = 0; n < FieldsProperties.size(); n++)
    if(FieldsProperties[n].Exist and FieldsProperties[n].Phase == PhaseIndex)
    {
        avgVol += FieldsProperties[n].Volume;
        count += 1.0;
    }

    fstream out(FileName.c_str(), ios::out | ios::app);
    if(count)
    {
        out << tStep << "\t" << avgVol/count << endl;
    }
    else
    {
        out << tStep << "\t" << 0.0 << endl;
    }
    out.close();
}

void PhaseField::MoveFrame(const int dx, const int dy, const int dz,
                           const BoundaryConditions& BC)
{
    MoveFrameSR(dx,dy,dz,BC);
    if(Grid.Resolution == Resolutions::Dual)
    {
        MoveFrameDR(dx,dy,dz,BC);
    }
    Finalize(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Frame moved");
}

void PhaseField::MoveFrameSR(const int dx, const int dy, const int dz,
                             const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    int xEnd = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    int yEnd = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    int zEnd = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Fields(i, j, k) = Fields(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        Fractions(i, j, k) = Fractions(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);

    }
    SetBoundaryConditionsSR(BC);
}

void PhaseField::MoveFrameDR(const int dx, const int dy, const int dz,
                             const BoundaryConditions& BC)
{
    for(int n = 0; n <= 1; n++)
    {
        int xBeg = (dx >= 0) + (dx < 0)*(FieldsDR.sizeX()) - 1;
        int xEnd = (dx >= 0)*(FieldsDR.sizeX()) + (dx < 0) - 1;
        int xInc = 1 - 2*(dx < 0);

        int yBeg = (dy >= 0) + (dy < 0)*(FieldsDR.sizeY()) - 1;
        int yEnd = (dy >= 0)*(FieldsDR.sizeY()) + (dy < 0) - 1;
        int yInc = 1 - 2*(dy < 0);

        int zBeg = (dz >= 0) + (dz < 0)*(FieldsDR.sizeZ()) - 1;
        int zEnd = (dz >= 0)*(FieldsDR.sizeZ()) + (dz < 0) - 1;
        int zInc = 1 - 2*(dz < 0);

        for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
        for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
        for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
        {
            FieldsDR(i, j, k) = FieldsDR(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        }
    }
    SetBoundaryConditionsDR(BC);
}

void PhaseField::ConsumePlane(const int dx, const int dy, const int dz,
                              const int x, const int y, const int z,
                              const BoundaryConditions& BC)
{
    ConsumePlaneSR(dx,dy,dz,x,y,z,BC);
    if(Grid.Resolution == Resolutions::Dual)
    {
        ConsumePlaneDR(dx,dy,dz,x,y,z,BC);
    }
    Finalize(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Plane consumed");
}

void PhaseField::ConsumePlaneSR(const int dx, const int dy, const int dz,
                                const int x, const int y, const int z,
                                const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    int xEnd = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    int yEnd = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    int zEnd = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz >= 0)
    {
        Fields(i, j, k) = Fields(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        Fractions(i, j, k) = Fractions(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
    }
    xBeg = (dx >= 0)*(Grid.Nx) + (dx < 0) - 1;
    xEnd = (dx >= 0) + (dx < 0)*(Grid.Nx) - 1;
    xInc = 2*(dx < 0) - 1;

    yBeg = (dy >= 0)*(Grid.Ny) + (dy < 0) - 1;
    yEnd = (dy >= 0) + (dy < 0)*(Grid.Ny) - 1;
    yInc = 2*(dy < 0) - 1;

    zBeg = (dz >= 0)*(Grid.Nz) + (dz < 0) - 1;
    zEnd = (dz >= 0) + (dz < 0)*(Grid.Nz) - 1;
    zInc = 2*(dz < 0) - 1;

    for(int i = xBeg; ((dx >= 0) and (i >= xEnd)) or ((dx < 0) and (i <= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j >= yEnd)) or ((dy < 0) and (j <= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k >= zEnd)) or ((dz < 0) and (k <= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz < 0)
    {
        Fields(i, j, k) = Fields(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
        Fractions(i, j, k) = Fractions(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
    }
    SetBoundaryConditions(BC);
}

void PhaseField::ConsumePlaneDR(const int dx, const int dy, const int dz,
                                const int x, const int y, const int z,
                                const BoundaryConditions& BC)
{
    for(int n = 0; n <= 1; n++)
    {
        int xBeg = (dx >= 0) + (dx < 0)*(FieldsDR.sizeX()) - 1;
        int xEnd = (dx >= 0)*(FieldsDR.sizeX()) + (dx < 0) - 1;
        int xInc = 1 - 2*(dx < 0);

        int yBeg = (dy >= 0) + (dy < 0)*(FieldsDR.sizeY()) - 1;
        int yEnd = (dy >= 0)*(FieldsDR.sizeY()) + (dy < 0) - 1;
        int yInc = 1 - 2*(dy < 0);

        int zBeg = (dz >= 0) + (dz < 0)*(FieldsDR.sizeZ()) - 1;
        int zEnd = (dz >= 0)*(FieldsDR.sizeZ()) + (dz < 0) - 1;
        int zInc = 1 - 2*(dz < 0);

        for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
        for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
        for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
        if((i-x)*dx+(j-y)*dy+(k-z)*dz >= 0)
        {
            FieldsDR(i, j, k) = FieldsDR(i + dx*Grid.dNx, j + dy*Grid.dNy, k + dz*Grid.dNz);
        }
        xBeg = (dx >= 0)*(FieldsDR.sizeX()) + (dx < 0) - 1;
        xEnd = (dx >= 0) + (dx < 0)*(FieldsDR.sizeX()) - 1;
        xInc = 2*(dx < 0) - 1;

        yBeg = (dy >= 0)*(FieldsDR.sizeY()) + (dy < 0) - 1;
        yEnd = (dy >= 0) + (dy < 0)*(FieldsDR.sizeY()) - 1;
        yInc = 2*(dy < 0) - 1;

        zBeg = (dz >= 0)*(FieldsDR.sizeZ()) + (dz < 0) - 1;
        zEnd = (dz >= 0) + (dz < 0)*(FieldsDR.sizeZ()) - 1;
        zInc = 2*(dz < 0) - 1;

        for(int i = xBeg; ((dx >= 0) and (i >= xEnd)) or ((dx < 0) and (i <= xEnd)); i += xInc)
        for(int j = yBeg; ((dy >= 0) and (j >= yEnd)) or ((dy < 0) and (j <= yEnd)); j += yInc)
        for(int k = zBeg; ((dz >= 0) and (k >= zEnd)) or ((dz < 0) and (k <= zEnd)); k += zInc)
        if((i-x)*dx+(j-y)*dy+(k-z)*dz < 0)
        {
            FieldsDR(i, j, k) = FieldsDR(i - dx*Grid.dNx, j - dy*Grid.dNy, k - dz*Grid.dNz);
        }
    }
    SetBoundaryConditionsDR(BC);
}

void PhaseField::PrintPFVolumes() const
{
    ConsoleOutput::WriteLineInsert("Phase-field volumes","=");

    for(unsigned int idx = 0; idx < FieldsProperties.size(); idx++)
    {
        if(FieldsProperties[idx].Volume > 0.0)
        {
            ConsoleOutput::WriteStandardNarrow("PF",      idx);
            ConsoleOutput::WriteStandardNarrow("Variant", FieldsProperties[idx].Variant);
            switch(FieldsProperties[idx].Stage)
            {
                case GrainStages::Seed:
                {
                    ConsoleOutput::WriteStandardNarrow("Stage", "Seed");
                    break;
                }
                case GrainStages::Nucleus:
                {
                    ConsoleOutput::WriteStandardNarrow("Stage", "Nucleus");
                    break;
                }
                case GrainStages::Stable:
                {
                    ConsoleOutput::WriteStandardNarrow("Stage", "Stable");
                    break;
                }
            }
            ConsoleOutput::WriteStandardNarrow("Volume",  FieldsProperties[idx].Volume);
            ConsoleOutput::WriteLine("-");
        }
    }
    ConsoleOutput::WriteLine("=");
}

void PhaseField::PrintVolumeFractions()
{
    ConsoleOutput::WriteSimple("Phase fractions:");
    for(size_t n = 0; n < Nphases; n++)
    {
        ConsoleOutput::WriteStandardNarrow("Phase " + std::to_string(n) + " (" + PhaseNames[n] + ") [%]", 100.0*FractionsTotal[n]);
    }
}

void PhaseField::WriteVolumePercentages(const std::string& filename, double time, char separator) const
{
    /** This function will create tabulated data on the volume percent of each
    thermodynamic phase. Each time this function is called, a new row will be
    written in the specified file name for the current time step. If the file is
    not present it will be created, if the phase already exists, the new data
    will be appended!*/

    // NOTE USED BY OPSTUDIO GUI !!!!!!

#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            // Write data header
            std::ofstream file(filename, ios::out);
            file << "time";
            for(size_t idx = 0; idx < Nphases; idx++)
            {
                file << separator << PhaseNames[idx];
            }
            file << "\n";
            file.close();
        }

        // Write data
        std::ofstream file(filename, ios::app);
        file << std::scientific << time;
        for(size_t idx = 0; idx < Nphases; idx++)
        {
            file << std::scientific << separator << 100.0*FractionsTotal[idx];
        }
        file << "\n";
        file.close();
    }
}

void PhaseField::WriteGrainsVolume(int time_step, double time, std::string filename)
{
    /** This function will write a list of all grain volumes in a file, each
    time this function is called in a new row. The first column is the time,
    the second row beginning with a "#" marks the thermodynamic index for each
    grain.*/

string separator = " ";

    if (!std::filesystem::exists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# Volume of each grain in the order of their field-index! "
             << "First row contains the simulation time. "
             << "The next line contains the thermodynamic phase field index "
             << "of each grain for identification\n";

        file << "TimeStep " << separator << "Time";
        for(size_t n = 0; n < FieldsProperties.size(); n++)
        {
            file << separator << FieldsProperties[n].Phase;
        }
        file << "\n";
        file.close();
    }

    ofstream file(filename, ios::app);
    file.precision(10);
    file << time_step << separator << time << separator;

    for (size_t i = 0; i < FieldsProperties.size(); i++)
    {
        file << FieldsProperties[i].Volume << separator;
    }
    file << "\n";
    file.close();
}

void PhaseField::Remesh(int newNx, int newNy, int newNz,
                        const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            Fields.Remesh(Grid.Nx, Grid.Ny, Grid.Nz);
            FieldsDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
            {
                Fields(i,j,k).flag = 2*(Fields(i,j,k).size() > 1);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            Fractions.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            break;
        }
        case Resolutions::Dual:
        {
            Fields.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            FieldsDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

            FieldsDR.Remesh(2*Grid.Nx, 2*Grid.Ny, 2*Grid.Nz);
            FieldsDotDR.Reallocate(2*Grid.Nx, 2*Grid.Ny, 2*Grid.Nz);

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,0,)
            {
                FieldsDR(i,j,k).flag = 2*(FieldsDR(i,j,k).size() > 1);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            Fractions.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            break;
        }
    }

    Finalize(BC);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

size_t PhaseField::PlantGrainNucleus(size_t PhaseIndex, int x, int y, int z)
{
    size_t locIndex = FieldsProperties.add_grain(PhaseIndex);

    FieldsProperties[locIndex].Stage  = GrainStages::Seed;

    FieldsProperties[locIndex].Rcm[0] = x;
    FieldsProperties[locIndex].Rcm[1] = y;
    FieldsProperties[locIndex].Rcm[2] = z;
    FieldsProperties[locIndex].RefVolume = CalculateReferenceVolume(Grid.iWidth);
    FieldsProperties[locIndex].VolumeRatio = 0.0;
    FieldsProperties[locIndex].State = PhaseAggregateStates[PhaseIndex];
    FieldsProperties[locIndex].Density = PhaseEquilibriumDensities[PhaseIndex];

    NucleationPresent = true;

    if (x - Grid.OffsetX >= 0 && x - Grid.OffsetX < Grid.Nx and
        y - Grid.OffsetY >= 0 && y - Grid.OffsetY < Grid.Ny and
        z - Grid.OffsetZ >= 0 && z - Grid.OffsetZ < Grid.Nz)
    {
        Fields(x - Grid.OffsetX, y - Grid.OffsetY, z - Grid.OffsetZ).set_value(locIndex, 0.0);

        long int fx = 1 + Grid.dNx;
        long int fy = 1 + Grid.dNy;
        long int fz = 1 + Grid.dNz;

        if(Grid.Resolution == Resolutions::Dual)
        {
            for(int di = -Grid.dNx; di <= Grid.dNx; di+=2)
            for(int dj = -Grid.dNy; dj <= Grid.dNy; dj+=2)
            for(int dk = -Grid.dNz; dk <= Grid.dNz; dk+=2)
            {
                FieldsDR(fx*(x - Grid.OffsetX)+(di+1)/2,fy*(y - Grid.OffsetY)+(dj+1)/2,fz*(z - Grid.OffsetZ)+(dk+1)/2).set_value(locIndex, 0.0);
            }
        }
    }

    return locIndex;
}

size_t PhaseField::AddGrainInfo(size_t PhaseIndex)
{
    size_t locIndex = FieldsProperties.add_grain(PhaseIndex);
    FieldsProperties[locIndex].RefVolume = CalculateReferenceVolume(Grid.iWidth);
    FieldsProperties[locIndex].VolumeRatio = 0.0;
    FieldsProperties[locIndex].State = PhaseAggregateStates[PhaseIndex];
    FieldsProperties[locIndex].Density = PhaseEquilibriumDensities[PhaseIndex];

    return locIndex;
}

PhaseField& PhaseField::operator= (const PhaseField& rhs)
{
    // protect against self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "PhaseField")
    {
        thisclassname = rhs.thisclassname;
        NucleationPresent = rhs.NucleationPresent;

        Grid = rhs.Grid;

        Nphases = rhs.Nphases;
        PhaseNames = rhs.PhaseNames;
        PhaseAggregateStates = rhs.PhaseAggregateStates;
        PhaseEquilibriumDensities = rhs.PhaseEquilibriumDensities;

        FractionsTotal = rhs.FractionsTotal;

        ConsiderNucleusVolume = rhs.ConsiderNucleusVolume;

        PhaseFieldLaplacianStencil = rhs.PhaseFieldLaplacianStencil;
        PhaseFieldGradientStencil = rhs.PhaseFieldGradientStencil;

        LStencil = rhs.LStencil;
        GStencil = rhs.GStencil;

        FieldsProperties = rhs.FieldsProperties;

        PairwiseGrowthFactors = rhs.PairwiseGrowthFactors;

        if (Fields.IsNotAllocated())
        {
            Fields.Allocate(Grid, rhs.Fields.Bcells());
            FieldsDot.Allocate(Grid, rhs.FieldsDot.Bcells());
            Fractions.Allocate(Grid, {Nphases},rhs.Fractions.Bcells());
        }
        else if (not Fields.IsSize(Grid.Nx, Grid.Ny, Grid.Nz))
        {
            Fields.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            FieldsDot.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
            Fractions.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        {
            Fields(i,j,k) = rhs.Fields(i,j,k);
            FieldsDot(i,j,k) = rhs.FieldsDot(i,j,k);

        }
        OMP_PARALLEL_STORAGE_LOOP_END
        CalculateFractions();

        if(Grid.Resolution == Resolutions::Dual)
        {
            GridParameters DoubleDimensions = Grid.DoubleResolution();

            if (FieldsDR.IsNotAllocated())
            {
                FieldsDR.Allocate(DoubleDimensions, rhs.Fields.Bcells());
                //FieldsDerivativesDR.Allocate(DoubleDimensions, rhs.Fields.Bcells());
                FieldsDotDR.Allocate(DoubleDimensions, rhs.FieldsDot.Bcells());
            }
            else if (not FieldsDR.IsSize(DoubleDimensions.Nx, DoubleDimensions.Ny, DoubleDimensions.Nz))
            {
                FieldsDR.Reallocate(DoubleDimensions.Nx, DoubleDimensions.Ny, DoubleDimensions.Nz);
                //FieldsDerivativesDR.Reallocate(DoubleDimensions.Nx, DoubleDimensions.Ny, DoubleDimensions.Nz);
                FieldsDotDR.Reallocate(DoubleDimensions.Nx, DoubleDimensions.Ny, DoubleDimensions.Nz);
            }

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
            {
                FieldsDR(i,j,k) = rhs.FieldsDR(i,j,k);
                FieldsDotDR(i,j,k) = rhs.FieldsDotDR(i,j,k);
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            //CopyDerivativesDR();
        }
    }
    return *this;
}

std::vector<int> PhaseField::GetPresentPhaseFields() const
{
    std::vector<int> indeces;
    for(unsigned int n = 0; n < FieldsProperties.size(); n++)
    if(FieldsProperties[n].Exist)
    {
        indeces.push_back(n);
    }

    return indeces;
}

std::vector<int> PhaseField::VicinityPhaseFields(const int i,
                                           const int j, const int k) const
{
    // Analyze vicinity
    std::vector<int> tempPFindex;
    tempPFindex.push_back (Fields(i,j,k).front().index);
    for(int ii = -Grid.dNx; ii <= Grid.dNx; ii++)
    for(int jj = -Grid.dNy; jj <= Grid.dNy; jj++)
    for(int kk = -Grid.dNz; kk <= Grid.dNz; kk++)
    {
        for(auto alpha  = Fields(i+ii,j+jj,k+kk).cbegin();
                 alpha != Fields(i+ii,j+jj,k+kk).cend(); ++alpha)
        {
            if(alpha->index != Fields(i,j,k).front().index)
            {
                tempPFindex.push_back (alpha->index);
            }
        }
    }

    // Erase double entries
    std::sort(tempPFindex.begin(), tempPFindex.end() );
    tempPFindex.erase( std::unique( tempPFindex.begin(),
                tempPFindex.end() ), tempPFindex.end() );

    return tempPFindex;
}

void PhaseField::CombinePhaseFields(void)
{
    for(size_t p = 0; p < Nphases; p++)
    if(Combine[p])
    {
        size_t MajorityIndex  = 0;
        double MajorityVolume = 0.0;
        bool   MajorityFound  = false;
        int    Nfields        = 0;
        for(size_t n = 0; n < FieldsProperties.size(); n++)
        {
            if(FieldsProperties[n].Exist and
               FieldsProperties[n].Phase == p and
               FieldsProperties[n].Stage == GrainStages::Stable)
            {
                Nfields ++;

                if(FieldsProperties[n].Volume > MajorityVolume)
                {
                    MajorityIndex  = n;
                    MajorityVolume = FieldsProperties[n].Volume;
                    MajorityFound  = true;
                }
            }
        }

        if(MajorityFound and Nfields > 1)
        switch(Grid.Resolution)
        {
            case Resolutions::Single:
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
                {
                    double value_to_combine = 0.0;

                    for (auto it  = Fields(i,j,k).begin();
                              it != Fields(i,j,k).end(); ++it)
                    {
                        if(FieldsProperties[it->index].Phase == p and
                           FieldsProperties[it->index].Stage == GrainStages::Stable and
                           it->index != MajorityIndex)
                        {
                            value_to_combine += it->value;
                            it->value = 0.0;
                        }
                    }

                    if(value_to_combine > DBL_EPSILON)
                    {
                        Fields(i,j,k).add_value(MajorityIndex, value_to_combine);
                    }
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                break;
            }
            case Resolutions::Dual:
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
                {
                    double value_to_combine = 0.0;

                    for (auto it  = FieldsDR(i,j,k).begin();
                              it != FieldsDR(i,j,k).end(); ++it)
                    {
                        if(FieldsProperties[it->index].Phase == p and
                           FieldsProperties[it->index].Stage == GrainStages::Stable and
                           it->index != MajorityIndex)
                        {
                            value_to_combine += it->value;
                            it->value = 0.0;
                        }
                    }

                    if(value_to_combine > DBL_EPSILON)
                    {
                        FieldsDR(i,j,k).add_value(MajorityIndex, value_to_combine);
                    }
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                break;
            }
        }
    }
}

void PhaseField::SelectiveCombinePhaseFields(const BoundaryConditions& BC,
                                             const size_t TargetPFIndex,
                                             const size_t SourcePFIndex)
{
    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            SelectiveCombinePhaseFieldsSR(BC,TargetPFIndex,SourcePFIndex);
            break;
        }
        case Resolutions::Dual:
        {
            SelectiveCombinePhaseFieldsDR(BC,TargetPFIndex,SourcePFIndex);
            break;
        }
    }
}

void PhaseField::SelectiveCombinePhaseFieldsSR(const BoundaryConditions& BC,
                                               const size_t TargetPFIndex,
                                               const size_t SourcePFIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        double loc_value = 0.0;
        bool ChangedPF = false;

        for (auto it  = Fields(i,j,k).begin();
                  it != Fields(i,j,k).end(); ++it)
        {
            if(it->index == SourcePFIndex)
            {
                loc_value += it->value;
                it->value = 0.0;
                ChangedPF = true;
            }
        }
        if(ChangedPF == true)
        {
            Fields(i,j,k).add_value(TargetPFIndex, loc_value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

void PhaseField::SelectiveCombinePhaseFieldsDR(const BoundaryConditions& BC,
                                               const size_t TargetPFIndex,
                                               const size_t SourcePFIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDR,FieldsDR.Bcells(),)
    {
        double loc_value = 0.0;
        bool ChangedPF = false;

        for (auto it  = FieldsDR(i,j,k).begin();
                  it != FieldsDR(i,j,k).end(); ++it)
        {
            if(it->index == SourcePFIndex)
            {
                loc_value += it->value;
                it->value = 0.0;
                ChangedPF = true;
            }
        }
        if(ChangedPF == true)
        {
            FieldsDR(i,j,k).add_value(TargetPFIndex, loc_value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

pair<NodeA<double>, NodeA<double>> PhaseField::PrincipalCurvatures(const int i,
                                                   const int j,
                                                   const int k) const
{
    NodeA<double> kappa1;
    NodeA<double> kappa2;

    kappa1.clear();
    kappa2.clear();

    const double GWeights[3] = {-0.5/Grid.dx, 0.0, 0.5/Grid.dx};

    // Check if neighbour cell are also in the interface
    if (Fields(i,j,k).interface())
    {
        // Calculate gradients of phase normal fields (Jacobian matrix)
        Tensor<NodeA<double>, 2 > NormalGradient;
        NormalGradient.Allocate({3,3});
        NormalGradient.set_to_zero();

        // Calculate gradients of normals
        for (int ii = -1; ii <= +1; ii += 2)
        {
            const double GWeight = GWeights[ii+1];

            if (Grid.dNx)
            {
                NodeA<dVector3> locNormalsX = NormalsPhase(i+ii,j,k);
                for (auto it = locNormalsX.cbegin(); it < locNormalsX.cend(); ++it)
                {
                    NormalGradient({0,0}).add_value(it->index, GWeight*(it->value[0]));
                    NormalGradient({1,0}).add_value(it->index, GWeight*(it->value[1]));
                    NormalGradient({2,0}).add_value(it->index, GWeight*(it->value[2]));
                }
            }
            if (Grid.dNy)
            {
                NodeA<dVector3> locNormalsY = NormalsPhase(i,j+ii,k);
                for (auto it = locNormalsY.cbegin(); it < locNormalsY.cend(); ++it)
                {
                    NormalGradient({0,1}).add_value(it->index, GWeight*(it->value[0]));
                    NormalGradient({1,1}).add_value(it->index, GWeight*(it->value[1]));
                    NormalGradient({2,1}).add_value(it->index, GWeight*(it->value[2]));
                }
            }
            if (Grid.dNz)
            {
                NodeA<dVector3> locNormalsZ = NormalsPhase(i,j,k+ii);
                for (auto it = locNormalsZ.cbegin(); it < locNormalsZ.cend(); ++it)
                {
                    NormalGradient({0,2}).add_value(it->index, GWeight*(it->value[0]));
                    NormalGradient({1,2}).add_value(it->index, GWeight*(it->value[1]));
                    NormalGradient({2,2}).add_value(it->index, GWeight*(it->value[2]));
                }
            }
        }

        // Calculate basis of tangent space at (i,j,k)
        // the basis of a spherical coordinate system at (i,j,k) will be
        // used here. The phase normal vector is in this case equal to the
        // radial normal vector of the spherical coordinate and the
        // remaining basis vectors will be the basis of the tangent space

        NodeA<dVector3> locNormals = NormalsPhase(i,j,k);
        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            double phi   = atan2(it->value[1],it->value[0]);
            double theta = acos(it->value[2]);

            // Calculate basis vectors of tangent space
            dVector3 e_n;
            e_n.set_to_zero();
            e_n[0] = it->value[0];
            e_n[1] = it->value[1];
            e_n[2] = it->value[2];

            dVector3 e_theta;
            e_theta.set_to_zero();
            e_theta[0] = cos(theta) * cos(phi);
            e_theta[1] = cos(theta) * sin(phi);
            e_theta[2] =            - sin(theta);

            dVector3 e_phi;
            e_phi.set_to_zero();
            e_phi[0] = - sin(phi);
            e_phi[1] =   cos(phi);

            // Calculate projection and inclusion matrices
            double Projection [2][3];
            double Inclusion  [3][2];
            for (size_t m = 0; m < 3; ++m)
            {
                Projection[1][m] = e_phi  [m];
                Projection[0][m] = e_theta[m];

                Inclusion[m][1] = Projection[1][m];
                Inclusion[m][0] = Projection[0][m];
            }

            // Calculate local Weingarten map W
            complex<double> W[2][2];
            for (size_t l = 0; l < 2; ++l)
            for (size_t m = 0; m < 2; ++m)
            {
                W[l][m] = 0.0;
                for (size_t n = 0; n < 3; ++n)
                for (size_t o = 0; o < 3; ++o)
                {
                    W[l][m] -= Projection[l][n] *
                        NormalGradient({n,o}).get_value(it->index)
                        * Inclusion[o][m];
                }
            }

            // Calculate eigenvalues of local Weingarten map
            const double part1 = real((W[0][0] + W[1][1])/2.0);
            const double W2    = real((W[0][0] - W[1][1]) * (W[0][0] - W[1][1]));
            const double part2 = real(0.5 * sqrt(W2) + 4.0*W[0][1] * W[1][0]);

            const double locKappa1 = real(part1 - part2);
            const double locKappa2 = real(part1 + part2);

            // Add eigenvalues to Node storage
            kappa1.add_value(it->index, locKappa1);
            kappa2.add_value(it->index, locKappa2);
        }
    }

    return make_pair(kappa1, kappa2);
}

std::array<double,2> PhaseField::PrincipalCurvatures(const int i, const int j, const int k, const size_t phase) const
{
    double phaseKappa1 = 0; // principle curvature phase
    double phaseKappa2 = 0; // principle curvature phase

    if (Fields(i,j,k).interface())
    {
          NodeA<double> locKappa1; // principle curvature of each phase field
          NodeA<double> locKappa2; // principle curvature of each phase field
          // Calculate curvature of all phase fields
          tie(locKappa1, locKappa2) = PrincipalCurvatures(i,j,k);

          // Calculate curvature of desired phase
          for (auto it  = locKappa1.cbegin();
                    it != locKappa1.cend(); it++)
          {
              size_t pIndex = FieldsProperties[it->index].Phase;
              if (pIndex == phase) phaseKappa1 += it->value;
          }
          for (auto it  = locKappa2.cbegin();
                    it != locKappa2.cend(); it++)
          {
              size_t pIndex = FieldsProperties[it->index].Phase;
              if (pIndex == phase) phaseKappa2 += it->value;
          }
    }

    std::array<double,2> result({phaseKappa1,phaseKappa2});
    return result;
}

double PhaseField::CurvaturePhase(const int i, const int j, const int k, const size_t phase) const
{
    double locKappa = 0.0;
    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        const int ii = gs->di;
        const int jj = gs->dj;
        const int kk = gs->dk;
        const NodeA<dVector3> locNormals = NormalsPhase(i + ii, j + jj, k + kk);
        for (auto it = locNormals.cbegin(); it != locNormals.cend(); ++it)
        if (it->index == phase)
        {
            const double vx = gs->weightX * it->value[0];
            const double vy = gs->weightY * it->value[1];
            const double vz = gs->weightZ * it->value[2];
            locKappa += vx+vy+vz;
        }
    }
    return locKappa;
}

double PhaseField::CurvaturePhaseDR(const int i, const int j, const int k, const size_t phase) const
{
    double locKappa = 0.0;
    for (auto gs = GStencil.cbegin(); gs != GStencil.cend(); ++gs)
    {
        const int ii = gs->di;
        const int jj = gs->dj;
        const int kk = gs->dk;
        const NodeA<dVector3> locNormals = NormalsPhaseDR(i + ii, j + jj, k + kk);
        for (auto it = locNormals.cbegin(); it != locNormals.cend(); ++it)
        if (it->index == phase)
        {
            const double vx = gs->weightX * it->value[0];
            const double vy = gs->weightY * it->value[1];
            const double vz = gs->weightZ * it->value[2];
            locKappa += vx+vy+vz;
        }
    }
    return locKappa;
}

bool PhaseField::PhaseFieldPresent(const int i, const int j, const int k,
                                                        const size_t Index) const
{
    if(!Fields(i,j,k).interface())
    {
        if (Fields(i,j,k).front().index == Index) return true;
    }
    else
    {
        for (auto alpha  = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend(); ++alpha)
        {
            if (alpha->index == Index) return true;
        }
    }
    return false;
}

bool PhaseField::ThermodynamicPhasePresent(size_t alpha)
{
    bool result = false;

    for(size_t beta = 0; beta < FieldsProperties.size(); beta++)
    {
        if(FieldsProperties[beta].Phase == alpha) result = true;
    }

    return result;
}

bool PhaseField::ThermodynamicPhasePairPresent(size_t alpha, size_t beta)
{
    bool result = false;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Fields,0,reduction(or:result))
    if (Fields(x,y,z).interface())
    {
        bool alphapresent = false;
        bool betapresent  = false;
        for (auto it  = Fields(x,y,z).cbegin();
                  it != Fields(x,y,z).cend(); ++it)
        {
            if(it->index == alpha)
            {
                alphapresent = true;
            }

            if(it->index == beta)
            {
                betapresent = true;
            }
        }
        if((alphapresent) and (betapresent)) result = true;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return result;
}

std::vector<int> PhaseField::GetMaxPhaseFieldOverlap(const size_t thPhase1,
        const size_t thPhase2)
{
    std::vector<int> maxOverlap;
    Matrix<int> Overlap;
    size_t numberOfGrains = FieldsProperties.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Fields(i,j,k).interface())
        {
            for (auto it  = Fields(i,j,k).cbegin();
                      it != Fields(i,j,k).cend(); ++it)
            {
                size_t locIndex1 = it->index;
                size_t thPhaseIndex1 = FieldsProperties[locIndex1].Phase;
                if (thPhaseIndex1 == thPhase1)
                {
                    for (auto jt = it+1; jt != Fields(i,j,k).cend(); ++jt)
                    {
                        size_t locIndex2 = jt->index;
                        size_t thPhaseIndex2 = FieldsProperties[locIndex2].Phase;
                        if (thPhaseIndex2 == thPhase2)
                        {
                            Overlap.add(locIndex1, locIndex2, 1);
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    int maxNx = -1;
    int maxNy = -1;
    for (size_t it = 0; it < numberOfGrains-1; it++)
    for (size_t jt = it+1; jt < numberOfGrains; jt++)
    {
        if (maxNx < Overlap.get(it, jt))
        {
            maxNx = it;
            maxNy = jt;
        }
    }

    if (maxNx == -1 or maxNy == -1)
    {
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        return maxOverlap;
    }
    maxOverlap.push_back(maxNx);
    maxOverlap.push_back(maxNy);
    maxOverlap.push_back(Overlap.get(maxNx, maxNy));
    return maxOverlap;
}

Matrix<int> PhaseField::GetPhaseFieldOverlap(const size_t thPhase1,
        const size_t thPhase2)
{
    Matrix<int> Overlap;
    int numberOfGrains = FieldsProperties.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Fields(i,j,k).interface())
        for (auto it  = Fields(i,j,k).cbegin();
                  it != Fields(i,j,k).cend(); ++it)
        {
            size_t phaseIndex = it->index;
            size_t thPhaseIndex = FieldsProperties[phaseIndex].Phase;
            if (thPhaseIndex == thPhase1)
            {
                for (auto jt = it+1;
                          jt < Fields(i,j,k).cend(); ++jt)
                {
                    size_t phaseIndex2 = jt->index;
                    size_t thPhaseIndex2 = FieldsProperties[phaseIndex2].Phase;

                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    if (thPhaseIndex2 == thPhase2)
                    {
                        Overlap.add(phaseIndex, phaseIndex2, 1);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Overlap;
}

void PhaseField::Advect(AdvectionHR& Adv, const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt, const double tStep,
        const bool finalize)
{
    if(FieldsAdvectionDot   .IsNotAllocated()) FieldsAdvectionDot   .Allocate(Grid, Fields.Bcells());
    if(FieldsAdvectionBackup.IsNotAllocated()) FieldsAdvectionBackup.Allocate(Grid, Fields.Bcells());

    if(not FieldsAdvectionDot   .IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionDot   .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    if(not FieldsAdvectionBackup.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionBackup.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    Adv.AdvectPhaseField(*this, Vel, BC, dt, tStep, finalize);
    CalculateFractions();
    //ConsoleOutput::WriteStandard(thisclassname, "Advected");
}

void PhaseField::AdvectALE(AdvectionHR& Adv, const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt, const double tStep,
        const bool finalize)
{
    if(FieldsAdvectionDot   .IsNotAllocated()) FieldsAdvectionDot   .Allocate(Grid, Fields.Bcells());
    if(FieldsAdvectionBackup.IsNotAllocated()) FieldsAdvectionBackup.Allocate(Grid, Fields.Bcells());

    if(not FieldsAdvectionDot   .IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionDot   .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    if(not FieldsAdvectionBackup.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionBackup.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    Adv.AdvectPhaseFieldALE(*this, Vel, BC, dt, tStep, finalize);
    Finalize(BC);
    ConsoleOutput::WriteStandard(thisclassname, "Advected");
}

void PhaseField::Advect(AdvectionHR& Adv, const Velocities& Vel,
        const BoundaryConditions& BC, FlowSolverLBM& LBM,
        const double dt, const double tStep,
        const bool finalize)
{
    if(FieldsAdvectionDot   .IsNotAllocated()) FieldsAdvectionDot   .Allocate(Grid, Fields.Bcells());
    if(FieldsAdvectionBackup.IsNotAllocated()) FieldsAdvectionBackup.Allocate(Grid, Fields.Bcells());

    if(not FieldsAdvectionDot   .IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionDot   .Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);
    if(not FieldsAdvectionBackup.IsSize(Grid.Nx, Grid.Ny, Grid.Nz)) FieldsAdvectionBackup.Reallocate(Grid.Nx, Grid.Ny, Grid.Nz);

    LBM.DetectObstacles(*this);
    Adv.AdvectPhaseField(*this, Vel, BC, dt, tStep, finalize);
    LBM.DetectObstaclesAdvection(*this,Vel,BC);
    CalculateFractions();
    //ConsoleOutput::WriteStandard(thisclassname, "Advected");
}

Tensor<double,1> PhaseField::NewFractions(const int x, const int y, const int z, double locdt) const
{
    Tensor<double,1> newFractions = Fractions(x,y,z);
    for(auto it = FieldsDot(x,y,z).cbegin();
             it != FieldsDot(x,y,z).cend(); ++it)
    {
        size_t pIndexA = FieldsProperties[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsProperties[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta

        if(pIndexA != pIndexB)
        {
            double factor = 1.0;
            if(FieldsProperties[it->indexA].GrowthConstraintsViolation != GrowthConstraintsViolations::None and
               FieldsProperties[it->indexB].GrowthConstraintsViolation != GrowthConstraintsViolations::None)
            {
                factor *= PairwiseGrowthFactors.get_sym1(it->indexA,it->indexB);
            }
        
            newFractions({pIndexA}) += factor*it->value1*locdt;
            newFractions({pIndexB}) -= factor*it->value1*locdt;
        }
    }
    return newFractions;
}

Tensor<double,2> PhaseField::FractionsPsi(const int x, const int y, const int z, double locdt) const
{
    Tensor<double,2> Psi({Nphases,Nphases});
    Psi.set_to_zero();
    for(auto it = FieldsDot(x,y,z).cbegin();
             it != FieldsDot(x,y,z).cend(); ++it)
    {
        size_t pIndexA = FieldsProperties[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsProperties[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta

        if(pIndexA != pIndexB)
        {
            Psi({pIndexA,pIndexB}) += it->value1*locdt;
            Psi({pIndexB,pIndexA}) -= it->value1*locdt;
        }
    }
    return Psi;
}

Tensor<double,1> PhaseField::CalculateNewFractions(Tensor<double,1> oldFractions,
                                                   NodeAB<double,double>& FieldsDot,
                                                   double locdt)
{
    Tensor<double,1> newFractions = oldFractions;
    for(auto it = FieldsDot.begin(); it != FieldsDot.end(); ++it)
    {
        size_t pIndexA = FieldsProperties[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsProperties[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta
        if(pIndexA != pIndexB)
        {
            newFractions({pIndexA}) += it->value1*locdt;
            newFractions({pIndexB}) -= it->value1*locdt;
        }
    }
    return newFractions;
}

Tensor<double,2> PhaseField::CalculatePsi(NodeAB<double,double>& FieldsDot, double locdt)
{
    Tensor<double,2> Psi({Nphases,Nphases});
    Psi.set_to_zero();
    for(auto it = FieldsDot.begin(); it != FieldsDot.end(); ++it)
    {
        size_t pIndexA = FieldsProperties[it->indexA].Phase;                    // Index of thermodynamic phase alpha in the phase pair alpha-beta
        size_t pIndexB = FieldsProperties[it->indexB].Phase;                    // Index of thermodynamic phase beta in the phase pair alpha-beta
        if(pIndexA != pIndexB)
        {
            Psi({pIndexA,pIndexB}) += it->value1*locdt;
            Psi({pIndexB,pIndexA}) -= it->value1*locdt;
        }
    }
    return Psi;
}

Storage3D<double,1> PhaseField::NewFractions(double dt)
{
    Storage3D<double,1> result(Grid,{Nphases},0);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,result,0,)
    if(Fields(x,y,z).wide_interface())
    {
        result(x,y,z) = CalculateNewFractions(Fractions(x,y,z),
                                              FieldsDot(x,y,z),dt);
    }
    else
    {
        result(x,y,z) = Fractions(x,y,z);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return result;
}

double PhaseField::Interfaces(const int i, const int j, const int k) const
{
    double sum = 0.5;
    if(Fields(i,j,k).interface())
    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    for (auto  beta  = alpha + 1;
               beta != Fields(i,j,k).cend(); ++beta)
    {
        sum -= alpha->value*beta->value;
    }
    return 1.0/(sum*2.0);
}

double PhaseField::InterfacesDR(const int i, const int j, const int k) const
{
    double sum = 0.5;
    if(FieldsDR(i,j,k).interface())
    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    for (auto  beta  = alpha + 1;
               beta != FieldsDR(i,j,k).cend(); ++beta)
    {
        sum -= alpha->value*beta->value;
    }
    return 1.0/(sum*2.0);
}

double PhaseField::Variants(const int i, const int j, const int k) const
{
    int locVariant = 0;
    double locValue = 0.0;
    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locVariant = FieldsProperties[alpha->index].Variant;
        }
    }
    return locVariant;
}

double PhaseField::VariantsDR(const int i, const int j, const int k) const
{
    int locVariant = 0;
    double locValue = 0.0;
    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            locVariant = FieldsProperties[alpha->index].Variant;
        }
    }
    return locVariant;
}

double PhaseField::ParentGrain(const int i, const int j, const int k) const
{
    int ParentGrain = 0;
    double locValue = 0.0;
    for (auto alpha  = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            ParentGrain = FieldsProperties[alpha->index].Parent;
        }
    }
    return ParentGrain;
}

double PhaseField::ParentGrainDR(const int i, const int j, const int k) const
{
    int ParentGrain = 0;
    double locValue = 0.0;
    for (auto alpha  = FieldsDR(i,j,k).cbegin();
              alpha != FieldsDR(i,j,k).cend(); ++alpha)
    {
        if(alpha->value > locValue)
        {
            locValue = alpha->value;
            ParentGrain = FieldsProperties[alpha->index].Parent;
        }
    }
    return ParentGrain;
}

double PhaseField::SolidPhaseFraction(const int i, const int j, const int k) const
{
    double SolidPhaseFraction = 0.0;
    for (auto alpha = Fields(i,j,k).cbegin(); alpha != Fields(i,j,k).cend(); alpha++)
    {
       if (FieldsProperties[alpha->index].is_solid())
       {
           SolidPhaseFraction += alpha->value;
       }
    }
    return SolidPhaseFraction;
}

double PhaseField::SolidMassDensity(const int i, const int j, const int k) const
{
    double SolidMassDensity = 0.0;
    for (auto alpha = Fields(i,j,k).cbegin(); alpha != Fields(i,j,k).cend(); alpha++)
    if (FieldsProperties[alpha->index].is_solid())
    {
        SolidMassDensity += alpha->value*FieldsProperties[alpha->value].Density;
    }
    return SolidMassDensity;
}

double PhaseField::CalculateSolidMass() const
{
    double locSolidMass = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,reduction(+:locSolidMass))
    {
        locSolidMass += SolidMassDensity(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    double SolidMass = locSolidMass;
    #ifdef MPI_PARALLEL
    OP_MPI_Allreduce(&locSolidMass, &(SolidMass), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    #endif

    SolidMass *= Grid.CellVolume();
    return SolidMass;
}
}// namespace openphase
