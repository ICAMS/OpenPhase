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

 *   File created :
 *   Main contributors :
 *
 */

#include "Includes.h"
#include "Stabilization.h"
#include "GrainsProperties.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "VTK.h"

namespace openphase
{

class BoundaryConditions;
class InterfaceProperties;
class PhaseField;
class Settings;
class DrivingForce;

class StencilDirection
{
 public:
    double alpha;
    double beta;
    double scal_prod;
    double phi_tl;
    double phi_tu;

    double stencil_weight;

    int d_x, d_y, d_z;
};

class StabilizationImpl : public OPObject
{
 public:
    StabilizationImpl();
    StabilizationImpl(Settings& locSettings);

    void Initialize(Settings& locSettings);
    void Initialize(Settings& locSettings, LaplacianStencil _DStencil);
    void ReadInput(const std::string FileName) override;
    void ReadInput(std::stringstream& inp) override;
    void CalculateCurvature(PhaseField& Phase);
    void CalculateCurvature_SPF(PhaseField& Phase);
    void CalculateAverageCurvature(PhaseField& Phase);
    void CalculateAverageCurvatureRegularization(PhaseField& Phase);
    void CalculateAverageCurvatureNormal(PhaseField& Phase);
    void CalculateAverageCurvatureNormalSnap(PhaseField& Phase);
    void CalculateAverageCurvatureSimple(PhaseField& Phase);
    void ComputeVelocity(PhaseField& Phase,
                                InterfaceProperties& IP, DrivingForce& dG);
    void ComputeUpdate(PhaseField& Phase,
                                DrivingForce& dG, double dt);
    double GetMaxDT();
    void Smooth(PhaseField& Phase);
    void WriteVTK(const Settings& locSettings, const int tStep,
                            const size_t indexA, const size_t indexB,
                            const int precision) const;
    Storage3D<NodePF, 0> AverageCurvature;
    Storage3D<NodePF, 0> AverageCurvatureWeight;
    Storage3D<NodePF, 0> AverageCurvatureAux;
    Storage3D<NodePF, 0> AverageCurvatureAuxWeight;
    Storage3D<NodePF, 0> Curvature;
    Storage3D<NodePF, 0> AverageInterface;
    Storage3D<NodePF, 0> AverageInterfaceWeight;
    Storage3D<NodePF, 0> AverageInterfaceAux;
    Storage3D<NodePF, 0> AverageInterfaceAuxWeight;
    Storage3D<NodePF, 0> Interface;
    Storage3D<NodeAB<double,double>, 0> Velocity;
    double maxvelocity;
    double maxdt;
    GridParameters Grid;                                                        ///< Simulation grid parameters

 protected:
 private:
    LaplacianStencil DStencil;
    double exponent;
    double range;
    double STAB;
    std::vector<double> n_vector;
    std::vector<StencilDirection>  StencilDirections;
};

Stabilization::Stabilization() : impl_(new StabilizationImpl) {}

Stabilization::Stabilization(Settings& locSettings) : impl_(new StabilizationImpl) {Initialize(locSettings);}

Stabilization::~Stabilization() = default;

void Stabilization::Initialize(Settings& locSettings)
{
    impl_->Initialize(locSettings);
}
void Stabilization::Initialize(Settings& locSettings, LaplacianStencil _DStencil)
{
    impl_->Initialize(locSettings, _DStencil);
}

void Stabilization::ReadInput(const std::string FileName)
{
    impl_->ReadInput(FileName);
}
void Stabilization::ReadInput(std::stringstream& inp)
{
    impl_->ReadInput(inp);
}
void Stabilization::CalculateCurvature(PhaseField& Phase)
{
    impl_->CalculateCurvature(Phase);
}
void Stabilization::CalculateCurvature_SPF(PhaseField& Phase)
{
    impl_->CalculateCurvature_SPF(Phase);
}
void Stabilization::CalculateAverageCurvature(PhaseField& Phase)
{
    impl_->CalculateAverageCurvature(Phase);
}
void Stabilization::CalculateAverageCurvatureRegularization(PhaseField& Phase)
{
    impl_->CalculateAverageCurvatureRegularization(Phase);
}
void Stabilization::CalculateAverageCurvatureNormal(PhaseField& Phase)
{
    impl_->CalculateAverageCurvatureNormal(Phase);
}
void Stabilization::CalculateAverageCurvatureNormalSnap(PhaseField& Phase)
{
    impl_->CalculateAverageCurvatureNormalSnap(Phase);
}
void Stabilization::CalculateAverageCurvatureSimple(PhaseField& Phase)
{
    impl_->CalculateAverageCurvatureSimple(Phase);
}
void Stabilization::ComputeVelocity(PhaseField& Phase,
                            InterfaceProperties& IP, DrivingForce& dG)
{
    impl_->ComputeVelocity(Phase, IP, dG);
}
void Stabilization::ComputeUpdate(PhaseField& Phase,
                            DrivingForce& dG, double dt)
{
    impl_->ComputeUpdate(Phase, dG, dt);
}
double Stabilization::GetMaxDT()
{
    return impl_->GetMaxDT();
}
void Stabilization::Smooth(PhaseField& Phase)
{
    impl_->Smooth(Phase);
}
void Stabilization::WriteVTK(const Settings& locSettings, const int tStep,
                            const size_t indexA, const size_t indexB,
                            const int precision) const
{
    impl_->WriteVTK(locSettings, tStep, indexA, indexB, precision);
}

StabilizationImpl::StabilizationImpl()
{

}

StabilizationImpl::StabilizationImpl(Settings& locSettings)
{
    Initialize(locSettings);
}

void StabilizationImpl::Initialize(Settings& locSettings)
{
    thisclassname = "Stabilization";

    int Bcells = 0;

    if(Grid.Resolution == Resolutions::Dual)
    {
        Grid = locSettings.Grid.DoubleResolution();

        Bcells = 2 * std::max(locSettings.Grid.Bcells, int(Grid.iWidth/2));
    }
    else
    {
        Grid = locSettings.Grid;
        
        Bcells = std::max(locSettings.Grid.Bcells, int(Grid.iWidth)-1);
    }

    Curvature.Allocate(Grid, Bcells);
    AverageCurvature.Allocate(Grid, Bcells);
    AverageCurvatureWeight.Allocate(Grid, Bcells);
    AverageCurvatureAux.Allocate(Grid, Bcells);
    AverageCurvatureAuxWeight.Allocate(Grid, Bcells);

    Interface.Allocate(Grid, Bcells);
    AverageInterface.Allocate(Grid, Bcells);
    AverageInterfaceWeight.Allocate(Grid, Bcells);
    AverageInterfaceAux.Allocate(Grid, Bcells);
    AverageInterfaceAuxWeight.Allocate(Grid, Bcells);

    Velocity.Allocate(Grid, Bcells);

    initialized = true;

    maxdt = 1.;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void StabilizationImpl::Initialize(Settings& locSettings, LaplacianStencil _DStencil)
{
    Initialize(locSettings);

    DStencil = _DStencil;
    double dx = locSettings.Grid.dx;
    double Eta = locSettings.Grid.iWidth * dx;

    StencilDirections = std::vector<StencilDirection>();

    n_vector.resize(3);
    n_vector[0] = 1.0;
    n_vector[1] = 0.0;
    n_vector[2] = 0.0;

    double norm = std::sqrt(n_vector[0]*n_vector[0] + n_vector[1]*n_vector[1] + n_vector[2]*n_vector[2] );

    n_vector[0] /= norm ;
    n_vector[1] /= norm ;
    n_vector[2] /= norm ;
    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
    {
        double d_x = ds->di;
        double d_y = ds->dj;
        double d_z = ds->dk;

        StencilDirection* dirNew = new StencilDirection();
        dirNew->d_x = d_x;
        dirNew->d_y = d_y;
        dirNew->d_z = d_z;

        dirNew->stencil_weight = ds->weight;

        dirNew->scal_prod   = (d_x *n_vector[0] +  d_y*n_vector[1] +  d_z*n_vector[2])*dx ;

        dirNew->alpha        = std::cos(Pi*(dirNew->scal_prod)/Eta)   ;
        dirNew->beta         = std::sin(Pi*(dirNew->scal_prod)/Eta)   ;
        dirNew->phi_tl       = 0.5*std::cos(Pi*(1-(dirNew->scal_prod)/Eta )) + 0.5  ;
        dirNew->phi_tu       = 0.5*std::cos(Pi*(dirNew->scal_prod)/Eta) + 0.5  ;

        StencilDirections.push_back(*dirNew);
    }
}

void StabilizationImpl::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());
    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);
    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(EXIT_FAILURE);
    };

    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    ReadInput(inp);
}

void StabilizationImpl::ReadInput(std::stringstream& inp)
{
    ConsoleOutput::WriteLineInsert("Stabilization input");
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);
    
    range = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Range"), false, 2.5);
    exponent = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Exponent"), false, 4.);
    STAB = FileInterface::ReadParameterD(inp, moduleLocation, std::string("STAB"), false, 1./(10.*Pi*Pi));
    
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void StabilizationImpl::CalculateCurvature(PhaseField& Phase)
{
    //Phase.CalculateDerivativesSR();

    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        Curvature(i,j,k).clear();
        Interface(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    const double Prefactor = Pi*Pi/(Phase.Grid.Eta*Phase.Grid.Eta);

    size_t bcells = std::max((long int)0,PF.Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,bcells,)
    {
        if (PF(i,j,k).wide_interface())
        {
            NodePF& locPF = (Grid.Resolution == Resolutions::Single)?Phase.Fields(i,j,k):Phase.FieldsDR(i,j,k);
            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            {
                double dPhi_dt = (alpha->laplacian + Prefactor*(alpha->value-0.5));
                double norm_alpha = alpha->gradient.length();
                double curv = 0.;
                if (std::fabs(norm_alpha) > 1e-16*Phase.Grid.dx) curv = dPhi_dt/norm_alpha;
                Curvature(i,j,k).add_value(alpha->index, curv);
                Interface(i,j,k).add_value(alpha->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateCurvature_SPF(PhaseField& Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        Curvature(i,j,k).clear();
        Interface(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    size_t bcells = std::max((long int)0,PF.Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,bcells,)
    {
        if (PF(i,j,k).wide_interface())
        {
            NodePF& locPF = (Grid.Resolution == Resolutions::Single)?Phase.Fields(i,j,k):Phase.FieldsDR(i,j,k);
            for(auto alpha  = locPF.cbegin();
                     alpha != locPF.cend(); ++alpha)
            {
                double pot_term = 0;
                for(auto StDir = StencilDirections.begin() ; StDir != StencilDirections.end() ; StDir++)
                {
                    double fik_term = 0;
                    if(StDir->scal_prod < 0 && alpha->value > StDir->phi_tu )
                    {
                        fik_term = 1;
                    }
                    else if(StDir->scal_prod > 0 && alpha->value < StDir->phi_tl )
                    {
                        fik_term = 0;
                    }
                    else
                    {
                        fik_term = 0.5 + StDir->alpha*(alpha->value - 0.5 ) - StDir->beta *std::sqrt( alpha->value *(1.0 - alpha->value) );
                    }
                    pot_term += -StDir->stencil_weight * (fik_term - alpha->value);
                }
                double dPhi_dt = (alpha->laplacian + pot_term);
                double norm_alpha = alpha->gradient.length();
                double curv = 0.;
                if (std::fabs(norm_alpha) > 1e-16*Phase.Grid.dx) curv = dPhi_dt/norm_alpha;
                Curvature(i,j,k).add_value(alpha->index, curv);
                Interface(i,j,k).add_value(alpha->index, dPhi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateAverageCurvature(PhaseField &Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageCurvature(i,j,k).clear();
        AverageCurvatureWeight(i,j,k).clear();
        AverageCurvatureAux(i,j,k).clear();
        AverageCurvatureAuxWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface()  && Curvature(i+ii,j+jj,k+kk).get_value(alpha->index) != 0.0)
                    {
                        double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.0-locPhiAlphaValue),exponent)*(range-dist);
                        AverageCurvatureAux(i,j,k).add_value(alpha->index, weight*Curvature(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageCurvatureAuxWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureAuxWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvatureAux(i,j,k).set_value(alpha->index, AverageCurvatureAux(i,j,k).get_value(alpha->index)/AverageCurvatureAuxWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i + ii < PF.sizeX()+PF.Bcells() -1)
                if (j + jj >= -PF.Bcells()+1 && j + jj < PF.sizeY()+PF.Bcells() -1)
                if (k + kk >= -PF.Bcells()+1 && k + kk < PF.sizeZ()+PF.Bcells() -1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface() && AverageCurvatureAux(i+ii,j+jj,k+kk).get_value(alpha->index) != 0.0)
                    {
                        AverageCurvature(i,j,k).add_value(alpha->index,AverageCurvatureAux(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageCurvatureWeight(i,j,k).add_value(alpha->index, 1);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvature(i,j,k).set_value(alpha->index, AverageCurvature(i,j,k).get_value(alpha->index)/AverageCurvatureWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageInterface(i,j,k).clear();
        AverageInterfaceWeight(i,j,k).clear();
        AverageInterfaceAux(i,j,k).clear();
        AverageInterfaceAuxWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface())
                    {
                        double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.-locPhiAlphaValue),exponent)*(range-dist);
                        AverageInterfaceAux(i,j,k).add_value(alpha->index, weight*Interface(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageInterfaceAuxWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageInterfaceAuxWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterfaceAux(i,j,k).set_value(alpha->index, AverageInterfaceAux(i,j,k).get_value(alpha->index)/AverageInterfaceAuxWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha = PF(i,j,k).cbegin();
                    alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface())
                    {
                        AverageInterface(i,j,k).add_value(alpha->index,AverageInterfaceAux(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageInterfaceWeight(i,j,k).add_value(alpha->index, 1);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageInterfaceWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterface(i,j,k).set_value(alpha->index, AverageInterface(i,j,k).get_value(alpha->index)/AverageInterfaceWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateAverageCurvatureRegularization(PhaseField &Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageCurvature(i,j,k).clear();
        AverageCurvatureWeight(i,j,k).clear();
        AverageCurvatureAux(i,j,k).clear();
        AverageCurvatureAuxWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface()  && Curvature(i+ii,j+jj,k+kk).get_value(alpha->index) != 0.0)
                    {
                        double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.0-locPhiAlphaValue),0.5)*(range-dist);
                        AverageCurvatureAux(i,j,k).add_value(alpha->index, weight*Curvature(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageCurvatureAuxWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureAuxWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvatureAux(i,j,k).set_value(alpha->index, AverageCurvatureAux(i,j,k).get_value(alpha->index)/AverageCurvatureAuxWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i + ii < PF.sizeX()+PF.Bcells() -1)
                if (j + jj >= -PF.Bcells()+1 && j + jj < PF.sizeY()+PF.Bcells() -1)
                if (k + kk >= -PF.Bcells()+1 && k + kk < PF.sizeZ()+PF.Bcells() -1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface() && AverageCurvatureAux(i+ii,j+jj,k+kk).get_value(alpha->index) != 0.0)
                    {
                        //double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.0-locPhiAlphaValue),0.5);
                        AverageCurvature(i,j,k).add_value(alpha->index,weight*AverageCurvatureAux(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageCurvatureWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvature(i,j,k).set_value(alpha->index, AverageCurvature(i,j,k).get_value(alpha->index)/AverageCurvatureWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageInterface(i,j,k).clear();
        AverageInterfaceWeight(i,j,k).clear();
        AverageInterfaceAux(i,j,k).clear();
        AverageInterfaceAuxWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface())
                    {
                        double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.-locPhiAlphaValue),0.5)*(range-dist);
                        AverageInterfaceAux(i,j,k).add_value(alpha->index, weight*Interface(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageInterfaceAuxWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageInterfaceAuxWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterfaceAux(i,j,k).set_value(alpha->index, AverageInterfaceAux(i,j,k).get_value(alpha->index)/AverageInterfaceAuxWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -range*PF.dNx(); ii <= range*PF.dNx(); ++ii)
                for (int jj = -range*PF.dNy(); jj <= range*PF.dNy(); ++jj)
                for (int kk = -range*PF.dNz(); kk <= range*PF.dNz(); ++kk)
                if (ii*ii+jj*jj+kk*kk <= range*range)
                if (i + ii >= -PF.Bcells()+1 && i+ ii < PF.sizeX()+PF.Bcells()-1)
                if (j + jj >= -PF.Bcells()+1 && j+ jj < PF.sizeY()+PF.Bcells()-1)
                if (k + kk >= -PF.Bcells()+1 && k+ kk < PF.sizeZ()+PF.Bcells()-1)
                {
                    if (PF(i+ii,j+jj,k+kk).wide_interface())
                    {
                        //double dist = sqrt(ii*ii + jj*jj + kk*kk);
                        double locPhiAlphaValue = PF(i+ii,j+jj,k+kk).get_value(alpha->index);
                        double weight = std::pow(locPhiAlphaValue * (1.0-locPhiAlphaValue),0.5);
                        AverageInterface(i,j,k).add_value(alpha->index,weight*AverageInterfaceAux(i+ii,j+jj,k+kk).get_value(alpha->index));
                        AverageInterfaceWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageInterfaceWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterface(i,j,k).set_value(alpha->index, AverageInterface(i,j,k).get_value(alpha->index)/AverageInterfaceWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateAverageCurvatureNormal(PhaseField &Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageCurvature(i,j,k).clear();
        AverageCurvatureWeight(i,j,k).clear();
        AverageInterface(i,j,k).clear();
        AverageInterfaceWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells()-1,)
    {
        if (PF(i,j,k).wide_interface())
        {
            const NodePF& locPF = (Grid.Resolution == Resolutions::Single)?Phase.Fields(i,j,k):Phase.FieldsDR(i,j,k);
            for (auto alpha  = locPF.cbegin();
                      alpha != locPF.cend(); ++alpha)
            for (auto  beta  = locPF.cbegin();
                       beta != locPF.cend(); ++beta)
            if (alpha != beta)
            {
                dVector3 normal = Phase.Normal(alpha,beta);                
                for (int r = -2; r <= 2; ++r)
                {
                    double ii = i+0.5*r*range*normal[0];
                    double jj = j+0.5*r*range*normal[1];
                    double kk = k+0.5*r*range*normal[2];
                    if (ii >= -PF.Bcells()+1 && ii < PF.sizeX()+PF.Bcells()-1)
                    if (jj >= -PF.Bcells()+1 && jj < PF.sizeY()+PF.Bcells()-1)
                    if (kk >= -PF.Bcells()+1 && kk < PF.sizeZ()+PF.Bcells()-1)
                    {
                        double weight = std::pow(PF.at(ii,jj,kk).get_value(alpha->index)*PF.at(ii,jj,kk).get_value(beta->index), exponent);
                        AverageCurvature(i,j,k).add_value(alpha->index, weight*Curvature.at(ii,jj,kk).get_value(alpha->index));
                        AverageInterface(i,j,k).add_value(alpha->index, weight*Interface.at(ii,jj,kk).get_value(alpha->index));
                        AverageCurvatureWeight(i,j,k).add_value(alpha->index, weight);
                        AverageInterfaceWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }                
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvature(i,j,k).set_value(alpha->index, AverageCurvature(i,j,k).get_value(alpha->index)/AverageCurvatureWeight(i,j,k).get_value(alpha->index));
                }
                if (AverageInterfaceWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterface(i,j,k).set_value(alpha->index, AverageInterface(i,j,k).get_value(alpha->index)/AverageInterfaceWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateAverageCurvatureNormalSnap(PhaseField &Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageCurvature(i,j,k).clear();
        AverageCurvatureWeight(i,j,k).clear();
        AverageInterface(i,j,k).clear();
        AverageInterfaceWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells()-1,)
    {
        if (PF(i,j,k).wide_interface())
        {
            const NodePF& locPF = (Grid.Resolution == Resolutions::Single)?Phase.Fields(i,j,k):Phase.FieldsDR(i,j,k);
            for (auto alpha  = locPF.cbegin();
                      alpha != locPF.cend(); ++alpha)
            for (auto  beta  = locPF.cbegin();
                       beta != locPF.cend(); ++beta)
            if (alpha != beta)
            {
                dVector3 normal = Phase.Normal(alpha,beta);                
                for (int r = -2; r <= 2; ++r)
                {
                    int ii = round(i+0.5*r*range*normal[0]);
                    int jj = round(j+0.5*r*range*normal[1]);
                    int kk = round(k+0.5*r*range*normal[2]);
                    if (ii >= -PF.Bcells()+1 && ii < PF.sizeX()+PF.Bcells()-1)
                    if (jj >= -PF.Bcells()+1 && jj < PF.sizeY()+PF.Bcells()-1)
                    if (kk >= -PF.Bcells()+1 && kk < PF.sizeZ()+PF.Bcells()-1)
                    {
                        double weight = std::pow(PF.at(ii,jj,kk).get_value(alpha->index)*PF(ii,jj,kk).get_value(beta->index), exponent);
                        AverageCurvature(i,j,k).add_value(alpha->index, weight*Curvature(ii,jj,kk).get_value(alpha->index));
                        AverageInterface(i,j,k).add_value(alpha->index, weight*Interface(ii,jj,kk).get_value(alpha->index));
                        AverageCurvatureWeight(i,j,k).add_value(alpha->index, weight);
                        AverageInterfaceWeight(i,j,k).add_value(alpha->index, weight);
                    }
                }                
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                if (AverageCurvatureWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageCurvature(i,j,k).set_value(alpha->index, AverageCurvature(i,j,k).get_value(alpha->index)/AverageCurvatureWeight(i,j,k).get_value(alpha->index));
                }
                if (AverageInterfaceWeight(i,j,k).get_value(alpha->index) != 0.0)
                {
                    AverageInterface(i,j,k).set_value(alpha->index, AverageInterface(i,j,k).get_value(alpha->index)/AverageInterfaceWeight(i,j,k).get_value(alpha->index));
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void StabilizationImpl::CalculateAverageCurvatureSimple(PhaseField &Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        AverageCurvature(i,j,k).clear();
        AverageCurvatureWeight(i,j,k).clear();
        AverageInterface(i,j,k).clear();
        AverageInterfaceWeight(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,PF.Bcells(),)
    {
        if (PF(i,j,k).interface_halo())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                AverageCurvature(i,j,k).set_value(alpha->index, 0.0);
                AverageInterface(i,j,k).set_value(alpha->index, 0.0);
            }
        }
        if (PF(i,j,k).interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {

                AverageCurvature(i,j,k).set_value(alpha->index, Curvature(i,j,k).get_value(alpha->index));
                AverageInterface(i,j,k).set_value(alpha->index, Interface(i,j,k).get_value(alpha->index));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void StabilizationImpl::ComputeVelocity(PhaseField& Phase,
                                InterfaceProperties& IP, DrivingForce& dG)
{
    auto& PF  = (Grid.Resolution == Resolutions::Single) ? Phase.Fields : Phase.FieldsDR;
    auto& IPP = (Grid.Resolution == Resolutions::Single) ? IP.Properties : IP.PropertiesDR;
    
    double maxv = 0.;
    double maxdG = 0.;
    double maxSigma = 0.;
    int divisor = (Grid.Resolution == Resolutions::Single)?1:2;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,reduction(max:maxv,maxdG,maxSigma))
    {
        Velocity(i,j,k).clear();
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend() - 1; ++alpha)
            for(auto  beta  = alpha + 1;
                      beta != PF(i,j,k).cend(); ++beta)
            {
               
                double Sigma = IPP(i,j,k).get_energy(alpha->index, beta->index);
                double Mu = IPP(i,j,k).get_mobility(alpha->index, beta->index);
                double locdG = dG.Force(i/divisor,j/divisor,k/divisor).get_average(alpha->index, beta->index);
                double norm_1 = 1.0/PF(i,j,k).size();

                double avgCurvA = AverageCurvature(i, j, k).get_value(alpha->index);
                double avgCurvB = AverageCurvature(i, j, k).get_value(beta->index);
                double kappa_alpha = avgCurvA;
                double kappa_beta  = avgCurvB;
                double dPhi_dt = Mu*(norm_1*Sigma*(kappa_alpha - kappa_beta)  + locdG);
                for(auto gamma  = PF(i,j,k).cbegin();
                         gamma != PF(i,j,k).cend(); ++gamma)
                if (gamma != alpha && gamma != beta and Phase.FieldsProperties[gamma->index].Stage == GrainStages::Stable)
                {
                    double avgCurvC = AverageCurvature(i, j, k).get_value(gamma->index);
                    double kappa_gamma = avgCurvC;
                    dPhi_dt += norm_1*Mu*(IPP(i,j,k).get_energy(gamma->index,  beta->index) -
                                          IPP(i,j,k).get_energy(alpha->index, gamma->index))*kappa_gamma;
                }
                Velocity(i,j,k).add_asym1(alpha->index, beta->index,  dPhi_dt);
                maxv  = std::max(maxv, std::fabs(dPhi_dt));
                maxdG = std::max(maxdG,std::fabs(locdG));
                maxSigma = std::max(maxSigma, std::fabs(dPhi_dt-norm_1*Mu*locdG));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (maxv > 0)
    {
        maxdt = 16 * STAB*Phase.Grid.dx / maxv;
    }
    else
    {
        maxdt = 16 * STAB;
    }
    maxvelocity = maxv;
    if (maxdt < 1e-6*Phase.Grid.dx)
    {
        maxdt = 1e-6*Phase.Grid.dx;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,0,)
        {
            Velocity(i,j,k).clear();
            if (PF(i,j,k).wide_interface())
            {
                for(auto alpha  = PF(i,j,k).cbegin();
                         alpha != PF(i,j,k).cend() - 1; ++alpha)
                for(auto  beta  = alpha + 1;
                          beta != PF(i,j,k).cend(); ++beta)
                {
                    double vel =  Velocity(i,j,k).get_asym1(alpha->index, beta->index);
                    if (std::fabs(vel) > 16*STAB*Phase.Grid.dx/maxdt)
                    {
                        Velocity(i,j,k).set_asym1(alpha->index, beta->index,  ((vel>0)-(vel<0))*16*STAB*Phase.Grid.dx/maxdt);
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

double StabilizationImpl::GetMaxDT()
{
    return maxdt;
}

void StabilizationImpl::ComputeUpdate(PhaseField& Phase,
                                DrivingForce& dG, double dt)
{
    auto& PF  = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;
    auto& PFD = (Grid.Resolution == Resolutions::Single)?Phase.FieldsDot:Phase.FieldsDotDR;
    
    if (dt > 0)
    {
        double eta = Phase.Grid.Eta;
        double pi = Pi;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PF, 0, )
        {
            if (PF(i, j, k).wide_interface())
            {
                NodePF NodePFDR;
                if(Grid.Resolution == Resolutions::Dual)
                {
                     NodePFDR = Phase.FieldsDR(i,j,k);
                }
                NodePF& locPF = (Grid.Resolution == Resolutions::Single)?Phase.Fields(i,j,k):NodePFDR;
                for (auto alpha  = locPF.cbegin();
                          alpha != locPF.cend()-1; ++alpha)
                for (auto  beta  = alpha+1;
                           beta != locPF.cend(); ++beta)
                {
                    double velocity = Velocity(i, j, k).get_asym1(alpha->index, beta->index);
                    double locSTAB = STAB / dt;
                    double norm_1 = 1.0/PF(i,j,k).size();
                    double prefactor = pi / eta * sqrt(alpha->value*beta->value);
                    if (Phase.FieldsProperties[alpha->index].Stage == GrainStages::Seed or
                        Phase.FieldsProperties[ beta->index].Stage == GrainStages::Seed)
                    {
                        prefactor = std::max(prefactor, pi / eta * sqrt(1e-4));
                    }
                    double dPhi_dt1 = (norm_1*locSTAB*(eta*eta*Interface(i, j, k).get_value(alpha->index) - eta*eta*AverageInterface(i, j, k).get_value(alpha->index))
                                    -  norm_1*locSTAB*(eta*eta*Interface(i, j, k).get_value(beta->index)  - eta*eta*AverageInterface(i, j, k).get_value(beta->index)))
                                    *std::min((long double)1., std::min((long double)Phase.FieldsProperties[alpha->index].Volume / (long double)Phase.FieldsProperties[alpha->index].RefVolume,
                                    (long double)Phase.FieldsProperties[beta->index].Volume / (long double)Phase.FieldsProperties[beta->index].RefVolume));
                    double dPhi_dt2 = prefactor * velocity;
                    PFD(i, j, k).add_asym1(alpha->index, beta->index, dPhi_dt1+dPhi_dt2);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void StabilizationImpl::WriteVTK(const Settings& locSettings, const int tStep,
                                 const size_t indexA, const size_t indexB,
                                 const int precision) const
{
    std::vector<VTK::Field_t> ListOfFields;
    {
        std::stringstream converter;
        converter << indexA;
        std::string phase = converter.str();

        ListOfFields.push_back((VTK::Field_t) {"Interface(" + phase + ")", [indexA,this](int i,int j,int k){return Interface(i,j,k).get_value(indexA);}});
        ListOfFields.push_back((VTK::Field_t) {"AverageInterface(" + phase + ")", [indexA,this](int i,int j,int k){return AverageInterface(i,j,k).get_value(indexA);}});
        ListOfFields.push_back((VTK::Field_t) {"Curvature(" + phase + ")", [indexA,this](int i,int j,int k){return Curvature(i,j,k).get_value(indexA);}});
        ListOfFields.push_back((VTK::Field_t) {"AverageCurvature(" + phase + ")", [indexA,this](int i,int j,int k){return AverageCurvature(i,j,k).get_value(indexA);}});
    }
    {
    
        std::stringstream converter;
        converter << indexA << "," << indexB;
        std::string phases = converter.str();

        ListOfFields.push_back((VTK::Field_t) {"Velocity(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Velocity(i,j,k).get_asym1(indexA, indexB);}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "Stabilization_", tStep, ".vts");

    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void StabilizationImpl::Smooth(PhaseField& Phase)
{
    auto& PF = (Grid.Resolution == Resolutions::Single)?Phase.Fields:Phase.FieldsDR;

    size_t bcells = std::max(long(0),PF.Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,bcells,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = PF(i,j,k).cbegin();
                     alpha != PF(i,j,k).cend(); ++alpha)
            {
                for (int ii = -1; ii <= 1; ++ii)
                for (int jj = -1; jj <= 1; ++jj)
                for (int kk = -1; kk <= 1; ++kk)
                {
                    AverageInterface(i,j,k).add_value(alpha->index,PF(i+ii,j+jj,k+kk).get_value(alpha->index));
                    AverageInterfaceWeight(i,j,k).add_value(alpha->index, 1);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PF,bcells,)
    {
        if (PF(i,j,k).wide_interface())
        {
            for(auto alpha  = AverageInterface(i,j,k).cbegin();
                     alpha != AverageInterface(i,j,k).cend(); ++alpha)
            {
                PF(i,j,k).set_value(alpha->index, AverageInterface(i,j,k).get_value(alpha->index)/AverageInterfaceWeight(i,j,k).get_value(alpha->index));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

}// namespace openphase
