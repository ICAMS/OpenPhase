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

 *   File created :   2023
 *   Main contributors :   Oleg Shchyglo;
 *
 */

#include "InterfaceRegularization.h"
#include "GrainsProperties.h"
#include "Settings.h"
#include "VTK.h"
#include "PhaseField.h"
#include "InterfaceProperties.h"
#include "BoundaryConditions.h"

namespace openphase
{

using namespace std;
InterfaceRegularization::InterfaceRegularization(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void InterfaceRegularization::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "InterfaceRegularization";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;

    Averaging = true;

    PhiThreshold = DBL_EPSILON;

    Range = (Grid.iWidth + 1)/2;

    switch(Grid.Resolution)
    {
        case Resolutions::Single:
        {
            Curvature.Allocate(Grid, Range);
            break;
        }
        case Resolutions::Dual:
        {
            GridParameters DoubleDimensions = Grid.DoubleResolution();
            Curvature.Allocate(DoubleDimensions, Range);
            break;
        }
    }

    locSettings.AddForRemeshing(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void InterfaceRegularization::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("InterfaceRegularization input");
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

void InterfaceRegularization::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    Averaging    = FileInterface::ReadParameterB(inp, moduleLocation, string("Average"),   false, true);
    Range        = FileInterface::ReadParameterI(inp, moduleLocation, string("Range"),     false, Range);
    PhiThreshold = FileInterface::ReadParameterD(inp, moduleLocation, string("Threshold"), false, PhiThreshold);

    string tmp   = FileInterface::ReadParameterK(inp, moduleLocation, string("ClearingMode"), false, "AUTOMATIC");

    if(tmp == string("AUTOMATIC")) ClearingMode = ClearingModes::Automatic;
    if(tmp == string("MANUAL")) ClearingMode = ClearingModes::Manual;

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void InterfaceRegularization::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,Curvature.Bcells(),)
    {
        Curvature(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::SetBoundaryConditions(const BoundaryConditions& BC)
{
    if(Grid.dNx) BC.SetX(Curvature);
    if(Grid.dNy) BC.SetY(Curvature);
    if(Grid.dNz) BC.SetZ(Curvature);
}

void InterfaceRegularization::Average(const PhaseField& Phase, const BoundaryConditions& BC)
{
    if(Averaging)
    {
        if(Grid.Resolution == Resolutions::Single)
        {
            SetWeightsSR(Phase);
            SetBoundaryConditions(BC);
            CollectAverageSR(Phase);
            SetBoundaryConditions(BC);
            DistributeAverageSR(Phase);
        }
        else
        {
            SetWeightsDR(Phase);
            SetBoundaryConditions(BC);
            CollectAverageDR(Phase);
            SetBoundaryConditions(BC);
            DistributeAverageDR(Phase);
        }
    }
}

void InterfaceRegularization::SetWeightsSR(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get_value(it->indexB);

            it->weight = sqrt(PhiAlpha*PhiBeta);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::CollectAverageSR(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            if (it->weight > weight_threshold)
            {
                double value = 0.0;

                double SumWeights = 0.0;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                    if((Range - dist) > 0.0 and Phase.Fields(i+ii,j+jj,k+kk).interface())
                    {
                        double weight = Curvature(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB)*(Range - dist);

                        if (weight > DBL_EPSILON)
                        {
                            SumWeights += weight;
                            value += weight*Curvature(i+ii, j+jj, k+kk).get_raw(it->indexA, it->indexB);
                        }
                    }
                }
                if(SumWeights > DBL_EPSILON)
                {
                    it->tmp = value/SumWeights;
                }
            }
            else if(Phase.FieldsProperties[it->indexA].Stage == GrainStages::Seed or
                    Phase.FieldsProperties[it->indexB].Stage == GrainStages::Seed)
            {
                it->tmp = it->raw;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::DistributeAverageSR(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            if (it->weight > DBL_EPSILON)
            {
                int  counter = 0;
                double value = 0.0;
                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                    if((Range - dist) > 0.0 and Phase.Fields(i+ii,j+jj,k+kk).interface())
                    {
                        double weight = Curvature(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB);

                        if (weight > weight_threshold)
                        {
                            counter ++;
                            value += Curvature(i+ii, j+jj, k+kk).get_tmp(it->indexA, it->indexB);
                        }
                    }
                }
                if(counter != 0)
                {
                    it->average = value/counter;
                }
                else
                {
                    it->average = it->tmp;
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

void InterfaceRegularization::SetWeightsDR(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.FieldsDR(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.FieldsDR(i,j,k).get_value(it->indexA);
            double PhiBeta  = Phase.FieldsDR(i,j,k).get_value(it->indexB);

            it->weight = sqrt(PhiAlpha*PhiBeta);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::CollectAverageDR(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.FieldsDR(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            if (it->weight > weight_threshold)
            {
                double value = 0.0;

                double SumWeights = 0.0;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                    if((Range - dist) > 0.0 and Phase.FieldsDR(i+ii,j+jj,k+kk).interface())
                    {
                        double weight = Curvature(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB)*(Range - dist);

                        if (weight > DBL_EPSILON)
                        {
                            SumWeights += weight;
                            value += weight*Curvature(i+ii, j+jj, k+kk).get_raw(it->indexA, it->indexB);
                        }
                    }
                }
                if(SumWeights > DBL_EPSILON)
                {
                    it->tmp = value/SumWeights;
                }
            }
            else if(Phase.FieldsProperties[it->indexA].Stage == GrainStages::Seed or
                    Phase.FieldsProperties[it->indexB].Stage == GrainStages::Seed)
            {
                it->tmp = it->raw;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::DistributeAverageDR(const PhaseField& Phase)
{
    const int Xrange = min(Range, Grid.Nx-1)*Grid.dNx;
    const int Yrange = min(Range, Grid.Ny-1)*Grid.dNy;
    const int Zrange = min(Range, Grid.Nz-1)*Grid.dNz;

    const double weight_threshold = sqrt(PhiThreshold*(1.0 - PhiThreshold));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Curvature,0,)
    {
        if (Phase.FieldsDR(i,j,k).interface())
        for(auto it  = Curvature(i,j,k).begin();
                 it != Curvature(i,j,k).end(); ++it)
        {
            if (it->weight > DBL_EPSILON)
            {
                int counter = 0;
                double value = 0.0;
                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                    if((Range - dist) > 0.0 and Phase.FieldsDR(i+ii,j+jj,k+kk).interface())
                    {
                        double weight = Curvature(i+ii, j+jj, k+kk).get_weight(it->indexA, it->indexB);

                        if (weight > weight_threshold)
                        {
                            counter ++;
                            value += Curvature(i+ii, j+jj, k+kk).get_tmp(it->indexA, it->indexB);
                        }
                    }
                }
                if(counter != 0)
                {
                    it->average = value/counter;
                }
                else
                {
                    it->average = it->tmp;
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

void InterfaceRegularization::MergePhaseFieldIncrements(PhaseField& Phase, InterfaceProperties& IP)
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

    if(ClearingMode == ClearingModes::Automatic)
    {
        Clear();
    }
}

void InterfaceRegularization::MergePhaseFieldIncrementsSR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    const double Prefactor = Pi/Phase.Grid.Eta;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Curvature, 0, )
    {
        if (Phase.Fields(i,j,k).interface())
        {
            for(auto it  = Curvature(i,j,k).begin();
                     it != Curvature(i,j,k).end(); ++it)
            {
                double locCurvature = 0.0;
                if(Averaging)
                {
                    locCurvature = it->average;
                }
                else
                {
                    locCurvature = it->raw;
                }

                double norm = it->weight;

                double dPsi_dt = locCurvature*IP.Properties(i,j,k).get_mobility(it->indexA,it->indexB)*norm*Prefactor;

                Phase.FieldsDot(i,j,k).add_asym1(it->indexA,it->indexB,dPsi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::MergePhaseFieldIncrementsDR(PhaseField& Phase,
                                                     InterfaceProperties& IP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Curvature, 0, )
    {
        if (Phase.FieldsDR(i,j,k).interface())
        {
            double Prefactor = Pi/Phase.Grid.Eta;

            for(auto it  = Curvature(i,j,k).begin();
                     it != Curvature(i,j,k).end(); ++it)
            {
                double locCurvature = 0.0;
                if(Averaging)
                {
                    locCurvature = it->average;
                }
                else
                {
                    locCurvature = it->raw;
                }

                double norm = it->weight;

                double dPsi_dt = locCurvature*IP.PropertiesDR(i,j,k).get_mobility(it->indexA,it->indexB)*norm*Prefactor;

                Phase.FieldsDotDR(i,j,k).add_asym1(it->indexA,it->indexB,dPsi_dt);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceRegularization::PrintPointStatistics(const int x, const int y, const int z) const
{
    std::stringstream pointstat;
    pointstat << "Curvature Indices:\t";

    for (auto alpha  = Curvature(x,y,z).cbegin();
              alpha != Curvature(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->indexA << ", " << alpha->indexB << "\t\t";
    }
    pointstat << endl;
    pointstat << "Curvature Values:\t";

    for (auto alpha  = Curvature(x,y,z).cbegin();
              alpha != Curvature(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->average << "\t\t";
    }
    pointstat << endl;

    ConsoleOutput::WriteSimple(pointstat.str());
}

void InterfaceRegularization::WriteVTK(const Settings& locSettings,
                                       const int tStep,
                                       const size_t indexA,
                                       const size_t indexB,
                                       const int precision) const
{
    stringstream converter;
    converter << indexA << "," << indexB;
    string phases = converter.str();

    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Curvature_avg(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Curvature(i,j,k).get_average(indexA, indexB);}});
    ListOfFields.push_back((VTK::Field_t) {"Curvature_raw(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Curvature(i,j,k).get_raw(indexA, indexB);}});
    ListOfFields.push_back((VTK::Field_t) {"Curvature_tmp(" + phases + ")", [indexA,indexB,this](int i,int j,int k){return Curvature(i,j,k).get_tmp(indexA, indexB);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir, "InterfaceRegularization_", tStep, ".vts");

    if(Grid.Resolution == Resolutions::Single)
    {
        VTK::Write(Filename, locSettings, ListOfFields, precision);
    }
    else
    {
        VTK::Write(Filename, locSettings, ListOfFields, precision, 2);
    }
}

void InterfaceRegularization::Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC)
{
    if(Grid.Resolution == Resolutions::Dual)
    {
        Curvature.Reallocate((1+Grid.dNx)*newNx, (1+Grid.dNy)*newNy, (1+Grid.dNz)*newNz);
    }
    else
    {
        Curvature.Reallocate(newNx, newNy, newNz);
    }

    Grid.SetDimensions(newNx, newNy, newNz);

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

InterfaceRegularization& InterfaceRegularization::operator= (const InterfaceRegularization& rhs)
{
    // protect against self-assignment and copy of uninitialized object
    if (this != &rhs and rhs.thisclassname == "InterfaceRegularization")
    {
        thisclassname = rhs.thisclassname;

        Grid = rhs.Grid;

        Range        = rhs.Range;
        PhiThreshold = rhs.PhiThreshold;

        Averaging  = rhs.Averaging;

        Curvature  = rhs.Curvature;
    }
    return *this;
}

NodeDF InterfaceRegularization::at(const double x, const double y, const double z) const
{
#ifdef DEBUG
    if(x > Curvature.sizeX() + Curvature.BcellsX() - 1 or
       y > Curvature.sizeY() + Curvature.BcellsY() - 1 or
       z > Curvature.sizeZ() + Curvature.BcellsZ() - 1 or
       x < -Curvature.BcellsX() or
       y < -Curvature.BcellsY() or
       z < -Curvature.BcellsZ())
    {
        std::stringstream message;
        message << "ERROR: InterfaceCurvature::Curvature_at()\n"
                << "Access beyond storage range -> ("
                << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                << -Curvature.BcellsX() << ", " << -Curvature.BcellsY() << ", " << -Curvature.BcellsZ() << "] and ("
                << Curvature.sizeX() + Curvature.BcellsX() << ", "
                << Curvature.sizeY() + Curvature.BcellsY() << ", "
                << Curvature.sizeZ() + Curvature.BcellsZ() << ")\n"
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

    NodeDF locCurvature;
    double weight = 0.0;
    if(Curvature(x0, y0, z0).size())
    {
        locCurvature += (Curvature(x0, y0, z0)*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)));
        weight += ((1.0 - dx)*(1.0 - dy)*(1.0 - dz));
    }
    if(Curvature(x0+Grid.dNx, y0, z0).size())
    {
        locCurvature += (Curvature(x0+Grid.dNx, y0, z0)*(dx*(1.0 - dy)*(1.0 - dz)));
        weight += (dx*(1.0 - dy)*(1.0 - dz));
    }
    if(Curvature(x0, y0+Grid.dNy, z0).size())
    {
        locCurvature += (Curvature(x0, y0+Grid.dNy, z0)*((1.0 - dx)*dy*(1.0 - dz)));
        weight += ((1.0 - dx)*dy*(1.0 - dz));
    }
    if(Curvature(x0, y0, z0+Grid.dNz).size())
    {
        locCurvature += (Curvature(x0, y0, z0+Grid.dNz)*((1.0 - dx)*(1.0 - dy)*dz));
        weight += ((1.0 - dx)*(1.0 - dy)*dz);
    }
    if(Curvature(x0+Grid.dNx, y0+Grid.dNy, z0).size())
    {
        locCurvature += (Curvature(x0+Grid.dNx, y0+Grid.dNy, z0)*(dx*dy*(1.0 - dz)));
        weight += (dx*dy*(1.0 - dz));
    }
    if(Curvature(x0+Grid.dNx, y0, z0+Grid.dNz).size())
    {
        locCurvature += (Curvature(x0+Grid.dNx, y0, z0+Grid.dNz)*(dx*(1.0 - dy)*dz));
        weight += (dx*(1.0 - dy)*dz);
    }
    if(Curvature(x0, y0+Grid.dNy, z0+Grid.dNz).size())
    {
        locCurvature += (Curvature(x0, y0+Grid.dNy, z0+Grid.dNz)*((1.0 - dx)*dy*dz));
        weight += ((1.0 - dx)*dy*dz);
    }
    if(Curvature(x0+Grid.dNx, y0+Grid.dNy, z0+Grid.dNz).size())
    {
        locCurvature += (Curvature(x0+Grid.dNx, y0+Grid.dNy, z0+Grid.dNz)*(dx*dy*dz));
        weight += (dx*dy*dz);
    }
    if(weight >= DBL_EPSILON)
    {
        locCurvature *= 1.0/weight;
    }
    return locCurvature;
}

}// namespace openphase

