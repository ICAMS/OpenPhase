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
 *   Main contributors : Raphael Schiedung
 *
 */

#include "DoubleObstacle.h"
#include "GrandPotential/GrandPotentialSolver.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Tools/AnalysisSintering.h"

namespace openphase
{

void AnalysisSintering::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "AnalysisSintering";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
};

void AnalysisSintering::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName.c_str());

    std::fstream inpF(InputFileName.c_str(), std::ios::in | std::ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();
    ReadInput(inp);
}

void AnalysisSintering::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);
    ReferenceDensity     = FileInterface::ReadParameterD(inp, moduleLocation, "Rho0", false, 1.0);
    Do_LineIntersection  = FileInterface::ReadParameterB(inp, moduleLocation, "LI",   false,false);
    Do_LineIntersectionX = FileInterface::ReadParameterB(inp, moduleLocation, "LIX",  false,false);
    Do_LineIntersectionY = FileInterface::ReadParameterB(inp, moduleLocation, "LIY",  false,false);
    Do_LineIntersectionZ = FileInterface::ReadParameterB(inp, moduleLocation, "LIZ",  false,false);

    BoundingBoxSizeX = std::round(FileInterface::ReadParameterD(inp, moduleLocation, "BBSX")/Grid.dx);
    BoundingBoxSizeY = std::round(FileInterface::ReadParameterD(inp, moduleLocation, "BBSY")/Grid.dx);
    BoundingBoxSizeZ = std::round(FileInterface::ReadParameterD(inp, moduleLocation, "BBSZ")/Grid.dx);
    
    Do_BoundingBox = (BoundingBoxSizeX > 0) or (BoundingBoxSizeY > 0) or (BoundingBoxSizeZ > 0);
    if (Do_BoundingBox)
    {
        BoundingBoxSizeX = std::clamp(BoundingBoxSizeX,1l,(long)Grid.TotalNx);
        BoundingBoxSizeY = std::clamp(BoundingBoxSizeY,1l,(long)Grid.TotalNy);
        BoundingBoxSizeZ = std::clamp(BoundingBoxSizeZ,1l,(long)Grid.TotalNz);
        SetBoundingBoxCoordinates();
    }
}

void AnalysisSintering::SetBoundingBoxCoordinates()
{
    // Calculate global bounding box coordinates
    long BoundingBox_min_i = (Grid.TotalNx - BoundingBoxSizeX)/2;
    long BoundingBox_max_i = (Grid.TotalNx - BoundingBoxSizeX)/2 + BoundingBoxSizeX;
    long BoundingBox_min_j = (Grid.TotalNy - BoundingBoxSizeY)/2;
    long BoundingBox_max_j = (Grid.TotalNy - BoundingBoxSizeY)/2 + BoundingBoxSizeY;
    long BoundingBox_min_k = (Grid.TotalNz - BoundingBoxSizeZ)/2;
    long BoundingBox_max_k = (Grid.TotalNz - BoundingBoxSizeZ)/2 + BoundingBoxSizeZ;
    BoundingBox_min_i = std::clamp(BoundingBox_min_i,0l,(long)Grid.TotalNx);
    BoundingBox_max_i = std::clamp(BoundingBox_max_i,0l,(long)Grid.TotalNx);
    BoundingBox_min_j = std::clamp(BoundingBox_min_j,0l,(long)Grid.TotalNy);
    BoundingBox_max_j = std::clamp(BoundingBox_max_j,0l,(long)Grid.TotalNy);
    BoundingBox_min_k = std::clamp(BoundingBox_min_k,0l,(long)Grid.TotalNz);
    BoundingBox_max_k = std::clamp(BoundingBox_max_k,0l,(long)Grid.TotalNz);
    BoundingBoxSizeX = BoundingBox_max_i - BoundingBox_min_i;
    BoundingBoxSizeY = BoundingBox_max_j - BoundingBox_min_j;
    BoundingBoxSizeZ = BoundingBox_max_k - BoundingBox_min_k;

    // Calculate local bounding box coordinates (on MPI-thread)
    locBoundingBox_min_i = BoundingBox_min_i - Grid.OffsetX;
    locBoundingBox_max_i = BoundingBox_max_i - Grid.OffsetX;
    locBoundingBox_min_j = BoundingBox_min_j - Grid.OffsetY;
    locBoundingBox_max_j = BoundingBox_max_j - Grid.OffsetY;
    locBoundingBox_min_k = BoundingBox_min_k - Grid.OffsetZ;
    locBoundingBox_max_k = BoundingBox_max_k - Grid.OffsetZ;
    if (locBoundingBox_min_i < 0)  locBoundingBox_min_i = 0;
    if (locBoundingBox_min_j < 0)  locBoundingBox_min_j = 0;
    if (locBoundingBox_min_k < 0)  locBoundingBox_min_k = 0;
    if (locBoundingBox_max_i > Grid.Nx) locBoundingBox_max_i = Grid.Nx;
    if (locBoundingBox_max_j > Grid.Ny) locBoundingBox_max_j = Grid.Ny;
    if (locBoundingBox_max_k > Grid.Nz) locBoundingBox_max_k = Grid.Nz;

    // Check if Bounding box size is non zero (on MPI-thread)
    if (locBoundingBox_min_i < Grid.Nx)
    if (locBoundingBox_min_j < Grid.Ny)
    if (locBoundingBox_min_k < Grid.Nz)
    if (locBoundingBox_max_i >= 0)
    if (locBoundingBox_max_j >= 0)
    if (locBoundingBox_max_k >= 0)
    {
        Do_BoundingBox = true;
    }

    BoundingBoxVolume = BoundingBoxSizeX*Grid.dx*
                        BoundingBoxSizeY*Grid.dx*
                        BoundingBoxSizeZ*Grid.dx;
    assert(BoundingBoxVolume > 0.0);
    assert(BoundingBoxVolume < std::numeric_limits<double>::max());

}

double AnalysisSintering::CalculateDensity_BoundingBox(const Field3D& MassDensity) const
{
    if (Do_BoundingBox)
    {
        double locMass = 0;
        #pragma omp parallel for collapse(3) reduction(+:locMass)
        for (long i = locBoundingBox_min_i; i < locBoundingBox_max_i; ++i)
        for (long j = locBoundingBox_min_j; j < locBoundingBox_max_j; ++j)
        for (long k = locBoundingBox_min_k; k < locBoundingBox_max_k; ++k)
        {
            locMass += MassDensity(i,j,k);
        }
        locMass *= Grid.CellVolume(true);
        double Mass = locMass;
        #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(&locMass, &Mass, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif

        return Mass/BoundingBoxVolume;
    }
    return 0;
}

void AnalysisSintering::PrintSolidPhaseFractionInBoundingBox(const PhaseField& Phase, const GrandPotentialDensity& omega, const GrandPotentialSolver& GPS)
{
    auto SolidPhaseFraction = [&Phase](long i, long j, long k){return Phase.SolidPhaseFraction(i,j,k);};
    ConsoleOutput::Write("Solid phase fraction in Bounding Box [%]", CalculateDensity_BoundingBox(SolidPhaseFraction));
}

void AnalysisSintering::WriteSolidPhaseFractionInBoundingBox(const PhaseField& Phase, const GrandPotentialDensity& omega, const GrandPotentialSolver& GPS, std::string filename, double time, char separator)
{
    auto SolidPhaseFraction = [&Phase](long i, long j, long k){return Phase.SolidPhaseFraction(i,j,k);};
    double Density  = CalculateDensity_BoundingBox(SolidPhaseFraction);
#ifdef MPI_PARALLEL
    if (MPI_RANK == 0)
#endif
    {
        bool WriteHeader = !(std::ifstream(filename).good());

        std::ofstream outp(filename, std::ios::out | std::ios::app );
        if (!outp)
        {
            ConsoleOutput::WriteExit("File \""+filename+"\" could not be opened\n", "AnalysisSintering", "WriteDensityInBoundingBox");
            OP_Exit(EXIT_FAILURE);
        };
        if (WriteHeader)
        {
            outp << "Time[s]" << separator << "SolidPhaseFraction[%]\n";
        }
        outp << std::defaultfloat;
        outp << time;
        outp << std::scientific;
        outp << separator << Density;
        outp << "\n";
        outp.close();
    }
}

long AnalysisSintering::First_SolidXCoordinate(const PhaseField &Phase, size_t GasPhaseIdx, long j, long k)
{
    long loc_x_min = Phase.Grid.TotalNx;
    for (long n = 0; n < Phase.Fields.sizeX(); ++n)
    {
        if (Phase.Fractions(n,j,k,{GasPhaseIdx}) <= 0.05)
        {
            loc_x_min = n + Phase.Grid.OffsetX;
            break;
        }
    }
    long x_min = loc_x_min;
    #ifdef MPI_PARALLEL
    OP_MPI_Allreduce(&loc_x_min, &x_min, 1, OP_MPI_LONG, OP_MPI_MIN, OP_MPI_COMM_WORLD);
    #endif
    return x_min;
}
long AnalysisSintering::Last_SolidXCoordinate(const PhaseField &Phase, size_t GasPhaseIdx, long j, long k)
{
    long loc_x_max = 0;
    for (long n = Phase.Fields.sizeX()-1; n >= 0; --n)
    {
        if (Phase.Fractions(n,j,k,{GasPhaseIdx}) <= 0.05)
        {
            loc_x_max = n + Phase.Grid.OffsetX;
            break;
        }
    }
    long x_max = loc_x_max;
    #ifdef MPI_PARALLEL
    OP_MPI_Allreduce(&loc_x_max, &x_max, 1, OP_MPI_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    #endif
    return x_max;
}
double AnalysisSintering::CalculateDensity_LineIntersectionX(const PhaseField &Phase,
        const Field3D& MassDensity)
{
    if (Phase.Grid.Nx > 1)
    {
        double LocMass   = 0;
        double LocVolume = 0;
        #ifndef MPI_PARALLEL
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        #endif
        for (long j = 0; j < Phase.Grid.Ny; ++j)
        for (long k = 0; k < Phase.Grid.Nz; ++k)
        {
            long ii_min = First_SolidXCoordinate(Phase,0,j,k);
            long ii_max =  Last_SolidXCoordinate(Phase,0,j,k);
            long i_min = std::clamp(ii_min - long(Phase.Grid.OffsetX), 0l, long(Phase.Grid.Nx));
            long i_max = std::clamp(ii_max - long(Phase.Grid.OffsetX), 0l, long(Phase.Grid.Nx));

            for (long i = i_min; i < i_max; ++i)
            {
                LocMass   += MassDensity(i,j,k);
                LocVolume += 1.0;
            }
        }
        LocMass*=Phase.Grid.CellVolume();
        LocVolume*=Phase.Grid.CellVolume();
        double Mass   = LocMass;
        double Volume = LocVolume;
        #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(&LocMass,   &Mass,   1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&LocVolume, &Volume, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif
       if (Volume > 0.0) return Mass/Volume;
       else return 0.0;
    }
    return 0.0;
}

double AnalysisSintering::CalculateDensity_LineIntersectionY(
        const PhaseField &Phase, const Field3D& MassDensity)
{
    #ifdef MPI_PARALLEL
    if (MPI_3D_DECOMPOSITION)
    {
        ConsoleOutput::WriteExit("MPI_3D_DECOMPOSITION is not supported", "", "LineIntersectionY");
        OP_Exit(EXIT_FAILURE);
    }
    #endif

    if (Phase.Grid.Ny > 1)
    {
        double LocMass   = 0;
        double LocVolume = 0;
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long i = 0; i < Phase.Grid.Nx; ++i)
        for (long k = 0; k < Phase.Grid.Nz; ++k)
        {
            long j_min = Phase.Grid.Ny-1;
            for (int n = 0; n < Phase.Grid.Ny; ++n)
            if (Phase.Fractions(i,n,k,{0}) <= 0.05)
            {
                j_min = n;
                break;
            }
            long j_max = 0;
            for (int n = Phase.Grid.Ny-1; n > 0; --n)
            if (Phase.Fractions(i,n,k,{0}) <= 0.05)
            {
                j_max = n;
                break;
            }
            for (long j = j_min; j < j_max; ++j)
            {
                LocMass   += MassDensity(i,j,k);
                LocVolume += 1.0;
            }
        }
        LocMass*=Phase.Grid.CellVolume();
        LocVolume*=Phase.Grid.CellVolume();
        double Mass   = LocMass;
        double Volume = LocVolume;
        #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(&LocMass,   &Mass,   1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&LocVolume, &Volume, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif
        if (Volume > 0.0) return Mass/Volume;
        else return 0.0;
    }
    return 0.0;
}
double AnalysisSintering::CalculateDensity_LineIntersectionZ(
        const PhaseField &Phase, const Field3D& MassDensity)
{
    #ifdef MPI_PARALLEL
    if (MPI_3D_DECOMPOSITION)
    {
        ConsoleOutput::WriteExit("MPI_3D_DECOMPOSITION is not supported", "", "LineIntersectionZ");
        OP_Exit(EXIT_FAILURE);
    }
    #endif

    double LocMass   = 0;
    double LocVolume = 0;

    if (Phase.Grid.Nz > 1)
    {
        #pragma omp parallel for collapse(2) reduction(+:LocVolume,LocMass)
        for (long i = 0; i < Phase.Grid.Nx; ++i)
        for (long j = 0; j < Phase.Grid.Ny; ++j)
        {
            long k_min = Phase.Grid.Nz-1;
            for (int n = 0; n < Phase.Grid.Nz; ++n)
            if (Phase.Fractions(i,j,n,{0}) <= 0.05)
            {
                k_min = n;
                break;
            }
            long k_max = 0;
            for (int n = Phase.Grid.Nz-1; n > 0; --n)
            if (Phase.Fractions(i,j,n,{0}) <= 0.05)
            {
                k_max = n;
                break;
            }
            for (long k = k_min; k <= k_max; ++k)
            {
                LocMass   += MassDensity(i,j,k);
                LocVolume += 1.0;
            }
        }
        LocMass*=Phase.Grid.CellVolume();
        LocVolume*=Phase.Grid.CellVolume();
        double Mass   = LocMass;
        double Volume = LocVolume;
        #ifdef MPI_PARALLEL
        OP_MPI_Allreduce(&LocMass,   &Mass,   1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        OP_MPI_Allreduce(&LocVolume, &Volume, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        #endif
        if (Volume > 0.0) return Mass/Volume;
        else return 0.0;
    }
    return 0.0;
}

double AnalysisSintering::CalculateDensity_LineIntersection(
        const PhaseField &Phase, const Field3D& MassDensity)
{
    double dim = 0;
    if (Phase.Grid.Nx > 1) dim+=1.0;
    if (Phase.Grid.Ny > 1) dim+=1.0;
    if (Phase.Grid.Nz > 1) dim+=1.0;
    double DensityLIX = CalculateDensity_LineIntersectionX(Phase,MassDensity);
    double DensityLIY = CalculateDensity_LineIntersectionY(Phase,MassDensity);
    double DensityLIZ = CalculateDensity_LineIntersectionZ(Phase,MassDensity);
    double DensityLI  = 0.0;
    if (Phase.Grid.Nx > 1) DensityLI += DensityLIX;
    if (Phase.Grid.Ny > 1) DensityLI += DensityLIY;
    if (Phase.Grid.Nz > 1) DensityLI += DensityLIZ;
    DensityLI /= dim;
    return DensityLI;
}

void AnalysisSintering::DoDiagnostics(Settings& locSettings, PhaseField& Phase,
        DoubleObstacle& DO, InterfaceProperties& IP, 
        GrandPotentialDensity& omega, GrandPotentialSolver& GPS, RunTimeControl& RTC)
{
    std::fstream log(locSettings.TextDir + thisclassname+".csv", std::ios::app|std::ios::out);
    double Density    = 0.0;
    double DensityLI  = 0.0;
    double DensityLIX = 0.0;
    double DensityLIY = 0.0;
    double DensityLIZ = 0.0;

    auto locMassDensity = [&GPS,&Phase,&omega] (long i, long j, long k) {return  GPS.MassDensity(i,j,k,Phase,omega);};
    auto SolidPhaseFraction = [&Phase](long i, long j, long k){return Phase.SolidPhaseFraction(i,j,k);};
    double RhoRel = CalculateDensity_LineIntersection(Phase,SolidPhaseFraction);
    if (Do_BoundingBox) Density = CalculateDensity_BoundingBox(locMassDensity);
    if (Do_LineIntersectionX) DensityLIX = CalculateDensity_LineIntersectionX(Phase,locMassDensity);
    if (Do_LineIntersectionY) DensityLIY = CalculateDensity_LineIntersectionY(Phase,locMassDensity);
    if (Do_LineIntersectionZ) DensityLIZ = CalculateDensity_LineIntersectionZ(Phase,locMassDensity);
    if (Do_LineIntersection)  DensityLI  = CalculateDensity_LineIntersection (Phase,locMassDensity);

    std::array<std::stringstream,2> line;
    line[1] << std::scientific << std::setprecision(16);
    ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Time [s]", RTC.SimulationTime);
    ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Interface energy [J]"      , DO.Energy(Phase, IP));
    //TODO Readd Surface an Interface area again
    //ConsoleOutput::WriteWithLog(line , RTC.tStep , "Surface energy [J]"        , DO.Energy(Phase, IP,0,1));
    //ConsoleOutput::WriteWithLog(line , RTC.tStep , "Grain boundary energy [J]" , DO.Energy(Phase, IP,1,1));
    for (size_t comp = 0; comp < locSettings.Ncomp; comp++)
    {
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Amount of "             + locSettings.ElementNames[comp]+" [mol]"    , GPS.TotalAmountOfComponent(comp));
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Change of Amount of "   + locSettings.ElementNames[comp]+" [%]  "    , GPS.TotalAmountOfComponentChange(comp)*100);
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Vapor concentration of "+ locSettings.ElementNames[comp]+" [mol/m^3]", GPS.Concentration(0,0,0,comp,Phase,omega));
    }
    ConsoleOutput::Write("-");
    if (Do_BoundingBox) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Density BB [kg/m^3]", Density);
    if (Do_LineIntersection) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Density LI [kg/m^3]", DensityLI);
    if (Do_LineIntersectionX and locSettings.Grid.Nx > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Density LIX [kg/m^3]" , DensityLIX);
    if (Do_LineIntersectionY and locSettings.Grid.Ny > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Density LIY [kg/m^3]" , DensityLIY);
    if (Do_LineIntersectionX and locSettings.Grid.Nz > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Density LIZ [kg/m^3]" , DensityLIZ);
    ConsoleOutput::Write("-");
    if (ReferenceDensity != 1.0)
    {
        if (Do_BoundingBox) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Density BB"  , Density/ReferenceDensity);
        if (Do_LineIntersection) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Density LI"  , DensityLI/ReferenceDensity);
        if (Do_LineIntersectionX and locSettings.Grid.Nx > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Density LIX" , DensityLIX/ReferenceDensity);
        if (Do_LineIntersectionY and locSettings.Grid.Ny > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Density LIY" , DensityLIY/ReferenceDensity);
        if (Do_LineIntersectionX and locSettings.Grid.Nz > 1) ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Density LIZ" , DensityLIZ/ReferenceDensity);
        ConsoleOutput::Write("-");
    }

    unsigned int dim = 0;
    if (locSettings.Grid.Nx > 1) dim++;
    if (locSettings.Grid.Ny > 1) dim++;
    if (locSettings.Grid.Nz > 1) dim++;
    const double dx = locSettings.Grid.dx;
    const double dV = std::pow(dx,dim);
    ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Relative Solid Phase Fraction LI", RhoRel);
    ConsoleOutput::Write("-");
    for (size_t n = 1; n < locSettings.Nphases; n++)
    {
        size_t NGrains      = Phase.FieldsProperties.NumberOfGrains();
        double MeanVolume   = Phase.FieldsProperties.MeanGrainVolume(n)*dV;
        double StddevVolume = Phase.FieldsProperties.StandardDeviationOfGrainVolumes(n)*dV;
        //double MeanRadius   = Phase.FieldsProperties.MeanGrainRadius(n,dim)*dx;
        //double StddevRadius = Phase.FieldsProperties.StandardDeviationOfGrainRadii(n,dim)*dx;
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Number of grains of "       +locSettings.PhaseNames[n]         , NGrains     );
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Mean grain volume of "      +locSettings.PhaseNames[n]+" [m^3]", MeanVolume  );
        ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Stddev of grain volumes of "+locSettings.PhaseNames[n]+" [m^3]", StddevVolume);
        //ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Mean grain radius of "      +locSettings.PhaseNames[n]+" [m]"  , MeanRadius  );
        //ConsoleOutput::WriteWithLog(line , RTC.TimeStep , "Stddev of grain radii of "  +locSettings.PhaseNames[n]+" [m]"  , StddevRadius);
    }
    ConsoleOutput::WriteLineToLogfile(log, line, RTC.TimeStep);
}
}
