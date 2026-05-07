/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2025
 *  Main contributors :   Reza Namdar
 *
 */ 
#include "ReactiveFlows/Energy.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <limits>

#include <vector>
#include <random>
#include <iostream>

using namespace std;
using namespace openphase;

void Energy::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("Energy");
    ConsoleOutput::WriteStandard("Source", InputFile);
    std::fstream inp(InputFile, std::ios::in);
    if (!inp)
    {
        std::stringstream message;
        message << "File \"" << InputFile << "\" could not be opened";
        throw std::runtime_error(message.str());
    };
    std::stringstream inp_data;
    inp_data << inp.rdbuf();
    inp.close();
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "Energy");
    Pr					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Pr"), false, 0.71);
    Cp					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Cp"), false, 1005.0);
    Mw					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Mw"), false, 0.02896);
    REACTION          = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("REACTION"), false, false);
    ENERGY            = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("ENERGY"), false, false);
    CONJUGATE		 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("CONJUGATE"), false, false);
    Ts0					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Ts0"), false, 300.0);
    Cps					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Cps"), false, 900.0);
    Rhos				 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Rhos"), false, 5000.0);
    Ks					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Ks"), false, 237);
    T0					 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("T0"), false, 300.0);
    ATstar  			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ATstar"), false, 100.0);
    TempSolid	    	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempSolid"), false, 300.0);
    TempSolidCold     	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempSolidCold"), false, 300.0);
    HeatFlux	    	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("HeatFlux"), false, 300.0);
    TMu0              	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TMu0"), false, 273.0);
    SMu0     	         = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("SMu0"), false, 110.5);
    Mu0     	         = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Mu0"), false, 1.68e-05);
    BCOrder              = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("BCOrder"), false, 0);
    std::string schemeString = FileInterface::ReadParameterK(inp_data, moduleLocation, "AdvScheme", false, "UPWIND");

    if (schemeString == "UPWIND")
    {
        AdvScheme = AdvectionSchemes::Upwind;
    }
    else if (schemeString == "CENTRAL")
    {
        AdvScheme = AdvectionSchemes::Central;
    }
    else if (schemeString == "VANLEER")
    {
        AdvScheme = AdvectionSchemes::VanLeer;
    }
    else
    {
        ConsoleOutput::WriteWarning("Wrong advection scheme specified!\nThe default \"Upwind\" scheme is used!", "Energy", "ReadInput()");
    }
}

void Energy::WriteVTKHHR(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"omega_T (W/m^3)", [this](int i,int j,int k){return HRR(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"HHR_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void Energy::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    size_t Bcells = locSettings.Grid.Bcells;
    T.Allocate     (locSettings.Grid, Bcells);
    TOld.Allocate  (locSettings.Grid, Bcells);
    CpMix.Allocate (locSettings.Grid, Bcells);
    KMix.Allocate  (locSettings.Grid, Bcells);

    if(CONJUGATE)
    {
        Ts.Allocate            (locSettings.Grid, Bcells);
        TsOld.Allocate         (locSettings.Grid, Bcells);
        CpSol.Allocate      (locSettings.Grid, Bcells);
        KSol.Allocate       (locSettings.Grid, Bcells);
        DensitySol.Allocate (locSettings.Grid, Bcells);
    }
    if(REACTION)
    {
        HRR.Allocate (locSettings.Grid, Bcells);
    }
}

void Energy::UpdateFields()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,T,T.Bcells(),)
    {
        TOld(i,j,k)=T(i,j,k);
        if(CONJUGATE)
        {
            TsOld(i,j,k)=Ts(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}



void Energy::SetOpenBoundary(const dVector3& A, const dVector3& B, int n)
{
    if (n != 1 && n != -1) return;

    const int ia = A.getX(), ja = A.getY(), ka = A.getZ();
    const int ib = B.getX(), jb = B.getY(), kb = B.getZ();

    const bool xConst = (ia == ib);
    const bool yConst = (ja == jb);
    const bool zConst = (ka == kb);

    // Must be exactly one axis-aligned face
    if ((int)xConst + (int)yConst + (int)zConst != 1) return;

    const int Nx = Grid.Nx, Ny = Grid.Ny, Nz = Grid.Nz;
    const int G  = T.Bcells();               // number of ghost layers for T

    auto clamp = [](int v, int lo, int hi){ return std::max(lo, std::min(hi, v)); };

    auto owns_face = [&](int axis)->bool
    {
        if (axis == 0) return ( (n > 0 && Grid.OffsetX == 0) ||
                                (n < 0 && Grid.OffsetX + Nx == Grid.TotalNx) );
        if (axis == 1) return ( (n > 0 && Grid.OffsetY == 0) ||
                                (n < 0 && Grid.OffsetY + Ny == Grid.TotalNy) );
        return            ( (n > 0 && Grid.OffsetZ == 0) ||
                            (n < 0 && Grid.OffsetZ + Nz == Grid.TotalNz) );
    };

    const int axis = xConst ? 0 : (yConst ? 1 : 2);
    if (!owns_face(axis)) return;

    // Local face index (interior boundary plane index)
    const int faceLocal =
        (axis == 0) ? (ia - Grid.OffsetX) :
        (axis == 1) ? (ja - Grid.OffsetY) :
                      (ka - Grid.OffsetZ);

    // Must correspond to actual local boundary plane
    if (axis == 0)
    {
        if (faceLocal < 0 || faceLocal >= Nx) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nx-1)) return;
    }
    else if (axis == 1)
    {
        if (faceLocal < 0 || faceLocal >= Ny) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Ny-1)) return;
    }
    else
    {
        if (faceLocal < 0 || faceLocal >= Nz) return;
        if ((n > 0 && faceLocal != 0) || (n < 0 && faceLocal != Nz-1)) return;
    }

    // Patch bounds (LOCAL coordinates for tangential directions)
    const int la_i = ia - Grid.OffsetX, lb_i = ib - Grid.OffsetX;
    const int la_j = ja - Grid.OffsetY, lb_j = jb - Grid.OffsetY;
    const int la_k = ka - Grid.OffsetZ, lb_k = kb - Grid.OffsetZ;

    int iMin=0,iMax=0,jMin=0,jMax=0,kMin=0,kMax=0;

    if (axis == 0)
    {
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else if (axis == 1)
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        kMin = clamp(std::min(la_k, lb_k), 0, Nz-1);
        kMax = clamp(std::max(la_k, lb_k), 0, Nz-1);
    }
    else
    {
        iMin = clamp(std::min(la_i, lb_i), 0, Nx-1);
        iMax = clamp(std::max(la_i, lb_i), 0, Nx-1);
        jMin = clamp(std::min(la_j, lb_j), 0, Ny-1);
        jMax = clamp(std::max(la_j, lb_j), 0, Ny-1);
    }

    // ---- compute extrapolated temperature at one ghost node along normal ----
    auto extrapolate_T = [&](int i, int j, int k, int d)->double
    {
        // d = 1..G is the ghost-layer index (distance from boundary plane)
        // sample positions along normal inside the domain:
        // s1 = faceLocal + 0*n  (boundary interior)
        // s2 = faceLocal + 1*n  (1 cell inside)
        // ...
        // we use points: faceLocal + (d)*n, (d+1)*n, ...
        // But note: for ghost nodes we write at faceLocal - d*n (outside),
        // and read from faceLocal + (d-1)*n, faceLocal + d*n, ...
        // Easiest: define base = faceLocal + (d-1)*n (first interior used)
        int ii=i, jj=j, kk=k;

        auto T_at = [&](int off)->double
        {
            int ti=ii, tj=jj, tk=kk;
            if (axis == 0) ti = faceLocal + off*n;
            if (axis == 1) tj = faceLocal + off*n;
            if (axis == 2) tk = faceLocal + off*n;
            return T(ti, tj, tk);
        };

        // Need enough interior points for chosen BCOrder
        // Points we will use (offsets from face plane into domain): 0,1,2,3,...
        // BCOrder 0: uses offset 0
        // BCOrder 1: uses offsets 0,1
        // BCOrder 2: uses offsets 0,1,2
        // BCOrder 3: uses offsets 0..3 (4 points)
        // Check availability (local domain indices)
        auto inside_ok = [&](int off)->bool
        {
            int ti=0, tj=0, tk=0;
            if (axis == 0) { ti = faceLocal + off*n; tj = j; tk = k; }
            if (axis == 1) { tj = faceLocal + off*n; ti = i; tk = k; }
            if (axis == 2) { tk = faceLocal + off*n; ti = i; tj = j; }
            return (ti >= 0 && ti < Nx && tj >= 0 && tj < Ny && tk >= 0 && tk < Nz);
        };

        if (BCOrder == 0)
        {
            if (!inside_ok(0)) return T0;
            return T_at(0);
        }
        else if (BCOrder == 1)
        {
            if (!inside_ok(0) || !inside_ok(1)) return inside_ok(0) ? T_at(0) : T0;
            // linear extrapolation outward: T(-1) = 2*T(0) - T(1)
            return 2.0*T_at(0) - T_at(1);
        }
        else if (BCOrder == 2)
        {
            if (!inside_ok(0) || !inside_ok(1) || !inside_ok(2)) return inside_ok(0) ? T_at(0) : T0;
            // your quadratic form
            const double T0i = T_at(0);
            const double T1i = T_at(1);
            const double T2i = T_at(2);
            return T0i + (T0i - T1i) + 0.5*(T0i - 2.0*T1i + T2i);
        }
        else // BCOrder == 3
        {
            if (!inside_ok(0) || !inside_ok(1) || !inside_ok(2) || !inside_ok(3))
                return inside_ok(0) ? T_at(0) : T0;

            // Your Lagrange polynomial expects values at x=1..4 and evaluates at x=0.
            // Map: T(faceLocal+0*n) -> x=4? or x=1? In your old code you used:
            // Lagrange_polynomial4(0,1,2,3,4, T(iEnd-4),T(iEnd-3),T(iEnd-2),T(iEnd-1))
            // which corresponds to points at x=1..4 and extrapolate to x=0.
            // For consistency:
            // off=3 -> "x=1", off=2 -> "x=2", off=1 -> "x=3", off=0 -> "x=4"
            const double v1 = T_at(3);
            const double v2 = T_at(2);
            const double v3 = T_at(1);
            const double v4 = T_at(0);

            return Lagrange_polynomial4(0.0, 1.0, 2.0, 3.0, 4.0, v1, v2, v3, v4);
        }
    };

    // ---- write ghost layers outside the face: index = faceLocal - d*n ----
    auto write_ghost = [&](int d, int i, int j, int k)
    {
        int gi=i, gj=j, gk=k;
        if (axis == 0) gi = faceLocal - d*n;
        if (axis == 1) gj = faceLocal - d*n;
        if (axis == 2) gk = faceLocal - d*n;

        T(gi, gj, gk) = extrapolate_T(i, j, k, d);
    };

    // Loop patch and ghost depth
    if (axis == 0)
    {
        for (int d = 1; d <= G; ++d)
        for (int k = kMin; k <= kMax; ++k)
        for (int j = jMin; j <= jMax; ++j)
            write_ghost(d, /*i unused*/0, j, k);
        return;
    }

    if (axis == 1)
    {
        for (int d = 1; d <= G; ++d)
        for (int k = kMin; k <= kMax; ++k)
        for (int i = iMin; i <= iMax; ++i)
            write_ghost(d, i, /*j unused*/0, k);
        return;
    }

    // axis == 2
    for (int d = 1; d <= G; ++d)
    for (int j = jMin; j <= jMax; ++j)
    for (int i = iMin; i <= iMax; ++i)
        write_ghost(d, i, j, /*k unused*/0);
}

void Energy::SetBoundaryConditions(BoundaryConditions& BC)
{
    if (Grid.dNx > 0) BC.SetX(T);
    if (Grid.dNy > 0) BC.SetY(T);
    if (Grid.dNz > 0) BC.SetZ(T);

    if(CONJUGATE)
    {
        if (Grid.dNx > 0) BC.SetX(Ts);
        if (Grid.dNy > 0) BC.SetY(Ts);
        if (Grid.dNz > 0) BC.SetZ(Ts);
    }
}

double Energy::InterfaceTemperature(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB,
                                           const int i, const int j, const int k)
{
    constexpr double kTiny = 1e-30;
    const double D  = SB.nDist;
    const double dx = Grid.dx; 

    // ---- Find which solid phase this (solid) cell belongs to ----
    size_t solidIdx = 1; // your default
    {
        const auto& cellFields = Phase.Fields(i, j, k);
        for (auto it = cellFields.cbegin(); it != cellFields.cend(); ++it)
        {
            if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
            {
                solidIdx = it->index;
                break;
            }
        }
    }

    // ---- Find matching thermal surface condition for this solid phase ----
    const ThermalCondition* surfPtr = nullptr;
    for (const auto& surf : ThermalSurfaceCondition)
    {
        if (surf.solidIdx == solidIdx) { surfPtr = &surf; break; }
    }

    const ThermalCondition& surf = *surfPtr;

    // Constant temperature BC
    if (surf.type == ThermalCondition::Type::ConstantTemp)
        return surf.value;

    // Phase masks (consistent with your FL.Obstacle convention)
    auto keepFluidCorner = [&](int ii, int jj, int kk) -> bool { return !FL.Obstacle(ii, jj, kk); };
    auto keepSolidCorner = [&](int ii, int jj, int kk) -> bool { return  FL.Obstacle(ii, jj, kk); };

    // ---- Geometry: project solid cell center to interface, then sample +/- D along normal ----
    const dVector3 X  { double(i), double(j), double(k) };
    const dVector3 n0 = SB.NormField(i, j, k);          // assumed unit normal pointing solid->fluid
    const dVector3 Xb = X  - n0 * SB.DistanceField(i, j, k);   // interface point
    const dVector3 Xf = Xb + n0 * D;    // fluid-side sample point at distance D
    const dVector3 Xs = Xb - n0 * D;    // solid-side sample point at distance D

    const dVector3 Xfb = Xb + n0 * SB.nFallBack ;    
    const double   Tfb  = T(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
    const double   Kfb  = KMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

    const double Tf = SB.MaskedTrilinearScalar(Xf.getX(), Xf.getY(), Xf.getZ(),
        [&](int ii, int jj, int kk) { return T(ii, jj, kk); }, keepFluidCorner,Tfb);

    const double Kf = SB.MaskedTrilinearScalar(Xf.getX(), Xf.getY(), Xf.getZ(),
        [&](int ii, int jj, int kk) { return KMix(ii, jj, kk); }, keepFluidCorner,Kfb);

    // Constant heat flux BC at interface: Tw = Tf - D * q"/Kf
    if (surf.type == ThermalCondition::Type::ConstantFlux)
    {
        const double Qb = surf.value;
        return Tf - D * dx * Qb / std::max(Kf, kTiny);
    }

    // Conjugate heat transfer: sample solid-side and compute interface temperature by flux continuity
    if (surf.type == ThermalCondition::Type::Conjugate)
    {
        const double   Tsfb  = Ts(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
        const double   Ksfb  = KSol(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
        const double Ts1 = SB.MaskedTrilinearScalar(
            Xs.getX(), Xs.getY(), Xs.getZ(),
            [&](int ii, int jj, int kk) { return Ts(ii, jj, kk); },
            keepSolidCorner,Tsfb);

        const double Ks1 = SB.MaskedTrilinearScalar(
            Xs.getX(), Xs.getY(), Xs.getZ(),
            [&](int ii, int jj, int kk) { return KSol(ii, jj, kk); },
            keepSolidCorner,Ksfb);

        const double denom = std::max(Kf + Ks1, kTiny);
        return (Kf * Tf + Ks1 * Ts1) / denom;
    }
    
    return surf.value; // fallback, should not reach here if all types handled
}

double Energy::SurfaceTemperature(PhaseField& Phase, FlowSolverLBM& FL, SolidBody& SB,
                                            const dVector3& Xb, const dVector3& n, const string time)
{
    constexpr double kTiny = 1e-30;
    const double D  = SB.nDist;
    const double dx = Grid.dx;

    // --- Helper: find which solid phase the solid cell belongs to ---
    auto detectSolidPhaseIndex = [&](int i, int j, int k) -> size_t
    {
        size_t solidIdx = 1; // your original default
        const auto& cellFields = Phase.Fields(i, j, k);
        for (auto it = cellFields.cbegin(); it != cellFields.cend(); ++it)
        {
            if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
            {
                solidIdx = it->index;
                break;
            }
        }
        return solidIdx;
    };
    const dVector3 Xs = Xb - n * 1.6;  // solid-side sample

    // --- Compute solid phase index from Xs cell ---
    const int is = (int)Xs.getX();
    const int js = (int)Xs.getY();
    const int ks = (int)Xs.getZ();
    const size_t solidIdx = detectSolidPhaseIndex(is, js, ks);
    // --- Find matching thermal surface condition for this solid phase ---
    const ThermalCondition* surfPtr = nullptr;
    for (const auto& surf : ThermalSurfaceCondition)
    {
        if (surf.solidIdx == solidIdx) { surfPtr = &surf; break; }
    }
    const ThermalCondition& surf = *surfPtr;

    // Constant temperature BC
    if (surf.type == ThermalCondition::Type::ConstantTemp)
        return surf.value;

    // --- Interface point and sampling points ---
    const dVector3 Xf1 = Xb + n * D;  // fluid-side sample (n points solid->fluid)
    const dVector3 Xs1 = Xb - n * D;  // solid-side sample
    
    const dVector3 Xfb = Xb + n * SB.nFallBack ;    
    const double   TfbO  = TOld(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
    const double   Tfb  = T(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
    const double   Kfb  = KMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

    // Masks: adjust if your FL.Obstacle convention differs
    auto keepFluidCorner = [&](int i, int j, int k) -> bool 
    {
        return !FL.Obstacle(i, j, k);
    };
    auto keepSolidCorner = [&](int i, int j, int k) -> bool 
    {
        return FL.Obstacle(i, j, k);
    };

    // --- Fluid-side properties ---
    double Tf1 = surf.value;
    if(time=="old")
    {
        Tf1 = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(), 
                                       [&](int i, int j, int k) { return TOld(i, j, k); },
                                       keepFluidCorner,TfbO);
    }
    else
    {
        Tf1 = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(), 
                                       [&](int i, int j, int k) { return T(i,j,k); },
                                       keepFluidCorner,Tfb);
    }

    const double Kf1 = SB.MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                                       [&](int i, int j, int k) { return KMix(i, j, k); },
                                       keepFluidCorner,Kfb);

    // Constant heat flux BC: Tw = Tf1 - D * q"/k
    if (surf.type == ThermalCondition::Type::ConstantFlux)
    {
        const double Qb = surf.value;
        return Tf1 - D * dx * Qb / std::max(Kf1, kTiny);
    }

    // Conjugate: sample solid-side Ts, Ks and compute interface temperature
    if (surf.type == ThermalCondition::Type::Conjugate)
    {
        const double   TsfbO  = TsOld(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
        const double   Tsfb  = Ts(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 
        const double   Ksfb  = KSol(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

        // Fallback: solid cell center (Xs) (or use TsOld(is,js,ks) if guaranteed valid)
        double Ts1 = surf.value;
        if(time=="old")
        {
            Ts1 = SB.MaskedTrilinearScalar(Xs1.getX(), Xs1.getY(), Xs1.getZ(),
                                           [&](int i, int j, int k) { return TsOld(i, j, k); },
                                           keepSolidCorner,TsfbO);
        }
        else
        {
            Ts1 = SB.MaskedTrilinearScalar(Xs1.getX(), Xs1.getY(), Xs1.getZ(),
                                           [&](int i, int j, int k) { return Ts(i, j, k); },
                                           keepSolidCorner,Tsfb);
        }
        const double Ks1 = SB.MaskedTrilinearScalar(Xs1.getX(), Xs1.getY(), Xs1.getZ(),
                                           [&](int i, int j, int k) { return KSol(i, j, k); },
                                           keepSolidCorner,Ksfb);

        const double denom = std::max(Kf1 + Ks1, kTiny);
        return (Kf1 * Tf1 + Ks1 * Ts1) / denom;
    }
    // Unknown type: safe fallback
    //return Tf1;
    return surf.value; 
}

double Energy::CalculatePhiVanLeer(double r)
{
    return (r+fabs(r))/(1.0+fabs(r));
}

void Energy::CalculateThermodynamicPressure(FlowSolverLBM& FL, double P0, double T0)
{
    double SumTinv=0.0;
    double fcells=0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,T,0,)
    {
        if(!FL.Obstacle(i,j,k))
        {
            SumTinv += 1.0/T(i,j,k);
            fcells++; 
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    #ifdef MPI_PARALLEL
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &SumTinv, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &fcells, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#endif
    FL.Pth = P0/T0*fcells/SumTinv;    
}

void Energy::WriteVTKTemperature(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Temperature", [this](int i,int j,int k){return T(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"Temperature_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void Energy::WriteVTKSolidTemperature(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Temperature", [this](int i,int j,int k){return Ts(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"TemperatureSolid_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

double Energy::Lagrange_polynomial4(double T1, double T2, double T3, double T4, double T5,
                                          double f1, double f2, double f3, double f4)
{
    return ((T5 - T2) * (T5 - T3) * (T5 - T4) * f1) / (-6.0) +
           ((T5 - T1) * (T5 - T3) * (T5 - T4) * f2) / (2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T4) * f3) / (-2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T3) * f4) / (6.0);
}

void Energy::SetSolidPhaseTemp(PhaseField& Phase, bool DI)
{
    double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,T,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t PhaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.FieldsProperties[PhaseIdx].State == AggregateStates::Solid)
            {
                if(Phase.Fractions(i,j,k)[PhaseIdx]>=CRI)
                {
                    size_t PFIdx = it->index;
                    for (size_t is = 0; is < ThermalSurfaceCondition.size(); is++)
                    {
                        if(ThermalSurfaceCondition.at(is).solidIdx==PFIdx)
                        {
                            if(ThermalSurfaceCondition.at(is).type==ThermalCondition::Type::ConstantTemp)
                            {
                                T(i,j,k) = ThermalSurfaceCondition.at(is).value;
                            }
                            else if(ThermalSurfaceCondition.at(is).type==ThermalCondition::Type::ConstantFlux)
                            {
                                T(i,j,k) = ThermalSurfaceCondition.at(is).visualTemp;
                            }
                            else if(ThermalSurfaceCondition.at(is).type==ThermalCondition::Type::Conjugate)
                            {
                                T(i,j,k) = Ts(i,j,k);
                            }
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Energy::SetInitialNeighboringFluidPoints(PhaseField& Phase, SolidBody &SB, bool DI)
{
    //double CRI = (DI)?(1.0-DBL_EPSILON):(0.5-DBL_EPSILON);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,T,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
        {
            size_t PhaseIdx =  Phase.FieldsProperties[it->index].Phase;
            if(Phase.FieldsProperties[PhaseIdx].State == AggregateStates::Gas)
            {
                auto keepCorner = [&](int ii, int jj, int kk) -> bool { return true; };
                // ---- Geometry: project solid cell center to interface, then sample +/- D along normal ----
                const dVector3 X  { double(i), double(j), double(k) };
                const dVector3 n0  = SB.NormField(i, j, k);          // assumed unit normal pointing solid->fluid
                const dVector3 Xb  = X  - n0 * SB.DistanceField(i, j, k);   // interface point
                const dVector3 Xsi = Xb - n0 * Grid.iWidth;    // solid-side sample point at distance D
                const dVector3 Xfb = Xb + n0 * SB.nFallBack ;    
                const double   Tfb = T(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

                const double Tsi = SB.MaskedTrilinearScalar(Xsi.getX(), Xsi.getY(), Xsi.getZ(),
                    [&](int ii, int jj, int kk) { return T(ii, jj, kk); }, keepCorner,Tfb);

                T(i,j,k) = (T(i,j,k) + Tsi)/2.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Energy::SetInitial(PhaseField& Phase, FlowSolverLBM& FL, MixtureFlow &MF, SolidBody &SB)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, T, T.Bcells(),)
    {
        T(i,j,k) = T0;
        CpMix(i,j,k) = Cp;
        if(CONJUGATE)
        {
            Ts(i,j,k) = T0;
            CpSol(i,j,k)=Cps;
            KSol(i,j,k)=Ks;
            DensitySol(i,j,k)=Rhos;
            if(FL.Obstacle(i,j,k))
            {
                Ts(i,j,k) = Ts0;
            }
        }  
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double Energy::CalculateSutherlandViscosity(double Temp)
{
    double mu=Mu0*(Temp/TMu0)*sqrt(Temp/TMu0)*(TMu0+SMu0)/(Temp+SMu0);
    return mu;
}

double Energy::CalculateIdealGasDensity(double p, double Rm, double Temp)
{
    double rho= p/(Rm * Temp);
    return rho;
}

void Energy::CalculateProperties(PhaseField &Phase, FlowSolverLBM& FL, SolidBody &SB, MixtureFlow &MF, BoundaryConditions &BC) 
{
    double R = 8.314462618; // J/(mol*K)
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, MF.DensityMix, MF.DensityMix.Bcells(),)
    {
        if(!FL.Obstacle(i,j,k))
        {
            MF.DensityMix(i,j,k)   = CalculateIdealGasDensity(FL.Pth, R/Mw, T(i,j,k));
            MF.ViscosityMix(i,j,k) = CalculateSutherlandViscosity(T(i,j,k))/MF.DensityMix(i,j,k);
            KMix(i,j,k) = MF.DensityMix(i,j,k) * MF.ViscosityMix(i,j,k) * CpMix(i,j,k) / Pr;
        }
        else if(Phase.Fields(i,j,k).interface())
        {
            double Tint = InterfaceTemperature(Phase,FL,SB,i,j,k);
            MF.DensityMix(i,j,k)=CalculateIdealGasDensity(FL.Pth, R/Mw, Tint);
            MF.ViscosityMix(i,j,k)=CalculateSutherlandViscosity(Tint)/MF.DensityMix(i,j,k);
            KMix(i,j,k) = MF.DensityMix(i,j,k) * MF.ViscosityMix(i,j,k) * CpMix(i,j,k) / Pr;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    if (Grid.dNx) BC.SetX(MF.DensityMix);
    if (Grid.dNy) BC.SetY(MF.DensityMix);
    if (Grid.dNz) BC.SetZ(MF.DensityMix);

	if (Grid.dNx) BC.SetX(MF.ViscosityMix);
    if (Grid.dNy) BC.SetY(MF.ViscosityMix);
    if (Grid.dNz) BC.SetZ(MF.ViscosityMix);

	if (Grid.dNx) BC.SetX(CpMix);
    if (Grid.dNy) BC.SetY(CpMix);
    if (Grid.dNz) BC.SetZ(CpMix);

	if (Grid.dNx) BC.SetX(KMix);
    if (Grid.dNy) BC.SetY(KMix);
    if (Grid.dNz) BC.SetZ(KMix);
}
