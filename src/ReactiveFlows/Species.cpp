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
 
#include "ReactiveFlows/Species.h"
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

void Species::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("Species");
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
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "Species");
    BCOrder    	         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("BCOrder"), false, 0);
    TempBurntGas		 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("TempBurntGas"), false, 0.0);
    X0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0BurntZone"), false, -10);
    XNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNBurntZone"), false, -10);
    Y0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0BurntZone"), false, -10);
    YNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNBurntZone"), false, -10);
    Z0BurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0BurntZone"), false, -10);
    ZNBurntZone          = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNBurntZone"), false, -10);
    X0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0FuelZone"), false, -10);
    XNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNFuelZone"), false, -10);
    Y0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0FuelZone"), false, -10);
    YNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNFuelZone"), false, -10);
    Z0FuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0FuelZone"), false, -10);
    ZNFuelZone           = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNFuelZone"), false, -10);
    nCoeffs    	         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("nCoeffs"), false, 15);
}

void Species::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    size_t Bcells = locSettings.Grid.Bcells;

    nSpecies = locSettings.Ncomp;
    CpSp.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    DSp.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    WSp.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    hSp.Allocate (locSettings.Grid,{nSpecies}, Bcells);
    MassFrac.Allocate(locSettings.Grid,{nSpecies}, Bcells);
    MassFracOld.Allocate(locSettings.Grid,{nSpecies}, Bcells);    
    MwSp.Allocate({nSpecies});
    PolyNomCoeffs.Allocate({nSpecies, nCoeffs});
}

void Species::UpdateFields()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN (i,j,k,MassFrac,MassFrac.Bcells(),)
    {
        MassFracOld(i,j,k) = MassFrac(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Species::CalculateCpSpAndhSp(Energy &EN, FlowSolverLBM& FL)
{
	double R = 8.314510;   // J/mol.K
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CpSp,CpSp.Bcells(),)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			double T = EN.T(i,j,k);
			for (size_t ic = 0 ; ic < nSpecies; ic++)
    		{
				if(T < PolyNomCoeffs({ic,0}))
				{
					CpSp(i,j,k,{ic}) = (PolyNomCoeffs({ic,8}) + PolyNomCoeffs({ic,9})*T + PolyNomCoeffs({ic,10})*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,11})*pow(T,3.0) + PolyNomCoeffs({ic,12})*pow(T,4.0) ) 
                                            * R / MwSp({ic});
					hSp(i,j,k,{ic})  = (PolyNomCoeffs({ic,8}) + PolyNomCoeffs({ic,9})/2.0*T + PolyNomCoeffs({ic,10})/3.0*pow(T,2.0)
                                            + PolyNomCoeffs({ic,11})/4.0*pow(T,3.0) + PolyNomCoeffs({ic,12})/5.0*pow(T,4.0) + PolyNomCoeffs({ic,13})/T) 
                                            * R * T / MwSp({ic});
				}
				else
				{
					CpSp(i,j,k,{ic}) = (PolyNomCoeffs({ic,1}) + PolyNomCoeffs({ic,2})*T + PolyNomCoeffs({ic,3})*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,4})*pow(T,3.0) + PolyNomCoeffs({ic,5})*pow(T,4.0) ) 
                                            * R / MwSp({ic});
					hSp(i,j,k,{ic})  = (PolyNomCoeffs({ic,1}) + PolyNomCoeffs({ic,2})/2.0*T + PolyNomCoeffs({ic,3})/3.0*pow(T,2.0) 
                                            + PolyNomCoeffs({ic,4})/4.0*pow(T,3.0) + PolyNomCoeffs({ic,5})/5.0*pow(T,4.0) + PolyNomCoeffs({ic,6})/T)
                                            * R * T / MwSp({ic});
				}
				
    		}
		}
	}
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Species::SetOpenBoundary(const dVector3& A, const dVector3& B, int n)              
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
    const int G = MassFrac.Bcells();

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

    // ---- access helpers for "interior samples" at offsets 0..3 along the normal ----
    auto inside_ok = [&](int i, int j, int k, int off)->bool
    {
        int ti=i, tj=j, tk=k;
        if (axis == 0) ti = faceLocal + off*n;
        if (axis == 1) tj = faceLocal + off*n;
        if (axis == 2) tk = faceLocal + off*n;
        return (ti >= 0 && ti < Nx && tj >= 0 && tj < Ny && tk >= 0 && tk < Nz);
    };

    auto Y_at = [&](int i, int j, int k, size_t ic, int off)->double
    {
        int ti=i, tj=j, tk=k;
        if (axis == 0) ti = faceLocal + off*n;
        if (axis == 1) tj = faceLocal + off*n;
        if (axis == 2) tk = faceLocal + off*n;
        return MassFrac(ti, tj, tk, {ic});
    };

    // ---- extrapolate one scalar using your BCOrder rules ----
    auto extrap0_1_2 = [&](double v0, double v1, double v2)->double
    {
        // BCOrder==2 formula you used: v0 + (v0-v1) + 0.5*(v0 - 2v1 + v2)
        return v0 + (v0 - v1) + 0.5*(v0 - 2.0*v1 + v2);
    };

    // ---- set one ghost node at depth d (d=1..G) for given tangential (i,j,k) ----
    auto write_ghost = [&](int d, int i, int j, int k)
    {
        // ghost coordinates (outside domain): faceLocal - d*n
        int gi=i, gj=j, gk=k;
        if (axis == 0) gi = faceLocal - d*n;
        if (axis == 1) gj = faceLocal - d*n;
        if (axis == 2) gk = faceLocal - d*n;

        // Need enough interior points for chosen order
        const int need = (BCOrder == 0 ? 1 : (BCOrder == 1 ? 2 : (BCOrder == 2 ? 3 : 4)));
        for (int off = 0; off < need; ++off)
        {
            if (!inside_ok(i,j,k,off))
            {
                // fallback: copy nearest interior if possible, else EN.T0 and do nothing exotic
                for (size_t ic = 0; ic < nSpecies; ++ic)
                    MassFrac(gi,gj,gk,{ic}) = inside_ok(i,j,k,0) ? Y_at(i,j,k,ic,0) : 0.0;
                return;
            }
        }

        // ---- Species extrapolation ----
        for (size_t ic = 0; ic < nSpecies; ++ic)
        {
            double Yg = 0.0;

            if (BCOrder == 0)
            {
                Yg = Y_at(i,j,k,ic,0);
            }
            else if (BCOrder == 1)
            {
                Yg = 2.0*Y_at(i,j,k,ic,0) - Y_at(i,j,k,ic,1);
            }
            else if (BCOrder == 2)
            {
                const double y0 = Y_at(i,j,k,ic,0);
                const double y1 = Y_at(i,j,k,ic,1);
                const double y2 = Y_at(i,j,k,ic,2);
                Yg = extrap0_1_2(y0, y1, y2);
            }
            else // BCOrder == 3
            {
                const double v1 = Y_at(i,j,k,ic,3);
                const double v2 = Y_at(i,j,k,ic,2);
                const double v3 = Y_at(i,j,k,ic,1);
                const double v4 = Y_at(i,j,k,ic,0);
                Yg = Lagrange_polynomial4(0.0, 1.0, 2.0, 3.0, 4.0, v1, v2, v3, v4);
            }

            MassFrac(gi,gj,gk,{ic}) = Yg;
        }
    };

    // ---- Loop patch + ghost depth ----
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

double Species::Lagrange_polynomial4(double T1, double T2, double T3, double T4, double T5,
                                          double f1, double f2, double f3, double f4)
{
    return ((T5 - T2) * (T5 - T3) * (T5 - T4) * f1) / (-6.0) +
           ((T5 - T1) * (T5 - T3) * (T5 - T4) * f2) / (2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T4) * f3) / (-2.0) +
           ((T5 - T1) * (T5 - T2) * (T5 - T3) * f4) / (6.0);
}

void Species::SetBoundaryConditions(BoundaryConditions& BC)
{
    if (Grid.dNx > 0) BC.SetX(MassFrac);
    if (Grid.dNy > 0) BC.SetY(MassFrac);
    if (Grid.dNz > 0) BC.SetZ(MassFrac);
}

double Species::CalculateMeanMolarMass(const int i, const int j, const int k, string whichtime)
{
    double MMWi=0.0;
    for(size_t n = 0; n < nSpecies; n++)
    {
        if (whichtime == "old")
        {
            MMWi += MassFracOld(i,j,k,{n}) / MwSp({n});
        }
        if (whichtime == "new")
        {
            MMWi += MassFrac(i,j,k,{n}) / MwSp({n});
        }
    }
    double MMW=1.0/MMWi;
    return MMW;
}

double Species::MoleFraction( const int i, const int j, const int k, double MMW, size_t kc)
{
    return ( MMW/MwSp({kc}) *  MassFracOld(i,j,k,{kc}) );
}

Tensor<double, 1> Species::InterfaceSpecies(FlowSolverLBM& FL, SolidBody& SB, const int i, const int j, const int k)
{
    //constexpr double kTiny = 1e-30;
    const double D = SB.nDist;

    // ---- Mask: fluid corners only ----
    auto keepFluidCorner = [&](int ii, int jj, int kk) -> bool
    {
        return !FL.Obstacle(ii, jj, kk);
    };

    // ---- Geometry: project to interface then move +D along normal to fluid side ----
    const dVector3 X  { double(i), double(j), double(k) };
    const dVector3 n0 = SB.NormField(i, j, k);       // assumed unit normal pointing solid->fluid

    const dVector3 Xb  = X  - n0 * SB.DistanceField(i, j, k);  // interface point
    const dVector3 Xf = Xb + n0 * D;    // fluid-side sample point at distance D

    const dVector3 Xfb = Xb + n0 * SB.nFallBack ;    
    const Tensor<double ,1> Yfb = MassFrac(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

    return SB.MaskedTrilinearTensor(Xf.getX(), Xf.getY(), Xf.getZ(),
                                    [&](int ii, int jj, int kk) { return MassFrac(ii, jj, kk); },
                                    keepFluidCorner, nSpecies, Yfb);
}

Tensor<double, 1> Species::SurfaceSpecies(FlowSolverLBM& FL, SolidBody& SB,
                                                    const dVector3& Xb, const dVector3& n, const string time)
{
    //constexpr double kTiny = 1e-30;
    const double D  = SB.nDist;
    //const double dx = Grid.dx;

    const dVector3 Xf1 = Xb + n * D;               // fluid-side normal sample

    const dVector3 Xfb = Xb + n * SB.nFallBack ;    
    const Tensor<double ,1> Yfb = MassFrac(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ())); 

    auto keepFluidCorner = [&](int i, int j, int k) -> bool 
    {
        return !FL.Obstacle(i, j, k);
    };

    // Fluid-side Y at Xf1
    Tensor<double,1> Yf1({nSpecies});
    if(time=="old")
    {
        Yf1 = SB.MaskedTrilinearTensor(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                                                    [&](int i,int j,int k){ return MassFracOld(i,j,k); },
                                                    keepFluidCorner, nSpecies, Yfb);
    }
    else
    {
        Yf1 = SB.MaskedTrilinearTensor(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                                                    [&](int i,int j,int k){ return MassFrac(i,j,k); },
                                                    keepFluidCorner, nSpecies, Yfb);
    }

    return Yf1;
}

void Species::WriteVTK(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    for(size_t comp = 0; comp < nSpecies; comp++)
    {
        const std::string nameComp1 = "MassFraction_" + locSettings.ElementNames[comp];
        ListOfFields.push_back((VTK::Field_t) {nameComp1, [this, comp](int i,int j,int k){return MassFrac(i,j,k,{comp});}});
    }
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"MassFractions_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void Species::SetInitial(PhaseField& Phase, FlowSolverLBM& FL, Energy &EN,
                            const Tensor<double, 1>& BM, const Tensor<double, 1>& FM , const Tensor<double, 1>& AM, 
                            double TB, double TF, double TA)
{
    //Simple composition set up, uniform component fractions per phase.
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, MassFrac, MassFrac.Bcells(),)
    {
        EN.HRR(i,j,k) = 0.0;
        if( (i+Phase.Grid.OffsetX)>=X0FuelZone and (i+Phase.Grid.OffsetX)<XNFuelZone and
                 (j+Phase.Grid.OffsetY)>=Y0FuelZone and (j+Phase.Grid.OffsetY)<YNFuelZone and
                 (k+Phase.Grid.OffsetZ)>=Z0FuelZone and (k+Phase.Grid.OffsetZ)<ZNFuelZone)
        {
            EN.T(i,j,k) = TF;
            MassFrac(i,j,k) = FM;
        }
        else if( (i+Phase.Grid.OffsetX)>=X0BurntZone and (i+Phase.Grid.OffsetX)<XNBurntZone and
            (j+Phase.Grid.OffsetY)>=Y0BurntZone and (j+Phase.Grid.OffsetY)<YNBurntZone and
            (k+Phase.Grid.OffsetZ)>=Z0BurntZone and (k+Phase.Grid.OffsetZ)<ZNBurntZone)
        {
            EN.T(i,j,k) = TB;
            MassFrac(i,j,k) = BM;
        }
        else 
        {
            EN.T(i,j,k) = TA;
            MassFrac(i,j,k) = AM;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
