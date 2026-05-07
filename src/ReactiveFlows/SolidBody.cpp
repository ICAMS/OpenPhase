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
 
#include "ReactiveFlows/SolidBody.h"
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

void SolidBody::CalculateDistancePeriodic(dVector3 A, dVector3 B,
                dVector3 &dist, double Nx, double Ny, double Nz)
{
    dist[0] = A[0]-B[0];
    if (std::abs(dist[0]) > std::abs(dist[0]+Nx)) dist[0] += Nx;
    else if (std::abs(dist[0]) > std::abs(dist[0]-Nx)) dist[0] -= Nx;
    dist[1] = A[1]-B[1];
    if (std::abs(dist[1]) > std::abs(dist[1]+Ny)) dist[1] += Ny;
    else if (std::abs(dist[1]) > std::abs(dist[1]-Ny)) dist[1] -= Ny;
    dist[2] = A[2]-B[2];
    if (std::abs(dist[2]) > std::abs(dist[2]+Nz)) dist[2] += Nz;
    else if (std::abs(dist[2]) > std::abs(dist[2]-Nz)) dist[2] -= Nz;
}

void SolidBody::ReadInput(string InputFile)
{
    ConsoleOutput::WriteBlankLine();
    ConsoleOutput::WriteLineInsert("SolidBody");
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
    int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "SolidBody");
    nParticles		     = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("nParticles"), false, 0);
    nRows		         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("nRows"), false, 0);
    FirstRow		     = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("FirstRow"), false, 0);
    nCols		         = FileInterface::ReadParameterI(inp_data, moduleLocation, std::string("nCols"), false, 0);
    PartDiameter    	 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("PartDiameter"), false, 0.01);
    Porosity			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Porosity"), false, 0.5);
    Clearance			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Clearance"), false, 5);
    nDist    			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("nDist"), false, 1);
    nFallBack    	     = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("nFallBack"), false, 1);

    Do_randDist    	 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Random_Dist"), false, false);
    X0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0DistZone"), false, -10);
    XNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNDistZone"), false, -10);
    Y0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0DistZone"), false, -10);
    YNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNDistZone"), false, -10);
    Z0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0DistZone"), false, -10);
    ZNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNDistZone"), false, -10);

    MinRadious       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MinRadious"), false,  10);
    MaxRadious       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MaxRadious"), false,  20);
    IntProf          = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("IntProf"),false, "SIN");

    for (size_t i = 0; i < nParticles; ++i)
    {
        rand_Bodies.push_back({0.0, 0.0, 0.0, 0.0});
    }
}

void SolidBody::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
    size_t Bcells = locSettings.Grid.Bcells;
    DistanceField.Allocate(Grid, Bcells);
    NormField.Allocate(Grid, Bcells);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DistanceField,DistanceField.Bcells(),)
    {
        DistanceField(i,j,k) = 0.0;
        NormField(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void SolidBody::CalculateDistanceField(PhaseField& Phase)
{
    const double dx = Phase.Grid.dx;

    // -------------------------------------------------------------------------
    // 1) Gas phase-field index
    // -------------------------------------------------------------------------
    size_t Gasidx = 0;
    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); ++idx)
    {
        if (Phase.FieldsProperties[idx].State == AggregateStates::Gas)
        {
            Gasidx = idx;
            break;
        }
    }

    // -------------------------------------------------------------------------
    // 2) Grid + constants
    // -------------------------------------------------------------------------
    const int Nx = Phase.Grid.Nx;
    const int Ny = Phase.Grid.Ny;
    const int Nz = Phase.Grid.Nz;

    const double phi_c    = 0.5;
    const double eta      = Phase.Grid.iWidth;
    const double maxD     = eta / 2.0 + 3.0;   // physical length
    const double maxDg    = maxD;         // grid units (DistanceField units)

    const double grad_eps = 1e-12;
    const double hess_eps = 1e-14;

    // Edge band: treat these as "problem points" and fix by shifted interpolation
    const double phi_min = 0.00000000005;
    const double phi_max = 0.99999999995;

    //const double d_max = eta / 2.0 - 0.1;

    // Shift distances in grid cells (you suggested 2 or 2.5)
    const double sstep1 = 2.4;
    //const double sstep2 = 2.5;

    auto clampi = [](int a, int lo, int hi) -> int
    {
        return (a < lo) ? lo : (a > hi ? hi : a);
    };

    auto sgn = [](double x) -> double
    {
        return (x > 0.0) ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
    };

    auto inCore = [&](double phi) -> bool
    {
        return (phi > phi_min && phi < phi_max);
    };

    //auto inCore = [&](double d) -> bool
    //{
    //    return (abs(d) < d_max);
    //};

    const double INF = 1e300;

    // We'll mark the edge points that need the shifted correction
    std::vector<unsigned char> needShift((size_t)Nx * (size_t)Ny * (size_t)Nz, 0);

    auto idx3D = [Ny, Nz](int i, int j, int k) -> size_t 
    {
        return (size_t)i * (size_t)Ny * (size_t)Nz
             + (size_t)j * (size_t)Nz
             + (size_t)k;
    };

    // -------------------------------------------------------------------------
    // 3) Pass A: compute accurate distances on wide_interface core only
    // -------------------------------------------------------------------------
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DistanceField, 0,)
    {
        const NodePF& node = Phase.Fields(i, j, k);
        const double  phi0 = node.get_value(Gasidx);
        

        // Default: bulk gets saturated sign distance
        if (!node.wide_interface())
        {
            DistanceField(i, j, k) = sgn(phi0 - phi_c) * maxDg;
        }
        else
        {
            // wide_interface but near pure phase -> postpone to shifted pass
            //if (node.interface())
            if (!inCore(phi0))
            //if (!inCore(DistanceField(i,j,k)))
            {
                DistanceField(i, j, k) = INF;
                needShift[idx3D(i, j, k)] = 1;
            }
            else
            {
                const dVector3 g     = node.get_gradient(Gasidx);
                const double   phi_n = g.length();

                if (phi_n < grad_eps)
                {
                    DistanceField(i, j, k) = 0.0;
                }
                else
                {
                    dVector3 n = g;
                    n.normalize();

                    double dphin_dx = 0.0, dphin_dy = 0.0, dphin_dz = 0.0;

                    for (auto gs = Phase.GStencil.cbegin(); gs != Phase.GStencil.cend(); ++gs)
                    {
                        int ii = clampi(i + gs->di, 0, Nx - 1);
                        int jj = clampi(j + gs->dj, 0, Ny - 1);
                        int kk = clampi(k + gs->dk, 0, Nz - 1);

                        const NodePF& nb = Phase.Fields(ii, jj, kk);
                        const double pn_nb = nb.get_gradient(Gasidx).length();

                        dphin_dx += gs->weightX * pn_nb;
                        dphin_dy += gs->weightY * pn_nb;
                        dphin_dz += gs->weightZ * pn_nb;
                    }

                    const double phi_nn = n[0]*dphin_dx + n[1]*dphin_dy + n[2]*dphin_dz;

                    const double C = (phi0 - phi_c);
                    const double A = 0.5 * phi_nn;
                    const double B = phi_n;

                    const double d_lin = (phi_c - phi0) / phi_n;
                    double d = d_lin;

                    const double disc = B*B - 4.0*A*C;

                    if (std::fabs(A) > hess_eps && disc >= 0.0)
                    {
                        const double sqrt_disc = std::sqrt(disc);
                        const double d1 = (-B + sqrt_disc) / (2.0 * A);
                        const double d2 = (-B - sqrt_disc) / (2.0 * A);

                        auto same_sign_or_zero = [&](double r) -> bool {
                            return (d_lin == 0.0) ? true : (r * d_lin >= 0.0);
                        };

                        const bool ok1 = same_sign_or_zero(d1);
                        const bool ok2 = same_sign_or_zero(d2);

                        if (ok1 && ok2)      d = (std::fabs(d1) < std::fabs(d2)) ? d1 : d2;
                        else if (ok1)        d = d1;
                        else if (ok2)        d = d2;
                    }

                    if (d >  maxD) d =  maxD;
                    if (d < -maxD) d = -maxD;

                    DistanceField(i, j, k) = -d / dx; // grid units
                   
                }
            }
        }

        // ---------------- NormField (unchanged) ----------------
        dVector3 nout = dVector3::ZeroVector();
        for (auto it = node.cbegin(); it != node.cend(); ++it)
        {
            if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
            {
                const size_t Solididx = it->index;
                const NodeAB<dVector3, dVector3> Normals = Phase.Normals(i, j, k);
                nout = Normals.get_asym1(Solididx, Gasidx);
                break;
            }
        }
        NormField(i, j, k) = nout;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DistanceField, 0,)
    {
        if (!needShift[idx3D(i, j, k)]) continue;
        
        const double phi0 = Phase.Fields(i, j, k).get_value(Gasidx);
        auto keepCorner = [&](int ii, int jj, int kk) -> bool
        {
            return true;
        };
        const dVector3 X  { double(i), double(j), double(k) };
        const dVector3 n0 = NormField(i, j, k);          // assumed unit normal pointing solid->fluid
        const dVector3 Xn = X  -  n0 * sstep1 * sgn(phi0 - phi_c);   // interface point
        const double Xnfb = sgn(phi0 - phi_c) * maxDg;

        double Dnew = MaskedTrilinearScalar(Xn.getX(), Xn.getY(), Xn.getZ(),
        [&](int ii, int jj, int kk) { return DistanceField(ii, jj, kk); },
        keepCorner,Xnfb);

        Dnew += sgn(phi0 - phi_c) * sstep1;

        if (abs(Dnew) > Grid.iWidth) Dnew = sgn(phi0 - phi_c) * maxDg;

        // Fallback: if still failing, use saturated sign distance

        DistanceField(i, j, k) = Dnew;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    // -------------------------------------------------------------------------
    // 4) Fill halos by clamping to nearest physical cell
    // -------------------------------------------------------------------------
    
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DistanceField, 0,)
    {
        const int ic = clampi((int)i, 0, Nx - 1);
        const int jc = clampi((int)j, 0, Ny - 1);
        const int kc = clampi((int)k, 0, Nz - 1);

        DistanceField(i, j, k) = DistanceField(ic, jc, kc);
        NormField(i, j, k)     = NormField(ic, jc, kc);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

}

double SolidBody::BilinearInterPolation(double xp, double zp, double qx1z1, double qx1z2, double qx2z1, double qx2z2)
{
    double x1 = floor(xp);
    double x2 = x1 + 1.0;
    double z1 = floor(zp);
    double z2 = z1 + 1.0;

    double qm = (x2 - xp) / (x2 - x1) * ((z2 - zp) / (z2 - z1) * qx1z1 + (zp - z1) / (z2 - z1) * qx1z2)
              + (xp - x1) / (x2 - x1) * ((z2 - zp) / (z2 - z1) * qx2z1 + (zp - z1) / (z2 - z1) * qx2z2);
    return qm;
}

double SolidBody::TrilinearInterPolation(double ic, double jc, double kc,  
                                                double q000, double q100, double q010, double q110, 
                                                double q001, double q101, double q011, double q111)
{
    double i0   = floor(ic);
    double i1   = i0 + 1.0;
    double j0   = floor(jc);
    double j1   = j0 + 1.0;
    double k0   = floor(kc);
    double k1   = k0 + 1.0;

    double CX = (ic-i0)/(i1-i0);
    double CY = (jc-j0)/(j1-j0);
    double CZ = (kc-k0)/(k1-k0);

    double q  = (1.0-CX)*(1.0-CY)*(1.0-CZ) * q000 + CX*(1.0-CY)*(1.0-CZ) * q100 + (1.0-CX)*CY*(1.0-CZ) * q010 + CX*CY*(1.0-CZ) * q110
               +(1.0-CX)*(1.0-CY)*CZ       * q001 + CX*(1.0-CY)*CZ       * q101 + (1.0-CX)*CY*CZ       * q011 + CX*CY*CZ       * q111;
    return q;
}

// Find 3 nearest neighbors (indices) for point i
std::array<int, 3> SolidBody::find_3_nearest(const std::vector<dVector3>& pts, size_t idx,BoundaryConditions& BC)
{
    std::array<double, 3> best_dist = {std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::max(),
                                       std::numeric_limits<double>::max()};
    std::array<int, 3> best_idx = {-1, -1, -1};

    for (size_t j = 0; j < pts.size(); ++j) {
        if (j == idx) continue;
        dVector3 dist = Tools::Distance<dVector3>(pts[idx], pts[j], Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, BC);
        double d = dist.abs();

        for (size_t k = 0; k < 3; ++k) {
            if (d < best_dist[k]) {
                for (size_t m = 2; m > k; --m) {
                    best_dist[m] = best_dist[m-1];
                    best_idx[m] = best_idx[m-1];
                }
                best_dist[k] = d;
                best_idx[k] = j;
                break;
            }
        }
    }
    return best_idx;
}

// Calculate area of quad as two triangles
double SolidBody::quad_area(const dVector3& a, const dVector3& b, const dVector3& c, const dVector3& d) 
{
    // Area = 0.5 * |(b-a) x (c-a)| + 0.5 * |(c-a) x (d-a)|
    dVector3 ab = b - a, ac = c - a, ad = d - a;
    dVector3 n1 = ab.cross(ac);
    dVector3 n2 = ac.cross(ad);
    return (n1.length() * 0.5 + n2.length()* 0.5);
}

// Centroid of 4 points
dVector3 SolidBody::quad_centroid(const dVector3& a, const dVector3& b, const dVector3& c, const dVector3& d) 
{
    return (a + b + c + d) / 4.0;
}

bool SolidBody::is_close(const dVector3& a, const dVector3& b, double tol) 
{
    for (int i = 0; i < 3; ++i)
        if (std::fabs(a[i] - b[i]) > tol) return false;
    return true;
}


double SolidBody::dist3(const dVector3& a, const dVector3& b) 
{
    return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2));
}

// Triangle area from 3 points
double SolidBody::triangle_area(const dVector3& a, const dVector3& b, const dVector3& c) 
{
    double ux = b[0]-a[0], uy = b[1]-a[1], uz = b[2]-a[2];
    double vx = c[0]-a[0], vy = c[1]-a[1], vz = c[2]-a[2];
    double cx = uy*vz - uz*vy;
    double cy = uz*vx - ux*vz;
    double cz = ux*vy - uy*vx;
    return 0.5*sqrt(cx*cx + cy*cy + cz*cz);
}

// For each point, estimate its local surface area
std::vector<double> SolidBody::surface_elements_3d(const std::vector<dVector3>& points, size_t kneigh)
{
    size_t N = points.size();
    std::vector<double> elem(N, 0.0);

    for (size_t i = 0; i < N; ++i)
    {
        const auto& pi = points[i];

        // Find k nearest neighbors
        std::vector<std::pair<double,int>> dists;
        for (size_t j = 0; j < N; ++j)
        {
            if (j == i) continue;
            double d = dist3(pi, points[j]);
            dists.push_back({d, j});
        }
        std::sort(dists.begin(), dists.end());
        std::vector<int> knn;
        for (size_t k = 0; k < kneigh && k < dists.size(); ++k)
            knn.push_back(dists[k].second);

        // For each neighbor pair, form a triangle
        double area = 0.0;
        for (size_t m = 0; m < knn.size(); ++m)
        {
            for (size_t n = m+1; n < knn.size(); ++n)
            {
                double triA = triangle_area(pi, points[knn[m]], points[knn[n]]);
                area += triA / 3.0; // assign 1/3 of triangle area to center point
            }
        }
        elem[i] = area;
    }
    return elem;
}

double SolidBody::curve_element_2d(const vector<dVector3>& points, size_t i, const BoundaryConditions& BC)
{
    const auto& pi = points[i];
    double min_dist1 = 1e30;
    size_t idx1 = 1e10;
    
    for(size_t j=0; j<points.size(); ++j)
    {
        if(j==i) continue;
        dVector3 dist = Tools::Distance<dVector3>(pi, points[j], Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, BC);
        double d = dist.abs();
        if(d < min_dist1) 
        {
            min_dist1 = d;
            idx1=j;
        } 
    }

    double min_dist2 = 5;
    size_t idx2 = 1e10;

    for(size_t k=0; k<points.size(); ++k)
    {
        if(k==i || k==idx1) continue;
        dVector3 v1 = Tools::Distance<dVector3>(points[k], points[idx1], Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, BC);
        dVector3 v2 = Tools::Distance<dVector3>(pi, points[idx1], Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, BC);
        double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        double norm2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
        if (norm2 < 1e-15) continue; // skip degenerate
        double t = dot / norm2;
        // t should be between 0 and 1 (not exactly 0 or 1)
        if (t > 0 && t < 1) 
        {
            dVector3 dist = Tools::Distance<dVector3>(pi, points[k], Grid.TotalNx, Grid.TotalNy, Grid.TotalNz, BC);
            double d = dist.abs();
            if(d < min_dist2) 
            {
                min_dist2 = d;
                idx2 = k;
            }
        }
    }
    double ds=0.0; 

    if (idx1 <=points.size() && idx2 <=points.size()) 
    {
        ds = 0.5 * (min_dist1 + min_dist2);
    } 
    else if (idx1 <=points.size()) 
    {
        ds = min_dist1;
    } 
    else
    {
        ds=0.0;
    }
    return ds;
}

void SolidBody::CalculateDragCoeff(PhaseField& Phase,
                                   MixtureFlow& MF,
                                   FlowSolverLBM& FL,
                                   double rhoRef,
                                   double Uref,
                                   double pRef,
                                   const dVector3& dragDir, // eD: flow / drag direction (unit)
                                   double dn,
                                   double Time,
                                   std::string dir,
                                   std::string fname,
                                   BoundaryConditions& BC)
{
    const double dx = Phase.Grid.dx;

    int rank = 0;
#ifdef MPI_PARALLEL
    rank = MPI_RANK;
#endif

    // Sampling distance along normal (lattice units)
    const double D = 1.5;

    // Unit drag/flow direction
    dVector3 eD = dragDir;
    {
        const double m2 = eD * eD;
        if (m2 > 0.0) eD = eD * (1.0 / std::sqrt(m2));
        else          eD = {1.0, 0.0, 0.0};
    }

    const double dynPress = 0.5 * rhoRef * Uref * Uref;
    size_t fc = 0;

    // Active directions (your convention)
    const int ax = (Phase.Grid.dNx == 1) ? 1 : 0;
    const int ay = (Phase.Grid.dNy == 1) ? 1 : 0;
    const int az = (Phase.Grid.dNz == 1) ? 1 : 0;
    const int nActive = ax + ay + az; // typically 2 or 3

    // Geometric element represented by one accepted sample:
    // 2D: dℓ  ~ dx    / max(|n| on active axes)
    // 3D: dA  ~ dx^2  / max(|n| on active axes)
    auto dS_element = [&](const dVector3& nUnit) -> double
    {
        const double eps = 1e-12;

        const double nx = ax ? std::abs(nUnit.getX()) : 0.0;
        const double ny = ay ? std::abs(nUnit.getY()) : 0.0;
        const double nz = az ? std::abs(nUnit.getZ()) : 0.0;

        double denom = 0.0;
        if (nActive == 2)
        {
            if (ax && ay) denom = std::max(nx, ny);
            else if (ax && az) denom = std::max(nx, nz);
            else               denom = std::max(ny, nz);
        }
        else // nActive==3 (or fallback)
        {
            denom = std::max(nx, std::max(ny, nz));
        }

        const double base = (nActive == 2) ? dx : (dx * dx);
        return base / (denom + eps);
    };

    // Reference projected measure per solid idx for this eD:
    // Aproj_ref = Σ ( max(0,-n·eD) * dS )
    static std::unordered_map<size_t, double>   AprojRef;
    static std::unordered_map<size_t, dVector3> eDRef;

    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); ++idx)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        std::vector<dVector3> bps;      // boundary points (for Cp output)
        std::vector<double>   CpLoc;    // local Cp

        // Integrated (local) contributions
        double sumFd_local    = 0.0; // Σ ( (t·eD) dS )
        double sumAproj_local = 0.0; // Σ ( w dS ), w=max(0,-n·eD)

        int nLocal = 0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, MF.VelocityMix, 0, )
        {
            const double X = DistanceField(i, j, k);

            if (X < dn && X > 0.0)
            {
                for (auto it = Phase.Fields(i, j, k).cbegin();
                          it != Phase.Fields(i, j, k).cend(); ++it)
                {
                    if (Phase.FieldsProperties[it->index].State == AggregateStates::Solid &&
                        it->index == idx)
                    {
                        auto keepFluidCorner = [&](int ii, int jj, int kk) -> bool
                        { return !FL.Obstacle(ii, jj, kk); };

                        const dVector3 Xi{double(i), double(j), double(k)};

                        // Normal (solid->fluid). Renormalize for safety.
                        dVector3 n0 = NormField(i, j, k);
                        {
                            const double n2 = n0 * n0;
                            if (n2 > 0.0) n0 = n0 * (1.0 / std::sqrt(n2));
                            else          n0 = {1.0, 0.0, 0.0};
                        }

                        // Interface + fluid sampling points
                        const dVector3 Xb  = Xi - n0 * X;
                        const dVector3 Xf1 = Xb + n0 * D;
                        const dVector3 Xf2 = Xb + n0 * (2.0 * D);

                        // Fallback point for masked interpolation
                        const dVector3 Xfb   = Xb + n0 * nFallBack;
                        const dVector3 ufb   = MF.VelocityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()));
                        const double   mufb  = MF.ViscosityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()));
                        const double   rhofb = MF.DensityMix(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()));
                        const double   pfb   = FL.HydroPressure(int(Xfb.getX()),int(Xfb.getY()),int(Xfb.getZ()),{fc});

                        const double s = D * dx;

                        const dVector3 u1 = MaskedTrilinearVector3(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return MF.VelocityMix(ii, jj, kk); },
                            keepFluidCorner, ufb);

                        const dVector3 u2 = MaskedTrilinearVector3(Xf2.getX(), Xf2.getY(), Xf2.getZ(),
                            [&](int ii, int jj, int kk) { return MF.VelocityMix(ii, jj, kk); },
                            keepFluidCorner, ufb);

                        // Keep your original definition: mu = ViscosityMix * DensityMix
                        const double mu1 = MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return MF.ViscosityMix(ii,jj,kk) * MF.DensityMix(ii,jj,kk); },
                            keepFluidCorner, mufb*rhofb);

                        const double mu2 = MaskedTrilinearScalar(Xf2.getX(), Xf2.getY(), Xf2.getZ(),
                            [&](int ii, int jj, int kk) { return MF.ViscosityMix(ii,jj,kk) * MF.DensityMix(ii,jj,kk); },
                            keepFluidCorner, mufb*rhofb);

                        const double p1 = MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return FL.HydroPressure(ii, jj, kk,{fc}); },
                            keepFluidCorner, pfb);

                        const double p2 = MaskedTrilinearScalar(Xf2.getX(), Xf2.getY(), Xf2.getZ(),
                            [&](int ii, int jj, int kk) { return FL.HydroPressure(ii, jj, kk,{fc}); },
                            keepFluidCorner, pfb);

                        // Wall velocity (stationary)
                        const dVector3 uw = {0.0, 0.0, 0.0};
                        const double   uwn = uw * n0;
                        const dVector3 uwt = uw - n0 * uwn;

                        // Tangential velocity at samples
                        const double   un1 = u1 * n0;
                        const dVector3 ut1 = u1 - n0 * un1;
                        const double   un2 = u2 * n0;
                        const dVector3 ut2 = u2 - n0 * un2;

                        // Your original one-sided stencil (kept)
                        const dVector3 dut = uwt * (-3.0) + ut1 * (4.0) + ut2 * (-1.0);
                        const double   dut_mag = dut.abs();

                        dVector3 t_visc = {0.0, 0.0, 0.0};
                        if (dut_mag > 0.0 && s > 1e-20)
                        {
                            const double dUtdn = dut_mag / (2.0 * s);
                            const double tau_w = 0.5 * (mu1 + mu2) * dUtdn;
                            const dVector3 t_hat = dut * (1.0 / dut_mag);
                            t_visc = t_hat * tau_w;
                        }

                        // Pressure traction (kept)
                        const double   pAvg    = 0.5 * (p1 + p2);
                        const dVector3 t_press = n0 * (-pAvg);

                        // Total traction projection along eD (force/area)
                        const double tD = (t_visc + t_press) * eD;

                        // Cp (kept)
                        const double Cp = (dynPress > 0.0) ? ((pAvg - pRef) / dynPress) : 0.0;

                        // Projected-area weight w = max(0, -n·eD) (kept)
                        const double ndot = n0 * eD;
                        const double w    = (ndot < 0.0) ? -ndot : 0.0;

                        // NEW: geometric element per sample (dℓ in 2D, dA in 3D)
                        const double dS = dS_element(n0);

                        // Integrate
                        sumFd_local    += tD * dS;
                        sumAproj_local += w  * dS;

                        // Store Cp samples (as you did)
                        bps.push_back(Xb);
                        CpLoc.push_back(Cp);

                        ++nLocal;

                        break; // done for this cell once idx matched
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // --- MPI: global sums (everyone participates)
        double sumFd    = sumFd_local;
        double sumAproj = sumAproj_local;
        int    nGlobal  = nLocal;

#ifdef MPI_PARALLEL
        {
            double gFd = 0.0, gAp = 0.0;
            int    gN  = 0;

            OP_MPI_Allreduce(&sumFd,    &gFd, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(&sumAproj, &gAp, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(&nGlobal,  &gN,  1, OP_MPI_INT,    OP_MPI_SUM, OP_MPI_COMM_WORLD);

            sumFd    = gFd;
            sumAproj = gAp;
            nGlobal  = gN;
        }
#endif

        // If no rank encountered this body, skip
        if (nGlobal <= 0) continue;

        // --- Choose ONE writer rank for this idx: the smallest rank with nLocal>0
        int writer = std::numeric_limits<int>::max();
#ifdef MPI_PARALLEL
        {
            const int candidate = (nLocal > 0) ? rank : std::numeric_limits<int>::max();
            OP_MPI_Allreduce(&candidate, &writer, 1, OP_MPI_INT, OP_MPI_MIN, OP_MPI_COMM_WORLD);
        }
#else
        writer = (nLocal > 0) ? rank : std::numeric_limits<int>::max();
#endif
        if (writer == std::numeric_limits<int>::max()) continue; // should not happen if nGlobal>0

        // --- Reference projected measure (frozen per idx & eD)
        bool needRef = (AprojRef.find(idx) == AprojRef.end());
        if (!needRef)
        {
            const dVector3 old = eDRef[idx];
            if ((old * eD) < 0.999999) needRef = true; // direction changed
        }
        if (needRef)
        {
            // IMPORTANT: store the GLOBAL projected measure
            AprojRef[idx] = sumAproj;
            eDRef[idx]    = eD;
        }
        const double Aproj_ref = AprojRef[idx];

        if (dynPress > 0.0 && Aproj_ref > 0.0)
        {
            const double Cd_global = sumFd / (dynPress * Aproj_ref);

            const std::string PFidx   = std::to_string(idx);
            const std::string Coreidx = std::to_string(rank);

            // Local Cp per rank (keep as before)
            std::ofstream floc(dir + "/PF_" + PFidx + "_R_" + Coreidx + "_Cp_local" + fname);
            for (size_t m = 0; m < bps.size(); ++m)
            {
                floc << bps[m].getX() << " "
                     << bps[m].getY() << " "
                     << bps[m].getZ() << " "
                     << CpLoc[m]      << "\n";
            }

            // Mean Cd: ONE file per solid body idx (no rank in filename)
            // Only the elected writer writes it.
            if (rank == writer)
            {
                std::ofstream fMean(dir + "/PF_" + PFidx + "_MeanCd" + fname,
                                    std::ios_base::app);

                // time, Cd_global, N_global, Aproj_ref, Aproj_current(global)
                fMean << Time << " " << Cd_global << " " << nGlobal
                      << " " << Aproj_ref << " " << sumAproj
                      << std::endl;
            }
        }
    }
}

void SolidBody::CalculateLocalNusseltNumber(PhaseField& Phase,
                                            Energy& EN,
                                            MixtureFlow& MF,
                                            FlowSolverLBM& FL,
                                            double Length,
                                            double Tin,
                                            double dn,
                                            double Time,
                                            std::string dir,
                                            std::string fname,
                                            BoundaryConditions& BC)
{
    const double R  = 8.314462618; // J/(mol*K)
    const double ep = 1e-15;
    const double dx = Phase.Grid.dx;

    // sampling distance in lattice units
    const double D  = 1.5;

    int rank = 0;
#ifdef MPI_PARALLEL
    rank = MPI_RANK;
#endif

    // Active directions (your convention)
    const int ax = (Phase.Grid.dNx == 1) ? 1 : 0;
    const int ay = (Phase.Grid.dNy == 1) ? 1 : 0;
    const int az = (Phase.Grid.dNz == 1) ? 1 : 0;
    const int nActive = ax + ay + az; // typically 2 or 3

    // Geometric element represented by one accepted sample:
    // 2D: dℓ  ~ dx    / max(|n| on active axes)
    // 3D: dA  ~ dx^2  / max(|n| on active axes)
    auto dS_element = [&](const dVector3& nUnit) -> double
    {
        const double eps = 1e-12;

        const double nx = ax ? std::abs(nUnit.getX()) : 0.0;
        const double ny = ay ? std::abs(nUnit.getY()) : 0.0;
        const double nz = az ? std::abs(nUnit.getZ()) : 0.0;

        double denom = 0.0;
        if (nActive == 2)
        {
            if (ax && ay) denom = std::max(nx, ny);
            else if (ax && az) denom = std::max(nx, nz);
            else               denom = std::max(ny, nz);
        }
        else
        {
            denom = std::max(nx, std::max(ny, nz));
        }

        const double base = (nActive == 2) ? dx : (dx * dx);
        return base / (denom + eps);
    };

    for (size_t idx = 0; idx < Phase.FieldsProperties.size(); ++idx)
    if (Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        std::vector<dVector3> bps;
        std::vector<double>   locNusselt;

        // For correct mean in 2D/3D: area/length-weighted average
        double sumNu_dS_local = 0.0;  // Σ (Nu_loc * dS)
        double sumdS_local    = 0.0;  // Σ dS
        int    nLocal         = 0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EN.T, 0, )
        {
            const double X = DistanceField(i, j, k);
            if (X < dn && X > dn / 4.0)
            {
                for (auto it = Phase.Fields(i, j, k).cbegin();
                          it != Phase.Fields(i, j, k).cend(); ++it)
                {
                    if (Phase.FieldsProperties[it->index].State != AggregateStates::Solid) continue;

                    const size_t locSolididx = it->index;
                    if (locSolididx != idx) continue;

                    for (size_t is = 0; is < EN.ThermalSurfaceCondition.size(); ++is)
                    {
                        if (EN.ThermalSurfaceCondition.at(is).solidIdx != locSolididx) continue;

                        const double Tb = EN.ThermalSurfaceCondition.at(is).value;

                        const double T0   = 0.5 * (Tb + Tin);
                        const double Rhob = EN.CalculateIdealGasDensity(FL.Pth0, R / EN.Mw, Tb);
                        const double nub  = EN.CalculateSutherlandViscosity(Tb) / Rhob;
                        const double Kb   = Rhob * nub * EN.Cp / EN.Pr;

                        const double Rho0 = EN.CalculateIdealGasDensity(FL.Pth0, R / EN.Mw, T0);
                        const double nu0  = EN.CalculateSutherlandViscosity(T0) / Rho0;
                        const double K0   = Rho0 * nu0 * EN.Cp / EN.Pr;

                        auto keepFluidCorner = [&](int ii, int jj, int kk) -> bool { return !FL.Obstacle(ii, jj, kk); };

                        const dVector3 Xi{double(i), double(j), double(k)};

                        // Normal (solid->fluid). Renormalize for safety.
                        dVector3 n0 = NormField(i, j, k);
                        {
                            const double n2 = n0 * n0;
                            if (n2 > 0.0) n0 = n0 * (1.0 / std::sqrt(n2));
                            else          n0 = {1.0, 0.0, 0.0};
                        }

                        // Interface point and fluid-side points
                        const dVector3 Xb  = Xi - n0 * X;
                        const dVector3 Xf1 = Xb + n0 * D;
                        const dVector3 Xf2 = Xb + n0 * (2.0 * D); // <-- FIX: was D in your code

                        // Fallback temperature for masked interpolation
                        const dVector3 Xfb = Xb + n0 * nFallBack;
                        const double   Tfb = EN.T(int(Xfb.getX()), int(Xfb.getY()), int(Xfb.getZ()));

                        const double s = D * dx;

                        const double Tf1 = MaskedTrilinearScalar(Xf1.getX(), Xf1.getY(), Xf1.getZ(),
                            [&](int ii, int jj, int kk) { return EN.T(ii, jj, kk); },
                            keepFluidCorner, Tfb);

                        const double Tf2 = MaskedTrilinearScalar(Xf2.getX(), Xf2.getY(), Xf2.getZ(),
                            [&](int ii, int jj, int kk) { return EN.T(ii, jj, kk); },
                            keepFluidCorner, Tfb);

                        // Your original one-sided stencil (kept)
                        const double dTdn  = std::abs(4.0 * Tf1 - 3.0 * Tb - Tf2) / (2.0 * s);

                        // Local Nusselt (kept)
                        const double locNu = (Kb / K0) * dTdn * Length / std::abs(Tin - Tb + ep);

                        // NEW: geometric element for correct averaging in 2D/3D
                        const double dS = dS_element(n0);

                        locNusselt.push_back(locNu);
                        bps.push_back(Xb);

                        sumNu_dS_local += locNu * dS;
                        sumdS_local    += dS;
                        ++nLocal;

                        break; // found thermal condition for this solid idx; stop scanning conditions
                    }

                    break; // matched idx for this cell; stop scanning Phase.Fields(i,j,k)
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // ---- Global reductions
        double sumNu_dS = sumNu_dS_local;
        double sumdS    = sumdS_local;
        int    nGlobal  = nLocal;

#ifdef MPI_PARALLEL
        {
            double gNuS = 0.0, gS = 0.0;
            int    gN   = 0;

            OP_MPI_Allreduce(&sumNu_dS, &gNuS, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(&sumdS,    &gS,   1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(&nGlobal,  &gN,   1, OP_MPI_INT,    OP_MPI_SUM, OP_MPI_COMM_WORLD);

            sumNu_dS = gNuS;
            sumdS    = gS;
            nGlobal  = gN;
        }
#endif

        // If nobody saw this solid body, skip
        if (nGlobal <= 0) continue;

        // ---- Elect one writer rank for this idx: smallest rank with nLocal>0
        int writer = std::numeric_limits<int>::max();
#ifdef MPI_PARALLEL
        {
            const int candidate = (nLocal > 0) ? rank : std::numeric_limits<int>::max();
            OP_MPI_Allreduce(&candidate, &writer, 1, OP_MPI_INT, OP_MPI_MIN, OP_MPI_COMM_WORLD);
        }
#else
        writer = (nLocal > 0) ? rank : std::numeric_limits<int>::max();
#endif
        if (writer == std::numeric_limits<int>::max()) continue;

        // ---- Write local (per-rank) file if this rank has samples (same post-processing strategy as drag)
        if (!bps.empty())
        {
            const std::string PFidx   = std::to_string(idx);
            const std::string Coreidx = std::to_string(rank);

            std::ofstream flocout(dir + "/PF_" + PFidx + "_R_" + Coreidx + "_loc" + fname);
            for (size_t m = 0; m < bps.size(); ++m)
            {
                flocout << bps[m].getX() << " "
                        << bps[m].getY() << " "
                        << bps[m].getZ() << " "
                        << locNusselt[m] << "\n";
            }
        }

        // ---- Write mean (ONE file per solid body idx), by elected writer rank
        // MeanNu = (Σ Nu*dS) / (Σ dS)  (valid in 2D and 3D; reduces to simple average if dS constant)
        if (rank == writer && sumdS > 0.0)
        {
            const double MeanNu = sumNu_dS / sumdS;

            const std::string PFidx = std::to_string(idx);
            std::ofstream fMeanout(dir + "/PF_" + PFidx + "_Mean" + fname, std::ios_base::app);

            // time, MeanNu, Nglobal, sumdS (optional diagnostic)
            fMeanout << Time << " " << MeanNu << " " << nGlobal << " " << sumdS << std::endl;
        }
    }
}

bool SolidBody::isOverlapping(const Body& a, const Body& b, int nDim, int DistPlaneAxis) const
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;

    double dist2 = 0.0;

    if (nDim == 3)
    {
        dist2 = dx*dx + dy*dy + dz*dz;
    }
    else // nDim == 2
    {
        // plane normal = DistPlaneAxis, so ignore that component
        switch (DistPlaneAxis)
        {
            case 0: dist2 = dy*dy + dz*dz; break; // yz plane
            case 1: dist2 = dx*dx + dz*dz; break; // xz plane
            default: dist2 = dx*dx + dy*dy; break; // xy plane
        }
    }

    const double rs = a.r + b.r + Clearance;
    return dist2 < rs*rs;
}


std::vector<Body> SolidBody::generateNonOverlappingBodies(
    int nDim,
    double Lx, double Ly, double Lz,      // physical box lengths
    int N,
    double minRadius, double maxRadius,
    double dx,
    int DistPlaneAxis,                   // 0: yz, 1: xz, 2: xy (only if nDim==2)
    double DistPlaneW0)                  // fixed normal coordinate (only if nDim==2)
{
    std::vector<Body> bodies;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> rDist(minRadius/dx, maxRadius/dx);

    const int maxAttempts = 100000;
    int attempts = 0;

    // Box in lattice units
    const double Nx = Lx;
    const double Ny = Ly;
    const double Nz = Lz;

    while ((int)bodies.size() < N && attempts < maxAttempts)
    {
        const double r = rDist(gen);

        Body cand{};
        cand.r = r;

        if (nDim == 3)
        {
            if (Nx <= 2*r || Ny <= 2*r || Nz <= 2*r) break;

            std::uniform_real_distribution<> xDist(r, Nx - r);
            std::uniform_real_distribution<> yDist(r, Ny - r);
            std::uniform_real_distribution<> zDist(r, Nz - r);

            cand.x = xDist(gen);
            cand.y = yDist(gen);
            cand.z = zDist(gen);
        }
        else // nDim == 2
        {
            // Sample only in-plane coords; set normal coord fixed.
            // DistPlaneAxis is the normal axis.
            switch (DistPlaneAxis)
            {
                case 0: // yz plane => x fixed
                {
                    if (Ny <= 2*r || Nz <= 2*r) { attempts = maxAttempts; break; }
                    std::uniform_real_distribution<> uDist(r, Ny - r); // y
                    std::uniform_real_distribution<> vDist(r, Nz - r); // z
                    cand.x = DistPlaneW0;
                    cand.y = uDist(gen);
                    cand.z = vDist(gen);
                    break;
                }
                case 1: // xz plane => y fixed
                {
                    if (Nx <= 2*r || Nz <= 2*r) { attempts = maxAttempts; break; }
                    std::uniform_real_distribution<> uDist(r, Nx - r); // x
                    std::uniform_real_distribution<> vDist(r, Nz - r); // z
                    cand.x = uDist(gen);
                    cand.y = DistPlaneW0;
                    cand.z = vDist(gen);
                    break;
                }
                default: // 2: xy plane => z fixed
                {
                    if (Nx <= 2*r || Ny <= 2*r) { attempts = maxAttempts; break; }
                    std::uniform_real_distribution<> uDist(r, Nx - r); // x
                    std::uniform_real_distribution<> vDist(r, Ny - r); // y
                    cand.x = uDist(gen);
                    cand.y = vDist(gen);
                    cand.z = DistPlaneW0;
                    break;
                }
            }
        }

        bool overlap = false;
        for (const auto& b : bodies)
        {
            if (isOverlapping(cand, b, nDim, DistPlaneAxis))
            {
                overlap = true;
                break;
            }
        }

        if (!overlap)
            bodies.push_back(cand);

        ++attempts;
    }
    // Convert to physical units
    for (auto& b : bodies)
    {
        b.x *= dx; b.y *= dx; b.z *= dx; b.r *= dx;
    }

    if ((int)bodies.size() < N)
    {
        std::cerr << "Warning: Only " << bodies.size()
                  << " bodies placed after " << maxAttempts << " attempts.\n";
    }

    return bodies;
}



void SolidBody::saveBodiesToFile(const std::vector<Body>& bodies, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open file for writing: " << filename << "\n";
        return;
    }

    for (const auto& b : bodies)
        out << b.x << " " << b.y << " " << b.z << " " << b.r << "\n";
}


std::vector<Body> SolidBody::loadBodiesFromFile(const std::string& filename)
{
    std::ifstream in(filename);
    std::vector<Body> bodies;

    if (!in) {
        std::cerr << "Warning: Cannot open file for reading: " << filename << "\n";
        return bodies;
    }

    Body b;
    while (in >> b.x >> b.y >> b.z >> b.r)
        bodies.push_back(b);

    return bodies;
}



void SolidBody::DistributeSolidBodies(int DistPlaneAxis, double DistPlaneW0, const std::string& file)
{
    const bool generateHere = (Grid.OffsetX == 0 &&
                               Grid.OffsetY == 0 &&
                               Grid.OffsetZ == 0);

    int nDim = Grid.Active();
    double dx = Grid.dx;

    if (Do_randDist)
    {
        if (generateHere)
        {
            const double Lx = (XNDistZone - X0DistZone);
            const double Ly = (YNDistZone - Y0DistZone);
            const double Lz = (ZNDistZone - Z0DistZone);

            rand_Bodies = generateNonOverlappingBodies(nDim, Lx, Ly, Lz,
                                static_cast<int>(nParticles),
                                MinRadious, MaxRadious, dx, DistPlaneAxis, DistPlaneW0);
            saveBodiesToFile(rand_Bodies, file);
        }
        #ifdef MPI_PARALLEL
        const size_t nSync = std::min(rand_Bodies.size(), nParticles);
        for (size_t i = 0; i < nSync; ++i)
        {
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Bodies[i].x, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Bodies[i].y, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Bodies[i].z, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Bodies[i].r, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        }
        #endif
    }
    else
    {
        rand_Bodies = loadBodiesFromFile(file);
        if (rand_Bodies.size() < nParticles)
            std::cerr << "Warning: file has " << rand_Bodies.size()
                      << " bodies but nParticles=" << nParticles << "\n";
    }
}


void SolidBody::WriteVTKDistanceField(Settings& locSettings, const int tStep, const int precision)
{
    std::vector<VTK::Field_t> ListOfFields;
    ListOfFields.push_back((VTK::Field_t) {"Distance", [this](int i,int j,int k){return DistanceField(i,j,k);}});
    ListOfFields.push_back((VTK::Field_t){"NormalVector",  [this](int i,int j,int k){return NormField(i,j,k);}});
    std::string Filename = FileInterface::MakeFileName(locSettings.VTKDir,"DistanceFields_", tStep, ".vts");
    VTK::Write(Filename, locSettings, ListOfFields, precision);
}

void SolidBody::SetBoundaryConditions(BoundaryConditions& BC)
{
    if (Grid.dNx > 0) BC.SetX(DistanceField);
    if (Grid.dNy > 0) BC.SetY(DistanceField);
    if (Grid.dNz > 0) BC.SetZ(DistanceField);

    if (Grid.dNx > 0) BC.SetXVector(NormField);
    if (Grid.dNy > 0) BC.SetYVector(NormField);
    if (Grid.dNz > 0) BC.SetZVector(NormField);
}

dVector3 SolidBody::DistanceToInterface(const dVector3& X1, const dVector3& X2)
{
    constexpr double kTinyLen     = 1e-30;
    constexpr double kParallelEps = 1e-14;
    constexpr double kTinySum     = 1e-12;

    // If 3D result deviates from 1D by more than this (in lattice units), reject it.
    // For your case (link length ~ 2), 0.5 is a safe, strong filter.
    constexpr double kMaxDevFrom1D = 0.5;

    // --- Distances (grid units) ---
    const double d1 = DistanceField(X1.getX(), X1.getY(), X1.getZ());
    const double d2 = DistanceField(X2.getX(), X2.getY(), X2.getZ());

    // Choose fluid-side endpoint as the one with positive distance (flip if your sign convention is opposite)
    const bool x1_is_fluid = (d1 >= 0.0);

    const dVector3& Xf = x1_is_fluid ? X1 : X2;
    const dVector3& Xs = x1_is_fluid ? X2 : X1;

    const double X  = x1_is_fluid ? d1 : d2;              // signed
    const double Xp = std::fabs(x1_is_fluid ? d2 : d1);   // positive

    // --- Periodic link Xf -> Xs ---
    dVector3 link;
    CalculateDistancePeriodic(Xs, Xf, link, Grid.Nx, Grid.Ny, Grid.Nz);

    const double l = link.length();
    if (l <= kTinyLen) return dVector3{0.0, 0.0, 0.0};

    const dVector3 e = link * (1.0 / l);

    // --- Always compute robust 1D distance (this is your bounce-back logic) ---
    double q1D = 0.0;
    {
        const double sum = std::fabs(X) + Xp;
        q1D = (sum > kTinySum) ? (l * std::fabs(X) / sum) : 0.0;
        if (q1D < 0.0) q1D = 0.0;
        if (q1D > l)   q1D = l;
    }

    // --- Normals ---
    dVector3 n0 = NormField(Xf.getX(), Xf.getY(), Xf.getZ());
    dVector3 n1 = NormField(Xs.getX(), Xs.getY(), Xs.getZ());

    const double n0l = n0.length();
    const double n1l = n1.length();

    // If normals unreliable => use 1D
    if (n0l <= kTinyLen || n1l <= kTinyLen)
    {
        // return displacement from original X1
        if (x1_is_fluid) return e * q1D;
        const dVector3 Xint = X2 + e * q1D; // here e points X2->X1 if Xf==X2
        dVector3 disp;
        CalculateDistancePeriodic(Xint, X1, disp, Grid.Nx, Grid.Ny, Grid.Nz);
        return disp;
    }

    n0 *= (1.0 / n0l);
    n1 *= (1.0 / n1l);

    // Make normals consistent
    if ((n0 * n1) < 0.0) n1 *= -1.0;

    // --- Reconstruct boundary points (same as your working code) ---
    const dVector3 Xb  = Xf - n0 * X;    // uses signed X
    const dVector3 Xbp = Xs + n1 * Xp;   // uses + and abs

    // Boundary segment direction must also be periodic
    dVector3 S;
    CalculateDistancePeriodic(Xbp, Xb, S, Grid.Nx, Grid.Ny, Grid.Nz);

    const dVector3 w = Xf - Xb;

    // --- 3D candidate ---
    double q3D = q1D; // start with safe value
    bool   use3D = true;

    const dVector3 n_vec = e.cross(S);
    const double denom   = n_vec * n_vec;

    if (denom < kParallelEps)
    {
        use3D = false;
    }
    else
    {
        q3D = - (w.cross(S) * n_vec) / denom;

        // clamp
        if (q3D < 0.0) q3D = 0.0;
        if (q3D > l)   q3D = l;

        // HARD sanity gate: if 3D disagrees with 1D too much, reject 3D.
        if (std::fabs(q3D - q1D) > kMaxDevFrom1D)
            use3D = false;
    }

    const double q = use3D ? q3D : q1D;

    // Return displacement from ORIGINAL X1 to interface point
    if (x1_is_fluid)
    {
        // Xf == X1
        return e * q;
    }
    else
    {
        // Xf == X2, e points X2->X1
        const dVector3 Xint = X2 + e * q;
        dVector3 disp;
        CalculateDistancePeriodic(Xint, X1, disp, Grid.Nx, Grid.Ny, Grid.Nz);
        return disp;
    }
}
