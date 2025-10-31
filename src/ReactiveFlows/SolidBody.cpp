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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
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

    Do_randDist    	 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("Random_Dist"), false, false);
    X0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("X0DistZone"), false, -10);
    XNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("XNDistZone"), false, -10);
    Y0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Y0DistZone"), false, -10);
    YNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("YNDistZone"), false, -10);
    Z0DistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("Z0DistZone"), false, -10);
    ZNDistZone       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("ZNDistZone"), false, -10);

    MinRadious       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MinRadious"), false,  10);
    MaxRadious       = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("MaxRadious"), false,  20);

    for (size_t i = 0; i < nParticles; ++i)
    {
        rand_Circles.push_back({0.0, 0.0, 0.0});
        rand_Spheres.push_back({0.0, 0.0, 0.0, 0.0});
    }
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

void SolidBody::CalculateLocalNusseltNumber(PhaseField& Phase, EnergyTransport& ET, FlowMixture& FM, FlowSolverLBM& FL, double Length,
                                                    double Tin, double dn, double Time, string dir, string fname, BoundaryConditions& BC)
{
    size_t Gasidx   = 0;
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Gas)
    {
        Gasidx = idx;
    }

    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        std::vector<dVector3> fps;
        std::vector<dVector3> bps;
        std::vector<double> locNusselt;
        size_t Solididx = idx;
        double ep = 1e-15;
        double Tb = ET.TempSolid;
        double T0 = (Tb+Tin)/2.0;

        double Kb = 1.0;
        double K0 = 1.0;

        double dx = Grid.dx;
        double R  = 8.314462618; // J/(mol*K)

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ET.Tx,0,)
        {
            //if(Phase.Interfaces(i,j,k)>1)
            if(!FL.Obstacle(i,j,k))
            {
                //if(FL.Obstacle(i,j,k))
                dVector3 Norm{0.0,0.0,0.0};
                for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                {
                    if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
                    {
                        size_t locSolididx=it->index;
                        if(locSolididx==Solididx)
                        {
                            NodeAB<dVector3,dVector3> Normals=Phase.Normals(i,j,k);
                            Norm=Normals.get_asym1(Solididx, Gasidx);
                            //Norm=Normals.get_sym(Solididx, Gasidx);
                            Tb = ET.SurfaceTemp[Solididx];
                            T0 = (Tb+Tin)/2.0;
                            double Rhob = FM.CalculateIdealGasDensity(FL.Pth0, R/ET.Mw, Tb);
                            double nub  = FM.CalculateSutherlandViscosity(Tb)/Rhob;
                            Kb =Rhob*nub*ET.Cp/ET.Pr;

                            double Rho0 = FM.CalculateIdealGasDensity(FL.Pth0, R/ET.Mw, T0);
                            double nu0  = FM.CalculateSutherlandViscosity(T0)/Rho0;
                            K0 =Rho0*nu0*ET.Cp/ET.Pr;

                            double Phif = Phase.Fractions(i,j,k,{Gasidx});
                            double DX   = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
                            //double DX   = fabs(Phase.iWidth/Pi * acos(2.0*Phif-1));
                            dVector3 X  = Norm*DX;
                            if(FL.Obstacle(i,j,k)) X=X*-1;
                            dVector3 Xp = {double(i+Grid.OffsetX), double(j+Grid.OffsetY), double(k+Grid.OffsetZ)};
                            dVector3 Xb = Xp - X;
                            dVector3 Xc = Xb + Norm*dn;
                            fps.push_back(Xc);
                            bps.push_back(Xb);
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        std::vector<dVector3> nfps;
        std::vector<dVector3> nbps;

        for (size_t i = 0; i < bps.size(); ++i)
        {
            const auto& p = bps[i];
            bool found = false;
            for (const auto& u : nbps)
            {
                if (is_close(p, u))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                nbps.push_back(p);
                nfps.push_back(fps[i]);
            }
        }

        double MeanNu=0.0;
        double nlocNu=0.0;
        for (size_t i = 0; i < nfps.size(); ++i)
        {
            double ic    = nfps[i].getX()-Grid.OffsetX;
            double jc    = nfps[i].getY()-Grid.OffsetY;
            double kc    = nfps[i].getZ()-Grid.OffsetZ;
        
            double Tf    = CalculateInterpolatedTemperature(ET,ic,jc,kc);
            double dTdn  = abs(Tf - Tb)/(dn*dx); 
            double locNu = Kb/K0 * dTdn * Length / abs(Tin-Tb+ep);
            //double locNu = dTdn * Length / abs(Tin-Tb);
            locNusselt.push_back(locNu);
            MeanNu += locNu;
            nlocNu++;
        }

        if(MeanNu!=0.0)
        {
            int rank=0;
            #ifdef MPI_PARALLEL
                rank=MPI_RANK;
            #endif
            string PFidx   = std::to_string(Solididx);
            string Coreidx = std::to_string(rank);

            // ---- Write centers to file ----

            std::ofstream flocout(dir+"/PF_"+PFidx+"_R_"+Coreidx+"_loc"+fname);
            if (!flocout) 
            {
                std::cerr << "Error: Could not open output file.\n";
            }
            for (size_t i = 0; i < nbps.size(); ++i)
            {
                flocout << nbps[i].getX() << " " << nbps[i].getY() << " " << nbps[i].getZ() << " " << locNusselt[i] << "\n";
            }
            flocout.close();

            std::ofstream fMeanout(dir+"/PF_"+PFidx+"_R_"+Coreidx+"_Mean"+fname, std::ios_base::app);
            if (!fMeanout) 
            {
                std::cerr << "Error: Could not open output file.\n";
            }
            fMeanout << Time << " " << MeanNu/nlocNu << " "<< nlocNu<< endl;
            fMeanout.close();
            //std::cout << "rank "<<rank<<", nbps.size(): "<<nbps.size()<<endl;
            std::cout << "boundary points and nusselt values written to "<<dir<<"PF_"<<PFidx<<"_R_"<<Coreidx<<"_"<<fname<<endl;
        }
    }
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

void SolidBody::CalculateDragCoeff(PhaseField& Phase, FlowSolverLBM& FL, FlowMixture& FM, double A_proj, double Rho0,
                                                 double U0, double dn, double Time, size_t axis,
                                                 string dir, string fname, BoundaryConditions& BC)
{
    size_t Gasidx   = 0;
    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Gas)
    {
        Gasidx = idx;
    }
    dVector3 nf  ({0,0,0});
    
    for(size_t id=0;id<3;id++)
    {
        if(id==axis)
        {
            nf[id]=1;
        }
    }
    double cte = 0.5*Rho0*U0*U0*A_proj;
    //double ep  = 1e-8;
    double dx  = Grid.dx;

    for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
    if(Phase.FieldsProperties[idx].State == AggregateStates::Solid)
    {
        std::vector<dVector3> fps;
        std::vector<dVector3> bps;
        //std::vector<double> locCp;
        size_t Solididx = idx;
        
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FL.DensityWetting,0,)
        {
            if(Phase.Interfaces(i,j,k)>1)
            {
                //if(!FL.Obstacle(i,j,k))
                for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
                {
                    if(Phase.FieldsProperties[it->index].State == AggregateStates::Solid)
                    {
                        size_t locSolididx=it->index;
                        dVector3 Norm{0.0,0.0,0.0};
                        if(locSolididx==Solididx)
                        {
                            NodeAB<dVector3,dVector3> Normals=Phase.Normals(i,j,k);
                            Norm=Normals.get_asym1(Solididx, Gasidx);
                            //Norm=Normals.get_sym(Solididx, Gasidx);
                            double Phif = Phase.Fractions(i,j,k,{Gasidx});
                            double DX   = fabs(Phase.Grid.iWidth/Pi * asin(1.0-2.0*Phif));
                            //double DX   = fabs(Phase.iWidth/Pi * acos(2.0*Phif-1));
                            dVector3 X  = Norm*DX;
                            if(FL.Obstacle(i,j,k)) X=X*-1;
                            dVector3 Xp = {double(i+Grid.OffsetX), double(j+Grid.OffsetY), double(k+Grid.OffsetZ)};
                            dVector3 Xb = Xp - X;
                            dVector3 Xc = Xb + Norm*dn;
                            fps.push_back(Xc);
                            bps.push_back(Xb);
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        std::vector<dVector3> nfps;
        std::vector<dVector3> nbps;
        std::vector<double> ds;

        for (size_t i = 0; i < bps.size(); ++i)
        {
            const auto& p = bps[i];
            bool found = false;
            for (const auto& u : nbps)
            {
                if (is_close(p, u))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                nbps.push_back(p);
                nfps.push_back(fps[i]);
            }
        }

        double Cd=0.0;
        double Stot=0.0;
        //double nds=0.0;
        
        //auto surf_neighbors = surface_elements_2d(nbps);
        //std::vector<double> ds(nbps.size());
        //for (int i = 0; i < nbps.size(); ++i) 
        //{
        //    int n1 = surf_neighbors[i].first;
        //    int n2 = surf_neighbors[i].second;
        //    ds[i] = 0.5 * (dist3(nbps[i], nbps[n1]) + dist3(nbps[i], nbps[n2]));
        //}
        //std::vector<double> surf_elem = surface_elements_3d(nbps);
        
        for (size_t i = 0; i < nfps.size(); ++i)
        {
            double locds  = curve_element_2d(nbps, i,BC);
            ds.push_back(locds);
            
            double ic     = nfps[i].getX()-Grid.OffsetX;
            double jc     = nfps[i].getY()-Grid.OffsetY;
            double kc     = nfps[i].getZ()-Grid.OffsetZ;
            dVector3 locn = (nfps[i] - nbps[i]).normalize();
            double Pf     = CalculateInterpolatedPressure(FL,ic,jc,kc);
            double muf    = CalculateInterpolatedViscosity(FL,ic,jc,kc);

            dVector3 Vf   = CalculateInterpolatedVelocity(FM,ic,jc,kc);
            dVector3 Vt   = Vf  - locn * (Vf  * locn); 
            dVector3 loct = Vt/Vt.abs();
            double tau_w  = muf * Vt.abs()/ (dn*dx);
            double Ft     = nf * loct * (tau_w * ds[i]*dx);
            double Fp     = nf * locn * (-1*Pf * ds[i]*dx);
            Cd += (Ft + Fp)/cte;
            //locCp.push_back((Ft + Fp)/cte);
            Stot += ds[i]*dx; 
            //nds++; 
        }

        if(Cd!=0.0)
        {
            int rank=0;
            #ifdef MPI_PARALLEL
                rank=MPI_RANK;
            #endif
            string PFidx   = std::to_string(Solididx);
            string Coreidx = std::to_string(rank);

            // ---- Write centers to file ----
            /*
            std::ofstream flocout(dir+"/PF_"+PFidx+"_R_"+Coreidx+"_loc"+fname);
            if (!flocout) 
            {
                std::cerr << "Error: Could not open output file.\n";
            }
            for (int i = 0; i < nbps.size(); ++i) 
            {
                flocout << nbps[i].getX() << " " << nbps[i].getY() << " " << nbps[i].getZ() << " " << locCp[i] << " " << ds[i] << "\n";
            }
            flocout.close();
            */

            std::ofstream fMeanout(dir+"/PF_"+PFidx+"_R_"+Coreidx+"_"+fname, std::ios_base::app);
            if (!fMeanout) 
            {
                std::cerr << "Error: Could not open output file.\n";
            }
            fMeanout << Time << " " << Cd << endl;
            fMeanout.close();
            std::cout << "rank "<<rank<<", Numerical Total Area: "<<Stot<<endl;
            std::cout << "boundary points and cd values written to "<<dir<<"PF_"<<PFidx<<"_R_"<<Coreidx<<"_"<<fname<<endl;
        }
    }
}

double SolidBody::CalculateInterpolatedTemperature(EnergyTransport& ET, double ic, double jc, double kc)
{
    double i0   = floor(ic);
    double i1   = i0 + 1.0;
    double j0   = floor(jc);
    double j1   = j0 + 1.0;
    double k0   = floor(kc);
    double k1   = k0 + 1.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    double Tf   =  TrilinearInterPolation(ic,jc,kc,
                                          ET.Tx(i0*idi,j0*jdi,k0*kdi),ET.Tx(i1*idi,j0*jdi,k0*kdi),
                                          ET.Tx(i0*idi,j1*jdi,k0*kdi),ET.Tx(i1*idi,j1*jdi,k0*kdi),
                                          ET.Tx(i0*idi,j0*jdi,k1*kdi),ET.Tx(i1*idi,j0*jdi,k1*kdi),
                                          ET.Tx(i0*idi,j1*jdi,k1*kdi),ET.Tx(i1*idi,j1*jdi,k1*kdi));
    return Tf;
}

dVector3 SolidBody::CalculateInterpolatedVelocity(FlowMixture& FM,double ic, double jc, double kc)
{
    dVector3 Vf({0.0,0.0,0.0});
    double i0   = floor(ic);
    double i1   = i0 + 1.0;
    double j0   = floor(jc);
    double j1   = j0 + 1.0;
    double k0   = floor(kc);
    double k1   = k0 + 1.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    for(size_t iv=0;iv<3;iv++)
    {
        Vf[iv] = TrilinearInterPolation(ic,jc,kc,
                                        FM.V_Mixture(i0*idi,j0*jdi,k0*kdi)[iv],FM.V_Mixture(i1*idi,j0*jdi,k0*kdi)[iv],
                                        FM.V_Mixture(i0*idi,j1*jdi,k0*kdi)[iv],FM.V_Mixture(i1*idi,j1*jdi,k0*kdi)[iv],
                                        FM.V_Mixture(i0*idi,j0*jdi,k1*kdi)[iv],FM.V_Mixture(i1*idi,j0*jdi,k1*kdi)[iv],
                                        FM.V_Mixture(i0*idi,j1*jdi,k1*kdi)[iv],FM.V_Mixture(i1*idi,j1*jdi,k1*kdi)[iv]);
    }
    
    return Vf;
}

double SolidBody::CalculateInterpolatedViscosity(FlowSolverLBM& FL, double ic, double jc, double kc)
{
    size_t n=0;
    double i0   = floor(ic);
    double i1   = i0 + 1.0;
    double j0   = floor(jc);
    double j1   = j0 + 1.0;
    double k0   = floor(kc);
    double k1   = k0 + 1.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    double muf   =  TrilinearInterPolation(ic,jc,kc,
                                          FL.DensityWetting(i0*idi,j0*jdi,k0*kdi,{n})*FL.nut(i0*idi,j0*jdi,k0*kdi,{n}),FL.DensityWetting(i1*idi,j0*jdi,k0*kdi,{n})*FL.nut(i1*idi,j0*jdi,k0*kdi,{n}),
                                          FL.DensityWetting(i0*idi,j1*jdi,k0*kdi,{n})*FL.nut(i0*idi,j1*jdi,k0*kdi,{n}),FL.DensityWetting(i1*idi,j1*jdi,k0*kdi,{n})*FL.nut(i1*idi,j1*jdi,k0*kdi,{n}),
                                          FL.DensityWetting(i0*idi,j0*jdi,k1*kdi,{n})*FL.nut(i0*idi,j0*jdi,k1*kdi,{n}),FL.DensityWetting(i1*idi,j0*jdi,k1*kdi,{n})*FL.nut(i1*idi,j0*jdi,k1*kdi,{n}),
                                          FL.DensityWetting(i0*idi,j1*jdi,k1*kdi,{n})*FL.nut(i0*idi,j1*jdi,k1*kdi,{n}),FL.DensityWetting(i1*idi,j1*jdi,k1*kdi,{n})*FL.nut(i1*idi,j1*jdi,k1*kdi,{n}));
    return muf;
}

double SolidBody::CalculateInterpolatedPressure(FlowSolverLBM& FL, double ic, double jc, double kc)
{
    size_t n=0;
    double i0   = floor(ic);
    double i1   = i0 + 1.0;
    double j0   = floor(jc);
    double j1   = j0 + 1.0;
    double k0   = floor(kc);
    double k1   = k0 + 1.0;

    int idi=Grid.dNx;
    int jdi=Grid.dNy;
    int kdi=Grid.dNz;

    double Pf   =  TrilinearInterPolation(ic,jc,kc,
                                          FL.HydroPressure(i0*idi,j0*jdi,k0*kdi,{n}),FL.HydroPressure(i1*idi,j0*jdi,k0*kdi,{n}),
                                          FL.HydroPressure(i0*idi,j1*jdi,k0*kdi,{n}),FL.HydroPressure(i1*idi,j1*jdi,k0*kdi,{n}),
                                          FL.HydroPressure(i0*idi,j0*jdi,k1*kdi,{n}),FL.HydroPressure(i1*idi,j0*jdi,k1*kdi,{n}),
                                          FL.HydroPressure(i0*idi,j1*jdi,k1*kdi,{n}),FL.HydroPressure(i1*idi,j1*jdi,k1*kdi,{n}));
    return Pf;
}

bool SolidBody::isOverlapping(const Circle& c1, const Circle& c2) 
{
    double dx = c1.x - c2.x;
    double dz = c1.z - c2.z;
    double distanceSquared = dx * dx + dz * dz;
    double radiusSum = c1.r + c2.r + Clearance;
    return distanceSquared < radiusSum * radiusSum;
}

vector<Circle> SolidBody::generateNonOverlappingCircles(double width, double height, int N, double minRadius, double maxRadius, double dx)
{
    std::vector<Circle> circles;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> radiusDist(minRadius/dx, maxRadius/dx);
    std::uniform_real_distribution<> xDist, yDist;

    int maxAttempts = 100000;
    int attempts = 0;

    while (circles.size() < static_cast<size_t>(N) && attempts < maxAttempts) 
    {
        double r = radiusDist(gen);
        xDist = std::uniform_real_distribution<>(r, width - r);
        yDist = std::uniform_real_distribution<>(r, height - r);

        Circle newCircle{xDist(gen), yDist(gen), r};

        bool overlap = false;
        for (const auto& c : circles) 
        {
            if (isOverlapping(newCircle, c)) 
            {
                overlap = true;
                break;
            }
        }

        if (!overlap) 
        {
            Circle phynewCircle {newCircle.x, newCircle.z, newCircle.r};
            circles.push_back(phynewCircle);
        }
        attempts++;
    }

    for (size_t ic = 0; ic < circles.size() ; ic++)
    {
        circles[ic].x *= dx; 
        circles[ic].z *= dx; 
        circles[ic].r *= dx; 
    }

    if (circles.size() < static_cast<size_t>(N)) 
    {
        std::cerr << "Warning: Only " << circles.size() << " circles placed after "
                  << maxAttempts << " attempts.\n";
    }
    return circles;
}

// Save circles to a file
void SolidBody::saveCirclesToFile(const std::vector<Circle>& circles, const std::string& filename, double dx) 
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open file for writing: " << filename << "\n";
        return;
    }

    for (const auto& c : circles) {
        out << c.x << " " << c.z << " " << c.r << "\n";
    }
    out.close();
}

// Load circles from a file
std::vector<Circle> SolidBody::loadCirclesFromFile(const std::string& filename) 
{
    std::ifstream in(filename);
    std::vector<Circle> circles;

    if (!in) {
        std::cerr << "Warning: Cannot open file for reading: " << filename << "\n";
        return circles;
    }

    Circle c;
    while (in >> c.x >> c.z >> c.r) 
    {
        circles.push_back(c);
    }
    in.close();
    return circles;
}

bool SolidBody::isOverlapping(const Sphere& s1, const Sphere& s2) 
{
    double dx = s1.x - s2.x;
    double dy = s1.y - s2.y;
    double dz = s1.z - s2.z;
    double distanceSquared = dx * dx + dy * dy + dz * dz;
    double radiusSum = s1.r + s2.r + Clearance;
    return distanceSquared < radiusSum * radiusSum;
}

std::vector<Sphere> SolidBody::generateNonOverlappingSpheres(double Nx, double Ny, double Nz, int N,double minRadius, double maxRadius, double dx) 
{
    std::vector<Sphere> spheres;
    std::random_device rd; 
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> radiusDist(minRadius/dx, maxRadius/dx);
    std::uniform_real_distribution<> xDist, yDist, zDist;

    int maxAttempts = 10000;
    int attempts = 0;

    while (spheres.size() < static_cast<size_t>(N) && attempts < maxAttempts) 
    {
        double r = radiusDist(gen);

        // Update bounds to keep the sphere entirely inside the box
        xDist = std::uniform_real_distribution<>(r, Nx - r);
        yDist = std::uniform_real_distribution<>(r, Ny - r);
        zDist = std::uniform_real_distribution<>(r, Nz - r);

        Sphere newSphere{ xDist(gen), yDist(gen), zDist(gen), r };

        bool overlap = false;
        for (const auto& s : spheres) 
        {
            if (isOverlapping(newSphere, s)) 
            {
                overlap = true;
                break;
            }
        }

        if (!overlap) 
        {
            Sphere phynewSphere{newSphere.x, newSphere.y, newSphere.z, newSphere.r};
            spheres.push_back(phynewSphere);
        }
        attempts++;
    }

    for (size_t is = 0; is < spheres.size() ; is++)
    {
        spheres[is].x *= dx; 
        spheres[is].y *= dx; 
        spheres[is].z *= dx; 
        spheres[is].r *= dx; 
    }

    if (spheres.size() < static_cast<size_t>(N)) 
    {
        std::cerr << "Warning: Only " << spheres.size() << " spheres were placed after "
                  << maxAttempts << " attempts.\n";
    }

    return spheres;
}

// Save spheres to a file
void SolidBody::saveSpheresToFile(const std::vector<Sphere>& spheres, const std::string& filename, double dx) 
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open file for writing: " << filename << "\n";
        return;
    }

    for (const auto& s : spheres) {
        out << s.x << " " << s.y << " " << s.z << " " << s.r << "\n";
    }
    out.close();
}

// Load spheres from a file
std::vector<Sphere> SolidBody::loadSpheresFromFile(const std::string& filename) 
{
    std::ifstream in(filename);
    std::vector<Sphere> spheres;

    if (!in) 
    {
        std::cerr << "Warning: Cannot open file for reading: " << filename << "\n";
        return spheres;
    }

    Sphere s;
    while (in >> s.x >> s.y >> s.z >> s.r) {
        spheres.push_back(s);
    }
    in.close();
    return spheres;
}

void SolidBody::DistributeRandomSolidBodies(PhaseField& Phase, Settings& locSettings, string dir)
{
    int nDim=locSettings.Grid.Active();
    double dx =locSettings.Grid.dx;
    if(Do_randDist)
    {
        if(Phase.Grid.OffsetX==0 and Phase.Grid.OffsetZ==0)
        {
            if(nDim==2)
            {
                rand_Circles = generateNonOverlappingCircles(XNDistZone-X0DistZone, ZNDistZone-Z0DistZone, 
                                                            nParticles, MinRadious, MaxRadious, dx);
                saveCirclesToFile(rand_Circles,dir+"cylindersData.txt",dx);
            }
            else if(nDim==3)
            {
                rand_Spheres = generateNonOverlappingSpheres(XNDistZone-X0DistZone, YNDistZone-Y0DistZone,  ZNDistZone-Z0DistZone,
                                                             nParticles, MinRadious, MaxRadious, dx);
                saveSpheresToFile(rand_Spheres,dir+"spheresData.txt",dx);
            }
        }
    
        if(nDim==2)
        {
            for (size_t i = 0; i < nParticles; ++i) 
            {
                #ifdef MPI_PARALLEL
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Circles[i].r, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Circles[i].x, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Circles[i].z, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                #endif
            }
        }
        else if(nDim==3)
        {
            for (size_t i = 0; i < nParticles; ++i) 
            {
                #ifdef MPI_PARALLEL
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Spheres[i].r, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Spheres[i].x, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Spheres[i].y, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &rand_Spheres[i].z, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
                #endif
            }
        }
    }

    if(!Do_randDist)
    {
        if(nDim==2)
        {
            rand_Circles = loadCirclesFromFile(dir+"cylindersData.txt");
        }
        else if(nDim==3)
        {
            rand_Spheres = loadSpheresFromFile(dir+"spheresData.txt");
        }
    }
}
