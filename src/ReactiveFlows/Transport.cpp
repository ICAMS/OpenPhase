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
 
#include "ReactiveFlows/Transport.h"
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

void Transport::Initialize(Settings& locSettings)
{
    Grid = locSettings.Grid;
}

double Transport::CalculatePhiVanLeer(double r)
{
    return (r+fabs(r))/(1.0+fabs(r));
}

/*
void  Transport::CalculateDiffusion(PhaseField& Phase, Energy& EN, FlowSolverLBM& FL, double dt)
{
    std::vector<int> dir(3);
    double dx = Grid.dx;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,SP.MassFrac,0,)
    {
        if(!FL.Obstacle(i,j,k))
        {
            double Tc     = EN.TOld(i,j,k);
            double Rhoc   = FL.DensityWetting(i,j,k,{MixtureComp});
            double cpc    = EN.CpMix(i,j,k);
            double ktc    = EN.KMix(i,j,k);
            for (int direction=0; direction<3; ++direction)
            {
                dir[0] = (direction == 0) ? 1*Grid.dNx : 0;
                dir[1] = (direction == 1) ? 1*Grid.dNy : 0;
                dir[2] = (direction == 2) ? 1*Grid.dNz : 0;
                int ip=i+dir[0];
                int jp=j+dir[1];
                int kp=k+dir[2];
                int im=i-dir[0];
                int jm=j-dir[1];
                int km=k-dir[2];
                double Tp   = EN.TOld(ip,jp,kp);
                double Tm   = EN.TOld(im,jm,km);
                double ktp  = EN.KMix(ip,jp,kp);
                double ktm  = EN.KMix(im,jm,km);

                EN.T(i, j, k) += dt/(Rhoc*cpc) * (ktc * (Tp+Tm-2.0*Tc)/dx/dx + (ktp-ktm)/dx/2.0 * (Tp-Tm)/dx/2.0);
                double Rhop    = FL.DensityWetting(ip,jp,kp,{MixtureComp});
                double Rhom    = FL.DensityWetting(im,jm,km,{MixtureComp});
                int ipp=i+2.0*dir[0];
                int jpp=j+2.0*dir[1];
                int kpp=k+2.0*dir[2];
                int imm=i-2.0*dir[0];
                int jmm=j-2.0*dir[1];
                int kmm=k-2.0*dir[2];
                */
                /*
                double YkVkc1[SP.nSpecies];
                double YkVkp1[SP.nSpecies];
                double YkVkm1[SP.nSpecies];

                double MMWc  = CalculateMeanMolarMass(i, j, k, "old");
                double MMWp  = CalculateMeanMolarMass(ip ,jp ,kp  ,"old");
                double MMWpp = CalculateMeanMolarMass(ipp,jpp,kpp ,"old");
                double MMWm  = CalculateMeanMolarMass(im ,jm ,km  ,"old");
                double MMWmm = CalculateMeanMolarMass(imm,jmm,kmm ,"old");

                double YkVcc=0.0;
                double YkVcp=0.0;
                double YkVcm=0.0;

                for(size_t ic =0; ic < SP.nSpecies; ic++)
                {
                    double Mk    = MolecularWeight({ic});
                    double Xc    = MoleFraction(i,j,k, MMWc ,ic);
                    double Xp    = MoleFraction(ip ,jp ,kp, MMWp,ic);
                    double Xpp   = MoleFraction(ipp,jpp,kpp, MMWpp,ic);
                    double Xm    = MoleFraction(im ,jm ,km, MMWm,ic);
                    double Xmm   = MoleFraction(imm,jmm,kmm, MMWmm,ic);

                    YkVkc1[ic]   = - SP.DSp(i,j,k,{ic})   *Mk/MMWc*(Xp-Xm) /dx/2.0;
                    YkVkp1[ic]   = - SP.DSp(ip,jp,kp,{ic})*Mk/MMWp*(Xpp-Xc)/dx/2.0;
                    YkVkm1[ic]   = - SP.DSp(im,jm,km,{ic})*Mk/MMWm*(Xc-Xmm)/dx/2.0;

                    YkVcc       +=  SP.MassFracOld(i,j,k)   ({ic})*SP.DSp(i,j,k,{ic})   * Mk/MMWc*(Xp-Xm) /dx/2.0;
                    YkVcp       +=  SP.MassFracOld(ip,jp,kp,{ic})*SP.DSp(ip,jp,kp,{ic})* Mk/MMWp*(Xpp-Xc)/dx/2.0;
                    YkVcm       +=  SP.MassFracOld(im,jm,km,{ic})*SP.DSp(im,jm,km,{ic})* Mk/MMWm*(Xc-Xmm)/dx/2.0;
                }

                for(size_t ic =0; ic < SP.nSpecies; ic++)
                {
                    EN.T(i, j, k) -= dt/cpc * (SP.CpSp(i,j,k,{ic}) * (YkVkc1[ic]+ YkVcc)) * (Tp-Tm)/dx/2.0;
                    SP.MassFracOld(i, j, k,{ic}) -= dt/Rhoc * (Rhop * (YkVkp1[ic]+ YkVcp) - Rhom*(YkVkm1[ic]+ YkVcm))/dx/2.0;
                }
                */
               /*
                std::vector<double> YkVkc(SP.nSpecies);
                std::vector<double> YkVkp(SP.nSpecies);
                std::vector<double> YkVkm(SP.nSpecies);

                for(size_t ic = 0; ic < SP.nSpecies; ic++)
                {
                    double MFc  = SP.MassFracOld(i,j,k,{ic});
                    double MFp  = SP.MassFracOld(ip,jp,kp,{ic});
                    double MFpp = SP.MassFracOld(ipp,jpp,kpp,{ic});
                    double MFm  = SP.MassFracOld(im,jm,km,{ic});
                    double MFmm = SP.MassFracOld(imm,jmm,kmm,{ic});

                    YkVkc[ic]   = - SP.DSp(i,j,k,{ic})   *(MFp-MFm)/dx/2.0;
                    YkVkp[ic]   = - SP.DSp(ip,jp,kp,{ic})*(MFpp-MFc)/dx/2.0;
                    YkVkm[ic]   = - SP.DSp(im,jm,km,{ic})*(MFc-MFmm)/dx/2.0;
                }

                for(size_t ic = 0; ic < SP.nSpecies; ic++)
                {
                    EN.T(i, j, k) -= dt/cpc * (SP.CpSp(i,j,k,{ic}) * (YkVkc[ic])) * (Tp-Tm)/dx/2.0;
                    SP.MassFrac(i, j, k,{ic}) -= dt/Rhoc * (Rhop * (YkVkp[ic]) - Rhom*(YkVkm[ic]))/dx/2.0;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
*/

void Transport::CalculateAdvectionDiffusion(PhaseField& Phase, Species &SP, Energy& EN, FlowSolverLBM& FL, MixtureFlow &MF,
                                                   SolidBody& SB, double dt)
{
    constexpr double eps = 1e-10;
    constexpr double alpha_min = 0.35;

    const double dx     = Grid.dx;
    const double invDx  = 1.0 / dx;
    const double inv2Dx = 0.5 * invDx;
    const double invDx2 = invDx * invDx;

    // Active directions only
    int activeDirs[3]; int nActive=0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto unit_or_fallback = [](dVector3 n, const dVector3& fallback) -> dVector3 {
        const double nl = n.length();
        if (nl > 1e-30) n *= (1.0/nl);
        else n = fallback;
        return n;
    };

    // Neighbor sampling (OLD) with obstacle reconstruction
    auto SampleNeighborOld = [&](int i,int j,int k, int in,int jn,int kn,
                                 const double Tc,const Tensor<double,1>& Yc,
                                 double& Tn,Tensor<double,1>& Yn,double& kn_eff)
    {
        if (!FL.Obstacle(in,jn,kn))
        {
            Tn = EN.TOld(in,jn,kn);
            Yn = SP.MassFracOld(in,jn,kn);
            kn_eff = EN.KMix(in,jn,kn);
            return;
        }

        double r = std::max({abs(in-i), abs(jn-j), abs(kn-k) });   

        dVector3 Xf{double(i),double(j),double(k)};
        dVector3 Xs{double(in),double(jn),double(kn)};
        dVector3 dist = SB.DistanceToInterface(Xf, Xs);

        const double alpha = std::clamp(dist.length(), alpha_min, 1.0);

        dVector3 nf = SB.NormField(i,j,k);
        dVector3 ns = SB.NormField(in,jn,kn);
        dVector3 n  = unit_or_fallback(nf + ns, nf);
        
        const double Tw = EN.SurfaceTemperature(Phase, FL, SB, Xf + dist, n, "old");
        Tn = r * (Tw - Tc) / alpha + Tc;

        // Diffusion: use center k as you already do
        kn_eff = EN.KMix(i,j,k);

        Tensor<double,1> Yw = SP.SurfaceSpecies(FL, SB, Xf + dist, n, "old");
        for (size_t s=0; s<SP.nSpecies; ++s)
            Yn[s] = Yc[s] + r * (Yw[s] - Yc[s]) / alpha;
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,SP.MassFrac,0,)
    {
        if (FL.Obstacle(i,j,k)) continue;

        const double Tc   = EN.TOld(i,j,k);
        const double Rhoc = MF.DensityMix(i,j,k);
        const double cpc  = EN.CpMix(i,j,k);
        const double ktc  = EN.KMix(i,j,k);

        Tensor<double,1> Yc = SP.MassFracOld(i,j,k);

        // reusable buffers
        Tensor<double,1> Yp({SP.nSpecies}), Ym({SP.nSpecies});
        Tensor<double,1> Ypp({SP.nSpecies}), Ymm({SP.nSpecies});
        double Tp=Tc, Tm=Tc, Tpp=Tc, Tmm=Tc;
        double ktp=ktc, ktm=ktc, dummyK=ktc;

        std::vector<double> YkVkc(SP.nSpecies), YkVkp(SP.nSpecies), YkVkm(SP.nSpecies);

        for (int a=0; a<nActive; ++a)
        {
            const int dir = activeDirs[a];
            int di=0,dj=0,dk=0;
            if (dir==0) di=Grid.dNx;
            else if (dir==1) dj=Grid.dNy;
            else dk=Grid.dNz;

            const int ip=i+di, jp=j+dj, kp=k+dk;
            const int im=i-di, jm=j-dj, km=k-dk;
            const int ipp=i+2*di, jpp=j+2*dj, kpp=k+2*dk;
            const int imm=i-2*di, jmm=j-2*dj, kmm=k-2*dk;

            // Sample ONCE (old)
            SampleNeighborOld(i,j,k, ip,jp,kp, Tc, Yc, Tp,  Yp,  ktp);
            SampleNeighborOld(i,j,k, im,jm,km, Tc, Yc, Tm,  Ym,  ktm);
            SampleNeighborOld(i,j,k, ipp,jpp,kpp, Tc, Yc, Tpp, Ypp, dummyK);
            SampleNeighborOld(i,j,k, imm,jmm,kmm, Tc, Yc, Tmm, Ymm, dummyK);

            const double Vc = MF.VelocityMix(i,j,k)[dir];

            // ---------------- Advection ----------------

            if (EN.AdvScheme==AdvectionSchemes::Upwind)
            {
                if (Vc >= 0.0)
                {
                    EN.T(i,j,k) += -dt*invDx * Vc * (Tc - Tm);
                    for (size_t s=0; s<SP.nSpecies; ++s)
                        SP.MassFrac(i,j,k,{s}) += -dt*invDx * Vc * (Yc[s] - Ym[s]);
                }
                else
                {
                    EN.T(i,j,k) += -dt*invDx * Vc * (Tp - Tc);
                    for (size_t s=0; s<SP.nSpecies; ++s)
                        SP.MassFrac(i,j,k,{s}) += -dt*invDx * Vc * (Yp[s] - Yc[s]);
                }
            }
            else if (EN.AdvScheme==AdvectionSchemes::Central)
            {
                EN.T(i,j,k) += -0.5 * dt*invDx * Vc * (Tp - Tm);
                for (size_t s=0; s<SP.nSpecies; ++s)
                    SP.MassFrac(i,j,k,{s}) += -0.5 * dt*invDx * Vc * (Yp[s] - Ym[s]);
            }
            else if(EN.AdvScheme==AdvectionSchemes::VanLeer)
			{
                double Thp = (Tc+Tp)/2.0;
			    double Thm = (Tc+Tm)/2.0;
			    if(Vc>=0)
                {
                    double r    = (Tc-Tm)/(Tp-Tc+eps);
			        double phi  = (r+fabs(r))/(1.0+fabs(r));
			        Thp         = Tc+0.5*phi*(Tp-Tc);
			        double rm  = (Tm-Tmm)/(Tc-Tm+eps);
			        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			        Thm=Tm+0.5*phim*(Tc-Tm);
			        EN.T(i,j,k) -=  dt*invDx * Vc * (Thp - Thm);
                    for(size_t s =0; s < SP.nSpecies; s++)
                    {
                        double r    = (Yc[s]-Ym[s])/(Yp[s]-Yc[s]+eps);
			            double phi  = (r+fabs(r))/(1.0+fabs(r));
			            double Yhp  = Yc[s]+0.5*phi*(Yp[s]-Yc[s]);
			            double rm   = (Ym[s]-Ymm[s])/(Yc[s]-Ym[s]+eps);
			            double phim = (rm+fabs(rm))/(1.0+fabs(rm));
			            double Yhm  = Ym[s]+0.5*phim*(Yc[s]-Ym[s]);
                        SP.MassFrac(i,j,k,{s}) -= dt*invDx * Vc * (Yhp - Yhm);
                    }
    	        }
    	        else
                {
			        double r  = (Tpp-Tp)/(Tp-Tc+eps);
			        double phi=(r+fabs(r))/(1.0+fabs(r));
			        Thp=Tp-0.5*phi*(Tp-Tc);
			        double rm  = (Tp-Tc)/(Tc-Tm+eps);
			        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			        Thm=Tc-0.5*phim*(Tc-Tm);
			        EN.T(i,j,k) -= dt*Vc*(Thp-Thm)/dx;

                    for(size_t s =0; s < SP.nSpecies; s++)
                    {
                        double r    = (Ypp[s]-Yp[s])/(Yp[s]-Yc[s]+eps);
			            double phi  = (r+fabs(r))/(1.0+fabs(r));
			            double Yhp  = Yp[s]-0.5*phi*(Yp[s]-Yc[s]);
			            double rm   = (Yp[s]-Yc[s])/(Yc[s]-Ym[s]+eps);
			            double phim = (rm+fabs(rm))/(1.0+fabs(rm));
			            double Yhm  = Yc[s]-0.5*phim*(Yc[s]-Ym[s]);
			            SP.MassFrac(i, j, k,{s}) -=  dt*invDx  * Vc * (Yhp-Yhm);
                    }
                }
			}

            // ---------------- Diffusion ----------------
            EN.T(i,j,k) += dt/(Rhoc*cpc) * ( ktc*(Tp + Tm - 2.0*Tc)*invDx2
                               + (ktp - ktm)*inv2Dx * (Tp - Tm)*inv2Dx );

            const double Rhop = MF.DensityMix(ip,jp,kp);
            const double Rhom = MF.DensityMix(im,jm,km);

            for (size_t s=0; s<SP.nSpecies; ++s)
            {
                const double MFc = Yc[s];
                YkVkc[s] = -SP.DSp(i,j,k,{s})    * (Yp[s]  - Ym[s])  * inv2Dx;
                YkVkp[s] = -SP.DSp(ip,jp,kp,{s}) * (Ypp[s] - MFc)    * inv2Dx;
                YkVkm[s] = -SP.DSp(im,jm,km,{s}) * (MFc   - Ymm[s])  * inv2Dx;
            }

            const double dTdx = (Tp - Tm) * inv2Dx;

            for (size_t s=0; s<SP.nSpecies; ++s)
            {
                EN.T(i,j,k) -= dt/cpc * (SP.CpSp(i,j,k,{s}) * YkVkc[s]) * dTdx;
                SP.MassFrac(i,j,k,{s}) -= dt/Rhoc * (Rhop*YkVkp[s] - Rhom*YkVkm[s]) * inv2Dx;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Transport::CalculateAdvectionDiffusion(PhaseField& Phase, Energy &EN, FlowSolverLBM& FL, MixtureFlow &MF,
                                                   SolidBody& SB, double dt)
{
    constexpr double eps = 1e-10;
    constexpr double alpha_min = 0.35;

    const double dx     = Grid.dx;
    const double invDx  = 1.0 / dx;
    const double inv2Dx = 0.5 * invDx;
    const double invDx2 = invDx * invDx;

    // Active directions only
    int activeDirs[3]; int nActive=0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto unit_or_fallback = [](dVector3 n, const dVector3& fallback) -> dVector3 {
        const double nl = n.length();
        if (nl > 1e-30) n *= (1.0/nl);
        else n = fallback;
        return n;
    };

    // Neighbor sampling (OLD) with obstacle reconstruction
    auto NeighborTemp = [&](int i,int j,int k, int in,int jn,int kn,
                                 const double Tc,
                                 double& Tn,
                                 double& kn_eff)
    {
        if (!FL.Obstacle(in,jn,kn))
        {
            Tn = EN.TOld(in,jn,kn);
            kn_eff = EN.KMix(in,jn,kn);
            return;
        }

        double r = std::max({abs(in-i), abs(jn-j), abs(kn-k) });   

        dVector3 Xf{double(i),double(j),double(k)};
        dVector3 Xs{double(in),double(jn),double(kn)};
        dVector3 dist = SB.DistanceToInterface(Xf, Xs);

        const double alpha = std::clamp(dist.length(), alpha_min, 1.0);

        dVector3 nf = SB.NormField(i,j,k);
        dVector3 ns = SB.NormField(in,jn,kn);
        dVector3 n  = unit_or_fallback(nf + ns, nf);

        const double Tw = EN.SurfaceTemperature(Phase, FL, SB, Xf + dist, n, "old");
        Tn = r * (Tw - Tc) / alpha + Tc;

        // Diffusion: use center k as you already do
        kn_eff = EN.KMix(i,j,k);
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.T,0,)
    {
        if (FL.Obstacle(i,j,k)) continue;

        const double Tc   = EN.TOld(i,j,k);
        const double Rhoc = MF.DensityMix(i,j,k);
        const double cpc  = EN.CpMix(i,j,k);
        const double ktc  = EN.KMix(i,j,k);

        double Tp=Tc, Tm=Tc;
        double ktp=ktc, ktm=ktc, dummyK=ktc;

        for (int a=0; a<nActive; ++a)
        {
            const int dir = activeDirs[a];
            int di=0,dj=0,dk=0;
            if (dir==0) di=Grid.dNx;
            else if (dir==1) dj=Grid.dNy;
            else dk=Grid.dNz;

            const int ip=i+di, jp=j+dj, kp=k+dk;
            const int im=i-di, jm=j-dj, km=k-dk;

            // Sample ONCE (old)
            NeighborTemp(i,j,k, ip,jp,kp, Tc, Tp, ktp);
            NeighborTemp(i,j,k, im,jm,km, Tc, Tm, ktm);

            const double Vc = MF.VelocityMix(i,j,k)[dir];

            // ---------------- Advection ----------------
            if (EN.AdvScheme==AdvectionSchemes::Upwind)
            {
                if (Vc >= 0.0)
                {
                    EN.T(i,j,k) += -dt*invDx * Vc * (Tc - Tm);
                }
                else
                {
                    EN.T(i,j,k) += -dt*invDx * Vc * (Tp - Tc);
                }
            }
            else if (EN.AdvScheme==AdvectionSchemes::Central)
            {
                EN.T(i,j,k) += -0.5 * dt*invDx * Vc * (Tp - Tm);
            }
            else if(EN.AdvScheme==AdvectionSchemes::VanLeer)
			{
                double Tpp=Tc, Tmm=Tc;
                const int ipp=i+2*di, jpp=j+2*dj, kpp=k+2*dk;
                const int imm=i-2*di, jmm=j-2*dj, kmm=k-2*dk;

                NeighborTemp(i,j,k, ipp,jpp,kpp, Tc, Tpp, dummyK);
                NeighborTemp(i,j,k, imm,jmm,kmm, Tc, Tmm, dummyK);

                double Thp = (Tc+Tp)/2.0;
			    double Thm = (Tc+Tm)/2.0;
			    if(Vc>=0)
                {
                    double r    = (Tc-Tm)/(Tp-Tc+eps);
			        double phi  = (r+fabs(r))/(1.0+fabs(r));
			        Thp         = Tc+0.5*phi*(Tp-Tc);
			        double rm  = (Tm-Tmm)/(Tc-Tm+eps);
			        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			        Thm=Tm+0.5*phim*(Tc-Tm);
			        EN.T(i,j,k) -=  dt*invDx * Vc * (Thp - Thm);
    	        }
    	        else
                {
			        double r  = (Tpp-Tp)/(Tp-Tc+eps);
			        double phi=(r+fabs(r))/(1.0+fabs(r));
			        Thp=Tp-0.5*phi*(Tp-Tc);
			        double rm  = (Tp-Tc)/(Tc-Tm+eps);
			        double phim=(rm+fabs(rm))/(1.0+fabs(rm));
			        Thm=Tc-0.5*phim*(Tc-Tm);
			        EN.T(i,j,k) -= dt*Vc*(Thp-Thm)/dx;
                }
			}

            // ---------------- Diffusion ----------------
            EN.T(i,j,k) += dt/(Rhoc*cpc) * ( ktc*(Tp + Tm - 2.0*Tc)*invDx2
                               + (ktp - ktm)*inv2Dx * (Tp - Tm)*inv2Dx );
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Transport::CalculateSolidDiffusion(PhaseField& Phase, Energy &EN, FlowSolverLBM& FL, 
                                            SolidBody& SB, double dt)
{
    constexpr double alpha_min = 0.35;
    const double dx    = Grid.dx;
    const double invDx2 = 1.0 / (dx * dx);

    // Active directions only
    int activeDirs[3]; int nActive = 0;
    if (Grid.dNx) activeDirs[nActive++] = 0;
    if (Grid.dNy) activeDirs[nActive++] = 1;
    if (Grid.dNz) activeDirs[nActive++] = 2;

    auto unit_or_fallback = [](dVector3 n, const dVector3& fallback) -> dVector3 {
        const double nl = n.length();
        if (nl > 1e-30) n *= (1.0 / nl);
        else n = fallback;
        return n;
    };

    // Sample neighbor solid temperature for diffusion stencil.
    // If neighbor is solid -> use EN.TsOld(neighbor).
    // If neighbor is fluid -> reconstruct ghost from interface Tw with alpha (curved interface).
    auto NeighborTemp = [&](int i,int j,int k, int in,int jn,int kn,
                                      const double Tc,
                                      const dVector3& Xs0,
                                      const dVector3& ns_center) -> double
    {
        if (FL.Obstacle(in,jn,kn))
        {
            return EN.TsOld(in,jn,kn);
        }

        // neighbor is fluid: reconstruct using interface temperature at the wall point
        const dVector3 Xf0{double(in), double(jn), double(kn)};

        // dist is vector from solid node toward interface along segment (as returned by your SB)
        const dVector3 dist = SB.DistanceToInterface(Xs0, Xf0);
        const double alpha  = std::clamp(dist.length(), alpha_min, 1.0);

        const dVector3 nf = SB.NormField(in,jn,kn);
        const dVector3 n  = unit_or_fallback(nf + ns_center, ns_center);

        const double Tw = EN.SurfaceTemperature(Phase, FL, SB, (Xs0 + dist), n, "old");

        // ghost solid-side neighbor temperature
        return Tc + (Tw - Tc) / alpha;
    };

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.Ts,0,)
    {
        if (!FL.Obstacle(i,j,k)) continue;

        const double Tc   = EN.TsOld(i,j,k);
        const double Rhoc = EN.DensitySol(i,j,k);
        const double cpc  = EN.CpSol(i,j,k);
        const double ktc  = EN.KSol(i,j,k);

        const dVector3 Xs0{double(i), double(j), double(k)};
        const dVector3 ns_center = SB.NormField(i,j,k);

        // accumulate Laplacian contributions
        double lapT = 0.0;

        for (int a = 0; a < nActive; ++a)
        {
            const int direction = activeDirs[a];

            int di=0, dj=0, dk=0;
            if (direction == 0) di = Grid.dNx;
            else if (direction == 1) dj = Grid.dNy;
            else dk = Grid.dNz;

            const int ip = i + di, jp = j + dj, kp = k + dk;
            const int im = i - di, jm = j - dj, km = k - dk;

            const double Tp = NeighborTemp(i,j,k, ip,jp,kp, Tc, Xs0, ns_center);
            const double Tm = NeighborTemp(i,j,k, im,jm,km, Tc, Xs0, ns_center);

            lapT += (Tp + Tm - 2.0 * Tc);
        }

        // Write new Ts (safe even if Ts wasn't pre-initialized)
        EN.Ts(i,j,k) = Tc + dt * (ktc / (Rhoc * cpc)) * lapT * invDx2;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Transport::CalculatePhaseFieldCoupling(PhaseField& Phase, Energy &EN, MixtureFlow& MF, double hT, double dt)
{
    double Tsolid=EN.TempSolid;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.T,0,)
    {
        if (Phase.Fields(i,j,k).interface())
        {
            double SF=0;
            for (auto it = Phase.Fields(i,j,k).cbegin(); it != Phase.Fields(i,j,k).cend(); ++it)
            {
                size_t PhaseIdx =  Phase.FieldsProperties[it->index].Phase;
                if(Phase.FieldsProperties[PhaseIdx].State == AggregateStates::Solid)
                {
                    for (size_t is = 0; is < EN.ThermalSurfaceCondition.size(); is++)
                    {
                        if(EN.ThermalSurfaceCondition.at(is).solidIdx==PhaseIdx)
                        {
                            if(EN.ThermalSurfaceCondition.at(is).type==ThermalCondition::Type::ConstantTemp)
                            {
                                Tsolid = EN.ThermalSurfaceCondition.at(is).value;
                                SF += it->value;
                            }
                            else if(EN.ThermalSurfaceCondition.at(is).type==ThermalCondition::Type::ConstantFlux)
                            {
                                Tsolid = EN.TOld(i,j,k);
                                SF += it->value;
                            }
                        }
                    }
                }
            }
            double ST  = hT * EN.KMix(i,j,k)* pow(SF,3)*(1.0-SF) / pow(Grid.iWidth*Grid.dx,2.0) * (EN.TOld(i,j,k)-Tsolid);
            EN.T(i,j,k) -=   dt / MF.DensityMix(i,j,k) / EN.CpMix(i,j,k) * ST; 
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Transport::CalculateReaction(Species &SP, Energy& EN, MixtureFlow &MF, FlowSolverLBM& FL, 
                                    int nDim, double &FuelConsumptionRate, size_t GasFuelIndex, double dt)
{
	FuelConsumptionRate=0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,SP.MassFrac,0,)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			EN.T(i,j,k) += dt/(MF.DensityMix(i, j, k) * EN.CpMix(i,j,k)) * EN.HRR(i,j,k);
			for (size_t ic=0; ic<SP.nSpecies; ic++)
    		{
        		SP.MassFrac(i, j, k,{ic}) += dt / MF.DensityMix(i, j, k) * (SP.WSp(i, j, k,{ic}));
    		}
    		// Compute Fuel Consumption only if species is fuel
    		FuelConsumptionRate += abs(SP.WSp(i, j, k,{GasFuelIndex}) * pow(MF.Grid.dx, nDim));
		}
	}
    OMP_PARALLEL_STORAGE_LOOP_END
	#ifdef MPI_PARALLEL
    	OP_MPI_Allreduce(OP_MPI_IN_PLACE, &FuelConsumptionRate, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#endif
}
