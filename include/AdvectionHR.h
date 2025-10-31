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

 *   File created :   2015
 *   Main contributors :   Philipp Engels; Marvin Tegeler; Raphael Schiedung
 *
 */

#ifndef ADVECTIONHR_H
#define ADVECTIONHR_H

#include "BoundaryConditions.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Velocities.h"
#include "Containers/NodeVectorN.h"
#include <set>

namespace openphase
{
class OP_EXPORTS AdvectionHR : OPObject
{
 public:

    AdvectionHR(){};
    AdvectionHR(Settings& locSettings, std::string InputFileName = DefaultInputFileName);

    void Initialize(Settings& Settings, std::string ObjectNameSuffix = "") override;
    void ReadInput(const std::string InputFileName) override;
    void ReadInput(std::stringstream& inp) override;
    void AdvectPhaseField(
            PhaseField& Phase,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dt,
            const int tStep,
            const bool finalize = true);                                        ///< single-time steps advection method for phase-fields that advects solid phases incompressible and fluid phases compressible
    void AdvectPhaseFieldALE(
                PhaseField& Phase,
                const Velocities& Vel,
                const BoundaryConditions& BC,
                const double dt,
                const int tStep,
                const bool finalize = true);                                    ///< single-time steps advection method for phase-fields that advects solid phases compressible

    template<bool compressible = true, class T, size_t Rank>
    void AdvectField(
            Storage3D<T, Rank>& Field,
            Storage3D<T, Rank>& FieldDot,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dt) const;                                             ///< single-time step compressible advection method

   template<bool compressible = true, class T, size_t Rank>
   void CalculateAdvectionIncrements(
            const Storage3D<T, Rank>& Field,
            Storage3D<T, Rank>& FieldDot,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dt) const;                                             ///< Calculates FieldDot without adding it to Field (GrandPotential Method)

    template<bool compressible = true, class T, size_t Rank>
    void AdvectField(
            Storage3D<T, Rank>& Field,
            Storage3D<T, Rank>& FieldDot,
            Storage3D<T, Rank>& FieldBackup,
            const Velocities& Vel,
            const BoundaryConditions& BC,
            const double dt,
            const int tStep) const;                                             ///< multi-time steps advection method (non-standard method!)

 private:

    AdvectionSchemes Scheme;
    double (*Limiter)(const double,const double);

    static double Slope_Limiter_Minmod(const double a, const double b)
    {

        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(std::abs(a),std::abs(b));
    }

    static double Slope_Limiter_MC(const double a, const double b)
    {
        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(0.5*std::abs(a+b),std::min(2.0*std::abs(a),2.0*std::abs(b)));
    }

    static double Slope_Limiter_Superbee(const double a, const double b)
    {
        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::max(std::min(2.0*std::abs(a),std::abs(b)),std::min(std::abs(a),2.0*std::abs(b)));
    }

    static double Slope_Limiter_Upwind(const double a, const double b)
    {
        return 0.0;
    }
    /*inline double Slope_Limiter_LaxWendroff(const double a, const double b)
    {
        return 1.0;
    }*/

    template<typename...>
        struct dependent_false : std::false_type {};                            ///< False if evaluated at compile time

    template<bool compressible>
    static double AdvectionKernel(int& CFL, double v, double vp, double vm,
            double q, double qp, double qm, double qpp, double qmm,
            double dt, double dx,
            double (*Limiter)(const double,const double));                      ///< Fundamental advection algorithm

    template<bool compressible, class T, size_t Rank>
    static void CalculacteLocalAdvectionIncrements(
            const int i, const int j, const int k,
            const Storage3D<T, Rank>& Field,
            Storage3D<T, Rank>& FieldDot,
            const Storage3D<dVector3, 0>& vel,
            int& CFL,
            const size_t direction,
            const double dt,
            const double dx,
            double (*Limiter)(const double, const double));                                    ///<< Calculates local advection in one direction in calling of AdvectionKernel adequate to data type T

    template<bool compressible, class T, size_t Rank>
    static void CalculateAdvection(
            Storage3D<T, Rank>& Field,
            Storage3D<T, Rank>& FieldDot,
            Storage3D<T, Rank>& FieldBackup,

            const BoundaryConditions& BC,
            const Velocities& Vel,
            const int direction,
            double dt,
            double (*Limiter)(const double,const double),
            bool& warning);                                                     ///<< Calculates advection in one direction by repeated calling of CalculacteLocalAdvectionIncrements if needed  (multi time-step with adding increments)

    size_t GetFluidIndex(const PhaseField& Phase, int long tStep);              ///< Gets fuild with which all solids will interact when moving

    void CalculateLocalAdvectionPhaseField(
            const int i, const int j, const int k,
            const int ii, const int jj, const int kk,
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Velocities& Vel,
            int& CFL,
            const int direction,
            const double dt,
            double (*Limiter)(const double,const double),
            size_t FluidPhaseIdx);                                              ///<< Calculates local advection in one direction in calling of AdvedtionKernel adequate to data type T

    void CalculateLocalAdvectionPhaseFieldALE(
            const int i, const int j, const int k,
            const int ii, const int jj, const int kk,
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Velocities& Vel,
            int& CFL,
            const int direction,
            const double dt,
            double (*Limiter)(const double,const double));                      ///<< Calculates local advection in one direction in calling of AdvedtionKernel adequate to data type T

    void CalculateAdvectionPhaseField(
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Velocities& Vel,
            const int direction,
            double dt,
            double (*Limiter)(const double,const double),
            size_t FluidPhaseIdx,
            bool& warning);                                                     ///<< Calculates advection in one direction by repeated calling of CalculacteLocalAdvection if needed

    void CalculateAdvectionPhaseFieldALE(
            PhaseField& Phase,
            const BoundaryConditions& BC,
            const Velocities& Vel,
            const int direction,
            double dt,
            double (*Limiter)(const double,const double),
            bool& warning);                                                     ///<< Calculates advection in one direction by repeated calling of CalculacteLocalAdvection if needed
};

template<bool compressible>
inline double AdvectionHR::AdvectionKernel(int& CFL,
        double v, double vp, double vm,
        double q, double qp, double qm, double qpp, double qmm,
        double dt, double dx, double (*Limiter)(const double,const double))
{
    if constexpr (compressible) // General case
    {
        const auto L0  = 0.5*Limiter(q-qm,qp-q);
        const auto Lp1 = 0.5*Limiter(qp-q,qpp-qp);
        const auto Lm1 = 0.5*Limiter(qm-qmm,q-qm);
        const auto ap  = 0.5*(vp+v);
        const auto am  = 0.5*(vm+v);
        const auto Fp  = std::max(ap,0.0)*(q+L0)   + std::min(ap,0.0)*(qp-Lp1);
        const auto Fm  = std::max(am,0.0)*(qm+Lm1) + std::min(am,0.0)*(q-L0);

        if (std::abs(ap)*dt > 0.49*dx or std::abs(am)*dt > 0.49*dx) CFL = 1;

        return - (Fp-Fm)/dx;
    }
    else // Incompressible - specific case
    {
        const auto L0  = 0.5*Limiter(q-qm,qp-q);
        const auto Lp1 = 0.5*Limiter(qp-q,qpp-qp);
        const auto Lm1 = 0.5*Limiter(qm-qmm,q-qm);
        const auto F   = std::max(v,0.0)*(q+L0-qm-Lm1) + std::min(v,0.0)*(qp-Lp1-q+L0);

        if (std::abs(v)*dt > 0.49*dx) CFL = 1;

        return - F/dx;
    }
}

template<bool compressible, class T, size_t Rank>
void AdvectionHR::CalculacteLocalAdvectionIncrements(
        const int i, const int j, const int k,
        const Storage3D<T, Rank>& Field,
        Storage3D<T, Rank>& FieldDot,
        const Storage3D<dVector3, 0>& vel,
        int& CFL,
        const size_t direction,
        const double dt,
        const double dx,
        double (*Limiter)(const double, const double))
{
    const int ii = (direction == 0) ? 1 : 0;
    const int jj = (direction == 1) ? 1 : 0;
    const int kk = (direction == 2) ? 1 : 0;

    const double v  = vel(i   ,j    ,k  )[direction];
    const double vp = vel(i+ii,j+jj,k+kk)[direction];
    const double vm = vel(i-ii,j-jj,k-kk)[direction];

    const auto& q   = Field(i     ,j     ,k     );
    const auto& qp  = Field(i+  ii,j+  jj,k+  kk);
    const auto& qm  = Field(i-  ii,j-  jj,k-  kk);
    const auto& qpp = Field(i+2*ii,j+2*jj,k+2*kk);
    const auto& qmm = Field(i-2*ii,j-2*jj,k-2*kk);

    if constexpr (std::is_fundamental<T>::value && Rank == 0 && !is_node<T>::value)
    {
        FieldDot(i,j,k) += AdvectionKernel<compressible>(CFL,v,vp,vm,q,qp,qm,qpp,qmm,dt,dx,Limiter);
    }
    else if constexpr ((!std::is_fundamental<T>::value != (Rank != 0)) && !is_node<T>::value)
    {
        // Tensor typed T will be flattened or Storage3<T,Rank>(i,j,k) with
        // Rank > 1 will be flattened if T is a scalar
        for (size_t n = 0; n < Field(i,j,k).size(); ++n)
        {
            FieldDot(i,j,k)[n] += AdvectionKernel<compressible>(CFL,v,vp,vm,q[n],qp[n],qm[n],qpp[n],qmm[n],dt,dx,Limiter);
        }
    }
    else if constexpr (std::is_same<T,NodeVectorN>::value)
    {
        if constexpr (Rank == 0)
        {
            std::set<size_t> GrainIdxs;
            for (auto alpha = q.  cbegin(); alpha != q.  cend(); alpha++) GrainIdxs.insert(alpha->index);
//            for (auto alpha = qp. cbegin(); alpha != qp. cend(); alpha++) GrainIdxs.insert(alpha->index);
//            for (auto alpha = qm. cbegin(); alpha != qm. cend(); alpha++) GrainIdxs.insert(alpha->index);
//            for (auto alpha = qpp.cbegin(); alpha != qpp.cend(); alpha++) GrainIdxs.insert(alpha->index);
//            for (auto alpha = qmm.cbegin(); alpha != qmm.cend(); alpha++) GrainIdxs.insert(alpha->index);

            for (auto idx = GrainIdxs.cbegin(); idx != GrainIdxs.cend(); idx++)
            for (size_t ss = 0; ss < Field(i, j, k).size_of_vector(); ss++)
            {
                double q_value   = q  .get(*idx, ss);
                double qp_value  = qp .exist(*idx) ? qp .get(*idx, ss) : q_value;
                double qpp_value = qpp.exist(*idx) ? qpp.get(*idx, ss) : qp_value;
                double qm_value  = qm .exist(*idx) ? qm .get(*idx, ss) : q_value;
                double qmm_value = qmm.exist(*idx) ? qmm.get(*idx, ss) : qm_value;

                const double dphi_dt   = AdvectionKernel<compressible>(CFL, v, vp, vm, q_value, qp_value, qm_value, qpp_value, qmm_value, dt, dx, Limiter);
                FieldDot(i, j, k).add(*idx, ss, dphi_dt);
            }
        }
        else
        {
            for (size_t n = 0; n < Rank; n++)
            {
                std::set<size_t> GrainIdxs;
                for (auto alpha = q  [n].cbegin(); alpha != q  [n].cend(); alpha++) GrainIdxs.insert(alpha->index);
//                for (auto alpha = qp [n].cbegin(); alpha != qp [n].cend(); alpha++) GrainIdxs.insert(alpha->index);
//                for (auto alpha = qm [n].cbegin(); alpha != qm [n].cend(); alpha++) GrainIdxs.insert(alpha->index);
//                for (auto alpha = qpp[n].cbegin(); alpha != qpp[n].cend(); alpha++) GrainIdxs.insert(alpha->index);
//                for (auto alpha = qmm[n].cbegin(); alpha != qmm[n].cend(); alpha++) GrainIdxs.insert(alpha->index);

                for (auto idx = GrainIdxs.cbegin(); idx != GrainIdxs.cend(); idx++)
                for (size_t ss = 0; ss < Field(i, j, k)[n].size_of_vector(); ss++)
                {
                    double q_value   = q[n].get(*idx, ss);
                    double qp_value  = qp[n].exist(*idx)  ? qp[n].get(*idx, ss)  : q_value;
                    double qpp_value = qpp[n].exist(*idx) ? qpp[n].get(*idx, ss) : qp_value;
                    double qm_value  = qm[n].exist(*idx)  ? qm[n].get(*idx, ss)  : q_value;
                    double qmm_value = qmm[n].exist(*idx) ? qmm[n].get(*idx, ss) : qm_value;

                    const double dphi_dt = AdvectionKernel<compressible>(CFL, v, vp, vm, q_value, qp_value, qm_value, qpp_value, qmm_value, dt, dx, Limiter);
                    FieldDot(i, j, k)[n].add(*idx, ss, dphi_dt);
                }
            }
        }
    }
    else static_assert(dependent_false<T>::value, "Advected can not be deducted for Datatype!");
}

template<bool compressible, class T, size_t Rank>
void AdvectionHR::CalculateAdvectionIncrements(
        const Storage3D<T, Rank>& Field,
        Storage3D<T, Rank>& FieldDot,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt) const
{
    if (Field.Bcells() < 2)
    {
        ConsoleOutput::WriteExit("Number of Bcells needs to be 2 or higher.", "AdvectionHR", "AdvectField");
        OP_Exit(EXIT_FAILURE);
    }

    int CFL = 0;
    if (Vel.Grid.dNx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,reduction(+:CFL))
        {
            CalculacteLocalAdvectionIncrements<compressible>(i,j,k,Field,FieldDot,Vel.Average,CFL,0,dt,Vel.Grid.dx,Limiter);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if (Vel.Grid.dNy)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,reduction(+:CFL))
        {
            CalculacteLocalAdvectionIncrements<compressible>(i,j,k,Field,FieldDot,Vel.Average,CFL,1,dt,Vel.Grid.dx,Limiter);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if (Vel.Grid.dNz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,reduction(+:CFL))
        {
            CalculacteLocalAdvectionIncrements<compressible>(i,j,k,Field,FieldDot,Vel.Average,CFL,2,dt,Vel.Grid.dx,Limiter);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

template<bool compressible, class T, size_t Rank>
void AdvectionHR::AdvectField(
        Storage3D<T, Rank>& Field,
        Storage3D<T, Rank>& FieldDot,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt) const
{
    CalculateAdvectionIncrements<compressible>(Field,FieldDot,Vel,BC,dt);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
    {
        Field(i,j,k) += FieldDot(i,j,k)*dt;
        FieldDot(i,j,k) *= 0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

template<bool compressible, class T, size_t Rank>
void AdvectionHR::CalculateAdvection(
        Storage3D<T, Rank>& Field,
        Storage3D<T, Rank>& FieldDot,
        Storage3D<T, Rank>& FieldBackup,
        const BoundaryConditions& BC,
        const Velocities& Vel,
        const int direction,
        double dt,
        double (*Limiter)(const double,const double),
        bool& warning)
{
    FieldBackup = Field;

    unsigned int iterations = 1;
    int CFL;
    do
    {
        CFL = 0;
        for (unsigned int it = 0; it < iterations; ++it)
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,reduction(+:CFL))
            {
                CalculacteLocalAdvectionIncrements<compressible>(i,j,k,Field,FieldDot,Vel.Average,CFL,direction,dt,Vel.Grid.dx,Limiter);
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
            {
                Field(i,j,k) += FieldDot(i,j,k)*dt;
                FieldDot(i,j,k) *= 0;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            BC.Set(Field,direction);
        }
        if (CFL != 0)
        {
            Field = FieldBackup;

            dt *= 0.5;
            warning = true;
        }
        iterations *= 2;
    }
    while (CFL > 0);
}

template<bool compressible, class T, size_t Rank>
void AdvectionHR::AdvectField(
        Storage3D<T, Rank>& Field,
        Storage3D<T, Rank>& FieldDot,
        Storage3D<T, Rank>& FieldBackup,
        const Velocities& Vel,
        const BoundaryConditions& BC,
        const double dt,
        const int tStep) const
{
    if (Field.Bcells() < 2)
    {
        ConsoleOutput::WriteExit("Number of Bcells needs to be 2 or higher.", "AdvectionHR", "AdvectField");
        OP_Exit(EXIT_FAILURE);
    }

    bool warning = false;
    if (tStep % 2 == 0)
    {
        if (Vel.Grid.dNx) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 0, dt, Limiter, warning);
        if (Vel.Grid.dNy) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 1, dt, Limiter, warning);
        if (Vel.Grid.dNz) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 2, dt, Limiter, warning);
    }
    else
    {
        if (Vel.Grid.dNz) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 2, dt, Limiter, warning);
        if (Vel.Grid.dNy) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 1, dt, Limiter, warning);
        if (Vel.Grid.dNx) CalculateAdvection<compressible>(Field, FieldDot, FieldBackup, BC, Vel, 0, dt, Limiter, warning);
    }

    if (warning)
    {
        ConsoleOutput::WriteWarning("Warning reduce time step, cfl > 0.5!", "AdvectionHR", "AdvectField");
    }
};
}// namespace openphase
#endif
