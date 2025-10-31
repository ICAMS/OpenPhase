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

 *   File created :   2021
 *   Main contributors :   Raphael Schiedung
 *
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifndef ROOTFINDINGALGORITHMS_H
#define ROOTFINDINGALGORITHMS_H

namespace openphase
{
struct OP_EXPORTS RootFindingAlgorithms ///< Numerically find the solution x of f(x) = 0
{
    //const static inline std::runtime_error ZeroGradient  = std::runtime_error("ZeroGradient ");
    //const static inline std::runtime_error NotConverging = std::runtime_error("NotConverging");
    //const static inline std::runtime_error std::runtime_error(MaxIterations) = std::runtime_error("std::runtime_error(MaxIterations)");

    static constexpr auto ZeroGradient  = "ZeroGradient";
    static constexpr auto NotConverging = "NotConverging";
    static constexpr auto MaxIterations = "MaxIterations";

    template<typename UNARY, typename T>
    static void Bisection(UNARY func, T& x0, T xa, T xb, T accuracy, unsigned long iterations)
    {
        T fb = func(xb);
        T fa = func(xa);
        if ((xa > xb)  or (fa > 0 and fb > 0) or (fa < 0 and fb < 0))
        {
            std::stringstream message;
            message << "Unsuitable initial conditions!\n"
                    << "xa = " << xa << "\n"
                    << "xb = " << xb << "\n"
                    << "fa = " << fa << "\n"
                    << "fb = " << fb;
            throw std::runtime_error(message.str());
        }

        T f0 = func(x0);
        unsigned long itr = 0;
        while (std::abs(f0) > accuracy && std::abs((xb - xa)/xa) > std::numeric_limits<double>::epsilon())
        {
            if ((fa > 0 and f0 > 0) or (fa < 0 and f0 < 0))
            {
                xa = x0;
                fa = f0;
            }
            else
            {
                xb = x0;
                fb = f0;
            }

            itr++; if (itr >= iterations) throw std::runtime_error(MaxIterations);
            x0 = (xb + xa)/2.0;
            f0 = func(x0);
        }
    }

    template<typename UNARY1, typename UNARY2, typename T>
    static void Newton(UNARY1 func, UNARY2 dfunc, T& x0, T accuracy, unsigned long iterations)
    {
        T f0 = func(x0);
        unsigned long itr = 0;
        while (std::abs(f0) > accuracy)
        {
            T df0 = dfunc(x0); if (df0 == 0) throw std::runtime_error(ZeroGradient);
            //std::cout << itr << " " << x0 << " " << f0 <<  " " << df0 << "\n";
            T p0  = - f0/df0;
            x0 += p0;
            itr++; if (itr >= iterations) throw std::runtime_error(MaxIterations);
            f0 = func(x0);
        }
    }

    template<typename UNARY1, typename UNARY2, typename T>
    static void PositivNewton(UNARY1 func, UNARY2 dfunc, T& x0, T accuracy, unsigned long iterations)
    {
        T f0 = func(x0);
        unsigned long itr = 0;
        while (std::abs(f0) > accuracy)
        {
            T df0 = dfunc(x0); if (df0 == 0.0) throw std::runtime_error(ZeroGradient);
            T p0  = - f0/df0;
            while (x0 + p0 <= 0.0 ) p0/=2.0;
            x0 += p0;
            itr++; if (itr >= iterations and x0 > 0.0) throw std::runtime_error(MaxIterations);
            f0 = func(x0);
        }
    }

    template<typename UNARY, typename T>
    static void Secant(UNARY func, T& x0, T& x1, T accuracy, unsigned long iterations)
    {
        T f0 = func(x0);
        T f1 = func(x1);
        if (std::abs(f0) > accuracy)
        {
            unsigned long itr = 0;
            do
            {
                T dx = (x0 - x1);
                T df = (f0 - f1); if (df == 0) throw std::runtime_error(ZeroGradient);
                T p0 = - f0*dx/df;
                x1  = x0;
                x0 += p0;
                f1  = f0;
                f0  = func(x0);
                itr++; if (itr >= iterations) throw std::runtime_error(MaxIterations);
            }
            while(std::abs(f0) > accuracy);
        }
    }

    template<typename UNARY, typename T>
    static void Broyden(UNARY func, T& x0, T accuracy, unsigned long iterations)
    {
        auto gradient = [func] (T x0, T dx)
        {
            T x_p1 = x0 + dx;
            T x_m1 = x0 - dx;
            T f_m1 = func(x_m1);
            T f_p1 = func(x_p1);
            return (f_p1 - f_m1)/dx/2;
        };

        T x1 = x0;
        T f0 = func(x0);
        T f1 = f0;

        if (std::abs(f0) > accuracy)
        {
            unsigned long itr = 1;
            T eps = 1.0;
            T del = 0.5;
            T dx  = f0/1000;
            T J0  = gradient(x0,dx);
            if (J0 == 0) throw std::runtime_error(ZeroGradient);
            T Jinv0 = 1.0/J0;
            T p0 = - Jinv0*f0;
            x0  += eps*p0;
            f0   = func(x0);
            itr++;

            do
            {
                dx  = del*(x0 - x1);
                J0  = gradient(x0,dx);
                if (J0 == 0) throw std::runtime_error(ZeroGradient);
                Jinv0 = 1.0/J0;
                eps = 1.0;
                del = 0.5;
                x1  = x0;
                f1  = f0;
                while(true)
                {
                    T p0 = - Jinv0*f0;
                    x0 += eps*p0;
                    f0  = func(x0);
                    if (f0*f0 < f1*f1 or itr < 3) break;
                    T eps_new = 0.75*eps;
                    del = eps_new-eps;
                    eps = eps_new;
                    x0 = x1;
                    if (eps < accuracy)  throw std::runtime_error(NotConverging);
                }
                itr++; if (itr>=iterations)  throw std::runtime_error(MaxIterations);
            }
            while (std::abs(f0) > accuracy);
        }
    }

    template<typename function, typename vector, typename T>
    static void Broyden(function func, vector& x0, T accuracy, unsigned long iterations)
    {
        std::size_t size = x0.size();
        typedef std::vector<std::vector<T>> matrix;
        vector f_p1;
        auto calculate_jacobian = [&func,&size,&f_p1] (matrix& J, const vector& x0, const vector& f0, const vector& dx)
        {
            f_p1 = f0;
            for (std::size_t m = 0; m < size; m++)
            {
                vector x_p1 = x0;
                x_p1[m] += dx[m];
                vector f_p1 = func(x_p1);
                for (std::size_t n = 0; n < size; n++)
                {
                    J[n][m] = (f_p1[n] - f0[n])/dx[m];
                }
            }
        };

        auto calculate_inverted = [size](matrix& Jinv, const matrix& J)
        {
            if (size == 1)
            {
                if (J[0][0] == 0) throw std::runtime_error(ZeroGradient);
                Jinv[0][0] = 1.0/J[0][0];
            }
            else if (size == 2)
            {
                T det = J[0][0]*J[1][1]-J[1][0]*J[1][0];
                if (det == 0) throw std::runtime_error(ZeroGradient);
                Jinv[0][0] =  J[1][1]/det;
                Jinv[0][1] = -J[0][1]/det;
                Jinv[1][0] = -J[1][0]/det;
                Jinv[1][1] =  J[0][0]/det;
            }
            else if (size == 3)
            {
                const T& a = J[0][0];
                const T& b = J[0][1];
                const T& c = J[0][2];
                const T& d = J[1][0];
                const T& e = J[1][1];
                const T& f = J[1][2];
                const T& g = J[2][0];
                const T& h = J[2][1];
                const T& i = J[2][2];
                T A =  (e*i-f*h);
                T B = -(d*i-f*g);
                T C =  (d*h-e*g);
                T D = -(b*i-c*h);
                T E =  (a*i-c*g);
                T F = -(a*h-b*g);
                T G =  (b*f-c*e);
                T H = -(a*f-c*d);
                T I =  (a*e-b*d);
                T det = a*A+b*B+c*C;
                if (det == 0) throw std::runtime_error(ZeroGradient);
                Jinv[0][0] = A/det;
                Jinv[0][1] = B/det;
                Jinv[0][2] = C/det;
                Jinv[1][0] = D/det;
                Jinv[1][1] = E/det;
                Jinv[1][2] = F/det;
                Jinv[2][0] = G/det;
                Jinv[2][1] = H/det;
                Jinv[2][2] = I/det;
            }
            else throw std::invalid_argument("Matrix inversion for more that 3 components is not yet implemented");
        };

        matrix J0    (size,vector(size,0.0));
        matrix J1    (size,vector(size,0.0));
        matrix Jinv0 (size,vector(size,0.0));
        matrix Jinv1 (size,vector(size,0.0));

        vector x1 = x0;
        vector f0 = func(x0);
        vector f1 = f0;
        vector dx  (size,0.0);
        vector p0  (size,0.0);
        vector y0  (size,0.0);

        auto L2Norm = [size](vector& x)
        {
            T norm = 0.0;
            for (std::size_t n = 0; n < size; n++) norm += x[n]*x[n];
            return std::sqrt(norm);
        };
        T f_norm0 = L2Norm(f0);
        T f_norm1 = f_norm0;

        auto iterate = [size,accuracy] (const vector& f)
        {
            for (std::size_t n = 0; n < size; n++)
            {
                if (std::abs(f[n]) > accuracy) return true;
            }
            return false;
        };
        if (iterate(f0))
        {
            //bool converging = true;
            T eps = 1.0;
            T del = 1.0;
            std::size_t itr = 0;

            // calculate inverted Jacobian
            auto calculate_step_size = [size] (vector& dx, const vector& f0)
            {
                T dx_init = std::abs(f0[0]);
                for (std::size_t n = 1; n < size; n++)
                {
                    if (std::abs(f0[n]) > dx_init) dx_init = std::abs(f0[n]);
                }
                dx_init /= 10;
                for (std::size_t n = 0; n < size; n++)
                {
                    dx[n] = dx_init;
                }
            };
            calculate_step_size(dx,f0);
            calculate_jacobian(J0,x0,f0,dx);
            calculate_inverted(Jinv0,J0);

            auto update_mu = [&]()
            {
                for (std::size_t n = 0; n < size; n++) f1[n] = f0[n];
                f_norm1 = f_norm0;
                for (std::size_t n = 0; n < size; n++)
                {
                    x1[n] = x0[n];
                }
                del = 1.0;
                eps = 1.0;
                while(true)
                {
                    for (std::size_t n = 0; n < size; n++)
                    {
                        p0[n] = 0.0;
                    }
                    for (std::size_t n = 0; n < size; n++)
                    for (std::size_t m = 0; m < size; m++)
                    {
                        p0[n] -= Jinv0[n][m]*f0[m];
                    }
                    for (std::size_t n = 0; n < size; n++)
                    {
                        x0[n] += eps*p0[n];
                    }
                    f0      = func(x0);
                    f_norm0 = L2Norm(f0);

                    if (f_norm0 < f_norm1 or itr < 3) break;
                    T eps_new = 0.75*eps;
                    del = eps_new-eps;
                    eps = eps_new;
                    for (std::size_t n = 0; n < size; n++)
                    {
                        x0[n] = x1[n];
                    }
                    if (eps < accuracy)
                    {
                        //converging = false;
                        break;
                    }
                }
                itr++;
            };
            update_mu();

            while (iterate(f0))
            {
                // Recalculate Jacobian
                {
                    calculate_step_size(dx,f0);
                    calculate_jacobian(J0,x0,f0,dx);
                    calculate_inverted(Jinv0,J0);
                }

                // Update Jacobian
                /*{
                    if(del == 1)
                    {
                        for (std::size_t n = 0; n < size; n++)
                        {
                            y0[n] = (f0[n]-f1[n]);
                        }
                    }
                    else
                    {
                        for (std::size_t n = 0; n < size; n++)
                        {
                            mus[n] = mu0[n] - del*p0[n];
                        }
                        calculate_resiudal(fs,mu0,mus);
                        for (std::size_t n = 0; n < size; n++)
                        {
                            y0[n] = (f0[n]-fs[n]);
                        }
                    }
                    for (std::size_t n = 0; n < size; n++)
                    for (std::size_t m = 0; m < size; m++)
                    {
                        J1[n][m] = J0[n][m];
                    }
                    T norm = 0;
                    for (std::size_t n = 0; n < size; n++)
                    {
                        norm += p0[n]*p0[n];
                    }
                    norm *= del;
                    if (norm != 0.0)
                    {
                        for (std::size_t n = 0; n < size; n++)
                        {
                            T an = 0;
                            for (std::size_t m = 0; m < size; m++)
                            {
                                an += J1[n][m]*p0[m];
                            }
                            an *= del;
                            for (std::size_t m = 0; m < size; m++)
                            {
                                J0[n][m] += (y0[n] - an)*p0[m]/norm;
                            }
                        }
                    }
                    else calculate_jacobian(J0,mu0,f0);
                    calculate_inverted(Jinv0,J0);
                }*/

                // Update inverted Jacobian
                /*{
                    for (std::size_t n = 0; n < size; n++)
                    for (std::size_t m = 0; m < size; m++)
                    {
                        Jinv1[n][m] = Jinv0[n][m];
                    }

                    T norm = 0;
                    for (std::size_t n = 0; n < size; n++)
                    for (std::size_t m = 0; m < size; m++)
                    {
                        norm -= p0[m]*Jinv1[m][n]*y0[n];
                    }
                    if (norm != 0.0)
                    {
                        for (std::size_t n = 0; n < size; n++)
                        {
                            T an = 0;
                            for (std::size_t m = 0; m < size; m++)
                            {
                                an += del*p0[n]-Jinv1[n][m]*y0[m];
                            }
                            for (std::size_t m = 0; m < size; m++)
                            {
                                T bm = 0;
                                for (std::size_t o = 0; o < size; o++)
                                {
                                    bm -= p0[m]*Jinv1[m][o];
                                }
                                //Jinv0[n][m] += a[n]*b[m]/norm;
                                Jinv0[n][m] += an*bm/norm;
                            }
                        }
                    }
                    else
                    {
                        calculate_step_size(dx,f0);
                        calculate_jacobian(J0,x0,f0,dx);
                        calculate_inverted(Jinv0,J0);
                    }
                }*/

                //print(itr,mu0,f0,J0,Jinv0,f_norm0);

                update_mu();
                itr++; if (itr>=iterations)  throw std::runtime_error(MaxIterations);
            }
        }
    }
};
}
#endif
