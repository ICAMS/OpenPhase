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

 *   File created :   2020
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef SOLVELINEARSYSTEM_H
#define SOLVELINEARSYSTEM_H

#include <stddef.h>
#include <cmath>
#include <iostream>
#include <exception>
#include <cassert>

#include "ConsoleOutput.h"

namespace openphase::SystemOfLinearEquationsSolvers ///< Functions for solving linear system of equations
{

template<typename Vector>
Vector OP_EXPORTS operator+(const Vector & lhs, const Vector& rhs)
{
    assert(lhs.size() == rhs.size());

    Vector tmp(lhs.size());
    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        tmp[i] = lhs[i]+rhs[i];
    }
    return tmp;
}

template<typename Vector>
Vector OP_EXPORTS operator+=(Vector& lhs, const Vector& rhs)
{
    assert(lhs.size() == rhs.size());

    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        lhs[i] += rhs[i];
    }
    return lhs;
}

template<typename Vector>
Vector OP_EXPORTS operator-(const Vector& lhs, const Vector& rhs)
{
    assert(lhs.size() == rhs.size());

    Vector tmp(lhs.size());
    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        tmp[i] = lhs[i]-rhs[i];
    }
    return tmp;
}

template<typename Vector>
Vector OP_EXPORTS operator-=(Vector& lhs, const Vector& rhs)
{
    assert(lhs.size() == rhs.size());

    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        lhs[i] -= rhs[i];
    }
    return lhs;
}

template<typename Vector>
double OP_EXPORTS operator*(const Vector& lhs, const Vector& rhs)
{
    assert(lhs.size() == rhs.size());

    double tmp = 0.;
    #pragma omp parallel for reduction(+:tmp)
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        tmp += lhs[i]*rhs[i];
    }
    return tmp;
}

template<typename Vector>
Vector OP_EXPORTS operator*(const double lhs, const Vector& rhs)
{
    Vector tmp(rhs.size());
    #pragma omp parallel for
    for (size_t i = 0; i < rhs.size(); ++i)
    {
        tmp[i] = lhs*rhs[i];
    }
    return tmp;
}

template<typename Vector>
Vector OP_EXPORTS operator*(const Vector& lhs, const double rhs)
{
    Vector tmp(lhs.size());
    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        tmp[i] = lhs[i]*rhs;
    }
    return tmp;
}

template<typename Vector>
Vector OP_EXPORTS operator/(const Vector& lhs, const double rhs)
{
    Vector tmp(lhs.size());
    #pragma omp parallel for
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        tmp[i] = lhs[i]/rhs;
    }
    return tmp;
}

template<class Vector>
inline double OP_EXPORTS  max(const Vector& x)
{
    const size_t N = x.size();
    double MaxValue = 0.0;
    #pragma omp parallel for reduction(max:MaxValue)
    for (size_t i = 0; i < N; i++)
    {
         if (std::abs(x[i]) > MaxValue) MaxValue = std::abs(x[i]);
    }
    return MaxValue;
}

template<class Matrix, class Vector>
inline Vector OP_EXPORTS  residual(const Matrix& A, Vector& x, const Vector& b)
{
    return (b-A*x);
}

/// Swaps rows i an j for Ax=b
template<class Matrix, class Vector>
void OP_EXPORTS swap(size_t i, size_t j, Matrix& A, Vector& b)
{
    A.SwapRows(i,j);
    double tmp = b[j];
    b[j] = b[i];
    b[i] = tmp;
}

/// Pivoting swaps rows so that the largest elements are on the diagonal A_ij
template<class Matrix, class Vector>
void OP_EXPORTS Pivote(Matrix& A, Vector& b)
{
    size_t N = b.size();

    size_t h = 0;
    size_t k = 0;
    while (h < N and k < N)
    {
        size_t i_max = A.MaxRow(h,k);

        if (A(i_max,k) == 0.0)
        {
            k++;
        }
        else
        {
            swap(i_max,h,A,b);
        }
        h++;
        k++;
    }
}

/// Pivoting swaps rows so that nonzero are on the diagonal A_ij
template<class Matrix, class Vector>
void OP_EXPORTS FastPivote(Matrix& A, Vector& b)
{
    size_t N = b.size();

    size_t h = 0;
    size_t k = 0;
    while (h < N and k < N)
    {
        size_t i_max = A.NonZeroRow(h,k);

        if (A(i_max,k) == 0.0)
        {
            k++;
        }
        else
        {
            swap(i_max,h,A,b);
        }
        h++;
        k++;
    }
}

template<class Matrix>
double OP_EXPORTS SpectralNumberMAX(Matrix& A, size_t N)
{
    double AMAX = 0.0;
    double AMIN = std::numeric_limits<double>::max();
    #pragma omp parallel
    {
        #pragma omp for reduction(max:AMAX)
        for (size_t i = 0; i < N; i++)
        {
            if (std::abs(A(i,i)) > AMAX) AMAX = A(i,i);
        }
        #pragma omp for reduction(min:AMIN)
        for (size_t i = 0; i < N; i++)
        {
            if (std::abs(A(i,i)) < AMIN) AMIN = A(i,i);
        }
    }
    if (AMIN > std::numeric_limits<double>::epsilon()) return AMAX/AMIN;
    else return std::numeric_limits<double>::max();
}

template<class Matrix, class Vector>
void OP_EXPORTS PreconditionJacobi(Matrix& A, Vector& b)
{
    size_t N = b.size();
    #pragma omp parallel for
    for (size_t i = 0; i < N; i++)
    {
        const double norm = A(i,i);
        A.DivideRowBy(i,norm);
        b[i]/=norm;
    }
}

template<class Matrix, class Vector>
void OP_EXPORTS Gauss(Matrix& A, Vector& x, Vector& b)
{
    // Source (09.07.2020): https://en.wikipedia.org/wiki/Gaussian_elimination
    const size_t N = x.size();

    size_t h = 0;
    size_t k = 0;
    while (h < N and k < N)
    {
        #pragma omp parallel for //NOTE parallelism will slow algorithm
        for (size_t i = 0; i < N; i++)
        if (i!=h)
        {
            if (A(i,k) != 0.0)
            {
                double f = A(i,k)/A(h,k);
                A(i,k) = 0.0;

                for (size_t j = k+1; j < N; j++)
                {
                    A(i,j) -= A(h,j)*f;
                }
                b[i] -= b[h]*f;
            }
        }

        h++;
        k++;
#ifdef DEBUG
        std::cout << "Element: (" << h << "," << k << ")/(" << N << "," << N << ")\n";
#endif
    }

    for (size_t i = 0; i < N; i++)
    {
        x[i] = b[i]/A(i,i);
    }
};

template<class Matrix, class Vector>
size_t OP_EXPORTS Jacobi(const Matrix& A, Vector& x, const Vector& b, const double MaxResidual)
{
    // Source (08.07.2020): https://en.wikipedia.org/wiki/Jacobi_method
    const size_t N = x.size();
    Vector xOld(x);
    for (size_t k = 0; k < N*N*N; k++)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < N; i++)
        {
            x[i] = b[i];

            for (size_t j = 0; j < N; j++)
            if (i!=j)
            {
                x[i] -= A(i,j)*xOld[j];
            }
            x[i]/=A(i,i);
        }

        #pragma omp parallel for
        for (size_t i = 0; i < N; i++)
        {
            xOld[i] = x[i];
        }

        const Vector r = residual(A,x,b);
        //const double locResidual = max(r);
        const double locResidual = sqrt(r*r);
#ifdef DEBUG
        std::cout << "Iteration: " << k << "\tResidual: " << locResidual << "\n";
#endif
        if (locResidual < MaxResidual) return k;
    }
    ConsoleOutput::WriteWarning("Solution did not converge", "SolveLinearSystem", "Jacob");
    return N*N*N;
};

template<class Matrix, class Vector>
size_t OP_EXPORTS GaussSeidel(const Matrix& A, Vector& x, const Vector& b, const double MaxResidual)
{
    //source (13.07.2020): http://www.netlib.org/linalg/html_templates/node14.html#SECTION00722000000000000000
    const size_t N = x.size();
    Vector xOld(x);
    for (size_t k = 0; k < N*N*N; k++)
    {
        for (size_t i = 0; i < N; i++)
        {
            double sigma = 0.0;
            //#pragma omp parallel for //NOTE parallelism will slow algorithm
            for (size_t j = 0; j < i; j++)
            {
                sigma += A(i,j)*x[j];
            }

            //#pragma omp parallel for //NOTE parallelism will slow algorithm
            for (size_t j = i+1; j < N; j++)
            {
                sigma += A(i,j)*xOld[j];
            }

            x[i] = (b[i] - sigma)/A(i,i);
        }

        //#pragma omp parallel for //NOTE parallelism will slow algorithm
        for (size_t i = 0; i < N; i++)
        {
            xOld[i] = x[i];
        }

        const Vector r = residual(A,x,b);
        //const double locResidual = max(r);
        const double locResidual = sqrt(r*r);
#ifdef DEBUG
        std::cout << "Iteration: " << k << "\tResidual: " << locResidual << "\n";
#endif
        if (locResidual < MaxResidual) return k;
    }
    ConsoleOutput::WriteWarning("Solution did not converge", "SolveLinearSystem", "GaussSeidel");
    return N*N*N;
};

template<class Matrix, class Vector>
size_t OP_EXPORTS GradientDescent(const Matrix& A, Vector& x, const Vector& b, const double MaxResidual)
{
    const size_t N = x.size();
    for (size_t k = 0; k < N*N*N; k++)
    {
        Vector r = residual(A,x,b);
        //double locResidual = max(r);
        const double locResidual = sqrt(r*r);
#ifdef DEBUG
        std::cout << "Iteration: " << k << "\tResidual: " << locResidual << "\n";
#endif
        if (locResidual < MaxResidual) return k;

        double alpha = (r*b - 0.5*(x*(A*r)+r*(A*x)))/(r*(A*r));
        x += alpha*r;
    }
    ConsoleOutput::WriteWarning("Solution did not converge", "SolveLinearSystem", "GradientDescent");
    return N*N*N;
};

template<class Matrix, class Vector>
size_t OP_EXPORTS ConjugateGradient(const Matrix& A, Vector& x, const Vector& b, const double MaxResidual)
{
    //source: Fletcher, Roger. "Conjugate gradient methods for indefinite systems." Numerical analysis. Springer, Berlin, Heidelberg, 1976. 73-89.
    const size_t N = x.size();

    Vector r = residual(A,x,b);
    Vector p = r;

    if (sqrt(r*r) > MaxResidual)
    for (size_t k = 0; k < N; k++)
    {
        //const double alpha = p*b/(p*(A*p));
        const double alpha = r*r/(p*(A*p));
        const double rr = r*r;
        x += alpha*p;
        r = residual(A,x,b);
        //const double beta = (-1.0)*p*(A*r)/(p*(A*p));
        const double beta = r*r/rr;
        p = r + beta*p;

        //const double locResidual = max(r);
        const double locResidual = sqrt(r*r);
#ifdef DEBUG
        std::cout << "Iteration: " << k << "\tResidual: " << locResidual << "\n";
#endif
        if (locResidual < MaxResidual) return k;
    }
    ConsoleOutput::WriteWarning("Solution did not converge", "SolveLinearSystem", "ConjugateGradient");
    return N;
};

template<class Matrix, class Vector>
size_t OP_EXPORTS BiconjugateGradient(const Matrix& A, Vector& x, const Vector& b, const double MaxResidual)
{
    //source: Fletcher, Roger. "Conjugate gradient methods for indefinite systems." Numerical analysis. Springer, Berlin, Heidelberg, 1976. 73-89.
    const size_t N  = x.size();
    const Matrix AT = A.transposed();

    Vector xb = x;
    Vector ra = (b- A*x );
    Vector rb = (b-AT*xb);
    Vector pa = ra;
    Vector pb = rb;

    //if (sqrt(ra*ra) > MaxResidual)
    if (max(ra) > MaxResidual)
    {
        for (size_t k = 0; k < N; k++)
        {
            const double alpha = rb*ra/(pb*(A*pb));

            x += alpha*pa;
            //xb += alpha*pb;

            const double denom = rb*ra;

            ra -= alpha*(A*pa);
            rb -= alpha*(A*pb);

            const double beta = rb*ra/denom;

            pa = ra + beta*pa;
            pb = rb + beta*pb;

            const Vector r = residual(A,x,b);
            const double locResidual = max(r);
            //const double locResidual = sqrt(r*r);
            #ifdef DEBUG
            std::cout << "Iteration: " << k << "\tResidual: " << locResidual << "\n";
            #endif
            if (locResidual < MaxResidual) return k;
        }
    }
    else
    {
        #ifdef DEBUG
        ConsoleOutput::WriteWarning("No iteration performed!", "SolveLinearSystem", "BiconjugateGradient");
        #endif

        return 0;
    }
    ConsoleOutput::WriteWarning("Solution did not converge", "SolveLinearSystem", "BiconjugateGradient");
    return N;
}

};
#endif
