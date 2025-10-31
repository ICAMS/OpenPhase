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

 *   File created :   2019
 *   Main contributors :   Marvin Tegeler, Raphael Schiedung
 *
 */

#include "Containers/SparseMatrix.h"
#include <algorithm>
#include <cassert>
#include <fenv.h>
#include <float.h>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <time.h>
#include <vector>

namespace openphase
{

void SparseMatrix::Allocate(const size_t n)
{
    size = n;
    reserve = 0;
    storage.resize(n);
}

void SparseMatrix::Allocate(const size_t n, const size_t m)
{
    size = n;
    reserve = m;
    storage.resize(n);
    for (auto it : storage) it.reserve(m);
}

double& SparseMatrixCSR::operator()(const size_t i, const size_t j)
{
    assert(i < storage.size() && "SparseMatrix::operator(): Access beyond matrix size ");

    for (auto it = storage[i].begin(); it != storage[i].end(); it++)
    {
        if (it->index == j) return it->value;
    }

    storage[i].push_back(entry({j,0.0}));
    return storage[i].back().value;
}

double& SparseMatrixCSC::operator()(const size_t i, const size_t j)
{
    assert(i < storage.size() && "SparseMatrix::operator(): Access beyond matrix size ");

    for (auto it = storage[j].begin(); it != storage[j].end(); it++)
    {
        if (it->index == i) return it->value;
    }

    storage[j].push_back(entry({i,0.0}));
    return storage[j].back().value;
}

double SparseMatrixCSR::operator()(const size_t i, const size_t j) const
{
    assert(i < storage.size() && "SparseMatrix::operator(): Access beyond matrix size ");

    for (auto it : storage[i]) if (it.index == j) return it.value;
    return 0.0;
}

double SparseMatrixCSC::operator()(const size_t i, const size_t j) const
{
    assert(i < storage.size() && "SparseMatrix::operator(): Access beyond matrix size ");

    for (auto it : storage[j]) if (it.index == i) return it.value;
    return 0.0;
}

//SparseMatrix& SparseMatrix::operator=(const SparseMatrix& rhs)
//{
//    for (size_t i = 0; i < size; ++i)
//    {
//        Matrix[i].clear();
//    }
//
//    #pragma omp parallel for
//    for (size_t i = 0; i < size; ++i)
//    for (auto it : rhs.Matrix[i])
//    {
//        if (it.value != zero)
//        {
//            Matrix[i].push_back(it);
//        }
//    }
//    return *this;
//}

SparseMatrix SparseMatrix::operator*(const double& scalar) const
{
    SparseMatrix tmp(*this);
    #pragma omp parallel for
    for (size_t i = 0; i < tmp.storage.size(); ++i)
    {
        for (auto& it : tmp.storage[i]) it.value *= scalar;
    }
    return tmp;
}

SparseMatrix& SparseMatrix::operator*(const double& scalar)
{
    #pragma omp parallel for
    for (size_t i = 0; i < storage.size(); ++i)
    {
        for (auto& it : storage[i]) it.value *= scalar;
    }
    return *this;
}

SparseMatrixCSR SparseMatrixCSR::operator*(const SparseMatrixCSR& rMatrix) const
{
    SparseMatrixCSR tmp(storage.size(),reserve);
    #pragma omp parallel for
    for (size_t i = 0; i < storage.size(); ++i)
    for (size_t j = 0; j < storage.size(); ++j)
    {
        for (auto k = storage[i].cbegin(); k != storage[i].cend(); ++k)
        {
            if (k->value*rMatrix(k->index,j) != 0.)
            {
                tmp(i,j) += k->value*rMatrix(k->index,j);
            }
        }
    }
    return tmp;
}

SparseMatrixCSC operator*(const SparseMatrixCSC& lhs, const SparseMatrixCSC& rhs)
{
    SparseMatrixCSC tmp(rhs.storage.size(),rhs.reserve);
    #pragma omp parallel for
    for (size_t i = 0; i < rhs.storage.size(); ++i)
    for (size_t j = 0; j < rhs.storage.size(); ++j)
    {
        for (auto k = rhs.storage[j].cbegin(); k != rhs.storage[j].cend(); ++k)
        {
            if (lhs(i,k->index)*k->value != 0.)
            {
                tmp(i,j) += lhs(i,k->index)*k->value;
            }
        }
    }
    return tmp;
}

SparseMatrixCSR& SparseMatrixCSR::operator=(const SparseMatrixCSC& rhs)
{
    SparseMatrixCSR tmp(*this);
    #pragma omp parallel for
    for (size_t j = 0; j < rhs.storage.size(); ++j)
    for (auto i = rhs.storage[j].cbegin(); i != rhs.storage[j].cend(); i++)
    {
        //NOTE: In debug mode the assignment may throw an exception if the sizes mismatch.
        //for *this to be unchanged in case of an exception tmp is used!
        tmp(i->index,j) = i->value;
    }
    *this = tmp;
    return *this;
}

SparseMatrixCSC& SparseMatrixCSC::operator=(const SparseMatrixCSR& rhs)
{
    SparseMatrixCSC tmp(*this);
    #pragma omp parallel for
    for (size_t i = 0; i < rhs.storage.size(); ++i)
    for (auto j = rhs.storage[i].cbegin(); j != rhs.storage[i].cend(); j++)
    {
        //NOTE: In debug mode the assignment may throw an exception if the sizes mismatch.
        //for *this to be unchanged in case of an exception tmp is used!
        tmp(i,j->index) = j->value;
    }
    *this = tmp;
    return *this;
}

SparseMatrixCSR SparseMatrixCSR::transposed() const
{
    SparseMatrixCSR tmp(storage.size(),reserve);
    //#pragma omp parallel for
    for (size_t i = 0; i < storage.size(); ++i)
    {
        for (auto j = storage[i].cbegin(); j != storage[i].cend(); ++j)
        {
            tmp(j->index,i) = j->value;
        }
    }
    return tmp;
}

SparseMatrixCSC SparseMatrixCSC::transposed() const
{
    SparseMatrixCSC tmp(storage.size(),reserve);
    //#pragma omp parallel for
    for (size_t j = 0; j < storage.size(); ++j)
    {
        for (auto i = storage[j].cbegin(); i != storage[j].cend(); ++i)
        {
            tmp(j, i->index) = i->value;
        }
    }
    return tmp;
}

SparseMatrixCSR& SparseMatrixCSR::transpose()
{
    *this = (*this).transposed();
    return *this;
}

SparseMatrixCSC& SparseMatrixCSC::transpose()
{
    *this = (*this).transposed();
    return *this;
}

void SparseMatrixCSR::SwapRows(const size_t i, const size_t j)
{
     auto tmp   = storage[i];
     storage[i] = storage[j];
     storage[j] = tmp;
}

void SparseMatrixCSC::SwapRows(const size_t i, const size_t j)
{
    #pragma omp parallel for
    for (size_t k = 0; k < storage.size(); ++k)
    for (auto it = storage[k].begin(); it != storage[k].end(); ++it)
    {
        if (it->index == i) it->index = j;
        else if (it->index == j) it->index = i;
    }
}

void SparseMatrixCSR::DivideRowBy(size_t i, double scalar)
{
    for (auto it = storage[i].begin(); it != storage[i].end(); ++it)
    {
       it->value/=scalar;
    }
}

void SparseMatrixCSC::DivideRowBy(size_t i, double scalar)
{
    #pragma omp parallel for
    for (size_t k = 0; k < storage.size(); ++k)
    for (auto it = storage[k].begin(); it != storage[k].end(); ++it)
    if (it->index == i)
    {
        it->value/=scalar;
        break;
    }
}

void SparseMatrix::Sort()
{
    for (size_t i = 0; i < storage.size(); i++)
    {
        //auto dis  = [&i](size_t x){return (x > i) ? x - i : i - x;};
        //auto func = [&dis](entry& a, entry& b) {return dis(a.index) < dis(b.index);};
        auto func = [](entry& a, entry& b) {return a.index < b.index;};
        std::sort(storage[i].begin(), storage[i].end(), func);
    }
}

void SparseMatrix::clear()
{
    for (size_t i = 0; i < storage.size(); i++)
    {
        storage[i].clear();
    }
}

size_t SparseMatrixCSR::MaxRow(const size_t i, const size_t j)
{
    assert(i < storage.size() && "MaxRow::operator(): Access beyond matrix size");

    size_t max_i = i;
    double max_value = 0.0;

    #pragma omp parallel
    {
        size_t loc_max_i = max_i;
        double loc_max_value = max_value;

        #pragma omp for
        for (size_t k = i; k < storage.size(); k++)
        for (auto it : storage[k])
        if (it.index == j)
        {
            if (std::abs(it.value) > std::abs(loc_max_value))
            {
                loc_max_i = k;
                loc_max_value = it.value;
            }
            break;
        }

        #pragma omp critical
        {
            if (std::abs(loc_max_value) > std::abs(max_value))
            {
                max_i = loc_max_i;
                max_value = loc_max_value;
            }
        }
    }
    return max_i;
}

size_t SparseMatrixCSR::NonZeroRow(const size_t i, const size_t j)
{
    assert(i < storage.size() && "NonZeroRow::operator(): Access beyond matrix size");

    for (size_t k = i; k < storage.size(); k++)
    for (auto it : storage[k])
    if (it.index == j)
    if (it.value != 0.0) return it.index;

    return i;
}

size_t SparseMatrixCSC::MaxRow(const size_t i, const size_t j)
{
    assert(i < storage.size() && "MaxRow::operator(): Access beyond matrix size");

    size_t max_i = i;
    double max_value = 0.0;

    for (auto it = storage[j].cbegin(); it != storage[j].cend(); ++it)
    if (it->index >= i)
    {
        if (std::abs(it->value) > std::abs(max_value))
        {
            max_i = it->index;
            max_value = it->value;
        }
    }
    return max_i;
}

size_t SparseMatrixCSC::NonZeroRow(const size_t i, const size_t j)
{
    assert(i < storage.size() && "NonZeroRow::operator(): Access beyond matrix size");

    for (auto it = storage[j].cbegin(); it != storage[j].cend(); ++it)
    if (it->index >= i)
    {
        if (it->value != 0.0) return it->index;
    }
    return i;
}

bool operator==(const SparseMatrixCSR& lhs, const SparseMatrixCSC& rhs)
{
    bool eq = true;
    for (size_t i = 0; i < lhs.storage.size(); i++)
    for (auto j = lhs.storage[i].cbegin(); j != lhs.storage[i].cend(); j++)
    {
        eq = eq and (j->value == rhs(i,j->index));
    }

    for (size_t j = 0; j < rhs.storage.size(); j++)
    for (auto i = rhs.storage[j].cbegin(); i != rhs.storage[j].cend(); i++)
    {
        eq = eq and (lhs(i->index,j) == i->value);
    }

    return eq;
}

//void BiCGStab::Solve_Unpreconditioned(Vector& x, const Matrix& A, const Vector& b, const double Accuracy)
//{
//    //Source (08.07.2020): https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
//    //Contact: raphael.schiedung@rub.de
//    //
//
//    const size_t size = b.size();
//    std::vector<double> r_0 = b-A*x;
//    std::vector<double> r_i(r_0);
//    std::vector<double> x_i = x;
//
//    double alpha     = 1.0;
//    double omega_i   = 1.0;
//    double rho_i     = 1.0;
//
//    std::vector<double> v_i (size,0.0);
//    std::vector<double> p_i (size,0.0);
//
//    auto div = [] (const double& nom, const double& div){return nom/div;};//(div > DBL_EPSILON) ? nom/div : DBL_MAX;};
//
//    double res = std::sqrt(r_0*r_0)/size;
//    std::cout << "Iteration: 0\tResiduum: " << res  << "\n";
//    if (res > Accuracy)
//    for (size_t i = 1; i < size; i++)
//    {
//        const double rho_im1 = rho_i;
//        rho_i = r_0*r_i;
//        const double beta  = (rho_i/rho_im1)*(alpha/omega_i);
//        const std::vector<double> p_im1 (p_i);
//        p_i = r_i + beta*(p_im1 - omega_i*v_i);
//        v_i = A*p_i;
//        alpha = div(rho_i,(r_0*v_i));
//        const std::vector<double> h = x_i + alpha*p_i;
//        res = std::sqrt((b-A*h)*(b-A*h))/size;
//        std::cout << "Iteration: " << i-1 << ".5\tResiduum: " << res  << "\n";
//        if (res < Accuracy) {x_i = h; break;};
//        const std::vector<double> s = r_i - v_i*alpha;
//        const std::vector<double> t = A*s;
//        omega_i = div(t*s,t*t);
//        x_i = h + omega_i*s;
//        r_i = s - omega_i*t;
//        res = std::sqrt(r_i*r_i)/size;
//        std::cout << "Iteration: " << i << ".0\tResiduum: " << res  << "\n";
//        if (res < Accuracy) break;
//    }
//    x = x_i;
//}

//void BiCGStab::BiCGstab_Jacobi(const Matrix& A, const Vector& b, Vector& x, const size_t l, const double TOL)
//{
//    std::cout << "Starte BiCGstab(" << l << ")" << std::endl;
//    std::vector<Vector> u(l + 1);
//    std::vector<Vector> r(l + 1);
//    std::vector<double> sigma(l + 1);
//    std::vector<double> gamma(l + 1);
//    std::vector<double> gamma_(l + 1);
//    std::vector<double> gamma__(l + 1);
//    Matrix Kinv(x.size());
//    for (size_t i = 0; i<x.size();++i)
//    {
//        Kinv(i,i) = 1./A(i,i);
//    }
//
//    std::vector<std::vector<double> > tau(l + 1);
//    for (size_t i = 0; i < l + 1; ++i)
//    {
//        tau[i].resize(l + 1);
//    }
//    u[0].resize(x.size(),0.);
//    long int k = -l;
//    auto x0 = x;
//    auto KA = Kinv*A;
//    r[0] = Kinv*b - KA * x0;
//    Vector r0_(r[0].size(),1.);
//    double rho0 = 1.;
//    double rho1;
//    double alpha = 0.;
//    double omega = 1.;
//    double res = (r[0]* r[0]);
//    double res0 = (r[0]* r[0]);
//    while (res > TOL * TOL*res0)
//    {
//        k=k+l;
//        rho0 = -omega * rho0;
//        size_t ll = l > 1 ? l -1 : 0;
//        for (size_t j = 0; j <= ll; ++j)
//        {
//            rho1 = (r[j]* r0_);
//            double beta = alpha * rho1 / rho0;
//            rho0 = rho1;
//            for (size_t i = 0; i <= j; ++i)
//            {
//                u[i] = r[i] - beta * u[i];
//            }
//            u[j + 1] = KA * u[j];
//            double gamma = (u[j + 1]* r0_);
//            alpha = rho0 / gamma;
//            for (size_t i = 0; i <= j; ++i)
//            {
//                r[i] -= alpha * u[i + 1];
//            }
//            r[j + 1] = KA * r[j];
//            x0 += alpha * u[0];
//        }
//        for (size_t j = 1; j <= l; ++j)
//        {
//            for (size_t i = 1; i <= j - 1; ++i)
//            {
//                tau[i][j] = (1. / sigma[i]) * (r[j]* r[i]);
//                r[j] -= tau[i][j] * r[i];
//            }
//            sigma[j] = (r[j]* r[j]);
//            gamma_[j] = (1. / sigma[j]) * (r[0]* r[j]);
//        }
//        gamma[l] = gamma_[l];
//        omega = gamma[l];
//        for (size_t j = ll; j >= 1; --j)
//        {
//            double sum = 0;
//            for (size_t i = j + 1; i <= l; ++i)
//            {
//                sum += tau[j][i] * gamma[i];
//            }
//            gamma[j] = gamma_[j] - sum;
//        }
//        for (size_t j = 1; j <= ll; ++j)
//        {
//            double sum = 0;
//            for (size_t i = j + 1; i <= ll; ++i)
//            {
//                sum += tau[j][i] * gamma[i + 1];
//            }
//            gamma__[j] = gamma[j + 1] + sum;
//        }
//        x0 += gamma[1] * r[0];
//        r[0] -= gamma_[l] * r[l];
//        u[0] -= gamma[l] * u[l];
//        for (size_t j = 1; j <= ll; ++j)
//        {
//            u[0] -= gamma[j]*u[j];
//            x0 += gamma__[j] * r[j];
//            r[0] -= gamma_[j] * r[j];
//        }
//        x = x0;
//        res = (r[0]* r[0]);
//        std::cout << sqrt(res) << std::endl;
//    }
//    std::cout << "Iterations: " << k << "\n";
//}
}//namespace openphase
