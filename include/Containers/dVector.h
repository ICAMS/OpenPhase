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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich;
 *                         Dmitry Medvedev; Philipp Engels
 *
 */

#ifndef DVECTOR_H
#define DVECTOR_H

#include <cassert>

namespace openphase
{

template<int N>
class dVector
{
 public:

    dVector()
    {
        set_to_zero();
    }
    dVector(const dVector<N>& vecinit)
    {
        storage = vecinit.storage;
    }
    dVector(std::initializer_list<double> vecinit)
    {
        assert(vecinit.size() == N && "Initialization list size is not equal to storage range.");
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }
    double& operator[](const size_t i)
    {
        assert(i < N && "Access beyond storage range");
        return storage[i];
    };
    double const& operator[](const size_t i) const
    {
        assert(i < N && "Access beyond storage range");
        return storage[i];
    };
    void set_to_zero(void)
    {
        storage.fill(0.0);
    };
    void set_to_value(const double value)
    {
        for(int i = 0; i < N; i++) storage[i] = value;
    };
    dVector<N> operator*(const double m) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i]*m;
        }
        return tmp;
    };
    dVector<N>& operator*=(const double m)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] *= m;
        }
        return *this;
    };
    dVector<N> operator+(const double rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] + rhs;
        }
        return tmp;
    };
    dVector<N> operator+(const dVector<N>& rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] + rhs[i];
        }
        return tmp;
    };
    dVector<N>& operator+=(const dVector<N>& rhs)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] += rhs[i];
        }
        return *this;
    };
    dVector<N> operator-(const double rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] - rhs;
        }
        return tmp;
    };
    dVector<N> operator-(const dVector<N>& rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] - rhs[i];
        }
        return tmp;
    };
    dVector<N>& operator-=(const dVector<N>& rhs)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] -= rhs[i];
        }
        return *this;
    };
    dVector<N>& operator=(const dVector<N>& rhs)
    {
        storage = rhs.storage;
        return *this;
    };
    double sum_of_entries(void) const
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            sum += storage[i];
        }
        return sum;
    };
    double average_entries(void) const
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            sum += storage[i];
        }
        return sum/double(N);
    };
    /*double* data(void)
    {
        return storage.data();
    };
    const double* const_data(void) const
    {
        return storage;
    };*/
    std::string print(void) const
    {
        std::stringstream out;
        out << "||";
        for(int i = 0; i < N-1; i++)
        {
            out << storage[i] << ", ";
        }
        out << storage[N-1] << "||";
        return out.str();
    };

 protected:
 private:
    std::array<double, N> storage;
};

} // namespace openphase
#endif
