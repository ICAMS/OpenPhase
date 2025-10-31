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
 *   Main contributors :   Muhammad Adil Ali; Oleg Shchyglo
 *
 */

#ifndef DVECTORN_H
#define DVECTORN_H

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <numeric>
#include <vector>

namespace openphase
{

class dVectorN
{
 public:
    dVectorN()
    {
    };
    dVectorN(const size_t size): storage(size, 0.0)
    {
    };
    dVectorN(const size_t size, const double value): storage(size, value)
    {
    };
    dVectorN(const dVectorN& other): storage(other.storage)
    {
    };
    dVectorN(const std::vector<double>& vec) : storage(vec)
    {
    };
    dVectorN(std::initializer_list<double> vecinit)
    {
        storage.resize(vecinit.size());
        std::transform(vecinit.begin(), vecinit.end(), storage.begin(),
                       [](double val) { return val; });
    }
    void set_to_zero()
    {
        std::fill(storage.begin(), storage.end(), 0.0);
    };
    void set_to_value(const double value)
    {
        std::fill(storage.begin(), storage.end(), value);
    };

    void Allocate(const size_t size)
    {
        if(storage.size() != size)
        {
            storage.resize(size);
            set_to_zero();
        }
    };
    void push_back(const double value)
    {
        storage.push_back(value);
    }
    [[nodiscard]]
    double norm() const
    {
        return std::sqrt(std::inner_product(storage.begin(), storage.end(), storage.begin(), 0.0));
    }
    [[nodiscard]]
    double sum() const
    {
        return std::accumulate(storage.begin(), storage.end(), 0.0);
    }
    [[nodiscard]]
    double max() const
    {
       return *std::max_element(storage.begin(), storage.end());
    }
    [[nodiscard]]
    double absMax() const
    {
        double maxValue = 0;
        for(size_t i = 0; i < storage.size(); i++)
        {
            maxValue = std::max(maxValue, std::abs(storage[i]));
        }
        return maxValue;
    }
    [[nodiscard]]
    double Average() const
    {
        return storage.empty() ? 0.0 : sum() / storage.size();
    }
    void normalize()
    {
        double dnorm = norm();
        if (dnorm != 0)
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= dnorm;
        }
    }
    dVectorN& pow(double exponent) 
    {
        for (double& element : storage) {
            element = std::pow(element, exponent);
        }
        return *this;
    }
    dVectorN& sqrt() 
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::sqrt(val);
        });
        return *this;
    }
    [[nodiscard]]
    dVectorN sqrted() const 
    {
        dVectorN result(storage.size());
        std::transform(storage.begin(), storage.end(), result.begin(), [](double val) 
        {
            return std::sqrt(val);
        });
        return result;
    }
    dVectorN& tanh() 
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::tanh(val);
        });
        return *this;
    }
    dVectorN& fabs() 
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::fabs(val);
        });
        return *this;
    }
    [[nodiscard]]
    dVectorN fabsd()
    {
        dVectorN absStorage(size());
        std::transform (begin(), end(), absStorage.begin(), [](double i) { return std::fabs(i); });
        return absStorage;
    }
    [[nodiscard]]
    double dot(const dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "dVectorN::dot Sizes of the vectors are not equal");
        return std::inner_product(storage.begin(), storage.end(), rhs.storage.begin(), 0.0);
    }
    dVectorN& min(const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::min Sizes of the vectors are not equal");
        std::transform(begin(), end(), rhs.cbegin(), begin(), 
                       [](double a, double b) { return std::min(a, b); });
        return *this;
    }
    [[nodiscard]]
    dVectorN minimized(const dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "dVectorN::max Sizes of the vectors are not equal");
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), rhs.cbegin(), result.begin(), 
                       [](double a, double b) { return std::min(a, b); });
        return result;
    }
    dVectorN& max(const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::max Sizes of the vectors are not equal");
        std::transform(begin(), end(), rhs.cbegin(), begin(), 
                       [](double a, double b) { return std::max(a, b); });
        return *this;
    }
    [[nodiscard]]
    dVectorN maximized(const dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "dVectorN::maximized Sizes of the vectors are not equal");
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), rhs.cbegin(), result.begin(), 
                       [](double a, double b) { return std::max(a, b); });
        return result;
    }
    dVectorN& min(double scalar)
    {
        std::transform(begin(), end(), begin(), [scalar](double a) {
            return std::min(a, scalar);
        });
        return *this;
    }
    [[nodiscard]]
    dVectorN minimized(double scalar) const 
    {
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), result.begin(), [scalar](double a) {
            return std::min(a, scalar);
        });

        return result;
    }
    dVectorN& max(double scalar)
    {
        std::transform(begin(), end(), begin(), [scalar](double a) {
            return std::max(a, scalar);
        });
        return *this;
    }
    [[nodiscard]]
    dVectorN maximized(double scalar) const 
    {
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), result.begin(), [scalar](double a) {
            return std::max(a, scalar);
        });

        return result;
    }
    double& operator[](const size_t i)
    {
        assert(i < storage.size() && "dVectorN::[] Access beyond storage range");
        return storage[i];
    };
    double const& operator[](const size_t i) const
    {
        assert(i < storage.size() && "dVectorN::[] Access beyond storage range");
        return storage[i];
    };
    dVectorN& operator=(const dVectorN& rhs)
    {
        storage = rhs.storage;
        return *this;
    };
    [[nodiscard]]
    dVectorN operator*(const double m) const
    {
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i]*m;
        }
        return tmp;
    };
    [[nodiscard]]
    dVectorN operator/(const double m) const
    {
        assert(m != 0 && "dVectorN::/ divided by zero");
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i]/m;
        }
        return tmp;
    };
    dVectorN& operator*=(const double m)
    {
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] *= m;
        }
        return *this;
    };
    dVectorN& operator/=(const double m)
    {
        assert(m != 0 && "dVectorN::/= divided by zero");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= m;
        }
        return *this;
    };
    [[nodiscard]]
    dVectorN operator+(const dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "dVectorN::+ Sizes of the vectors are not equal");
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i] + rhs[i];
        }
        return tmp;
    };
    [[nodiscard]]
    dVectorN operator-(const dVectorN& rhs) const
    {
        assert(rhs.size() == storage.size() && "dVectorN::- Sizes of the vectors are not equal");
        dVectorN tmp(storage.size(),0.0);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] = storage[i] - rhs[i];
        }
        return tmp;
    };
    [[nodiscard]]
    dVectorN operator*(const dVectorN& rhs) const
    {
        dVectorN result = *this;
        result *= rhs;
        return result;
    }

    [[nodiscard]]
    dVectorN operator/(const dVectorN& rhs) const
    {
        dVectorN result = *this;
        result /= rhs;
        return result;
    }
    dVectorN& operator+=( const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::+= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] + rhs[i];
        }
        return *this;
    };
    dVectorN& operator-=(const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::-= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] - rhs[i];
        }
        return *this;
    };
    dVectorN& operator/=(const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::/= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] / rhs[i];
        }
        return *this;
    };
    dVectorN& operator*=(const dVectorN& rhs)
    {
        assert(rhs.size() == storage.size() && "dVectorN::*= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] = storage[i] * rhs[i];
        }
        return *this;
    };

    static dVectorN ZeroVector(size_t size) 
    {
        dVectorN unity(size);
        unity.set_to_value(0.0);
        return unity;
    }

    static dVectorN UnitVector(size_t size) 
    {
        dVectorN unity(size);
        unity.set_to_value(1.0);
        return unity;
    }

    std::string print(const char separator = ' ', const size_t precision = 6) const
    {
        std::stringstream out;
        for (const auto& elem : storage)
        {
            out << std::setprecision(precision) << elem << separator;
        }
        return out.str();
    }

    [[nodiscard]]
    size_t size(void) const
    {
        return storage.size();
    }
    double* data(void)
    {
        return storage.data();
    }
    const double* data(void) const
    {
        return storage.data();
    }

    void Read(std::istream& inp)
    {
        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
        storage.resize(size);
        inp.read(reinterpret_cast<char*>(storage.data()), size * sizeof(double));
    }
    void Write(std::ostream& outp) const
    {
        size_t size = storage.size();
        outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(storage.data()), size * sizeof(double));
    }

    typedef typename std::vector<double>::iterator iterator;                    ///< Iterator over the vector fields
    typedef typename std::vector<double>::const_iterator citerator;             ///< Constant iterator over the vector fields
    iterator   begin()       { return storage.begin();  }                       ///< Iterator to the begin of vector fields
    iterator     end()       { return storage.end();    }                       ///< Iterator to the end of vector fields
    citerator cbegin() const { return storage.cbegin(); }                       ///< Constant iterator to the begin of vector fields
    citerator   cend() const { return storage.cend();   }                       ///< Constant iterator to the end of vector fields
    //size_t      size()       { return storage.size();   }
 protected:
 private:
    std::vector<double> storage;
};

} // namespace openphase
#endif
