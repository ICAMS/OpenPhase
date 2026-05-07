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
    dVectorN()                                                                  ///< Default constructor. Creates empty vector.
    {
    }

    dVectorN(const size_t size_in): storage(size_in, 0.0)                       ///< Constructor. Creates the vector of the given size and sets its components to zero.
    {
    }

    dVectorN(const size_t size_in, const double value): storage(size_in, value) ///< Constructor. Creates the vector of the given size and sets its components to the given value.
    {
    }

    dVectorN(const dVectorN& other): storage(other.storage)                     ///< Copy constructor. Initializes the current vector with the copy of the other vector.
    {
    }

    dVectorN(const std::vector<double>& vec) : storage(vec)                     ///< Constructor. Initializes the current vector with the content of the std::vector<>.
    {
    }

    dVectorN(std::initializer_list<double> vecinit)                             ///< Constructor. Initializes the current vector with the content of the std::initializer_list<>.
    {
        storage.resize(vecinit.size());
        std::transform(vecinit.begin(), vecinit.end(), storage.begin(),
                       [](double val) { return val; });
    }

    void Allocate(const size_t size_in)                                         ///< Allocates the storage of the current vector to the given size.
    {
        if(storage.size() != size_in)
        {
            storage.resize(size_in);
            set_to_zero();
        }
    }

    void push_back(const double value)                                          ///< Append new vector component at the back of the current vector. Increases vector size by 1.
    {
        storage.push_back(value);
    }

    static dVectorN ZeroVector(size_t size_in)                                  ///< Returns vector with the given size with all its element set to zero.
    {
        dVectorN zero(size_in);
        zero.set_to_value(0.0);
        return zero;
    }

    static dVectorN UnitVector(size_t size_in)                                  ///< Returns vector with the given size with all its element set to one.
    {
        dVectorN unity(size_in);
        unity.set_to_value(1.0);
        return unity;
    }

    double& operator[](const size_t i)                                          ///< Random access operator. Returns the reference to the vector component pointer to by the index i.
    {
        assert(i < storage.size() && "dVectorN::[] Access beyond storage range");
        return storage[i];
    }

    double const& operator[](const size_t i) const                              ///< Random access operator. Returns const reference to the vector component pointer to by the index i.
    {
        assert(i < storage.size() && "dVectorN::[] Access beyond storage range");
        return storage[i];
    }

    double* data(void)                                                          ///< Returns pointer to the stored data
    {
        return storage.data();
    }

    const double* data(void) const                                              ///< Returns const pointer to the stored data
    {
        return storage.data();
    }

    [[nodiscard]]
    double norm() const                                                         ///< Returns the norm of the current vector.
    {
        return std::sqrt(std::inner_product(storage.begin(), storage.end(), storage.begin(), 0.0));
    }

    [[nodiscard]]
    double sum() const                                                          ///< Returns the sum of all elements of the current vector.
    {
        return std::accumulate(storage.begin(), storage.end(), 0.0);
    }

    [[nodiscard]]
    double max() const                                                          ///< Returns the maximum of all elements of the current vector.
    {
       return *std::max_element(storage.begin(), storage.end());
    }

    [[nodiscard]]
    double absMax() const                                                       ///< Returns the maximum (by the absolute value) of all elements of the current vector.
    {
        double maxValue = 0;
        for(size_t i = 0; i < storage.size(); i++)
        {
            maxValue = std::max(maxValue, std::abs(storage[i]));
        }
        return maxValue;
    }

    [[nodiscard]]
    double Average() const                                                      ///< Returns the average of all elements of the current vector.
    {
        return storage.empty() ? 0.0 : sum() / storage.size();
    }

    void normalize()                                                            ///< Normalizes the current vector.
    {
        double dnorm = norm();
        if (dnorm != 0)
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= dnorm;
        }
    }

    dVectorN& pow(double exponent)                                              ///< Rises all elements of the current vector to the given exponent.
    {
        for (double& element : storage)
        {
            element = std::pow(element, exponent);
        }
        return *this;
    }

    dVectorN& sqrt()                                                            ///< Takes the square roots of all current vector components.
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::sqrt(val);
        });
        return *this;
    }

    [[nodiscard]]
    dVectorN sqrted() const                                                     ///< Returns the copy of the current vector with the square roots of all its components.
    {
        dVectorN result(storage.size());

        std::transform(storage.begin(), storage.end(), result.begin(), [](double val) 
        {
            return std::sqrt(val);
        });
        return result;
    }

    dVectorN& tanh()                                                            ///< Takes hyperbolic tangent of all current vector components.
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::tanh(val);
        });
        return *this;
    }

    dVectorN& fabs()                                                            ///< Applies absolute value to all current vector components.
    {
        std::transform(storage.begin(), storage.end(), storage.begin(), [](double val) 
        {
            return std::fabs(val);
        });
        return *this;
    }

    [[nodiscard]]
    dVectorN fabsd()                                                            ///< Returns the copy of the current vector with the absolute values of all its components.
    {
        dVectorN absStorage(size());
        std::transform (begin(), end(), absStorage.begin(), [](double i) { return std::fabs(i); });
        return absStorage;
    }

    [[nodiscard]]
    double dot(const dVectorN& rhs) const                                       ///< Dot product of the two vectors.
    {
        assert(rhs.size() == storage.size() && "dVectorN::dot Sizes of the vectors are not equal");
        return std::inner_product(storage.begin(), storage.end(), rhs.storage.begin(), 0.0);
    }

    dVectorN& min(const dVectorN& rhs)                                          ///< Replaces the elements of the current vector with the minimum elements between the two vectors.
    {
        assert(rhs.size() == storage.size() && "dVectorN::min Sizes of the vectors are not equal");
        std::transform(begin(), end(), rhs.cbegin(), begin(), 
                       [](double a, double b) { return std::min(a, b); });
        return *this;
    }

    [[nodiscard]]
    dVectorN minimized(const dVectorN& rhs) const                               ///< Returns vector containing minimum elements between the two vectors
    {
        assert(rhs.size() == storage.size() && "dVectorN::max Sizes of the vectors are not equal");
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), rhs.cbegin(), result.begin(), 
                       [](double a, double b) { return std::min(a, b); });
        return result;
    }

    dVectorN& max(const dVectorN& rhs)                                          ///< Replaces the elements of the current vector with the maximum elements between the two vectors.
    {
        assert(rhs.size() == storage.size() && "dVectorN::max Sizes of the vectors are not equal");
        std::transform(begin(), end(), rhs.cbegin(), begin(), 
                       [](double a, double b) { return std::max(a, b); });
        return *this;
    }

    [[nodiscard]]
    dVectorN maximized(const dVectorN& rhs) const                               ///< Returns vector containing maximum elements between the two vectors
    {
        assert(rhs.size() == storage.size() && "dVectorN::maximized Sizes of the vectors are not equal");
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), rhs.cbegin(), result.begin(), 
                       [](double a, double b) { return std::max(a, b); });
        return result;
    }

    dVectorN& min(double scalar)                                                ///< Applies upper value limit, set by the scalar factor, to all elements of the current vector
    {
        std::transform(begin(), end(), begin(), [scalar](double a)
        {
            return std::min(a, scalar);
        });
        return *this;
    }

    [[nodiscard]]
    dVectorN minimized(double scalar) const                                     ///< Returns vector containing the elements of the current vector with upper value limit, set by the scalar factor, applied to all elements
    {
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), result.begin(), [scalar](double a)
        {
            return std::min(a, scalar);
        });

        return result;
    }

    dVectorN& max(double scalar)                                                ///< Applies lower value limit, set by the scalar factor, to all elements of the current vector
    {
        std::transform(begin(), end(), begin(), [scalar](double a)
        {
            return std::max(a, scalar);
        });
        return *this;
    }

    [[nodiscard]]
    dVectorN maximized(double scalar) const                                     ///< Returns vector containing the elements of the current vector with lower value limit, set by the scalar factor, applied to all elements
    {
        dVectorN result(storage.size());
        std::transform(cbegin(), cend(), result.begin(), [scalar](double a)
        {
            return std::max(a, scalar);
        });

        return result;
    }

    dVectorN& operator=(const dVectorN& rhs)                                    ///< Assigns the copy of the rhs to the current vector.
    {
        storage = rhs.storage;
        return *this;
    }

    [[nodiscard]]
    dVectorN operator*(const double m) const                                    ///< Returns the copy of the current vector with all its components multiplied by the value m.
    {
        dVectorN tmp(*this);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] *= m;
        }
        return tmp;
    }

    [[nodiscard]]
    dVectorN operator/(const double m) const                                    ///< Returns the copy of the current vector with all its components divided by the value m.
    {
        assert(m != 0 && "dVectorN::/ divided by zero");
        dVectorN tmp(*this);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] /= m;
        }
        return tmp;
    }

    dVectorN& operator*=(const double m)                                        ///< Multiplies all elements of the current vector by the value m.
    {
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] *= m;
        }
        return *this;
    }

    dVectorN& operator/=(const double m)                                        ///< Divides all elements of the current vector by the value m.
    {
        assert(m != 0 && "dVectorN::/= divided by zero");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= m;
        }
        return *this;
    }

    [[nodiscard]]
    dVectorN operator+(const dVectorN& rhs) const                               ///< Returns the sum of the current vector and the rhs.
    {
        assert(rhs.size() == storage.size() && "dVectorN::+ Sizes of the vectors are not equal");
        dVectorN tmp(*this);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] += rhs[i];
        }
        return tmp;
    }

    [[nodiscard]]
    dVectorN operator-(const dVectorN& rhs) const                               ///< Returns the difference of the current vector and the rhs.
    {
        assert(rhs.size() == storage.size() && "dVectorN::- Sizes of the vectors are not equal");
        dVectorN tmp(*this);
        for(size_t i = 0; i < storage.size(); i++)
        {
            tmp[i] -= rhs[i];
        }
        return tmp;
    }

//    [[nodiscard]]
//    double operator*(const dVectorN& rhs) const
//    {
//        assert(rhs.size() == storage.size() && "dVectorN::* Sizes of the vectors are not equal");
//
//        double result = 0.0;
//        for(size_t i = 0; i < storage.size(); i++)
//        {
//            result += storage[i] * rhs[i];
//        }
//        return result;
//    }

    [[nodiscard]]
    dVectorN operator*(const dVectorN& rhs) const                               ///< Returns a copy of the current vector component-wise multiplied with the rhs (non-standard operator)
    {
        dVectorN result(*this);
        result *= rhs;
        return result;
    }

    [[nodiscard]]
    dVectorN operator/(const dVectorN& rhs) const                               ///< Returns a copy of the current vector component-wise divided by the rhs (non-standard operator)
    {
        dVectorN result(*this);
        result /= rhs;
        return result;
    }

    dVectorN& operator+=( const dVectorN& rhs)                                  ///< Adds rhs vector to the current vector.
    {
        assert(rhs.size() == storage.size() && "dVectorN::+= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] += rhs[i];
        }
        return *this;
    }

    dVectorN& operator-=(const dVectorN& rhs)                                   ///< Subtracts rhs vector from the current vector.
    {
        assert(rhs.size() == storage.size() && "dVectorN::-= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] -= rhs[i];
        }
        return *this;
    }

    dVectorN& operator/=(const dVectorN& rhs)                                   ///< Returns component-wise division of the current vector and the rhs (non-standard operator).
    {
        assert(rhs.size() == storage.size() && "dVectorN::/= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] /= rhs[i];
        }
        return *this;
    }

    dVectorN& operator*=(const dVectorN& rhs)                                   ///< Returns component-wise product of the current vector and the rhs (non-standard operator).
    {
        assert(rhs.size() == storage.size() && "dVectorN::*= Sizes of the vectors are not equal");
        for(size_t i = 0; i < storage.size(); i++)
        {
            storage[i] *= rhs[i];
        }
        return *this;
    }

    void set_to_zero()                                                          ///< Sets all component of the current vector to the zero.
    {
        std::fill(storage.begin(), storage.end(), 0.0);
    }

    void set_to_value(const double value)                                       ///< Sets all component of the current vector to the given value.
    {
        std::fill(storage.begin(), storage.end(), value);
    }

    [[nodiscard]]
    size_t size(void) const                                                     ///< Returns size of the vector.
    {
        return storage.size();
    }

    std::string print(void) const                                               ///< Returns formatted string with vector content
    {
        std::stringstream out;
        out << "(";
        if(storage.size() != 0)
        {
            for (size_t n = 0; n < storage.size() - 1; n++)
            {
                out << std::setprecision(6) << storage[n] << ", ";
            }
            out << std::setprecision(6) << storage.back();
        }
        out << ")";
        return out.str();
    }

    void read_binary(std::istream& inp)                                         ///< Reads vector content from the binary input stream
    {
        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
        storage.resize(size);
        inp.read(reinterpret_cast<char*>(storage.data()), size * sizeof(double));
    }

    void write_binary(std::ostream& outp) const                                 ///< Writes vector content into the binary output stream
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

 protected:
 private:
    std::vector<double> storage;
};

} // namespace openphase
#endif
