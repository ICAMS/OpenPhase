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
 *   Main contributors :   Oleg Shchyglo, Philipp Engels
 *
 */

#ifndef NODEVN_H
#define NODEVN_H

#include "dVector.h"

namespace openphase
{

template<size_t N>
struct VectorVnEntry                                                            /// Individual entries. Used in the NodeVn class as a storage unit.
{
    double value[N];
    size_t index;
};

/***************************************************************/
template<size_t N>
class NodeVn                                                                    /// Basic class to store the vector valued quantities for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    NodeVn(void){};                                                             /// Empty constructor
    NodeVn(const NodeVn<N>& rhs);                                               /// Copy constructor
    void   set_comp(const size_t idx, const size_t ss, const double value);     /// Set component ss for phase field index idx to value.
    void   set     (const size_t idx, const dVector<N> value);                  /// Set dVector<N> for phase field index idx.
    void   set     (const size_t idx, const double value);                      /// Set all entries to single value for phase field index idx.
    void   set_to_zero(const size_t idx);                                       /// Set all entries to zero for phase field index idx.
    double get_comp(const size_t idx,  const size_t ss) const;                  /// Return value of component ss for phase field index idx.
    dVector<N> get (const size_t idx) const;                                    /// Return vector corresponding to phase field index idx.
    void   add    (const size_t n, const double value[N]);                      /// Increment value using two indices.
    void   change_index  (const size_t n, const size_t m);

    NodeVn&  operator=(const NodeVn<N>& n);
    NodeVn   operator+(const NodeVn<N>& n) const;                                  /// Plus operator. Takes as input another Node type entry.
    NodeVn   operator-(const NodeVn<N>& n) const;                                  /// Minus operator. Takes as input another Node type entry.
    NodeVn   operator*(const double n) const;                                   /// Multiply all fields by a number.
    NodeVn   operator+(const double n) const;                                   /// Multiply all fields by a number.
    NodeVn   operator-(const double n) const;                                   /// Multiply all fields by a number.

    NodeVn   add     (const NodeVn<N>& value) const;                               /// Add two nodes.

    void   clear(){VectorVnFields.clear();};                                    /// Empties the vector fields.
    size_t  size() const {return VectorVnFields.size();};                        /// Returns the size of vector fields.
    typedef typename std::vector<VectorVnEntry<N>>::iterator iterator;          /// Iterator over the vector fields
    typedef typename std::vector<VectorVnEntry<N>>::const_iterator citerator;   /// Constant iterator over the vector fields
    iterator begin() {return VectorVnFields.begin();};                          /// Iterator to the begin of vector fields
    iterator end()   {return VectorVnFields.end();};                            /// Iterator to the end of vector fields
    citerator cbegin() const {return VectorVnFields.begin();};                  /// Constant iterator to the begin of vector fields
    citerator cend()   const {return VectorVnFields.end();};                    /// Constant iterator to the end of vector fields
    void pack(std::vector<double>& buffer)
    {
        buffer.push_back(VectorVnFields.size());
        for (size_t i = 0; i < VectorVnFields.size(); ++i)
        {
            buffer.push_back(VectorVnFields[i].index);
            for (size_t n = 0; n < N; ++n)
            {
                buffer.push_back(VectorVnFields[i].value[n]);
            }
        }
    }
    void unpack(std::vector<double>& buffer, size_t& it)
    {
        VectorVnFields.resize(buffer[it]); ++it;
        for (size_t i = 0; i < VectorVnFields.size(); ++i)
        {
            VectorVnFields[i].index = buffer[it]; ++it;
            for (size_t n = 0; n < N; ++n)
            {
                VectorVnFields[i].value[n] = buffer[it]; ++it;
            }
        }
    }
 protected:
 private:
    std::vector<VectorVnEntry<N>> VectorVnFields;                               /// List of nonvanishing vector fields.
};

template<size_t N>
inline NodeVn<N>::NodeVn(const NodeVn<N>& n)
{
    VectorVnFields = n.VectorVnFields;
}

template<size_t N>
inline void NodeVn<N>::set_comp(const size_t idx, const size_t ss, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        i->value[ss] = value;
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(size_t ii = 0; ii < N; ii++) NewEntry.value[ii] = 0.0;
    NewEntry.value[ss] = value;

    VectorVnFields.push_back(NewEntry);
}

template<size_t N>
inline void NodeVn<N>::set(const size_t idx, const dVector<N> value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        for (size_t ss = 0; ss < N; ss++)
        {
            i->value[ss] = value[ss];
        }
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(size_t ii = 0; ii < N; ii++)
    {
        NewEntry.value[ii] = value[ii];
    }
    VectorVnFields.push_back(NewEntry);
}

template<size_t N>
inline void NodeVn<N>::set(const size_t idx, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->index == idx)
        {
            for (size_t ii = 0; ii < N; ++ii)
            {
                i->value[ii] = value;
            }
            return;
        }
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(size_t ii = 0; ii < N; ii++) NewEntry.value[ii] = value;

    VectorVnFields.push_back(NewEntry);
}

template<size_t N>
inline void NodeVn<N>::set_to_zero(const size_t idx)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->index == idx)
        {
            for (size_t ii = 0; ii < N; ++ii)
            {
                i->value[ii] = 0.0;
            }
            return;
        }
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(size_t ii = 0; ii < N; ii++) NewEntry.value[ii] = 0.0;

    VectorVnFields.push_back(NewEntry);
}

template<size_t N>
inline void NodeVn<N>::change_index(const size_t n, const size_t m)
{
    // Set values and index of m to index n
    // Warning: Values are averaged arithmetically!
    if (n != m)
    {
        double returndV[N];
        for (size_t ss = 0; ss < N; ss++) returndV[ss] = 0.0;

        bool found = false;
        for (auto i = cbegin(); i < cend(); ++i)
        {
            if (i->index == m)
            {
                for (size_t ss = 0; ss < N; ss++)
                {
                    returndV[ss] = (i->value)[ss];
                }
                found = true;
                break;
            }
        }
        if (found == true)
        {
            for (auto i = begin(); i < end(); ++i)
            if (i->index == m)
            {
                for (size_t ss = 0; ss < N; ss++)
                {
                    (i->value)[ss] = 0.5*((i->value)[ss]+returndV[ss]);
                }
                i->index = n;
                break;
            }
        }
    }
}

template<size_t N>
inline double NodeVn<N>::get_comp(const size_t idx, const size_t ss) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        return i->value[ss];
    }
    return 0.0;
}

template<size_t N>
inline dVector<N> NodeVn<N>::get(const size_t idx) const
{
    dVector<N> returndV; returndV.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    {
        if (i->index == idx)
        {
            for (size_t ss = 0; ss < N; ss++)
            {
                returndV[ss] = (i->value)[ss];
            }
            break;
        }
    }
    return returndV;
}

template<size_t N>
inline void NodeVn<N>::add(const size_t n, const double summand[N])
{
    for (auto i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        for(size_t ss = 0; ss < N; ++ss)
        {
            i->value[ss] += summand[ss];
        }
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = n;
    for (size_t ss = 0; ss < N; ss ++) NewEntry.value[ss] = summand[ss];

    VectorVnFields.push_back(NewEntry);
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::add(const NodeVn<N>& n) const
{
    NodeVn<N> result = n;

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<size_t N>
inline NodeVn<N>& NodeVn<N>::operator=(const NodeVn<N>& n)
{
    VectorVnFields = n.VectorVnFields;
    return *this;
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::operator+(const NodeVn<N>& n) const
{
    NodeVn<N> result = n;

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::operator-(const NodeVn<N>& n) const
{
    NodeVn<N> result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    for(size_t ss = 0; ss < N; ++ss)
    {
        i->value[ss] = -i->value[ss];
    }

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::operator*(const double n) const
{
    NodeVn<N> result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    for(size_t ss = 0; ss < N; ++ss)
    {
        i->value[ss] *= n;
    }

    return result;
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::operator+(const double n) const
{
    NodeVn<N> result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    for(size_t ss = 0; ss < N; ++ss)
    {
        i->value[ss] += n;
    }

    return result;
}

template<size_t N>
inline NodeVn<N> NodeVn<N>::operator-(const double n) const
{
    NodeVn<N> result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    for(size_t ss = 0; ss < N; ++ss)
    {
        i->value[ss] -= n;
    }

    return result;
}

} //namespace openphase
#endif
