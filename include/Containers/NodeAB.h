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

 *   File created :   2009
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich;
 *                         Reza Darvishi Kamachali; Dmitry Medvedev
 *
 */

#ifndef NODEAB_H
#define NODEAB_H

#include <vector>
#include <iostream>

namespace openphase
{
/**********************************************************/
template<class T1, class T2>
struct DoubleIndexFieldEntry                                                    ///< Structure for storing the field entry with two indices. Used in the NodeAB class as a storage unit.
{
    size_t indexA;                                                              ///< First index.
    size_t indexB;                                                              ///< Second index.
    T1     value1;                                                              ///< First value.
    T2     value2;                                                              ///< Second value.

    DoubleIndexFieldEntry() :                                                   ///< Default constructor.
        indexA(0),
        indexB(0),
        value1(),
        value2()
    {

    }
    DoubleIndexFieldEntry(const DoubleIndexFieldEntry<T1,T2>& entry) :              ///< Copy constructor.
        indexA(entry.indexA),
        indexB(entry.indexB),
        value1(entry.value1),
        value2(entry.value2)
    {

    }
    DoubleIndexFieldEntry<T1,T2>& operator=(const DoubleIndexFieldEntry<T1,T2>& entry)  ///< Copy operator.
    {
        indexA = entry.indexA;
        indexB = entry.indexB;
        value1 = entry.value1;
        value2 = entry.value2;
        return *this;
    }
};

/**********************************************************/
template<class T1, class T2>
class NodeAB                                                                    ///< Stores the fields with two indices and two values at a grid point. Provides access and manipulation methods for the stored field entries.
{
 public:
    NodeAB<T1,T2>()                                                             ///< Constructor used to allocate space for at least 3 fields.
    {
        Fields.reserve(3);
    }
    NodeAB<T1,T2>(const NodeAB<T1,T2>& n) :                                     ///< Copy constructor.
        Fields(n.Fields)
    {

    }
    NodeAB<T1,T2>  operator+(const NodeAB<T1,T2>& n) const;                     ///< Plus operator. Takes as input another Node type entry.
    NodeAB<T1,T2>  operator-(const NodeAB<T1,T2>& n) const;                     ///< Minus operator. Takes as input another Node type entry.
    NodeAB<T1,T2>  operator*(const double  n) const;                            ///< Multiplies all fields by a number.
    NodeAB<T1,T2>& operator=(const NodeAB<T1,T2>& n);                           ///< Assignment operator

    NodeAB<T1,T2>& operator+=(const NodeAB<T1,T2>& n);                          ///< Plus-equal operator. Takes as input another Node type entry.
    NodeAB<T1,T2>& operator-=(const NodeAB<T1,T2>& n);                          ///< Minus-equal operator. Takes as input another Node type entry.
    NodeAB<T1,T2>& operator*=(const double  n);                                 ///< Multiplies all fields by a number.

    void set1(const size_t n, const size_t m, const T1 value);                  ///< Sets value1 following strict index order.
    void set2(const size_t n, const size_t m, const T2 value);                  ///< Sets value2 following strict index order.
    T1    get1(const size_t n, const size_t m) const;                           ///< Returns value1 following strict index order.
    T2    get2(const size_t n, const size_t m) const;                           ///< Returns value2 following strict index order.
    void add1(const size_t n, const size_t m, const T1 value);                  ///< Increments value1 following strict index order.
    void add2(const size_t n, const size_t m, const T2 value);                  ///< Increments value2 following strict index order.

    void set_sym1(const size_t n, const size_t m, const T1 value);              ///< Sets value1 in symmetric case: f(n,m) = f(m,n).
    void set_sym2(const size_t n, const size_t m, const T2 value);              ///< Sets value2 in symmetric case: f(n,m) = f(m,n).
    T1    get_sym1(const size_t n, const size_t m) const;                       ///< Returns value1 in symmetric case: f(n,m) = f(m,n).
    T2    get_sym2(const size_t n, const size_t m) const;                       ///< Returns value2 in symmetric case: f(n,m) = f(m,n).
    void add_sym1(const size_t n, const size_t m, const T1 value);              ///< Increments value1 in symmetric case: f(n,m) = f(m,n).
    void add_sym2(const size_t n, const size_t m, const T2 value);              ///< Increments value2 in symmetric case: f(n,m) = f(m,n).

    void set_asym1(const size_t n, const size_t m, const T1 value);             ///< Sets value1 in antisymmetric case: f(n,m) = -f(m,n).
    void set_asym2(const size_t n, const size_t m, const T2 value);             ///< Sets value2 in antisymmetric case: f(n,m) = -f(m,n).
    T1    get_asym1(const size_t n, const size_t m) const;                      ///< Returns value1 in antisymmetric case: f(n,m) = -f(m,n).
    T2    get_asym2(const size_t n, const size_t m) const;                      ///< Returns value2 in antisymmetric case: f(n,m) = -f(m,n).
    void add_asym1(const size_t n, const size_t m, const T1 value);             ///< Increments value1 in antisymmetric case: f(n,m) = -f(m,n).
    void add_asym2(const size_t n, const size_t m, const T2 value);             ///< Increments value2 in antisymmetric case: f(n,m) = -f(m,n).

    void set_pair(const size_t n, const size_t m, const T1 value1,
                                                  const T2 value2);             ///< Sets value1 and value2 following strict index order.
    void add_pair(const size_t n, const size_t m, const T1 value1,
                                                  const T2 value2);             ///< Increments value1 and value2 following strict index order.
    std::pair<T1,T2> get_pair(const size_t n, const size_t m) const;            ///< Returns pair<value1,value2> following strict index order.

//    T1    get_asym_sum(const size_t n, const size_t m) const;                    ///< Returns sum of value1 and value2 in antisymmetric case: f(n,m) = -f(m,n).
//    T2    get_sym_sum(const size_t n, const size_t m) const;                     ///< Returns sum of value1 and value2 in symmetric case: f(n,m) = f(m,n).

    void set_sym_pair(const size_t n, const size_t m, const T1 value1,
                                                      const T2 value2);         ///< Sets value1 and value2 in symmetric case: f(n,m) = f(m,n).
    void add_sym_pair(const size_t n, const size_t m, const T1 value1,
                                                      const T2 value2);         ///< Increments value1 and value2 in symmetric case: f(n,m) = f(m,n).
    std::pair<T1,T2> get_sym_pair(const size_t n, const size_t m) const;        ///< Returns pair<value1,value2>  in symmetric case: f(n,m) = f(m,n).

    void set_asym_pair(const size_t n, const size_t m, const T1 value1,
                                                          const T2 value2);     ///< Sets value1 and value2 in antisymmetric case: f(n,m) = -f(m,n).
    void add_asym_pair(const size_t n, const size_t m, const T1 value1,
                                                          const T2 value2);     ///< Increments value1 and value2 in antisymmetric case: f(n,m) = -f(m,n).
    std::pair<T1,T2> get_asym_pair(const size_t n, const size_t m) const;       ///< Returns pair<value1,value2>  in antisymmetric case: f(n,m) = -f(m,n).

    void add1(const NodeAB<T1,T2>& value);                                      ///< Add value1 of two nodes following strict index order.
    void add2(const NodeAB<T1,T2>& value);                                      ///< Add value2 of two nodes following strict index order.

    void add_sym1(const NodeAB<T1,T2>& value);                                  ///< Add value1 of two nodes in symmetric case: f(n,m) = f(m,n).
    void add_sym2(const NodeAB<T1,T2>& value);                                  ///< Add value2 of two nodes in symmetric case: f(n,m) = f(m,n).

    void add_asym1(const NodeAB<T1,T2>& value);                                 ///< Add value1 of two nodes in antisymmetric case: f(n,m) = -f(m,n).
    void add_asym2(const NodeAB<T1,T2>& value);                                 ///< Add value2 of two nodes in antisymmetric case: f(n,m) = -f(m,n).

    void add_sym_pairs(const NodeAB<T1,T2>& value);                             ///< Add value1 and value2 of two nodes in symmetric case: f(n,m) = f(m,n).
    void add_asym_pairs(const NodeAB<T1,T2>& value);                            ///< Add value1 and value2 of two nodes in antisymmetric case: f(n,m) = -f(m,n).

    void add_sym1_exist(const NodeAB<T1,T2>& value);                            ///< Add only value1 existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).
    void add_sym2_exist(const NodeAB<T1,T2>& value);                            ///< Add only value2 existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).

    void add_asym1_exist(const NodeAB<T1,T2>& value);                           ///< Add only value1 existing in two nodes simultaneously in antisymmetric case: f(n,m) = -f(m,n).
    void add_asym2_exist(const NodeAB<T1,T2>& value);                           ///< Add only value2 existing in two nodes simultaneously in antisymmetric case: f(n,m) = -f(m,n).

    bool    present(const size_t idxA, const size_t idxB) const;                ///< Returns true if the field with a given indices is present in the node, false otherwise.

    void      clear() {Fields.clear();};                                        ///< Empties the field storage. Sets flag to 0.
    size_t    size() const {return Fields.size();};                             ///< Returns the size of storage.
    typedef typename std::vector<DoubleIndexFieldEntry<T1,T2>>::iterator iterator;  ///< Iterator over storage vector
    typedef typename std::vector<DoubleIndexFieldEntry<T1,T2>>::const_iterator citerator; ///< Constant iterator over storage vector
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector
    size_t    capacity() const {return Fields.capacity();};                     ///< Returns the size of storage.
    iterator  erase(iterator it) {return Fields.erase(it);};                    ///< Erase a single record pointed by iterator it
    DoubleIndexFieldEntry<T1,T2>& front(void) {return Fields.front();};         ///< Reference to the first FieldEntry.
    const DoubleIndexFieldEntry<T1,T2>& front(void) const {return Fields.front();}; ///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void read(std::istream& inp);
    void write(std::ostream& outp) const;

 protected:
 private:
    std::vector<DoubleIndexFieldEntry<T1,T2>> Fields;                           ///< Fields storage vector.
};

/********************************* Implementation *****************************/

template<class T1, class T2>
inline bool NodeAB<T1,T2>::present(const size_t indexA, const size_t indexB) const
{
    for(auto i = Fields.cbegin(); i != Fields.cend(); ++i)
    {
        if(i->indexA == indexA and i->indexB == indexB) return true;
    }
    return false;
}
template<class T1, class T2>
inline NodeAB<T1,T2> NodeAB<T1,T2>::operator+(const NodeAB<T1,T2>& n) const
{
    NodeAB<T1,T2> result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return result;
}
template<class T1, class T2>
inline NodeAB<T1,T2> NodeAB<T1,T2>::operator-(const NodeAB<T1,T2>& n) const
{
    NodeAB<T1,T2> result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value1 *= -1.0;
        i->value2 *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return result;
}
template<class T1, class T2>
inline NodeAB<T1,T2>& NodeAB<T1,T2>::operator+=(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
    return *this;
}
template<class T1, class T2>
inline NodeAB<T1,T2>& NodeAB<T1,T2>::operator-=(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_pair(i->indexA, i->indexB, i->value1*(-1.0), i->value2*(-1.0));
    }
    return *this;
}
template<class T1, class T2>
inline NodeAB<T1,T2> NodeAB<T1,T2>::operator*(const double n) const
{
    NodeAB<T1,T2> result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value1 *= n;
        i->value2 *= n;
    }
    return result;
}
template<class T1, class T2>
inline NodeAB<T1,T2>& NodeAB<T1,T2>::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->value1 *= n;
        i->value2 *= n;
    }
    return *this;
}
template<class T1, class T2>
inline NodeAB<T1,T2>& NodeAB<T1,T2>::operator=(const NodeAB<T1,T2>& n)
{
    Fields = n.Fields;
    return *this;
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set1(const size_t n, const size_t m, const T1 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 = value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = T2();

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set2(const size_t n, const size_t m, const T2 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value2 = value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = T1();
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline T1 NodeAB<T1,T2>::get1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return i->value1;
    }
    return T1();
}
template<class T1, class T2>
inline T2 NodeAB<T1,T2>::get2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return i->value2;
    }
    return T2();
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_sym1(const size_t n, const size_t m, const T1 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 = value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = T2();

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_sym2(const size_t n, const size_t m, const T2 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value2 = value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = T1();
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym1(const size_t n, const size_t m, const T1 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 += value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = T2();

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym2(const size_t n, const size_t m, const T2 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value2 += value;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = T1();
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline T1 NodeAB<T1,T2>::get_sym1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value1;
    }
    return T1();
}
template<class T1, class T2>
inline T2 NodeAB<T1,T2>::get_sym2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value2;
    }
    return T2();
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_asym1(const size_t n, const size_t m, const T1 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value1 = value;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value1 = value*(-1.0);
            return;
        }
    }
    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = T2();

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_asym2(const size_t n, const size_t m, const T2 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value2 = value;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value2 = value*(-1.0);
            return;
        }
    }
    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = T1();
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym1(const size_t n, const size_t m, const T1 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value1 += value;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value1 += value*(-1.0);
            return;
        }
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value;
    NewEntry.value2 = T2();

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym2(const size_t n, const size_t m, const T2 value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value2 += value;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value2 += value*(-1.0);
            return;
        }
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = T1();
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline T1 NodeAB<T1,T2>::get_asym1(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            return i->value1;
        }
        if(i->indexA == m and i->indexB == n)
        {
            return i->value1*(-1.0);
        }
    }
    return T1();
}
template<class T1, class T2>
inline T2 NodeAB<T1,T2>::get_asym2(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            return i->value2;
        }
        if(i->indexA == m and i->indexB == n)
        {
            return i->value2*(-1.0);
        }
    }
    return T2();
}
//template<class T1, class T2>
//inline T NodeAB<T>::get_asym_sum(const size_t n, const size_t m) const
//{
//    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
//    {
//        if (i->indexA == n and i->indexB == m)
//        {
//            return i->value1 + i->value2;
//        }
//        if(i->indexA == m and i->indexB == n)
//        {
//            return -(i->value1 + i->value2);
//        }
//    }
//    return T();
//}
//template<class T1, class T2>
//inline T NodeAB<T1,T2>::get_sym_sum(const size_t n, const size_t m) const
//{
//    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
//    {
//        if ((i->indexA == n and i->indexB == m) or
//            (i->indexA == m and i->indexB == n))
//        {
//            return i->value1 + i->value2;
//        }
//    }
//    return T();
//}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 = value1;
        i->value2 = value2;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        i->value1 += value1;
        i->value2 += value2;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline std::pair<T1,T2> NodeAB<T1,T2>::get_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m))
    {
        return std::pair<double,double>(i->value1, i->value2);
    }
    return std::pair<T1,T2>();
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_sym_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 = value1;
        i->value2 = value2;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value1 += value1;
        i->value2 += value2;
        return;
    }

    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline std::pair<T1,T2> NodeAB<T1,T2>::get_sym_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return std::pair<T1,T2>(i->value1, i->value2);
    }
    return std::pair<T1,T2>();
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::set_asym_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value1 = value1;
            i->value2 = value2;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value1 = value1*(-1.0);
            i->value2 = value2*(-1.0);
            return;
        }
    }
    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym_pair(const size_t n, const size_t m, const T1 value1, const T2 value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            i->value1 += value1;
            i->value2 += value2;
            return;
        }
        if(i->indexA == m and i->indexB == n)
        {
            i->value1 += value1*(-1.0);
            i->value2 += value2*(-1.0);
            return;
        }
    }
    DoubleIndexFieldEntry<T1,T2> NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value1 = value1;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}
template<class T1, class T2>
inline std::pair<T1,T2> NodeAB<T1,T2>::get_asym_pair(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->indexA == n and i->indexB == m)
        {
            return std::pair<T1,T2>(i->value1, i->value2);
        }
        if(i->indexA == m and i->indexB == n)
        {
            return std::pair<T1,T2>(i->value1*(-1.0), i->value2*(-1.0));
        }
    }
    return std::pair<T1,T2>();
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym1(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym1(i->indexA, i->indexB, i->value1);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym2(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym2(i->indexA, i->indexB, i->value2);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym1(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym1(i->indexA, i->indexB, i->value1);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym2(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym2(i->indexA, i->indexB, i->value2);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym_pairs(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_sym_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym_pairs(const NodeAB<T1,T2>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_asym_pair(i->indexA, i->indexB, i->value1, i->value2);
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym1_exist(const NodeAB<T1,T2>& n)
{
    for (auto i = begin(); i < end(); )
    {
        if(n.present(i->indexA,i->indexB) or
           n.present(i->indexB,i->indexA))
        {
            i->value1 += n.get_sym1(i->indexA, i->indexB);
            i++;
        }
        else
        {
            erase(i);
        }
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym1_exist(const NodeAB<T1,T2>& n)
{
    for (auto i = begin(); i < end(); )
    {
        if(n.present(i->indexA,i->indexB))
        {
            i->value1 += n.get_asym1(i->indexA, i->indexB);
            i++;
        }
        else
        {
            erase(i);
        }
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_sym2_exist(const NodeAB<T1,T2>& n)
{
    for (auto i = begin(); i < end(); )
    {
        if(n.present(i->indexA,i->indexB) or
           n.present(i->indexB,i->indexA))
        {
            i->value2 += n.get_sym2(i->indexA, i->indexB);
            i++;
        }
        else
        {
            erase(i);
        }
    }
}
template<class T1, class T2>
inline void NodeAB<T1,T2>::add_asym2_exist(const NodeAB<T1,T2>& n)
{
    for (auto i = begin(); i < end(); )
    {
        if(n.present(i->indexA,i->indexB))
        {
            i->value2 += n.get_asym2(i->indexA, i->indexB);
            i++;
        }
        else
        {
            erase(i);
        }
    }
}
template<>
inline void NodeAB<double,double>::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end();++it)
    {
        buffer.push_back(it->indexA);
        buffer.push_back(it->indexB);
        buffer.push_back(it->value1);
        buffer.push_back(it->value2);
    }
}
template<>
inline void NodeAB<double,double>::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].indexA = buffer[it]; ++it;
        Fields[i].indexB = buffer[it]; ++it;
        Fields[i].value1 = buffer[it]; ++it;
        Fields[i].value2 = buffer[it]; ++it;
    }
}
template<>
inline void NodeAB<double,double>::read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);

    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.indexA), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.indexB), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.value1), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.value2), sizeof(double));
    }
}
template<>
inline void NodeAB<double,double>::write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.indexA), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.indexB), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.value1), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.value2), sizeof(double));
    }
}

//template<class T>
//inline void NodeAB<T>::pack(std::vector<double>& buffer)
//{
//    buffer.push_back(Fields.size());
//    for(auto it = Fields.begin(); it != Fields.end();++it)
//    {
//        buffer.push_back(it->indexA);
//        buffer.push_back(it->indexB);
//        buffer.push_back(it->value1);
//        buffer.push_back(it->value2);
//    }
//}
//template<class T>
//inline void NodeAB<T>::unpack(std::vector<double>& buffer, size_t& it)
//{
//    clear();
//    size_t size = buffer[it]; ++it;
//    Fields.resize(size);
//    for(size_t i = 0; i < size; ++i)
//    {
//        Fields[i].indexA = buffer[it]; ++it;
//        Fields[i].indexB = buffer[it]; ++it;
//        Fields[i].value1 = buffer[it]; ++it;
//        Fields[i].value2 = buffer[it]; ++it;
//    }
//}
//template<class T>
//inline void NodeAB<T>::read(std::istream& inp)
//{
//    size_t size = 0;
//    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
//    Fields.resize(size);
//
//    for(auto &Field : Fields)
//    {
//        inp.read(reinterpret_cast<char*>(&Field.indexA), sizeof(size_t));
//        inp.read(reinterpret_cast<char*>(&Field.indexB), sizeof(size_t));
//        inp.read(reinterpret_cast<char*>(&Field.value1), sizeof(double));
//        inp.read(reinterpret_cast<char*>(&Field.value2), sizeof(double));
//    }
//}
//template<class T>
//inline void NodeAB<T>::write(std::ostream& outp) const
//{
//    size_t size = Fields.size();
//    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
//
//    for(auto &Field : Fields)
//    {
//        outp.write(reinterpret_cast<const char*>(&Field.indexA), sizeof(size_t));
//        outp.write(reinterpret_cast<const char*>(&Field.indexB), sizeof(size_t));
//        outp.write(reinterpret_cast<const char*>(&Field.value1), sizeof(double));
//        outp.write(reinterpret_cast<const char*>(&Field.value2), sizeof(double));
//    }
//}

}//namespace openphase
#endif
