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

#ifndef NODEA_H
#define NODEA_H

#include <vector>
#include <iostream>
#include <cmath>

namespace openphase
{
/********************************* Declaration *******************************/
template<class T = double >
struct SingleIndexFieldEntry                                                    ///< Structure for storing the single index entries. Used in the NodeA class as a storage unit.
{
    size_t index;                                                               ///< Field index.
    T      value;                                                               ///< Stored value.

    SingleIndexFieldEntry<T>() :                                                ///< Default constructor.
        index(0)
    {

    }
    SingleIndexFieldEntry<T>(const SingleIndexFieldEntry<T>& entry) :           ///< Copy constructor.
        index(entry.index),
        value(entry.value)
    {

    }
    SingleIndexFieldEntry<T>& operator=(const SingleIndexFieldEntry<T>& entry)  ///< Copy operator.
    {
        index = entry.index;
        value = entry.value;
        return *this;
    }
};
template<class T>
class NodeA                                                                     ///< Stores the single-valued fields at a grid point. Provides access and manipulation methods for the entries.
{
 public:
    NodeA<T>()                                                                  ///< Constructor used to allocate space for at least 3 fields
    {
        Fields.reserve(3);
    }

    NodeA<T>(const NodeA& n) :                                                  ///< Copy constructor
        Fields(n.Fields)
    {

    }

    void clear()                                                                ///< Empties the field storage and resets flag and index
    {
        Fields.clear();
    };

    NodeA<T>   operator+(const NodeA<T>& n) const;                              ///< Plus operator. Takes as input another NodeA type entry.
    NodeA<T>   operator-(const NodeA<T>& n) const;                              ///< Minus operator. Takes as input another NodeA type entry.
    NodeA<T>   operator*(const double n) const;                                 ///< Multiplies all fields by a number.
    NodeA<T>&  operator=(const NodeA<T>& n);                                    ///< Assignment operator.

    NodeA<T>&  operator+=(const NodeA<T>& n);                                   ///< Plus-equal operator. Takes as input another NodeA type entry.
    NodeA<T>&  operator-=(const NodeA<T>& n);                                   ///< Minus-equal operator. Takes as input another NodeA type entry.
    NodeA<T>&  operator*=(const double n);                                      ///< Multiply all fields by a number.

    bool    present(const size_t idx) const;                                    ///< Returns true if the field with a given index is present in the node, false otherwise.

    void    set_value(const size_t n, const T& value);                          ///< Sets value.
    void    add_value(const size_t n, const T& value);                          ///< Increments value.
    T       get_value(const size_t n) const;                                    ///< Returns value.

    void    add_values         (const NodeA<T>& value);                         ///< Adds values of two nodes.
    void    add_existing_values(const NodeA<T>& value);                         ///< Adds only values existing in two nodes simultaneously.

    typedef typename std::vector<SingleIndexFieldEntry<T>>::iterator iterator;  ///< Iterator over storage vector.
    typedef typename std::vector<SingleIndexFieldEntry<T>>::const_iterator citerator; ///< Constant iterator over storage vector.
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector.
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector.
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector.
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector.
    size_t    size() const {return Fields.size();};                             ///< Returns the size of storage.
    size_t    capacity() const {return Fields.capacity();};                     ///< Returns the size of storage.
    iterator  erase(iterator it) {return Fields.erase(it);};                    ///< Erase a single record pointed by iterator it.
    SingleIndexFieldEntry<T>& front(void) {return Fields.front();};             ///< Reference to the first FieldEntry.
    const SingleIndexFieldEntry<T>& front(void) const {return Fields.front();}; ///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void read(std::istream& inp);
    void write(std::ostream& outp) const;

 protected:
 private:
    std::vector<SingleIndexFieldEntry<T>> Fields;                               ///< Storage vector.
};

/******************************* Implementation ******************************/
template<class T>
inline bool NodeA<T>::present(const size_t index) const
{
    for(auto i = Fields.cbegin(); i != Fields.cend(); ++i)
    {
        if(i->index == index) return true;
    }
    return false;
}
template<class T>
inline void NodeA<T>::set_value(const size_t n, const T& value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value = value;
        return;
    }
    SingleIndexFieldEntry<T> NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    Fields.push_back(NewEntry);
}
template<class T>
inline void NodeA<T>::add_value(const size_t n, const T& value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value += value;
        return;
    }

    SingleIndexFieldEntry<T> NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    Fields.push_back(NewEntry);
}
template<class T>
inline T NodeA<T>::get_value(const size_t n) const
{
    T result = T();
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n)
        {
            result = i->value;
            break;
        }
    }
    return result;
}
template<class T>
inline NodeA<T> NodeA<T>::operator+(const NodeA<T>& n) const
{
    NodeA result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_value(i->index, i->value);
    }
    return result;
}
template<class T>
inline NodeA<T> NodeA<T>::operator-(const NodeA<T>& n) const
{
    NodeA result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_value(i->index, i->value);
    }
    return result;
}
template<class T>
inline NodeA<T> NodeA<T>::operator*(const double n) const
{
    NodeA result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= n;
    }
    return result;
}
template<class T>
inline NodeA<T>& NodeA<T>::operator+=(const NodeA<T>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_value(i->index, i->value);
    }
    return *this;
}

template<class T>
inline NodeA<T>& NodeA<T>::operator-=(const NodeA<T>& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_value(i->index, -i->value);
    }
    return *this;
}

template<class T>
inline NodeA<T>& NodeA<T>::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->value *= n;
    }
    return *this;
}

template<class T>
inline NodeA<T>& NodeA<T>::operator=(const NodeA<T>& n)
{
    Fields = n.Fields;
    return *this;
}

template<class T>
inline void NodeA<T>::add_existing_values(const NodeA<T>& n)
{
    for (auto i = begin(); i < end(); )
    {
        if(n.present(i->index))
        {
            i->value += n.get_value(i->index);
            i++;
        }
        else
        {
            i = erase(i);
        }
    }
}

template<>
inline void NodeA<double>::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end();++it)
    {
        buffer.push_back(it->index);
        buffer.push_back(it->value);
    }
}
template<>
inline void NodeA<double>::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].index = buffer[it]; ++it;
        Fields[i].value = buffer[it]; ++it;
    }
}
template<>
inline void NodeA<double>::read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);

    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.value), sizeof(double));
    }
}
template<>
inline void NodeA<double>::write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.value), sizeof(double));
    }
}

}//namespace openphase
#endif
