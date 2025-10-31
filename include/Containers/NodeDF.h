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

#ifndef NODEDF_H
#define NODEDF_H

#include <vector>
#include <iostream>

namespace openphase
{
/**********************************************************/

struct DrivingForceEntry                                                        ///< Structure for storing the driving force entry. Used in the NodeDF class as a storage unit.
{
    size_t indexA;                                                              ///< First index.
    size_t indexB;                                                              ///< Second index.
    double raw;                                                                 ///< Raw driving force value.
    double tmp;                                                                 ///< Intermediary driving force value.
    double average;                                                             ///< Average driving force value.
    double weight;                                                              ///< Driving force weight.

    DrivingForceEntry() :                                                       ///< Default constructor.
        indexA(0),
        indexB(0),
        raw(0.0),
        tmp(0.0),
        average(0.0),
        weight(0.0)
    {

    }
    DrivingForceEntry(const DrivingForceEntry& entry) :                         ///< Copy constructor.
        indexA(entry.indexA),
        indexB(entry.indexB),
        raw(entry.raw),
        tmp(entry.tmp),
        average(entry.average),
        weight(entry.weight)
    {

    }
    DrivingForceEntry& operator=(const DrivingForceEntry& entry)                ///< Assignment operator.
    {
        indexA  = entry.indexA;
        indexB  = entry.indexB;
        raw     = entry.raw;
        tmp     = entry.tmp;
        average = entry.average;
        weight  = entry.weight;
        return *this;
    }
};

/**********************************************************/
class NodeDF                                                                    ///< Stores driving forces at a grid point. Provides access and manipulation methods for the stored entries.
{
 public:
    NodeDF() :                                                                  ///< Constructor used to allocate space for at least 3 fields.
        Fields(0)
    {
        Fields.reserve(3);
    }
    NodeDF(const NodeDF& n) :                                                   ///< Copy constructor.
        Fields(n.Fields)
    {

    }
    NodeDF  operator+(const NodeDF& n) const;                                   ///< Plus operator. Takes as input another NodeDF type entry. Adds all entries.
    NodeDF  operator-(const NodeDF& n) const;                                   ///< Minus operator. Takes as input another NodeDF type entry. Subtracts all entries.
    NodeDF  operator*(const double  n) const;                                   ///< Multiplies all entries by a number.
    NodeDF& operator=(const NodeDF& n);                                         ///< Assignment operator

    NodeDF& operator+=(const NodeDF& n);                                        ///< Plus-equal operator. Takes as input another NodeDF type entry. Adds all entries.
    NodeDF& operator-=(const NodeDF& n);                                        ///< Minus-equal operator. Takes as input another NodeDF type entry. Subtracts all entries.
    NodeDF& operator*=(const double  n);                                        ///< Multiplies all entries by a number.

    void    set_raw    (const size_t n, const size_t m, const double value);    ///< Sets raw in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n)
    void    set_tmp    (const size_t n, const size_t m, const double value);    ///< Sets tmp in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n)
    void    set_average(const size_t n, const size_t m, const double value);    ///< Sets average in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    void    set_weight (const size_t n, const size_t m, const double value);    ///< Sets weight in symmetric case (does not change the sign when index order changes): f(n,m) = f(m,n).
    void    set_all    (const size_t n, const size_t m, const double raw_value,
                                                        const double tmp_value,
                                                        const double avg_value,
                                                        const double wgt_value);///< Sets all values

    double  get_raw    (const size_t n, const size_t m) const;                  ///< Returns raw in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    double  get_tmp    (const size_t n, const size_t m) const;                  ///< Returns tmp in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    double  get_average(const size_t n, const size_t m) const;                  ///< Returns average in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    double  get_weight (const size_t n, const size_t m) const;                  ///< Returns weight in symmetric case (does not change the sign when index order changes): f(n,m) = f(m,n).

    void    add_raw    (const size_t n, const size_t m, const double value);    ///< Increments raw in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    void    add_tmp    (const size_t n, const size_t m, const double value);    ///< Increments tmp in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    void    add_average(const size_t n, const size_t m, const double value);    ///< Increments average in antisymmetric case (changes sign when index order changes): f(n,m) = -f(m,n).
    void    add_weight (const size_t n, const size_t m, const double value);    ///< Increments weight in symmetric case (does not change the sign when index order changes): f(n,m) = f(m,n).
    void    add_all    (const size_t n, const size_t m, const double raw_value,
                                                        const double tmp_value,
                                                        const double avg_value,
                                                        const double wgt_value);///< Adds all values

    void    clear() {Fields.clear();};                                          ///< Empties the field storage.
    size_t  size() const {return Fields.size();};                               ///< Returns the size of storage.
    typedef std::vector<DrivingForceEntry>::iterator iterator;                  ///< Iterator over storage vector
    typedef std::vector<DrivingForceEntry>::const_iterator citerator;           ///< Constant iterator over storage vector
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector
    iterator erase(iterator it) {return Fields.erase(it);};                     ///< Erase a single record pointed by iterator it
    DrivingForceEntry& front(void) {return Fields.front();};                    ///< Reference to the first FieldEntry.
    const DrivingForceEntry& front(void) const {return Fields.front();};        ///< Constant reference to the first DrivingForceEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void Read(std::istream& inp);
    void Write(std::ostream& outp) const;

 protected:
 private:
    std::vector<DrivingForceEntry> Fields;                                      ///< DrivingForceEntry storage vector.
};

/***************************************************************/

inline NodeDF NodeDF::operator+(const NodeDF& n) const
{
    NodeDF result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_all(i->indexA, i->indexB, i->raw, i->tmp, i->average, i->weight);
    }
    return result;
}

inline NodeDF NodeDF::operator-(const NodeDF& n) const
{
    NodeDF result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->raw     *= -1.0;
        i->tmp     *= -1.0;
        i->average *= -1.0;
        i->weight  *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_all(i->indexA, i->indexB, i->raw, i->tmp, i->average, i->weight);
    }
    return result;
}

inline NodeDF& NodeDF::operator+=(const NodeDF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_all(i->indexA, i->indexB, i->raw, i->tmp, i->average, i->weight);
    }
    return *this;
}

inline NodeDF& NodeDF::operator-=(const NodeDF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_all(i->indexA, i->indexB, -i->raw, -i->tmp, -i->average, -i->weight);
    }
    return *this;
}

inline NodeDF NodeDF::operator*(const double n) const
{
    NodeDF result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->raw     *= n;
        i->tmp     *= n;
        i->average *= n;
        i->weight  *= n;
    }
    return result;
}

inline NodeDF& NodeDF::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->raw     *= n;
        i->tmp     *= n;
        i->average *= n;
        i->weight  *= n;
    }
    return *this;
}

inline NodeDF& NodeDF::operator=(const NodeDF& n)
{
    Fields = n.Fields;
    return *this;
}

inline void NodeDF::set_raw(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->raw = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->raw = -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = value;
    NewEntry.tmp     = 0.0;
    NewEntry.average = 0.0;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::set_tmp(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->tmp = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->tmp = -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = value;
    NewEntry.average = 0.0;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::set_average(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->average = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->average = -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = 0.0;
    NewEntry.average = value;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::set_weight(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n))
    {
        i->weight = value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = 0.0;
    NewEntry.average = 0.0;
    NewEntry.weight  = value;

    Fields.push_back(NewEntry);
}

inline void NodeDF::set_all(const size_t n, const size_t m, const double raw_value,
                                                            const double tmp_value,
                                                            const double avg_value,
                                                            const double wgt_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->raw     = raw_value;
        i->tmp     = tmp_value;
        i->average = avg_value;
        i->weight  = wgt_value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->raw     = -raw_value;
        i->tmp     = -tmp_value;
        i->average = -avg_value;
        i->weight  =  wgt_value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = raw_value;
    NewEntry.tmp     = tmp_value;
    NewEntry.average = avg_value;
    NewEntry.weight  = wgt_value;

    Fields.push_back(NewEntry);
}

inline double NodeDF::get_raw(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->raw;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->raw;
    }
    return 0.0;
}

inline double NodeDF::get_tmp(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->tmp;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->tmp;
    }
    return 0.0;
}

inline double NodeDF::get_average(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->average;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->average;
    }
    return 0.0;
}

inline double NodeDF::get_weight(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n))
    {
        return i->weight;
    }
    return 0.0;
}

inline void NodeDF::add_raw(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->raw += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->raw += -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = value;
    NewEntry.tmp     = 0.0;
    NewEntry.average = 0.0;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::add_tmp(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->tmp += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->tmp += -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = value;
    NewEntry.average = 0.0;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::add_average(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->average += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->average += -value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = 0.0;
    NewEntry.average = value;
    NewEntry.weight  = 0.0;

    Fields.push_back(NewEntry);
}

inline void NodeDF::add_weight(const size_t n, const size_t m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n))
    {
        i->weight += value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = 0.0;
    NewEntry.tmp     = 0.0;
    NewEntry.average = 0.0;
    NewEntry.weight  = value;

    Fields.push_back(NewEntry);
}

inline void NodeDF::add_all(const size_t n, const size_t m, const double raw_value,
                                                            const double tmp_value,
                                                            const double avg_value,
                                                            const double wgt_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->raw     += raw_value;
        i->tmp     += tmp_value;
        i->average += avg_value;
        i->weight  += wgt_value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->raw     += -raw_value;
        i->tmp     += -tmp_value;
        i->average += -avg_value;
        i->weight  +=  wgt_value;
        return;
    }

    DrivingForceEntry NewEntry;
    NewEntry.indexA  = n;
    NewEntry.indexB  = m;
    NewEntry.raw     = raw_value;
    NewEntry.tmp     = tmp_value;
    NewEntry.average = avg_value;
    NewEntry.weight  = wgt_value;

    Fields.push_back(NewEntry);
}

inline void NodeDF::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end();++it)
    {
        buffer.push_back(it->indexA );
        buffer.push_back(it->indexB );
        buffer.push_back(it->raw    );
        buffer.push_back(it->tmp    );
        buffer.push_back(it->average);
        buffer.push_back(it->weight );
    }
}

inline void NodeDF::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].indexA  = buffer[it]; ++it;
        Fields[i].indexB  = buffer[it]; ++it;
        Fields[i].raw     = buffer[it]; ++it;
        Fields[i].tmp     = buffer[it]; ++it;
        Fields[i].average = buffer[it]; ++it;
        Fields[i].weight  = buffer[it]; ++it;
    }
}

inline void NodeDF::Read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);

    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.indexA ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.indexB ), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.raw    ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.tmp    ), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.average), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.weight ), sizeof(double));
    }
}

inline void NodeDF::Write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.indexA ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.indexB ), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.raw    ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.tmp    ), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.average), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.weight ), sizeof(double));
    }
}

}//namespace openphase
#endif
