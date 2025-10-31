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

#ifndef NODEIP_H
#define NODEIP_H

#include <vector>
#include <iostream>

namespace openphase
{
/**********************************************************/

struct InterfacePropertiesFieldEntry                                            ///< Structure for storing the field entry with two indices. Used in the NodeIP class as a storage unit.
{
    size_t indexA;                                                              ///< First index.
    size_t indexB;                                                              ///< Second index.
    double energy;                                                              ///< Energy value.
    double stiffness;                                                           ///< Stiffness value.
    double mobility;                                                            ///< Mobility value.

    InterfacePropertiesFieldEntry() :                                           ///< Default constructor.
        indexA(0),
        indexB(0),
        energy(0.0),
        stiffness(0.0),
        mobility(0.0)
    {

    }
    InterfacePropertiesFieldEntry(const InterfacePropertiesFieldEntry& entry) : ///< Copy constructor.
        indexA(entry.indexA),
        indexB(entry.indexB),
        energy(entry.energy),
        stiffness(entry.stiffness),
        mobility(entry.mobility)
    {

    }
    InterfacePropertiesFieldEntry& operator=(const InterfacePropertiesFieldEntry& entry)///< Copy operator.
    {
        indexA    = entry.indexA;
        indexB    = entry.indexB;
        energy    = entry.energy;
        stiffness = entry.stiffness;
        mobility  = entry.mobility;
        return *this;
    }
};

/**********************************************************/
class NodeIP                                                                    ///< Stores the interface properties with two indices and two values at a grid point. Provides access and manipulation methods for the stored field entries.
{
 public:
    NodeIP()                                                                    ///< Constructor used to allocate space for at least 3 fields.
    {
        Fields.reserve(3);
    }
    NodeIP(const NodeIP& n) :                                                   ///< Copy constructor.
        Fields(n.Fields)
    {

    }
    NodeIP  operator+(const NodeIP& n) const;                                   ///< Plus operator. Takes as input another NodeIP type entry.
    NodeIP  operator-(const NodeIP& n) const;                                   ///< Minus operator. Takes as input another NodeIP type entry.
    NodeIP  operator*(const double  n) const;                                   ///< Multiplies all fields by a number.
    NodeIP& operator=(const NodeIP& n);                                         ///< Assignment operator

    NodeIP& operator+=(const NodeIP& n);                                        ///< Plus-equal operator. Takes as input another NodeIP type entry.
    NodeIP& operator-=(const NodeIP& n);                                        ///< Minus-equal operator. Takes as input another NodeIP type entry.
    NodeIP& operator*=(const double  n);                                        ///< Multiplies all fields by a number.

    void   set_energy(const size_t n, const size_t m, const double en_value);   ///< Sets interface energy value in symmetric case: f(n,m) = f(m,n).
    void   set_stiffness(const size_t n, const size_t m, const double st_value);///< Sets interface stiffness value in symmetric case: f(n,m) = f(m,n).
    void   set_mobility(const size_t n, const size_t m, const double mob_value);///< Sets interface mobility value in symmetric case: f(n,m) = f(m,n).

    double get_energy(const size_t n, const size_t m) const;                    ///< Returns interface energy value in symmetric case: f(n,m) = f(m,n).
    double get_stiffness(const size_t n, const size_t m) const;                 ///< Returns interface stiffness value in symmetric case: f(n,m) = f(m,n).
    double get_mobility(const size_t n, const size_t m) const;                  ///< Returns interface mobility value in symmetric case: f(n,m) = f(m,n).

    void   add_energy(const size_t n, const size_t m, const double en_value);   ///< Increments interface energy value in symmetric case: f(n,m) = f(m,n).
    void   add_stiffness(const size_t n, const size_t m, const double st_value);///< Increments interface stiffness value in symmetric case: f(n,m) = f(m,n).
    void   add_mobility(const size_t n, const size_t m, const double mob_value);///< Increments interface mobility value in symmetric case: f(n,m) = f(m,n).

    void   set_energy_and_mobility(const size_t n, const size_t m,
                                   const double en_value,
                                   const double mob_value);                     ///< Sets en_value and mob_value in symmetric case: f(n,m) = f(m,n).

    void   set_stiffness_and_mobility(const size_t n, const size_t m,
                                      const double st_value,
                                      const double mob_value);                  ///< Sets st_value and mob_value in symmetric case: f(n,m) = f(m,n).

    void   set_energy_and_stiffness(const size_t n, const size_t m,
                                    const double en_value,
                                    const double st_value);                     ///< Sets en_value and st_value in symmetric case: f(n,m) = f(m,n).

    void   set_energy_stiffness_and_mobility(const size_t n, const size_t m,
                                             const double st_value,
                                             const double en_value,
                                             const double mob_value);           ///< Sets en_value and mob_value in symmetric case: f(n,m) = f(m,n).

    void   add_energy_and_mobility(const size_t n, const size_t m,
                                   const double en_value,
                                   const double mob_value);                     ///< Increments en_value and mob_value in symmetric case: f(n,m) = f(m,n).

    void   add_stiffness_and_mobility(const size_t n, const size_t m,
                                      const double stn_value,
                                      const double mob_value);                  ///< Increments st_value and mob_value in symmetric case: f(n,m) = f(m,n).

    void   add_energy_and_stiffness(const size_t n, const size_t m,
                                    const double st_value,
                                    const double en_value);                     ///< Increments en_value and st_value in symmetric case: f(n,m) = f(m,n).

    void   add_energy_stiffness_and_mobility(const size_t n, const size_t m,
                                             const double st_value,
                                             const double en_value,
                                             const double mob_value);           ///< Increments en_value and mob_value in symmetric case: f(n,m) = f(m,n).

    std::pair<double,double> get_energy_and_mobility(const size_t n,
                                                     const size_t m) const;     ///< Returns pair<en_value,mob_value> in symmetric case: f(n,m) = f(m,n).
    std::pair<double,double> get_stiffness_and_mobility(const size_t n,
                                                        const size_t m) const;  ///< Returns pair<st_value,mob_value> in symmetric case: f(n,m) = f(m,n).
    std::pair<double,double> get_energy_and_stiffness(const size_t n,
                                                      const size_t m) const;    ///< Returns pair<en_value,st_value> in symmetric case: f(n,m) = f(m,n).

    InterfacePropertiesFieldEntry get_entry(const size_t n, const size_t m) const;   ///< Returns InterfacePropertiesFieldEntry in symmetric case: f(n,m) = f(m,n).

    void   add_energies(const NodeIP& en_values);                               ///< Add en_value of two nodes in symmetric case: f(n,m) = f(m,n).
    void   add_stiffnesses(const NodeIP& st_values);                            ///< Add st_value of two nodes in symmetric case: f(n,m) = f(m,n).
    void   add_mobilities(const NodeIP& mob_values);                            ///< Add mob_value of two nodes in symmetric case: f(n,m) = f(m,n).

    void   add_energies_and_mobilities(const NodeIP& values);                   ///< Add en_value and mob_value of two nodes in symmetric case: f(n,m) = f(m,n).
    void   add_energies_and_stiffnesses(const NodeIP& values);                  ///< Add en_value and st_value of two nodes in symmetric case: f(n,m) = f(m,n).
    void   add_energies_stiffnesses_and_mobilities(const NodeIP& values);       ///< Add en_value, st_value and mob_value of two nodes in symmetric case: f(n,m) = f(m,n).

    void   add_energies_exist(const NodeIP& en_values);                         ///< Add only en_value existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).
    void   add_stiffnesses_exist(const NodeIP& en_values);                      ///< Add only st_value existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).
    void   add_mobilities_exist(const NodeIP& mob_values);                      ///< Add only mob_value existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).

    void    clear() {Fields.clear();};                                          ///< Empties the field storage. Sets flag to 0.
    size_t  size() const {return Fields.size();};                               ///< Returns the size of storage.
    typedef std::vector<InterfacePropertiesFieldEntry>::iterator iterator;      ///< Iterator over storage vector
    typedef std::vector<InterfacePropertiesFieldEntry>::const_iterator citerator;///< Constant iterator over storage vector
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector
    iterator erase(iterator it) {return Fields.erase(it);};                     ///< Erase a single record pointed by iterator it
    InterfacePropertiesFieldEntry& front(void) {return Fields.front();};        ///< Reference to the first FieldEntry.
    const InterfacePropertiesFieldEntry& front(void) const {return Fields.front();};///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void read(std::istream& inp);
    void write(std::ostream& outp) const;

 protected:
 private:
    std::vector<InterfacePropertiesFieldEntry> Fields;                          ///< Fields storage vector.
};

/***************************************************************/

inline NodeIP NodeIP::operator+(const NodeIP& n) const
{
    NodeIP result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_energy_stiffness_and_mobility(i->indexA, i->indexB, i->energy, i->stiffness, i->mobility);
    }
    return result;
}

inline NodeIP NodeIP::operator-(const NodeIP& n) const
{
    NodeIP result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->energy *= -1.0;
        i->stiffness *= -1.0;
        i->mobility *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_energy_stiffness_and_mobility(i->indexA, i->indexB, i->energy, i->stiffness, i->mobility);
    }
    return result;
}

inline NodeIP& NodeIP::operator+=(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy_stiffness_and_mobility(i->indexA, i->indexB, i->energy, i->stiffness, i->mobility);
    }
    return *this;
}

inline NodeIP& NodeIP::operator-=(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy_stiffness_and_mobility(i->indexA, i->indexB, -i->energy, -i->stiffness, -i->mobility);
    }
    return *this;
}

inline NodeIP NodeIP::operator*(const double n) const
{
    NodeIP result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->energy *= n;
        i->stiffness *= n;
        i->mobility *= n;
    }
    return result;
}

inline NodeIP& NodeIP::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->energy *= n;
        i->stiffness *= n;
        i->mobility *= n;
    }
    return *this;
}

inline NodeIP& NodeIP::operator=(const NodeIP& n)
{
    Fields = n.Fields;
    return *this;
}

inline void NodeIP::set_energy(const size_t n, const size_t m,
                               const double en_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy = en_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::set_stiffness(const size_t n, const size_t m,
                                  const double st_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->stiffness = st_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.stiffness = st_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::set_mobility(const size_t n, const size_t m,
                                 const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->mobility = mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_energy(const size_t n, const size_t m,
                               const double en_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy += en_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_stiffness(const size_t n, const size_t m,
                                  const double st_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->stiffness += st_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.stiffness = st_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_mobility(const size_t n, const size_t m,
                                 const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->mobility += mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline double NodeIP::get_energy(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->energy;
    }
    return 0.0;
}

inline double NodeIP::get_stiffness(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->stiffness;
    }
    return 0.0;
}

inline double NodeIP::get_mobility(const size_t n, const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->mobility;
    }
    return 0.0;
}

inline void NodeIP::set_energy_and_mobility(const size_t n, const size_t m,
                                            const double en_value,
                                            const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy = en_value;
        i->mobility = mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::set_stiffness_and_mobility(const size_t n, const size_t m,
                                               const double st_value,
                                               const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->stiffness = st_value;
        i->mobility = mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.stiffness = st_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::set_energy_and_stiffness(const size_t n, const size_t m,
                                             const double en_value,
                                             const double st_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy = en_value;
        i->stiffness = st_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.stiffness = st_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::set_energy_stiffness_and_mobility(const size_t n, const size_t m,
                                                      const double en_value,
                                                      const double st_value,
                                                      const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy = en_value;
        i->stiffness = st_value;
        i->mobility = mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.stiffness = st_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_energy_and_mobility(const size_t n, const size_t m,
                                            const double en_value,
                                            const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy += en_value;
        i->mobility += mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_stiffness_and_mobility(const size_t n, const size_t m,
                                               const double st_value,
                                               const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->stiffness += st_value;
        i->mobility += mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.stiffness = st_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_energy_and_stiffness(const size_t n, const size_t m,
                                             const double en_value,
                                             const double st_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy += en_value;
        i->stiffness += st_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.stiffness = st_value;

    Fields.push_back(NewEntry);
}

inline void NodeIP::add_energy_stiffness_and_mobility(const size_t n, const size_t m,
                                                      const double en_value,
                                                      const double st_value,
                                                      const double mob_value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->energy += en_value;
        i->stiffness += st_value;
        i->mobility += mob_value;
        return;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.energy = en_value;
    NewEntry.stiffness = st_value;
    NewEntry.mobility = mob_value;

    Fields.push_back(NewEntry);
}

inline std::pair<double,double> NodeIP::get_energy_and_mobility(const size_t n,
                                                                const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return std::pair<double,double>(i->energy, i->mobility);
    }
    return std::pair<double,double>(0.0, 0.0);
}

inline std::pair<double,double> NodeIP::get_stiffness_and_mobility(const size_t n,
                                                                   const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return std::pair<double,double>(i->stiffness, i->mobility);
    }
    return std::pair<double,double>(0.0, 0.0);
}

inline std::pair<double,double> NodeIP::get_energy_and_stiffness(const size_t n,
                                                                 const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return std::pair<double,double>(i->energy, i->stiffness);
    }
    return std::pair<double,double>(0.0, 0.0);
}

inline InterfacePropertiesFieldEntry NodeIP::get_entry(const size_t n,
                                                              const size_t m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return *i;
    }

    InterfacePropertiesFieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;

    return NewEntry;
}

inline void NodeIP::add_energies(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy(i->indexA, i->indexB, i->energy);
    }
}

inline void NodeIP::add_stiffnesses(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_stiffness(i->indexA, i->indexB, i->stiffness);
    }
}

inline void NodeIP::add_mobilities(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_mobility(i->indexA, i->indexB, i->mobility);
    }
}

inline void NodeIP::add_energies_and_mobilities(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy_and_mobility(i->indexA, i->indexB, i->energy, i->mobility);
    }
}

inline void NodeIP::add_energies_and_stiffnesses(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy_and_stiffness(i->indexA, i->indexB, i->energy, i->stiffness);
    }
}
inline void NodeIP::add_energies_stiffnesses_and_mobilities(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_energy_stiffness_and_mobility(i->indexA, i->indexB, i->energy, i->stiffness, i->mobility);
    }
}

inline void NodeIP::add_energies_exist(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_energy(i->indexA, i->indexB);
        if(locVal != 0.0 and i->energy != 0.0)
        {
            set_energy(i->indexA, i->indexB, i->energy + locVal);
        }
    }
}

inline void NodeIP::add_stiffnesses_exist(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_stiffness(i->indexA, i->indexB);
        if(locVal != 0.0 and i->stiffness != 0.0)
        {
            set_stiffness(i->indexA, i->indexB, i->stiffness + locVal);
        }
    }
}

inline void NodeIP::add_mobilities_exist(const NodeIP& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        double locVal = get_mobility(i->indexA, i->indexB);
        if(locVal != 0.0 and i->mobility != 0.0)
        {
            set_mobility(i->indexA, i->indexB, i->mobility + locVal);
        }
    }
}

inline void NodeIP::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end();++it)
    {
        buffer.push_back(it->indexA);
        buffer.push_back(it->indexB);
        buffer.push_back(it->energy);
        buffer.push_back(it->stiffness);
        buffer.push_back(it->mobility);
    }
}

inline void NodeIP::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].indexA = buffer[it]; ++it;
        Fields[i].indexB = buffer[it]; ++it;
        Fields[i].energy = buffer[it]; ++it;
        Fields[i].stiffness = buffer[it]; ++it;
        Fields[i].mobility = buffer[it]; ++it;
    }
}

inline void NodeIP::read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);

    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.indexA), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.indexB), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.energy), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.stiffness), sizeof(double));
        inp.read(reinterpret_cast<char*>(&Field.mobility), sizeof(double));
    }
}

inline void NodeIP::write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.indexA), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.indexB), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.energy), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.stiffness), sizeof(double));
        outp.write(reinterpret_cast<const char*>(&Field.mobility), sizeof(double));
    }
}

}//namespace openphase
#endif
