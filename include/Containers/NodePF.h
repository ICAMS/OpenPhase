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

#ifndef NODEPF_H
#define NODEPF_H

#include <vector>
#include <utility>

#include "dVector3.h"

namespace openphase
{
/********************************* Declaration *******************************/

struct PhaseFieldEntry                                                          ///< Structure for storing individual phase-fields. Used in the NodePF class as a storage unit.
{
    size_t index;                                                               ///< Phase-field index.
    double value;                                                               ///< Phase-field value.
    double laplacian;                                                           ///< Phase-field Laplacian.
    dVector3 gradient;                                                          ///< Phase-field gradient.

    PhaseFieldEntry() :                                                         ///< Default constructor.
    index(0),
    value(0.0),
    laplacian(0.0),
    gradient(dVector3::ZeroVector())
    {

    }
    PhaseFieldEntry(const PhaseFieldEntry& entry) :                             ///< Copy constructor.
        index(entry.index),
        value(entry.value),
        laplacian(entry.laplacian),
        gradient(entry.gradient)
    {

    }
    PhaseFieldEntry& operator=(const PhaseFieldEntry& entry)                    ///< Copy operator.
    {
        index     = entry.index;
        value     = entry.value;
        laplacian = entry.laplacian;
        gradient  = entry.gradient;
        return *this;
    }
};

class NodePF                                                                    ///< Stores the phase-fields and their derivatives at a grid point. Provides access and manipulation methods for the phase-field entries.
{
 public:
    NodePF() :                                                                  ///< Constructor used to allocate space for at least 3 fields
        flag(0),
        Fields(0),
        tmpFields(0)
    {
        Fields.reserve(3);
        tmpFields.reserve(3);
    }

    NodePF(const NodePF& n) :                                                   ///< Copy constructor
        flag(n.flag),
        Fields(n.Fields),
        tmpFields(n.tmpFields)
    {

    }

    void clear()                                                                ///< Empties the field storage. Sets flag to 0.
    {
        Fields.clear();
        tmpFields.clear();
        flag  = 0;
    };

    int     finalize(void);                                                     ///< Adjusts phase-field values to (0, 1] interval.
    int     majority_index(void) const;                                         ///< Calculates and returns majority phase-field index

    NodePF  operator+(const NodePF& n) const;                                   ///< Plus operator. Takes as input another NodePF type entry.
    NodePF  operator-(const NodePF& n) const;                                   ///< Minus operator. Takes as input another NodePF type entry.
    NodePF  operator*(const double n) const;                                    ///< Multiplies all fields by a number.
    NodePF& operator=(const NodePF& n);                                         ///< Assignment operator.

    NodePF& operator+=(const NodePF& n);                                        ///< Plus-equal operator. Takes as input another NodePF type entry.
    NodePF& operator-=(const NodePF& n);                                        ///< Minus-equal operator. Takes as input another NodePF type entry.
    NodePF& operator*=(const double n);                                         ///< Multiply all fields by a number.

    bool    present(const size_t idx) const;                                    ///< Returns true if the phase-field with a given index is present in the node, false otherwise.

    void    set_value(const size_t n, const double value);                      ///< Sets value.
    void    add_value(const size_t n, const double value);                      ///< Increments value.
    double  get_value(const size_t n) const;                                    ///< Returns value.

    void    set_laplacian(const size_t n, const double laplacian);              ///< Sets Laplacian.
    void    add_laplacian(const size_t n, const double laplacian);              ///< Increments Laplacian.
    void    add_laplacian_tmp(const size_t n, const double laplacian);          ///< Increments temporary Laplacian.

    double  get_laplacian(const size_t n) const;                                ///< Returns Laplacian.

    void    set_gradient(const size_t n, const dVector3 gradient);              ///< Sets gradient.
    void    add_gradient(const size_t n, const dVector3 gradient);              ///< Increments gradient.
    void    add_gradient_tmp(const size_t n, const dVector3 gradient);          ///< Increments temporary gradient.

    dVector3 get_gradient(const size_t n) const;                                ///< Returns gradient.
    NodeA<dVector3> get_gradients() const;                                      ///< Returns gradients of all phase-fields.

    void    set_derivatives(const size_t n, const double laplacian, const dVector3 gradient);///< Sets gradient and Laplacian.
    void    add_derivatives(const size_t n, const double laplacian, const dVector3 gradient);///< Increments value and Laplacian.
    std::pair<double,dVector3> get_derivatives(const size_t n) const;           ///< Returns pair<Laplacian, gradient>.

    void    set_all(const size_t n, const double value, const double laplacian, const dVector3 gradient);///< Sets gradient and Laplacian.
    void    add_all(const size_t n, const double value, const double laplacian, const dVector3 gradient);///< Increments value and Laplacian.
    void    get_all(const size_t n, double& value, double& laplacian, dVector3& gradient) const;///< Increments value and Laplacian.

    void    add_values         (const NodePF& value);                           ///< Adds values of two nodes.
    void    add_laplacians     (const NodePF& value);                           ///< Adds Laplacians of two nodes.
    void    add_gradients      (const NodePF& value);                           ///< Adds gradients of two nodes.
    void    add_derivatives    (const NodePF& value);                           ///< Adds derivatives of two nodes.
    void    add_existing_values(const NodePF& value);                           ///< Adds only values existing in two nodes simultaneously.
    void    add_existing_all   (const NodePF& value);                           ///< Adds only phase-field entries existing in two nodes simultaneously.

    void    set_temporary(void);                                                ///< Sets temporary Fields copy
    void    copy_from_temporary(void);                                          ///< Restores Fields from tmpFields

    PhaseFieldEntry get_max(void) const;                                        ///< Returns FieldEntry with max value

    bool    interface(void) const;                                              ///< Returns true if flag == 2
    bool    interface_halo(void) const;                                         ///< Returns true if flag == 1
    bool    wide_interface(void) const;                                         ///< Returns true if flag != 0
    bool    bulk(void) const;                                                   ///< Returns true if flag != 2

    typedef std::vector<PhaseFieldEntry>::iterator iterator;                    ///< Iterator over storage vector.
    typedef std::vector<PhaseFieldEntry>::const_iterator citerator;             ///< Constant iterator over storage vector.

    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector.
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector.
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector.
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector.

    size_t    size() const {return Fields.size();};                             ///< Returns the size of storage.
    size_t    capacity() const {return Fields.capacity();};                     ///< Returns the capacity of storage.
    iterator  erase(iterator it) {return Fields.erase(it);};                    ///< Erase a single record pointed by iterator it.
    PhaseFieldEntry& front(void) {return Fields.front();};                      ///< Reference to the first FieldEntry.
    const PhaseFieldEntry& front(void) const {return Fields.front();};          ///< Constant reference to the first FieldEntry.

    void pack(std::vector<double>& buffer);
    void unpack(std::vector<double>& buffer, size_t& it);
    void read(std::istream& inp);                                               ///< Reads NodePF content from the input stream.
    void write(std::ostream& outp) const;                                       ///< Writes NodePF content to the output stream.

    int flag;                                                                   ///< Interface flag: 1 if node is near the interface, 2 if it is in the interface and 0 if it is in the bulk.
 protected:
 private:
    std::vector<PhaseFieldEntry> Fields;                                        ///< Stores phase fields and their derivatives.
    std::vector<PhaseFieldEntry> tmpFields;                                     ///< Stores temporary phase fields.
};

/******************************* Implementation ******************************/
inline void NodePF::set_temporary(void)
{
    tmpFields = Fields;
}

inline void NodePF::copy_from_temporary(void)
{
    Fields = tmpFields;
    tmpFields.clear();
}

inline int NodePF::majority_index(void) const
{
    if(Fields.size() == 1)
    {
        return Fields.front().index;
    }

    double max_value = 0.0;
    int max_index = 0;
    for(auto i = Fields.begin(); i != Fields.end(); ++i)
    {
        if(i->value > max_value)
        {
            max_value = i->value;
            max_index = i->index;
        }
    }
    return max_index;
}

inline int NodePF::finalize(void)
{
    double total = 0.0;
    for(auto i = Fields.begin(); i != Fields.end();)
    {
        bool erase = false;
        if (i->value >= 1.0 - DBL_EPSILON) i->value = 1.0;
        if (i->value <= 0.0 + DBL_EPSILON)
        {
            i->value = 0.0;
            erase = true;
        }
        total += i->value;
        if (erase)
        {
            i = Fields.erase(i);
        }
        else
        {
            ++i;
        }
    }
    if(Fields.size() > 1)
    {
        double norm = 1.0/total;
        for(auto i = Fields.begin(); i < Fields.end(); ++i)
        {
            i->value     *= norm;
            i->laplacian  = 0.0;
            i->gradient.set_to_zero();
        }
        flag = 2;
    }
    else
    {
        Fields.front().value      = 1.0;
        Fields.front().laplacian  = 0.0;
        Fields.front().gradient.set_to_zero();
        flag = 0;
    }
    return flag;
}

inline bool NodePF::present(const size_t index) const
{
    for(auto i = Fields.cbegin(); i != Fields.cend(); ++i)
    {
        if(i->index == index) return true;
    }
    return false;
}

inline void NodePF::set_value(const size_t n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value = value;
        return;
    }

    if(Fields.size())
    {
        flag = 2;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    NewEntry.laplacian  = 0.0;
    NewEntry.gradient.set_to_zero();
    Fields.push_back(NewEntry);
}

inline void NodePF::add_value(const size_t n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value += value;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    NewEntry.laplacian  = 0.0;
    NewEntry.gradient.set_to_zero();
    Fields.push_back(NewEntry);
}

inline double NodePF::get_value(const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value;
    }
    return 0.0;
}

inline void NodePF::set_laplacian(const size_t n, const double laplacian)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->laplacian = laplacian;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient.set_to_zero();
    Fields.push_back(NewEntry);
}

inline void NodePF::add_laplacian(const size_t n, const double laplacian)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->laplacian += laplacian;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient.set_to_zero();
    Fields.push_back(NewEntry);
}

inline void NodePF::add_laplacian_tmp(const size_t n, const double laplacian)
{
    for(auto i = tmpFields.begin(); i < tmpFields.end(); ++i)
    if(i->index == n)
    {
        i->laplacian += laplacian;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient.set_to_zero();
    tmpFields.push_back(NewEntry);
}

inline double NodePF::get_laplacian(const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->laplacian;
    }
    return 0.0;
}

inline void NodePF::set_gradient(const size_t n, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->gradient = gradient;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = 0.0;
    NewEntry.gradient = gradient;
    Fields.push_back(NewEntry);
}

inline void NodePF::add_gradient(const size_t n, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->gradient += gradient;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = 0.0;
    NewEntry.gradient = gradient;
    Fields.push_back(NewEntry);
}

inline void NodePF::add_gradient_tmp(const size_t n, const dVector3 gradient)
{
    for(auto i = tmpFields.begin(); i < tmpFields.end(); ++i)
    if(i->index == n)
    {
        i->gradient += gradient;
        return;
    }

    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = 0.0;
    NewEntry.gradient = gradient;
    tmpFields.push_back(NewEntry);
}

inline dVector3 NodePF::get_gradient(const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->gradient;
    }
    return dVector3::ZeroVector();
}

inline NodeA<dVector3> NodePF::get_gradients() const
{
    NodeA<dVector3> result;
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.set_value(i->index, i->gradient);
    }
    return result;
}

inline NodePF NodePF::operator+(const NodePF& n) const
{
    NodePF result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_all(i->index, i->value, i->laplacian, i->gradient);
    }
    result.flag = std::max(flag, n.flag);
    return result;
}

inline NodePF NodePF::operator-(const NodePF& n) const
{
    NodePF result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value      *= -1.0;
        i->laplacian  *= -1.0;
        i->gradient   *= -1.0;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_all(i->index, i->value, i->laplacian, i->gradient);
    }
    result.flag = std::max(flag, n.flag);
    return result;
}

inline NodePF NodePF::operator*(const double n) const
{
    NodePF result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value      *= n;
        i->laplacian  *= n;
        i->gradient   *= n;
    }
    return result;
}

inline NodePF& NodePF::operator+=(const NodePF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_all(i->index, i->value, i->laplacian, i->gradient);
    }
    flag = std::max(flag, n.flag);
    return *this;
}

inline NodePF& NodePF::operator-=(const NodePF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_all(i->index, -i->value, -i->laplacian, i->gradient*(-1.0));
    }
    flag = std::max(flag, n.flag);
    return *this;
}

inline NodePF& NodePF::operator*=(const double n)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    {
        i->value      *= n;
        i->laplacian  *= n;
        i->gradient   *= n;
    }
    return *this;
}

inline NodePF& NodePF::operator=(const NodePF& n)
{
    Fields = n.Fields;
    tmpFields = n.tmpFields;
    flag = n.flag;
    return *this;
}

inline void NodePF::set_all(const size_t n, const double value, const double laplacian, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value      = value;
        i->laplacian  = laplacian;
        i->gradient   = gradient;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient   = gradient;
    Fields.push_back(NewEntry);
}

inline void NodePF::add_all(const size_t n, const double value, const double laplacian, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value      += value;
        i->laplacian  += laplacian;
        i->gradient   += gradient;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = value;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient   = gradient;
    Fields.push_back(NewEntry);
}

inline void NodePF::get_all(const size_t n, double& value, double& laplacian, dVector3& gradient) const
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        value = i->value;
        laplacian = i->laplacian;
        gradient = i->gradient;
        return;
    }
}

inline void NodePF::set_derivatives(const size_t n, const double laplacian, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->laplacian  = laplacian;
        i->gradient   = gradient;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient   = gradient;
    Fields.push_back(NewEntry);
}

inline void NodePF::add_derivatives(const size_t n, const double laplacian, const dVector3 gradient)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->laplacian  += laplacian;
        i->gradient   += gradient;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    PhaseFieldEntry NewEntry;
    NewEntry.index      = n;
    NewEntry.value      = 0.0;
    NewEntry.laplacian  = laplacian;
    NewEntry.gradient   = gradient;
    Fields.push_back(NewEntry);
}

inline std::pair<double,dVector3> NodePF::get_derivatives(const size_t n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return std::pair<double,dVector3>(i->laplacian,i->gradient);
    }
    return std::pair<double,dVector3>(0.0, dVector3::ZeroVector());
}

inline void NodePF::add_values(const NodePF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_value(i->index, i->value);
    }
    flag = std::max(flag, n.flag);
}

inline void NodePF::add_laplacians(const NodePF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_laplacian(i->index, i->laplacian);
    }
    flag = std::max(flag, n.flag);
}

inline void NodePF::add_gradients(const NodePF& n)
{
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add_gradient(i->index, i->gradient);
    }
    flag = std::max(flag, n.flag);
}

inline void NodePF::add_existing_values(const NodePF& n)
{
    for (auto i = begin(); i < end(); )
    {
        double locVal = n.get_value(i->index);
        if(locVal != 0.0 and i->value != 0.0)
        {
            i->value += locVal;
            i++;
        }
        else
        {
            erase(i);
        }
    }
    flag = std::max(flag, n.flag);
}

inline void NodePF::add_existing_all(const NodePF& n)
{
    for (auto i = begin(); i < end(); )
    {
        double locValue = 0;
        double locLaplacian = 0;
        dVector3 locGradient;
        n.get_all(i->index, locValue, locLaplacian, locGradient);
        if(locValue != 0.0 and i->value != 0.0)
        {
            i->value += locValue;
            i->laplacian += locLaplacian;
            i->gradient += locGradient;
            i++;
        }
        else
        {
            erase(i);
        }
    }
    flag = std::max(flag, n.flag);
}

inline PhaseFieldEntry NodePF::get_max(void) const
{
    PhaseFieldEntry returnFieldEntry;

    returnFieldEntry.value = 0.0;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (fabs(i->value) > returnFieldEntry.value)
    {
        returnFieldEntry = *i;
    }
    return returnFieldEntry;
}

inline void NodePF::pack(std::vector<double>& buffer)
{
    buffer.push_back(Fields.size());
    for(auto it = Fields.begin(); it != Fields.end(); ++it)
    {
        buffer.push_back(it->index);
        buffer.push_back(it->value);
        buffer.push_back(it->laplacian);
        buffer.push_back(it->gradient[0]);
        buffer.push_back(it->gradient[1]);
        buffer.push_back(it->gradient[2]);
    }
    buffer.push_back(flag);
}

inline void NodePF::unpack(std::vector<double>& buffer, size_t& it)
{
    clear();
    size_t size = buffer[it]; ++it;
    Fields.resize(size);
    for(size_t i = 0; i < size; ++i)
    {
        Fields[i].index       = buffer[it]; ++it;
        Fields[i].value       = buffer[it]; ++it;
        Fields[i].laplacian   = buffer[it]; ++it;
        Fields[i].gradient[0] = buffer[it]; ++it;
        Fields[i].gradient[1] = buffer[it]; ++it;
        Fields[i].gradient[2] = buffer[it]; ++it;
    }
    flag  = buffer[it]; ++it;
}

inline void NodePF::read(std::istream& inp)
{
    size_t size = 0;
    flag = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    Fields.resize(size);
    if(size > 1) flag = 2;
    for(auto &Field : Fields)
    {
        inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&Field.value), sizeof(double));
    }
}

inline void NodePF::write(std::ostream& outp) const
{
    size_t size = Fields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));

    for(auto &Field : Fields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&Field.value), sizeof(double));
    }
}

inline bool NodePF::interface(void) const
{
    return flag == 2;
}

inline bool NodePF::interface_halo(void) const
{
    return flag == 1;
}

inline bool NodePF::wide_interface(void) const
{
    return flag != 0;
}

inline bool NodePF::bulk(void) const
{
    return flag != 2;
}
}//namespace openphase
#endif
