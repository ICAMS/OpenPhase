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

#ifndef NodeVN_H
#define NodeVN_H

#include "dVectorN.h"

namespace openphase
{

class VectorEntry                                                              	///< Individual vector entry. Used in the NodeVectorN class as a storage unit.
{
public:
    VectorEntry():
        Vector(0),
        index(0)
    {
    };
    VectorEntry(size_t size, size_t idx):
        Vector(size),
        index(idx)
    {
    };
    VectorEntry(const VectorEntry& other):
        Vector(other.Vector),
        index(other.index)
    {
    };
    VectorEntry& operator=(const VectorEntry& other)
    {
        Vector = other.Vector;
        index = other.index;
        return *this;
    };
    dVectorN Vector;
    size_t index;
};

/******************************************************************************/

class NodeVN                                                                    ///< Basic class to store the vector valued quantities for each phase field.
{
public:

    NodeVN(): SIZE_X(0){}                                                       ///< Default constructor.
    NodeVN(const size_t size): SIZE_X(size){}                                   ///< Constructor. Sets the stored dVectorN size

    NodeVN(const NodeVN& other):
        VectorFields(other.VectorFields),
        SIZE_X(other.SIZE_X)                                                    ///< Copy constructor.
    {
    }
    void set_vector_size(const size_t size)                                     ///< Assigns the stored dVectorN size
    {
        SIZE_X = size;
    }
    size_t get_vector_size() const                                              ///< Returns stored dVectorN size
    {
        return SIZE_X;
    }

    bool exist(const size_t idx) const;                                         ///< True if entry with phase field "idx" is present in the container, false otherwise.

    void set_to_value(const size_t idx, const int value);                       ///< Sets all entries for phase field idx to the given value.
    void set_to_zero(const size_t idx);                                         ///< Sets all entries for phase field idx to zero.
    void set(const size_t idx, const size_t ii, const double value);            ///< Sets vector component ii for phase field idx to the given value.
    void set(const size_t idx, dVectorN value);                                 ///< Sets dVectorN for phase field idx.
    double   get(const size_t idx, const size_t ii) const;                      ///< Returns value of component ii for phase field idx.
    dVectorN get(const size_t idx) const;                                       ///< Returns dVectorN for phase field idx.
    void         add(const size_t idx, dVectorN value);                         ///< Increments dVectorN by value for phase field idx.
    NodeVN  add(const NodeVN& n) const;                                         ///< Adds two containers.
    void         add(const size_t idx, const size_t ii, const double value);    ///< Adds value to dVectorN component ii for phase field idx.
    void    multiply(const size_t idx, dVectorN value);

    NodeVN& operator=(const NodeVN& n);                                         ///< Assignment operator. Assigns the content of n to the current container.
    NodeVN& operator+=(const NodeVN& n);                                        ///< Adds the content of n to the current container
    NodeVN& operator*=(const NodeVN& n);                                        ///< Multiplies the content of the current container by the content of n
    NodeVN& operator-=(const NodeVN& n);                                        ///< Subtracts the content of n from the current container
    NodeVN& operator*=(const double n);                                         ///< Multiplies all entries in the current container by the factor n
    NodeVN  operator+(const  NodeVN& n) const;                                  ///< Plus operator. Returns the sum of the two containers.
    NodeVN  operator-(const  NodeVN& n) const;                                  ///< Minus operator. Returns the difference of the two containers.
    NodeVN  operator*(const double n) const;                                    ///< Returns copy of the current container with all entries multiplied by the factor n.

    void   clear(){VectorFields.clear();};                                      ///< Empties the vector fields.
    size_t  size() const {return VectorFields.size();};                         ///< Returns the number of the stored entries.

    typedef typename std::vector<VectorEntry>::iterator iterator;               ///< Iterator over the stored entries
    typedef typename std::vector<VectorEntry>::const_iterator citerator;        ///< Constant iterator over the entries storage
    iterator  begin()        {return VectorFields.begin();};                    ///< Iterator to the begin of the entries storage
    iterator  end()          {return VectorFields.end();};                      ///< Iterator to the end of the entries storage
    citerator cbegin() const {return VectorFields.cbegin();};                   ///< Constant iterator to the begin of the entries storage
    citerator cend()   const {return VectorFields.cend();};                     ///< Constant iterator to the end of the entries storage
    iterator  erase(iterator it) {return VectorFields.erase(it);};              ///< Erases single entry pointed to by the iterator it.

    std::string print();                                                        ///< Returns formatted string with container content

    void pack(std::vector<double>& buffer);                                     ///< Writes NodeVectorN content into the buffer (used for MPI mode communication)
    void unpack(std::vector<double>& buffer, size_t& it);                       ///< Read NodeVectorN content from the buffer (used for MPI mode communication)

    void read(std::istream& inp);                                               ///< Reads NodeVectorN content from the input stream.
    void write(std::ostream& outp) const;                                       ///< Writes NodeVectorN content to the output stream.

protected:
private:
    std::vector<VectorEntry> VectorFields;
    size_t SIZE_X;
};

inline bool NodeVN::exist(const size_t idx) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if(i->index == idx)
    {
        return true;
    }
    return false;
}

inline void NodeVN::set_to_value(const size_t idx, const int value)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        i->Vector.set_to_value(value);
        return;
    }

    VectorEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector.set_to_value(value);
    VectorFields.push_back(NewEntry);
}

inline void NodeVN::set_to_zero(const size_t idx)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        i->Vector.set_to_zero();
        return;
    }

    VectorEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector.set_to_zero();
    VectorFields.push_back(NewEntry);
}

inline void NodeVN::set(const size_t idx, const size_t ii, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::set() ii is greater than dVectorN.size()");
        i->Vector[ii] = value;
        return;
    }

    VectorEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector.set_to_zero();
    NewEntry.Vector[ii] = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeVN::set(const size_t idx, dVectorN Vector)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::set() Vector size is not equal");
        i->Vector = Vector;
        return;
    }

    VectorEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline double NodeVN::get(const size_t idx, const size_t ii) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::get() ii is greater than dVectorN.size()");
        return i->Vector[ii];
    }
    return 0.0;
}

inline dVectorN NodeVN::get(const size_t idx) const
{
    dVectorN returndV(SIZE_X); returndV.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        returndV = i->Vector;
        break;
    }
    return returndV;
}

inline void NodeVN::add(const size_t idx, dVectorN Vector)
{
    for(auto i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::add Vector size is not equal");
        i->Vector += Vector;
        return;
    }

    VectorEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline void NodeVN::multiply(const size_t idx, dVectorN Vector)
{
    for(auto i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::multiply() Vector size is not equal");
        i->Vector *= Vector;
        return;
    }

    VectorEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline NodeVN NodeVN::add(const NodeVN& n) const
{
    NodeVN result = n;
    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->Vector);
    }
    return result;
}

inline void NodeVN::add(const size_t idx, const size_t ii, const double value)
{
    for(iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::add() idx is greater than dVectorN.size()");
        i->Vector[ii] += value;
        return;
    }

    VectorEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector.set_to_zero();
    NewEntry.Vector[ii] += value;
    VectorFields.push_back(NewEntry);
}

inline NodeVN& NodeVN::operator=(const NodeVN& n)
{
    SIZE_X       = n.SIZE_X;
    VectorFields = n.VectorFields;
    return *this;
}

inline NodeVN& NodeVN::operator+=(const NodeVN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add(i->index, i->Vector);
    }
    return *this;
}

inline NodeVN& NodeVN::operator*=(const NodeVN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        multiply(i->index, i->Vector);
    }
    return *this;
}

inline NodeVN& NodeVN::operator-=(const NodeVN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        dVectorN vec = i->Vector * -1.0;
        add(i->index, vec);
    }
    return *this;
}

inline NodeVN NodeVN::operator+(const NodeVN& n) const
{
    NodeVN result = n;
    for (auto i = VectorFields.cbegin(); i < VectorFields.cend(); ++i)
    {
        result.add(i->index, i->Vector);
    }
    return result;
}

inline NodeVN NodeVN::operator-(const NodeVN& n) const
{
    NodeVN result = n;
    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->Vector *= -1.0;
    }

    for (auto i = VectorFields.cbegin(); i < VectorFields.cend(); ++i)
    {
        result.add(i->index, i->Vector);
    }
    return result;
}

inline NodeVN NodeVN::operator*(const double n) const
{
    NodeVN result = *this;
    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->Vector *= n;
    }
    return result;
}

inline NodeVN& NodeVN::operator*=(const double n)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        i->Vector *= n;
    }
    return *this;
}

inline void NodeVN::pack(std::vector<double>& buffer)
{
    buffer.push_back(SIZE_X);
    buffer.push_back(VectorFields.size());
    for (size_t i = 0; i < VectorFields.size(); ++i)
    {
        buffer.push_back(VectorFields[i].index);
        for (size_t n = 0; n < SIZE_X; ++n)
        {
            buffer.push_back(VectorFields[i].Vector[n]);
        }
    }
}
inline void NodeVN::unpack(std::vector<double>& buffer, size_t& it)
{
    SIZE_X = buffer[it]; ++it;
    VectorFields.resize(buffer[it]); ++it;
    for (size_t i = 0; i < VectorFields.size(); ++i)
    {
        VectorFields[i].index = buffer[it]; ++it;
        VectorFields[i].Vector.Allocate(SIZE_X);
        for (size_t n = 0; n < SIZE_X; ++n)
        {
            VectorFields[i].Vector[n] = buffer[it]; ++it;
        }
    }
}

inline std::string NodeVN::print()
{
    std::stringstream out;
    for (iterator i = begin(); i < end(); ++i)
    {
        out << "index :  "     << i->index         << "\n";
        out << "Size_X : "     << SIZE_X           << "\n";
        out << "vec.size() : " << i->Vector.size() << "\n";
        out << "value : "      << "\n";
        out << i->Vector.print();
    }
    return out.str();
}

inline void NodeVN::read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    inp.read(reinterpret_cast<char*>(&SIZE_X), sizeof(size_t));
    VectorFields.resize(size);
    for(auto &Field : VectorFields)
    {
        inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
        Field.Vector.read_binary(inp);
    }
}

inline void NodeVN::write(std::ostream& outp) const
{
    size_t size = VectorFields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
    outp.write(reinterpret_cast<const char*>(&SIZE_X), sizeof(size_t));

    for(auto &Field : VectorFields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
        Field.Vector.write_binary(outp);
    }
}

} //namespace openphase
#endif
