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

 *   File created :   2018
 *   Main contributors :    Muhammad Adil Ali; Oleg Shchyglo
 *
 */

#ifndef NodeVectorN_H
#define NodeVectorN_H

#include "dVectorN.h"

namespace openphase
{

class VecEntry                                                              	///< Individual vector entry. Used in the NodeVectorN class as a storage unit.
{
public:
    VecEntry():
        Vector(0),
        index(0)
    {
    };
    VecEntry(size_t size, size_t idx):
        Vector(size),
        index(idx)
    {
    };
    VecEntry(const VecEntry& other):
        Vector(other.Vector),
        index(other.index)
    {
    };
    VecEntry& operator=(const VecEntry& other)
    {
        Vector = other.Vector;
        index = other.index;
        return *this;
    };
    dVectorN Vector;
    size_t index;
};

/***********************************************************************************/

class NodeVectorN                                                           	///< Basic class to store the vector valued quantities for each phase field.
{
public:

    NodeVectorN():
        SIZE_X(0)
    {
    };
    NodeVectorN(const size_t size):
        SIZE_X(size)
    {
    };
    NodeVectorN(const NodeVectorN& other):
        VectorFields(other.VectorFields),
        SIZE_X(other.SIZE_X)
    {
    };
    void Allocate(const size_t size)
    {
        SIZE_X = size;
    }
    size_t size_of_vector() const
    {
        return SIZE_X;
    }

    bool exist(const size_t idx) const;                                         ///< True if entry with index "idx" is present, false otherwise.

    void set_to_value(const size_t idx, const int value);                       ///< Set all entries to value for phase field index idx.
    void set_to_zero(const size_t idx);                                         ///< Set all entries to zero for phase field index idx.
    void set(const size_t idx, const size_t ii, const double value);            ///< Set component ii for phase field index idx to value.
    void set(const size_t idx, dVectorN value);                                 ///< Set dVectorN for phase field index idx.
    double   get(const size_t idx, const size_t ii) const;                      ///< Return value of component ii for phase field index idx.
    dVectorN get(const size_t idx) const;                                       ///< Return values for phase field index idx.
    void         add(const size_t idx, dVectorN value);                         ///< Increment value using two indeces.
    NodeVectorN  add(const NodeVectorN& n) const;                               ///< Add two nodes.
    void         add(const size_t idx, const size_t ii, const double value);    ///< Add in component ii for phase field index idx to value.
    void    multiply(const size_t idx, dVectorN value);

    NodeVectorN& operator=(const NodeVectorN& n);
    NodeVectorN& operator+=(const NodeVectorN& n);
    NodeVectorN& operator*=(const NodeVectorN& n);
    NodeVectorN& operator-=(const NodeVectorN& n);
    NodeVectorN& operator*=(const double n);
    NodeVectorN  operator+(const  NodeVectorN& n) const;                        ///< Plus operator. Takes as input another Node type entry.
    NodeVectorN  operator-(const  NodeVectorN& n) const;                        ///< Minus operator. Takes as input another Node type entry.
    NodeVectorN  operator*(const double n) const;                               ///< Multiply all fields by a number.

    void   clear(){VectorFields.clear();};                                      ///< Empties the vector fields.
    size_t  size() const {return VectorFields.size();};                         ///< Returns the size of vector fields.
    void read(std::istream& inp);                                               ///< Reads NodeVectorN content from the input stream.
    void write(std::ostream& outp) const;                                       ///< Writes NodeVectorN content to the output stream.

    typedef typename std::vector<VecEntry>::iterator iterator;                  ///< Iterator over the vector fields
    typedef typename std::vector<VecEntry>::const_iterator citerator;           ///< Constant iterator over the vector fields
    iterator  begin()        {return VectorFields.begin();};                    ///< Iterator to the begin of vector fields
    iterator  end()          {return VectorFields.end();};                      ///< Iterator to the end of vector fields
    citerator cbegin() const {return VectorFields.cbegin();};                   ///< Constant iterator to the begin of vector fields
    citerator cend()   const {return VectorFields.cend();};                     ///< Constant iterator to the end of vector fields
    iterator  erase(iterator it) {return VectorFields.erase(it);};              ///< Erase a single record pointed by iterator it.

    std::string print()
    {
        std::stringstream out;
        for (iterator i = begin(); i < end(); ++i)
        {
            out << "index :  "     << i->index        << "\n";
            out << "Size_X : "     << SIZE_X          << "\n";
            out << "vec.size() : " << i->Vector.size() << "\n";
            out << "value : "      << "\n";
            out << i->Vector.print(',');
        }
        return out.str();
    }

//    dVectorN& operator[](size_t idx) dangerous method because it will grow the container if idx is not found
//    {
//        for (auto& node : VectorFields)
//        {
//            if (node.index == idx)
//            {
//                return node.Vector;
//            }
//        }
//        VecEntry NewEntry(SIZE_X, idx);
//        NewEntry.Vector.set_to_zero();
//        VectorFields.push_back(NewEntry);
//        return VectorFields.back().Vector;
//    }

    void pack(std::vector<double>& buffer)
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
    void unpack(std::vector<double>& buffer, size_t& it)
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
    void Write(std::ostream& outp) const
    {
        size_t size = VectorFields.size();
        outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
        outp.write(reinterpret_cast<const char*>(&SIZE_X), sizeof(size_t));

        for(auto &Field : VectorFields)
        {
            outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
            Field.Vector.Write(outp);
        }
    }
    void Read(std::istream& inp)
    {
        size_t size = 0;
        inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
        inp.read(reinterpret_cast<char*>(&SIZE_X), sizeof(size_t));
        VectorFields.resize(size);

        for(auto &Field : VectorFields)
        {
            inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
            Field.Vector.Read(inp);
        }

    }
protected:
private:
    std::vector<VecEntry> VectorFields;
    size_t SIZE_X;
};

inline bool NodeVectorN::exist(const size_t idx) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if(i->index == idx)
    {
        return true;
    }
    return false;
}

inline void NodeVectorN::set_to_value(const size_t idx, const int value)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        i->Vector.set_to_value(value);
        return;
    }

    VecEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector.set_to_value(value);
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set_to_zero(const size_t idx)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        i->Vector.set_to_zero();
        return;
    }

    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector.set_to_zero();
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set(const size_t idx, const size_t ii, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::set() ii is greater than dVectorN.size()");
        i->Vector[ii] = value;
        return;
    }

    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector.set_to_zero();
    NewEntry.Vector[ii] = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::set(const size_t idx, dVectorN Vector)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::set() Vector size is not equal");
        i->Vector = Vector;
        return;
    }

    VecEntry NewEntry(SIZE_X,idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline double NodeVectorN::get(const size_t idx, const size_t ii) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::get() ii is greater than dVectorN.size()");
        return i->Vector[ii];
    }
    return 0.0;
}

inline dVectorN NodeVectorN::get(const size_t idx) const
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

inline void NodeVectorN::add(const size_t idx, dVectorN Vector)
{
    for(auto i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::add Vector size is not equal");
        i->Vector += Vector;
        return;
    }

    VecEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline void NodeVectorN::multiply(const size_t idx, dVectorN Vector)
{
    for(auto i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(i->Vector.size() == Vector.size() && "NodeVectorN::multiply() Vector size is not equal");
        i->Vector *= Vector;
        return;
    }

    VecEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector = Vector;
    VectorFields.push_back(NewEntry);
}

inline NodeVectorN NodeVectorN::add(const NodeVectorN& n) const
{
    NodeVectorN result = n;
    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->Vector);
    }
    return result;
}

inline void NodeVectorN::add(const size_t idx, const size_t ii, const double value)
{
    for(iterator i = begin(); i < end(); ++i)
    if(i->index == idx)
    {
        assert(ii < i->Vector.size() && "NodeVectorN::add() idx is greater than dVectorN.size()");
        i->Vector[ii] += value;
        return;
    }

    VecEntry NewEntry(SIZE_X, idx);
    NewEntry.Vector.set_to_zero();
    NewEntry.Vector[ii] += value;
    VectorFields.push_back(NewEntry);
}

inline NodeVectorN& NodeVectorN::operator=(const NodeVectorN& n)
{
    SIZE_X       = n.SIZE_X;
    VectorFields = n.VectorFields;
    return *this;
}

inline NodeVectorN& NodeVectorN::operator+=(const NodeVectorN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        add(i->index, i->Vector);
    }
    return *this;
}

inline NodeVectorN& NodeVectorN::operator*=(const NodeVectorN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        multiply(i->index, i->Vector);
    }
    return *this;
}

inline NodeVectorN& NodeVectorN::operator-=(const NodeVectorN& n)
{
    SIZE_X = n.SIZE_X;
    for (auto i = n.cbegin(); i < n.cend(); ++i)
    {
        dVectorN vec = i->Vector * -1.0;
        add(i->index, vec);
    }
    return *this;
}

inline NodeVectorN NodeVectorN::operator+(const NodeVectorN& n) const
{
    NodeVectorN result = n;
    for (auto i = VectorFields.cbegin(); i < VectorFields.cend(); ++i)
    {
        result.add(i->index, i->Vector);
    }
    return result;
}

inline NodeVectorN NodeVectorN::operator-(const NodeVectorN& n) const
{
    NodeVectorN result = n;
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

inline NodeVectorN NodeVectorN::operator*(const double n) const
{
    NodeVectorN result = *this;
    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->Vector *= n;
    }
    return result;
}

inline NodeVectorN& NodeVectorN::operator*=(const double n)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        i->Vector *= n;
    }
    return *this;
}

inline void NodeVectorN::read(std::istream& inp)
{
    size_t size = 0;
    inp.read(reinterpret_cast<char*>(&size), sizeof(size_t));
    inp.read(reinterpret_cast<char*>(&SIZE_X), sizeof(size_t));
    VectorFields.resize(size);
    for(auto &Field : VectorFields)
    {
        inp.read(reinterpret_cast<char*>(&Field.index), sizeof(size_t));
        Field.Vector.Read(inp);
    }
}

inline void NodeVectorN::write(std::ostream& outp) const
{
    size_t size = VectorFields.size();
    outp.write(reinterpret_cast<const char*>(&size), sizeof(size_t));
    outp.write(reinterpret_cast<const char*>(&SIZE_X), sizeof(size_t));

    for(auto &Field : VectorFields)
    {
        outp.write(reinterpret_cast<const char*>(&Field.index), sizeof(size_t));
        Field.Vector.Write(outp);
    }
}

} //namespace openphase
#endif
