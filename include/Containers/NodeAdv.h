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

 *   File created :   2011
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Dmitry Medvedev
 *
 */

#ifndef NODEADV_H
#define NODEADV_H

#include<cstddef>
#include<vector>

namespace openphase
{
/***************************************************************/

struct AdvEntry                                                                 /// Structure for storing the individual phase field derivatives entry. Used in the NodeA class as a storage unit.
{
    size_t index;                                                               ///  Phase field index
    double dXY;                                                                 ///  Second derivative w.r.t. X and Y
    double dXZ;                                                                 ///  Second derivative w.r.t. X and Z
    double dYZ;                                                                 ///  Second derivative w.r.t. Y and Z
    double dXYZ;                                                                ///  Third derivative w.r.t. X, Y and z
};

/***************************************************************/
class NodeAdv                                                                   /// Stores the phase fields derivatives at a grid point. Provide access and manipulation methods for the entries.
{
 public:
    double operator[](const size_t n);                                          /// Index operator for accessing the n's phase field derivatives

    void   set_dXY (const size_t n, const double value);                        /// Set second derivative w.r.t X and Y
    void   set_dXZ (const size_t n, const double value);                        /// Set second derivative w.r.t X and Z
    void   set_dYZ (const size_t n, const double value);                        /// Set second derivative w.r.t Y and Z
    void   set_dXYZ(const size_t n, const double value);                        /// Set third derivative w.r.t X, Y and Z

    void   set(const size_t n, const double dXY,
                               const double dXZ,
                               const double dYZ,
                               const double dXYZ);                              /// Set all derivatives simultaneously.

    void   add_dXY (const size_t n, const double value);                        /// Increment second derivative w.r.t X and Y
    void   add_dXZ (const size_t n, const double value);                        /// Increment second derivative w.r.t X and Z
    void   add_dYZ (const size_t n, const double value);                        /// Increment second derivative w.r.t Y and Z
    void   add_dXYZ(const size_t n, const double value);                        /// Increment third derivative w.r.t X, Y and Z

    double get_dXY (const size_t n);                                            /// Return second derivative w.r.t X and Y
    double get_dXZ (const size_t n);                                            /// Return second derivative w.r.t X and Z
    double get_dYZ (const size_t n);                                            /// Return second derivative w.r.t Y and Z
    double get_dXYZ(const size_t n);                                            /// Return second derivative w.r.t X, Y and Z

    void   clear()                                                              /// Emptys the storage.
           {Fields.clear();};
    size_t  size()                                                              /// Returns the size of storage.
           {return Fields.size();};
    typedef std::vector<AdvEntry>::iterator iterator;                           /// Iterator over phase field derivatives storage
    typedef std::vector<AdvEntry>::const_iterator citerator;                    /// Constant iterator over phase field derivatives storage
    iterator begin()
             {return Fields.begin();};                                          /// Iterator to the begin of phase field derivatives storage
    iterator end()                                                              /// Iterator to the end of phase field derivatives storage
             {return Fields.end();};
    citerator cbegin()
             {return Fields.begin();};                                          /// Constant iterator to the begin of phase field derivatives storage
    citerator cend()                                                            /// Constant iterator to the end of phase field derivatives storage
             {return Fields.end();};

 protected:
 private:
    std::vector<AdvEntry> Fields;                                               /// List of nonvanishing phase field derivatives.
};

inline void NodeAdv::set_dXY(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXY   = value;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::set_dXZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::set_dYZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dYZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dYZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::set_dXYZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXYZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXYZ  = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::set(const size_t n, const double dXY,
                             const double dXZ,
                             const double dYZ,
                             const double dXYZ)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY  = dXY;
        i->dXZ  = dXZ;
        i->dYZ  = dYZ;
        i->dXYZ = dXYZ;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index     = n;
    NewEntry.dXY   = dXY;
    NewEntry.dXZ   = dXZ;
    NewEntry.dYZ   = dYZ;
    NewEntry.dXYZ  = dXYZ;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::add_dXY(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index     = n;
    NewEntry.dXY   = value;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::add_dXZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::add_dYZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dYZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dYZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeAdv::add_dXYZ(const size_t n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXYZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXYZ  = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    Fields.push_back(NewEntry);
}

inline double NodeAdv::get_dXY(const size_t n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXY;
    }
    return 0.0;
}

inline double NodeAdv::get_dXZ(const size_t n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXZ;
    }
    return 0.0;
}

inline double NodeAdv::get_dYZ(const size_t n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dYZ;
    }
    return 0.0;
}

inline double NodeAdv::get_dXYZ(const size_t n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXYZ;
    }
    return 0;
}
}// namespace openphase
#endif
