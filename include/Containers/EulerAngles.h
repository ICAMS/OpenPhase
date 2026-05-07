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

#ifndef EULERANGLES_H
#define EULERANGLES_H

#include "Globals.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

namespace openphase
{
//================================ Declarations ================================

enum class AngleModes : int                                                     ///< Angle representation modes
{
    Radians,
    Degrees
};

const std::vector<std::string> RadiansStrings
     {"RADIANS", "RAD"};

const std::vector<std::string> DegreesStrings
     {"DEGREES", "DEG"};

const std::unordered_map<std::string, AngleModes> AngleRepresentationTable =
{
    {"RADIANS", AngleModes::Radians},
    {"RAD", AngleModes::Radians},

    {"DEGREES", AngleModes::Degrees},
    {"DEG", AngleModes::Degrees}
};

AngleModes AngleRepresentation(std::string locAngleRepresentationString);       ///< Converts angle representation mode string to the corresponding enum entry

enum class EulerConventions: int
    {XZX, XYX, YXY, YZY, ZXZ, ZYZ, /* proper Euler angles*/
     XYZ, YZX, ZXY, XZY, ZYX, YXZ, /* Tait–Bryan angles*/
     NNN /*default convention -> convention not set*/};

const std::vector<std::string> EulerConventionStrings
     {"XZX", "XYX", "YXY", "YZY", "ZXZ", "ZYZ", /* proper Euler angles*/
      "XYZ", "YZX", "ZXY", "XZY", "ZYX", "YXZ", /* Tait–Bryan angles*/
      "NNN" /*default convention -> convention not set*/};

const std::unordered_map<std::string, EulerConventions> EulerConventionTable =
{
    {"XZX", EulerConventions::XZX},
    {"XYX", EulerConventions::XYX},
    {"YXY", EulerConventions::YXY},
    {"YZY", EulerConventions::YZY},
    {"ZXZ", EulerConventions::ZXZ},
    {"ZYZ", EulerConventions::ZYZ},

    {"XYZ", EulerConventions::XYZ},
    {"YZX", EulerConventions::YZX},
    {"ZXY", EulerConventions::ZXY},
    {"XZY", EulerConventions::XZY},
    {"ZYX", EulerConventions::ZYX},
    {"YXZ", EulerConventions::YXZ}
};

EulerConventions EulerConvention(std::string locConventionString);              ///< Converts Euler angles' convention string to the corresponding enum entry

class Quaternion;
class dMatrix3x3;
class dVector3;

class OP_EXPORTS EulerAngles                                                    ///< Euler angles and their Cos() and Sin().
{
 public:
    std::array<double,3> Q;                                                     ///< Stores three Euler angles in radians

    EulerAngles();                                                              ///< Default constructor
    EulerAngles(std::initializer_list<double> Angles,
                std::string locConventionString);                               ///< Constructor, initializes the Euler angles using initializer list and convention string
    EulerAngles(std::initializer_list<double> Angles,
                EulerConventions locConvention);                                ///< Constructor, initializes the Euler angles using initializer list and convention
    EulerAngles(const EulerAngles& rhs);                                        ///< Copy constructor

    bool operator==(const EulerAngles& rhs);                                    ///< Comparison operator. Returns true if two EulerAngles are of the same convention and their values are within DBL_EPSILON from one another.

    void set(const double q1, const double q2, const double q3,
            std::string locConventionString);                                   ///< Sets Euler angles to the specified values following given convention string
    void set(const double q1, const double q2, const double q3,
             EulerConventions locConvention);                                   ///< Sets Euler angles to the specified values following given convention
    void set_convention(const std::string locConventionString);                 ///< Sets Euler angles convention using given convention string
    void set_convention(const EulerConventions locConvention);                  ///< Sets Euler angles convention
    void set(const dMatrix3x3& RotMatrix, std::string locConventionString);     ///< Sets Euler angles form the rotation matrix using give convention string
    void set(const dMatrix3x3& RotMatrix, EulerConventions locConvention);      ///< Sets Euler angles form the rotation matrix using give convention
    void set(Quaternion Quat, std::string locConventionString,
             const bool active = true);                                         ///< Sets Euler angles from the quaternion using specified convention and type of rotation (active/passive) string
    void set(Quaternion Quat, EulerConventions locConvention,
             const bool active = true);                                         ///< Sets Euler angles from the quaternion using specified convention and type of rotation (active/passive)
    //void set(dVector3 Axis, const double Angle, std::string locConventionString);///< Sets Euler angles from axis-angle entries using given convention string.
    //void set(dVector3 Axis, const double Angle, EulerConventions locConvention);///< Sets Euler angles from axis-angle entries using given convention.

    void set_trigonometric_functions(void);                                     ///< Sets internal sines and cosines of Euler angles for faster computations.
    void set_to_zero(void);                                                     ///< Sets Euler angles to zero.

    void add(const double q1, const double q2, const double q3);                ///< Adds specified values to the corresponding Euler angles

    EulerAngles& operator= (const EulerAngles& rhs);                            ///< Assignment operator.
    EulerAngles  operator+ (const EulerAngles& rhs) const;                      ///< Returns the sum of two Euler angles
    EulerAngles  operator- (const EulerAngles& rhs) const;                      ///< Returns the difference of two Euler angles
    EulerAngles& operator+=(const EulerAngles& rhs);                            ///< Adds rhs to the current Euler angles
    EulerAngles& operator-=(const EulerAngles& rhs);                            ///< Subtracts rhs to the current Euler angles
    EulerAngles  operator* (const double rhs) const;                            ///< Returns Euler angles multiplied by the specified factor
    EulerAngles& operator*=(const double rhs);                                  ///< Multiplies Euler angles by the specified factor

    std::string get_convention(void) const;                                     ///< Returns Euler angle convention
    dMatrix3x3  get_rotation_matrix(void) const;                                ///< Returns rotation matrix corresponding to Euler angles
    Quaternion  get_quaternion(const bool Active = true) const;                 ///< Returns quaternion corresponding to Euler angles considering rotation type (active/passive)
    void        get_axis_angle(dVector3& Axis, double& Angle) const;            ///< Returns axis-angle representation of the stored Euler angles

    void read_binary(std::istream& inp);                                        ///< Reads Euler angles from a file stream in binary format
    void read_ASCII(std::istream& inp);                                         ///< Reads Euler angles from a file stream in ASCII format
    void read_degrees(std::istream& inp);                                       ///< Reads Euler angles in degrees from a file stream in ASCII format

    void write_binary(std::ostream& out) const;                                 ///< Writes Euler angles to a file stream in binary format
    void write_ASCII(std::ostream& out,
                     const int precision = 16, const char sep = ' ') const;     ///< Writes Euler angles to a file stream in ASCII format
    void write_degrees(std::ostream& out,
                       const int precision = 16, const char sep = ' ') const;   ///< Writes Euler angles in degrees to a file stream in ASCII format

    std::string get_output_string(std::string ModeString,
            const int precision = std::numeric_limits<double>::digits10 + 1,
            const char sep = ' ') const;                                        ///< Returns a string of Euler angles values with given precision separated by "sep"

    std::string get_output_string(AngleModes Mode,
            const int precision = std::numeric_limits<double>::digits10 + 1,
            const char sep = ' ') const;                                        ///< Returns a string of Euler angles values with given precision separated by "sep"

    std::string print(void) const;                                              ///< Returns formatted string containing Euler angles and corresponding convention in degrees
    std::string print_entire(void) const;                                       ///< Returns formatted string containing full container content.

 protected:
 private:
    bool is_set;                                                                ///< Indicates whether sin/cos are set.
    EulerConventions Convention;                                                ///< Stores the convention for the Euler angles
    std::array<double,3> CosQ;                                                  ///< Stores cosines of the stored Euler angles
    std::array<double,3> SinQ;                                                  ///< Stores sines of the stored Euler angles
};

}// namespace openphase
#endif
