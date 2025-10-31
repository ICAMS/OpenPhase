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

 *   File created :   2013
 *   Main contributors :   Philipp Engels; Raphael Schiedung
 *
 */

#ifndef CONSOLEOUTPUT_H
#define CONSOLEOUTPUT_H

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Globals.h"
#include "Containers/TypeTraits.h"

namespace openphase
{

class RunTimeControl;

enum class OP_EXPORTS VerbosityLevels : int                                                ///< Console output verbosity levels
{
    Silent  = 0,                                                                ///< No console output except for exits.
    Warning = 1,                                                                ///< No console output except for warnings and exits.
    Normal  = 2,                                                                ///< All outputs are enabled
    Debug   = 3                                                                 ///< All outputs plus debug outputs
};

enum class OP_EXPORTS floatfield                                                ///< Console output format for floating point numbers
{
    Default,                                                                    ///< Use default floating-point notation for output
    Scientific,                                                                 ///< Use scientific floating-point notation for output
    Fixed                                                                       ///< Use fixed floating-point notation for output
};

class OP_EXPORTS ConsoleOutput                                                  ///< A collection of methods to format console output in a standardized way

{
 public:

    static constexpr auto   StandardNotation    = floatfield::Default;
    static constexpr int    StandardPrecision   = 6;
    static constexpr size_t StandardColumnWidth = 40;
    static constexpr size_t LineLength          = 80;

    inline static VerbosityLevels OutputVerbosity = VerbosityLevels::Normal;    ///< Sets the level of console output verbosity

    /// Return white space for positive values to align positive and negative
    /// numbers vertically
    template <typename T>
    static inline std::string sgn(const T value)
    {
         return (value >= 0)?" ":"";
    }

    /// Writes a Message, Blank Line, or Line of char to screen
    static void Write(const std::string message)
    {
        if (message.size() > 1)
        {
            WriteSimple(message);
        }
        else if (message.size() == 1)
        {
            WriteLine(message);
        }
        else
        {
            WriteBlankLine();
        }
    }

    /// Writes name and value to screen
    ///
    /// @param name of the value printed to screen
    /// @param value printed to screen
    /// @param precision of value printed to screen
    template <typename T>
    static void Write(const std::string name, const T value,
            const int precision=4)
    {
        WriteStandard(name, value, precision);
    }

    /// Writes name and value to screen and to line buffer for later file write
    ///
    /// @param line is a buffer which can be written to log file
    /// @param tStep current time step
    /// @param name of the value printed to screen
    /// @param value printed to screen
    /// @param precision of value printed to screen
    /// @param sep is separator char used in the log file
    template <typename T>
    static void WriteWithLog(std::array<std::stringstream,2> &line,
            const int tStep, const std::string name, const  T value,
            const int precision=4, const char sep = ',')
    {
        if (tStep == 0)
        {
            line[0] << name << sep;
        }
        line[1] << value << sep;
        WriteStandard(name, value, precision);
    }

    /// Writes line buffer to log file
    ///
    /// @param log is the log file where line will be written to
    /// @param line is a buffer which can be written to log file
    /// @param tStep current time step
    /// @param precision of value printed to screen
    /// @param sep is separator char used in the log file
    static void WriteLineToLogfile(std::fstream &log,
            const std::array<std::stringstream,2> &line,
            const int tStep, const char sep = ',')
    {
        #ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
        {
        #endif
            if (tStep == 0)
            {
                log << tStep << sep << line[0].str() << std::endl;
            }
            log << tStep << sep << line[1].str() << std::endl;
        #ifdef MPI_PARALLEL
        }
        #endif
    }

    template <typename TargetType>
    static void SetFloationgPointOutputFormat(TargetType& Target,
            const int precision, floatfield notation)
    {
        Target << std::setprecision(precision);
        switch(notation)
        {
            case floatfield::Default:
            {
                Target << std::defaultfloat;
                break;
            }
            case floatfield::Scientific:
            {
                Target << std::scientific;
                break;
            }
            case floatfield::Fixed:
            {
                Target << std::fixed;
                break;
            }
        }
    };

    /// Returns string of left column of ColumnWidth
    static std::string StandardLHS(std::string Left, const size_t ColumnWidth)
    {
        if (Left.size() < ColumnWidth) Left.resize(ColumnWidth, ' ');
        Left.append(": ");
        return Left;
    };

    /// Returns string of matrix in standard output format
    template <typename T>
    static std::string StandardMatrix(const T value,
            const int precision = StandardPrecision,
            const floatfield notation = StandardNotation,
            const size_t ColumnWidth = StandardColumnWidth/5)
    {
        std::stringstream ss;
        SetFloationgPointOutputFormat(ss, precision, notation);
        for (size_t i = 0; i < value.sizeX(); i++)
        {
            ss << "\n|| ";
            for (size_t j = 0; j < value.sizeY(); j++)
            {
                ss << sgn(value(i,j)) << std::setw(ColumnWidth) << value(i,j) << " ";
            }
            ss << "||";
        }
        ss << "\n";
        return ss.str();
    };

    /// Returns string of vector in standard output format
    template <typename T>
    static std::string StandardVector(const T value,
            const int precision = StandardPrecision,
            const floatfield notation = StandardNotation,
            const size_t ColumnWidth = StandardColumnWidth/5)
    {
        std::stringstream ss;
        SetFloationgPointOutputFormat(ss, precision, notation);
        ss << "( ";
        for (size_t i = 0; i < value.size(); i++)
        {
            ss << sgn(value[i]) << std::setw(ColumnWidth)  << value[i] << " ";
        }
        ss << ")\n";
        return ss.str();
    };

    /// Returns string of floating point number in standard output format
    template <typename T>
    static std::string StandardFloat(const T value,
            const int precision = StandardPrecision,
            const floatfield notation = StandardNotation)
    {
        std::stringstream ss;
        SetFloationgPointOutputFormat(ss, precision, notation);
        ss << sgn(value) << value << "\n";
        return ss.str();
    };

    /// Returns string of integral number in standard output format
    template <typename T>
    static std::string StandardIntegral(const T value,
            const int precision = StandardPrecision,
            const floatfield notation = StandardNotation)
    {
        std::stringstream ss;
        ss << sgn(value) << value << "\n";
        return ss.str();
    };

    template <typename T>
    static std::string GetStandard(const std::string& Left, const T Right,
            [[maybe_unused]] const int precision = StandardPrecision,
            [[maybe_unused]] const floatfield notation = StandardNotation,
            const size_t ColumnWidth = StandardColumnWidth)
    {
        std::string out = StandardLHS(Left,ColumnWidth);
        if constexpr (is_matrix<T>::value)
        {
            out.append(StandardMatrix(Right,precision,notation,ColumnWidth/5));
        }
        else if constexpr (is_vector<T>::value)
        {
            out.append(StandardVector(Right,precision,notation,ColumnWidth/5));
        }
        else if constexpr (std::is_floating_point<T>::value)
        {
            out.append(StandardFloat(Right,precision,notation));
        }
        else if constexpr (std::is_signed<T>::value)
        {
            out.append(StandardIntegral(Right,precision,notation));
        }
        else
        {
            std::stringstream ss;
            ss << " " << Right << "\n";
            out.append(ss.str());
        }
        return out;
    }

    template <typename T>
    static std::string GetStandardNarrow(const std::string& Left, const T Right,
            const int precision = StandardPrecision,
            const floatfield notation = StandardNotation)
    {
        return GetStandard(Left, Right, precision, notation, StandardColumnWidth/2);
    }

    /// Writes GetStandard result to std::cout if MPI-Rank 0
    template <typename T>
    static void WriteStandard(const std::string& Left, const T Right,
            [[maybe_unused]] const int precision = StandardPrecision,
            [[maybe_unused]] floatfield notation = StandardNotation,
            const size_t ColumnWidth = StandardColumnWidth)
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
#endif
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            std::cout << GetStandard(Left,Right,precision,notation,ColumnWidth);
        }
    }

    /// Writes GetStandardNarrow result to std::cout if MPI-Rank 0
    template <typename T>
    static void WriteStandardNarrow(const std::string& Left, const T Right,
            [[maybe_unused]] const int precision = StandardPrecision,
            [[maybe_unused]] floatfield notation = StandardNotation)
    {
#ifdef MPI_PARALLEL
        if (MPI_RANK == 0)
#endif
        if(OutputVerbosity >= VerbosityLevels::Normal)
        {
            std::cout << GetStandard(Left,Right,precision,notation,StandardColumnWidth/2);
        }
    }

    // Static console output methods
    static std::string get_time();
    static void WriteLine(const std::string LineType = "-");
    static void WriteLineInsert(const std::string Insert,
            const std::string LineType = "-");
    static void WriteBlankLine(void);
    static void WriteTimeStep(const RunTimeControl& RTC,
            const std::string Message = "",
            const size_t ColumnWidth = StandardColumnWidth);
    static void WriteTimeStep(const int tStep, const int nSteps,
            const std::string Message = "",
            const size_t ColumnWidth = StandardColumnWidth);
    static void WriteTimeStep(const int tScreenWrite, const int tStep,
            const int nSteps, const std::string Message = "",
            const size_t ColumnWidth = StandardColumnWidth);
    static void WriteSimple(const std::string Message);
    static void WriteCoordinate(const int x, const int y,
            const int z, const double dx);
    static void WriteWithinMethod(const std::string Message,
            const std::string Instance = "", const std::string Method = "");
    static void WriteWarning(const std::string Message,
            const std::string Instance = "", const std::string Method = "");
    static void WriteExit(const std::string Message,
            const std::string Instance = "", const std::string Method = "");
    static void WriteStartScreen(void);
    static void PressEnterToContinue(void);

    static void StartProgressIndicator(std::string IndicatorTitle);
    static void AdvanceProgressIndicator(double start_value, double end_value, double current_value, int& position);
    static void EndProgressIndicator(void);

    template <typename T>
    static std::string to_string_with_precision(const T a_value, const int n = 6)
    {
        std::ostringstream out;
        out << std::setprecision(n) << a_value;
        return out.str();
    }
};
}// namespace openphase
#endif // CONSOLEOUTPUT_H
