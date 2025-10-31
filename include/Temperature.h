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

 *   File created :   2012
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "Includes.h"
#include "H5Interface.h"

namespace openphase
{
class AdvectionHR;
class BoundaryConditions;
class Composition;
class FlowSolverLBM;
class GrainsProperties;
class H5Interface;
class PhaseField;
class Settings;
class Velocities;

enum class ControlModes                                                         ///< Modes of temperature control
{
    Linear,                                                                     ///< Linear cooling/heating
    Newton,                                                                     ///< Newton's law of cooling/heating
    None                                                                        ///< No cooling/heating is applied
};

enum class LatentHeatModes                                                      ///< Latent heat consideration modes
{
    Off,                                                                        ///< No latent heat release is considered
    Global,                                                                     ///< Integrated latent heat effect on temperature is considered
    Local                                                                       ///< Latent heat release is considered locally in every grid cell
};

struct TemperatureParameters                                                    ///< Temperature control parameters
{
    double ActivationTime;                                                      ///< Time of activation
    double LinearRate;                                                          ///< Rate for "Linear" cooling/heating mode [K/s]
    double NewtonRate;                                                          ///< Rate for "Newton" cooling/heating mode [1/s]
    double HoldingTemperature;                                                  ///< Holding temperature
    ControlModes Mode;                                                          ///< Temperature control mode
};

class Temperature1Dextension                                                    ///< 1D temperature field extension storage
{
 public:
    void Initialize(size_t size, iVector3 direction)
    {
        if(size <= 0) // Checks for sufficient extension size
        {
            std::stringstream message;
            message << "1D temperature extension size = " << size << " is incorrect!\n"
                    << "It should be in the range [1, inf) for correct operation\n";
            ConsoleOutput::WriteExit(message.str(), "Temperature1Dextension", "Initialize()");
            OP_Exit(EXIT_FAILURE);
        }

        Data.Allocate(size, 1);
        DataTMP.Allocate(size, 1);
        Direction = direction;
        Qdot = 0.0;
    }
    Temperature1Dextension& operator=(const Temperature1Dextension& RHS)
    {
        Qdot      = RHS.Qdot;
        Data      = RHS.Data;
        DataTMP   = RHS.DataTMP;
        Direction = RHS.Direction;

        return *this;
    }
    void setBC(BoundaryConditionTypes extBC)                                    ///< Sets boundary condition at the far end of the extension
    {
        switch(extBC)
        {
            case BoundaryConditionTypes::Periodic:
            {
                std::string message = std::string("If 1D temperature extension is active, periodic boundary conditions\n")
                                    + std::string("along the extension dimension are not permitted.\n")
                                    + std::string("Permitted boundary conditions: NoFlux, Free and Fixed.");
                ConsoleOutput::WriteWarning(message, "Temperature1Dextension", "setBC()");
                break;
            }
            case BoundaryConditionTypes::NoFlux:
            case BoundaryConditionTypes::Mirror:
            {
                Data(size()) = Data(size()-1);
                break;
            }
            case BoundaryConditionTypes::Free:
            {
                Data(size()) = Data(size()-1)*2.0 - Data(size()-2);
                break;
            }
            default:
            case BoundaryConditionTypes::Fixed:
            {
                break;
            }
        }
    }
    void setFixedBC(double value)                                               ///< Sets fixed temperature at the far end of the extension to the specified value
    {
        Data(size()) = value;
    }
    void moveFrame(const int dx, const BoundaryConditionTypes extBC)            ///< Moves the data according to the moving frame motion
    {
        if(dx > 0)
        for(long int x = 0; x < (long int) Data.size(); x++)
        {
            Data(x) = Data(x+dx);
        }
        if(dx < 0)
        for(long int x = (long int) Data.size() - 1; x >= 0; x--)
        {
            Data(x) = Data(x+dx);
        }
        setBC(extBC);
    }
    bool isActive() const                                                       ///< Returns true if extension was activated, false otherwise
    {
        return Data.size() > 0;
    }
    size_t size() const
    {
        return Data.size();
    }
    void store_temporary(void)
    {
        DataTMP = Data;
    }
    void read(std::ifstream& out)
    {
        out.read(reinterpret_cast<char*>(Data.data()),Data.total_size()*sizeof(double));
    }
    void write(std::ofstream& out) const
    {
        out.write(reinterpret_cast<const char*>(Data.data()),Data.total_size()*sizeof(double));
    }

    double Qdot;                                                                ///< Heat source at the far end of the extension
    Storage1D<double,0> Data;                                                   ///< Data storage array
    Storage1D<double,0> DataTMP;                                                ///< Temporary data storage for iterative heat diffusion solver
    iVector3 Direction;                                                         ///< Selects the extension's direction: 0 -> direction inactive, 1 -> upper boundary extension, -1 -> lower boundary extension

 protected:
 private:

};

class OP_EXPORTS Temperature : public OPObject                                  ///< Storage for the temperature
{
 public:
    Temperature(){};
    Temperature(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates internal storages, initializes internal parameters
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file.
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream.
    void Remesh(const int newNx, const int newNy, const int newNz,
                                        const BoundaryConditions& BC) override; ///< Changes system size while keeping the data
    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC) override;                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondigly.
    void ConsumePlane(const int dx, const int dy, const int dz,
                   const int x, const int y, const int z,
                   const BoundaryConditions& BC);                               ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.

    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets boundary conditions for the temperature by setting the appropriate values in ghost nodes
    bool Write(const Settings& locSettings, const int tStep = -1) const override; ///< Writes the raw temperature into a file
    void WriteH5(H5Interface& H5, const int tStep);
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep = -1) override;     ///< Reads the raw temperature from a file
    bool ReadH5(H5Interface& H5, const int tStep);
    void WriteVTK(Settings& locSettings, const int tStep) const;                ///< Writes temperature in the VTK format into a file
    void WriteGradientVTK(Settings& locSettings, const int tStep) const;        ///< Writes temperature gradient in the VTK format into a file

    void WriteMinMaxAverage(int time_step, double time,
                            std::string filename = "Temperature.txt") const;    ///< Write min, max and average temperature over time into an ASCII file

    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints temperature at a given point (x, y, z) to screen
    void PrintStatistics() const;                                               ///< Prints min, max and average temperature to screen
    void Set(const BoundaryConditions& BC,
             const PhaseField& Phase,
             const double simulation_time,
             const double dt);                                                  ///< Sets temperature according to the cooling rate
    void CalculateMinMaxAvg();                                                  ///< Evaluates min, max and average temperatures in the simulation domain
    void CalculateInterfaceAverage(PhaseField& Phi);                            ///< Averages temperature in all interfaces of all phase pairs

    void SetInitial(const BoundaryConditions& BC);                              ///< Sets initial temperature according to the starting temperature and temperature gradient

    double& operator()(const int x, const int y, const int z)                   ///< Bidirectional access operator
    {
        return Tx(x, y, z);
    };
    double const& operator()(const int x, const int y, const int z) const       ///< Constant reference access operator
    {
        return Tx(x, y, z);
    };
    double at(const double x, const double y, const double z) const             ///< Interpolating access operator
    {
        return Tx.at(x, y, z);
    };

    Temperature& operator=(const Temperature& rhs);                             ///< Copy operator for Temperature class
    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt, const double tStep) override;                  ///< Advects temperature

    double T0;                                                                  ///< Starting temperature
    double TBC0X;                                                               ///< x0 Boundary starting temperature
    double TBCNX;                                                               ///< xN Boundary starting temperature
    double TBC0Y;                                                               ///< y0 Boundary starting temperature
    double TBCNY;                                                               ///< yN Boundary starting temperature
    double TBC0Z;                                                               ///< z0 Boundary starting temperature
    double TBCNZ;                                                               ///< zN Boundary starting temperature

    dVector3 dT_dr;                                                             ///< Initial temperature gradient
    dVector3 r0;                                                                ///< Initial position with T0 temperature (temperature gradient will be applied with respect to this position)

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases

    Storage3D<double, 0> Tx;                                                    ///< Array of temperature
    Storage3D<double, 0> TxDot;                                                 ///< Array of temperature increments needed for advection
    Storage3D<double, 0> TxOld;                                                 ///< Temporary temperature storage

    double Tmin;                                                                ///< Minimum temperature in the simulation domain
    double Tmax;                                                                ///< Maximum temperature in the simulation domain
    double Tavg;                                                                ///< Average temperature in the simulation domain
    Tensor<double, 2> Tiavg;                                                    ///< Average temperature in interfaces between different phases

    //ControlModes ControlMode;                                                   ///< Temperature control mode
    std::vector<TemperatureParameters> ControlParameters;                       ///< Cooling/heating/holding parameters

    LatentHeatModes LatentHeatMode;                                             ///< Indicates if and how latent heat release should be considered
    Storage<double> HeatCapacity;                                               ///< Volumetric heat capacity for all phases
    Tensor<double, 2> LatentHeat;                                               ///< Latent heat values for each phase pair

    bool ExtensionsActive;
    Temperature1Dextension ExtensionX0;                                         ///< 1D temperature field extension at the lower X boundary
    Temperature1Dextension ExtensionXN;                                         ///< 1D temperature field extension at the upper X boundary
    Temperature1Dextension ExtensionY0;                                         ///< 1D temperature field extension at the lower Y boundary
    Temperature1Dextension ExtensionYN;                                         ///< 1D temperature field extension at the upper Y boundary
    Temperature1Dextension ExtensionZ0;                                         ///< 1D temperature field extension at the lower Z boundary
    Temperature1Dextension ExtensionZN;                                         ///< 1D temperature field extension at the upper Z boundary

 protected:
 private:
    bool ReadFromFile;                                                          ///< True if temperature should be read from the input file
    std::vector<std::pair<double,double> > TemperatureProfile;                  ///< Temperature entries from the input file

    void SetInitial1Dextension(Temperature1Dextension& TxExt);                  ///< Sets initial temperature values in 1D extension
    double CalculateLatentHeatEffect(const PhaseField& Phase, const double dt); ///< Calculates temperature change due to release of latent heat

};

} // namespace openphase
#endif
