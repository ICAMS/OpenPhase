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
 *   Main contributors :   Oleg Shchyglo; Efim Borukhovich; Matthias Stratmann
 *
 */

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include "Includes.h"
#include "PhaseField.h"
#include "Thermodynamics/Element.h"
#include "Thermodynamics/ThermodynamicPhase.h"
#include "Thermodynamics/PeriodicTable.h"

namespace openphase
{
class PhaseField;
class BoundaryConditions;
class ElasticProperties;
class Velocities;
class PhaseField;
class Settings;
class PeriodicTable;
class Composition1Dextension;

class Composition1Dextension                                                    ///< 1D composition field extension storage
{
 public:
    void Initialize(const size_t size, const size_t phase_index,
                    const size_t n_comp, iVector3 direction)
    {
        if(size <= 0) // Checks for sufficient extension size
        {
            std::stringstream message;
            message << "1D composition extension size = " << size << " is incorrect!\n"
                    << "It should be in the range [1, inf) for correct operation\n";
            ConsoleOutput::WriteExit(message.str(), "Composition1Dextension", "Initialize()");
            OP_Exit(EXIT_FAILURE);
        }

        PhaseIndex = phase_index;
        Ncomp = n_comp;
        Data.Allocate(size, {Ncomp},1);
        DataDot.Allocate(size, {Ncomp},1);
        Direction = direction;
    }
    Composition1Dextension& operator=(const Composition1Dextension& RHS)
    {
        Data    = RHS.Data;
        DataDot = RHS.DataDot;
        Direction = RHS.Direction;
        PhaseIndex = RHS.PhaseIndex;
        Ncomp = RHS.Ncomp;

        return *this;
    }
    void setBC(const BoundaryConditionTypes extBC)                              ///< Sets no-flux boundary condition at the far end of the extension
    {
        switch(extBC)
        {
            case BoundaryConditionTypes::Periodic:
            {
                std::string message = std::string("If 1D composition extension is active, periodic boundary conditions\n")
                                    + std::string("along the extension dimension are not permitted.\n")
                                    + std::string("Permitted boundary conditions: NoFlux, Free and Fixed.");
                ConsoleOutput::WriteWarning(message, "Composition1Dextension", "setBC()");
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
    void setFixedBC(Tensor<double,1> value)                                     ///< Sets fixed composition at the far end of the extension to the specified value
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
    void read(std::ifstream& out)
    {
        out.read(reinterpret_cast<char*>(Data.data()),Data.total_size()*Ncomp*sizeof(double));
    }
    void write(std::ofstream& out)
    {
        out.write(reinterpret_cast<char*>(Data.data()),Data.total_size()*Ncomp*sizeof(double));
    }

    size_t PhaseIndex;                                                          ///< Index of thermodynamic phase for the extension
    size_t Ncomp;                                                               ///< Number of chemical components

    Storage1D<double,1> Data;                                                   ///< Data storage array
    Storage1D<double,1> DataDot;                                                ///< Diffusion increments storage for diffusion solver
    iVector3 Direction;                                                         ///< Selects the extension's direction: 0 -> direction inactive, 1 -> upper boundary extension, -1 -> lower boundary extension

 protected:
 private:

};

class OP_EXPORTS Composition : public OPObject                                  ///< Stores the composition as concentrations or molar densities or ... etc.
{
 public:
    Composition(){};
    Composition(Settings& locSettings,
                const std::string InputFileName = DefaultInputFileName)         ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream

    void Remesh(int newNx, int newNy, int newNz,
                const BoundaryConditions& BC) override;                         ///< Remesh the storage while keeping the data

    void SetInitialMoleFractions(PhaseField& Phi);

    void CalculateTotalMoleFractions(PhaseField& Phase);                        ///< Calculates total mole fractions from the phase mole fractions
    void WriteVTK(Settings& locSettings,
                  const int tStep,
                  const int precision=16) const;                                ///< Writes composition in VTK format (.vts file)

    void WriteDistortedVTK(Settings& locSettings,
                           const ElasticProperties& EP,
                           const int tStep);                                    ///< Writes composition in VTK format on a distorted grid (.vtk file)
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Writes raw composition into a file
    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep = -1) override;     ///< Reads raw composition from a file
    void WriteH5(H5Interface& H5, const int tStep);
    bool ReadH5(H5Interface& H5, const int tStep);

    void SetBoundaryConditions(const BoundaryConditions& BC) override;          ///< Sets the boundary conditions
    void SetLimitsBoundaryConditions(const BoundaryConditions& BC);             ///< Sets the boundary conditions for limits
    void SetTotalLimitsBoundaryConditions(const BoundaryConditions& BC);        ///< Sets the boundary conditions for limits

    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC) override;                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondingly.
    void ConsumePlane(const int dx, const int dy, const int dz,
                      const int x, const int y, const int z,
                      const BoundaryConditions& BC);                            ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.
    void PrintPointStatistics(int x, int y, int z);                             ///< Prints to screen composition at a given point (x, y, z)
    void WriteStatistics(const Settings& locSettings,
                         const int tStep, double dt);                           ///< Writes composition statistics into file, input: time step

    void Advect(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                const BoundaryConditions& BC,
                const double dt, const double tStep) override;                  ///< Advects composition

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components

    GridParameters Grid;                                                        ///< Simulation grid parameters

    double AtomicWeightMixture;                                                 ///< Atomic weight of the Mixture or Composition

    /* template specialization integer number stands for
     * the number of extra dimensions
     * order:
     * Rank = 2: phase, component
     * Rank = 1: phase or component
     */

    Storage3D<double, 1> MassFractionsTotal;                                    ///< Total mass fractions
    Storage3D<double, 1> MassFractionsTotalOld;                                 ///< Total mass fractions for previous timestep
    Tensor<   double, 1> MolecularWeight;                                       ///< Molecular weight of components

    Storage3D<double, 2> MoleFractions;                                         ///< Mole fractions for each phase
    Storage3D<double, 1> MoleFractionsTotal;                                    ///< Total mole fractions

    Storage3D<double, 2> MoleFractionsDotIn;                                    ///< Phase composition incoming increments storage
    Storage3D<double, 2> MoleFractionsDotOut;                                   ///< Phase composition outgoing increments storage
    Storage3D<double, 1> MoleFractionsTotalDot;                                 ///< Total composition increments storage

    Storage3D<   int, 0> Limiting;                                              ///< Storage for the limiting of the composition increments

    Storage3D<double, 2> NormIn;                                                ///< Storage for the normalization of incoming phase composition
    Storage3D<double, 2> NormOut;                                               ///< Storage for the normalization of outgoing phase composition
    Storage3D<double, 1> NormTotal;                                             ///< Storage for the normalization of total composition

    Tensor<double, 2> Initial;                                                  ///< Initial composition of components in all phases
    Tensor<double, 2> MoleFractionsAverage;                                     ///< Average composition of components in all phases
    Tensor<double, 1> MoleFractionsTotalAverage;                                ///< Average composition of components in the simulation domain
    Tensor<double, 3> MoleFractionsInterfaceAverage;                            ///< Average composition of components in the interfaces

    Tensor<   int, 2> Interface;                                                ///< Indicates if there is an interface between two phases
    Tensor<double, 1> Reference;                                                ///< Reference composition

    Tensor<double, 1> TotInitial;                                               ///< Initial amount of each component
    double TotalMolarVolume;                                                    ///< Total molar volume of the system

    Tensor<double, 1> WeightFractionsTotal(int x, int y, int z) const;          ///< Returns local weight fractions of all elements
    void CalculateMoleFractionsTotalAverage(void);
    void CalculateMoleFractionsAverage(PhaseField& Phase);
    void CalculateTotalMolarVolume(void);

    std::vector<Element> Component;                                             ///< List of all chemical components present in the phases
    std::vector<ThermodynamicPhase> Phase;                                      ///< Phases considered in the simulation

    size_t PhaseNumber(std::string name);                                       ///< Returns the index number of the thermodynamic phase
    //TODO:Make the following functions private:

    std::vector<double> TotalAverage;
    std::vector<size_t> CalculateSortedElementMatrix(void);
    std::vector<size_t> SortedElementMatrix;
    bool ElementsAreSorted;

    double MaxMoleFraction(size_t alpha, size_t comp)
    {
        return Phase[alpha].Component[comp].Max;
    };
    double MinMoleFraction(size_t alpha, size_t comp)
    {
        return Phase[alpha].Component[comp].Min;
    };

    PeriodicTable PT;                                                           ///< Periodic table of elements

    bool ExtensionsActive;                                                      ///< True if at least one 1D extension is active, false otherwise
    Composition1Dextension ExtensionX0;                                         ///< 1D composition field extension at the lower X boundary
    Composition1Dextension ExtensionXN;                                         ///< 1D composition field extension at the upper X boundary
    Composition1Dextension ExtensionY0;                                         ///< 1D composition field extension at the lower Y boundary
    Composition1Dextension ExtensionYN;                                         ///< 1D composition field extension at the upper Y boundary
    Composition1Dextension ExtensionZ0;                                         ///< 1D composition field extension at the lower Z boundary
    Composition1Dextension ExtensionZN;                                         ///< 1D composition field extension at the upper Z boundary

 protected:
    int AtStart;
 private:
    void PrintData(void);                                                       ///< Prints chemical data to screen
    void CalculateMolefractionLimits(void);                                     ///< Calculates mole fraction limits for each phase
    double GetMolarVolume(std::string Element);                                 ///< Returns the molar volume of a given element
    double GetAtomicWeight(std::string Element);                                ///< Returns the atomic weight of a given element
    void SetInitiaMoleFractions1Dextension(Composition1Dextension& CxExt);      ///< Sets initial mole fraction values in 1D extension
};

} // namespace openphase

#endif // COMPOSITION_H
