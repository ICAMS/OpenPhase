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

#ifndef SETTINGS_H
#define SETTINGS_H

#include "Includes.h"
#include "BoundaryConditions.h"
#ifndef WIN32
#include "MetaData.h"
#endif
namespace openphase
{
class BoundaryConditions;
class OPObject;
#ifndef WIN32
class MetaData;
#endif
class GridParameters;

class OP_EXPORTS Settings                                                       ///< System settings module. Reads and stores system settings
{
 public:
    std::string thisclassname;                                                  ///< Object's implementation class name
    std::string thisobjectname;                                                 ///< Object's name

    GridParameters Grid;                                                        ///< Simulation grid parameters

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical elements
    std::vector<size_t> Nvariants;                                              ///< Number of crystallographic (symmetry/translation/...) variants

    bool RemeshingAllowed;                                                      ///< Remeshing allowed (Yes or No)
    bool initialized = false;                                                   ///< Indicates if the settings have been initialized

    Settings()                                                                  ///< Default constructor
    {
        Initialize();
    }

    Settings(const std::string InputFileName)                                   ///< Constructor
    {
        Initialize();
        ReadInput(InputFileName);
    }
    void Initialize(std::string ObjectNameSuffix = "");                         ///< Initializes the class
    void ReadInput(const std::string InputFileName = DefaultInputFileName);     ///< Reads input from the specified input file
    void ReadInput(std::stringstream& data);                                    ///< Reads input from the specified input stream
    void ReadJSON(const std::string InputFileName);

    void AddForRemeshing(OPObject& obj);                                        ///< Adds object to be remeshed to the ObjectsToRemesh
    void RemeshAll(int newNx, int newNy, int newNz, int tStep,
                   const BoundaryConditions& BC);                               ///< Calls Remesh() on all objects in ObjectsToRemesh
    void AddForAdvection(OPObject& obj);                                        ///< Adds object to be advected to the ObjectsToAdvect
    void AdvectAll(AdvectionHR& Adv, const Velocities& Vel, PhaseField& Phi,
                   const BoundaryConditions& BC, double dt, int tStep);         ///< Calls Advect() on all objects in ObjectsToAdvect
    void Advect(std::string ObjectName, AdvectionHR& Adv, const Velocities& Vel,
                PhaseField& Phi, const BoundaryConditions& BC,
                double dt, int tStep);                                          ///< Calls Advect() on object(s) with the given name base
    void AddForReading(OPObject& obj);                                          ///< Adds object which can read raw data to the ObjectsToRead
    bool ReadAll(const BoundaryConditions& BC, const int tStep);                ///< Calls Read() on all objects in ObjectsToRead
    bool Read(std::string ObjectName, const BoundaryConditions& BC,
              const int tStep);                                                 ///< Calls Read() on object(s) with the given name base

    GridParameters GridHistoryParameters(int time_step) const;                  ///< Returns grid parameters for a given time step based on grid history records

    std::vector<OPObject*> ObjectsToRemesh;                                     ///< Stores pointers to objects which have remeshing capability
    std::vector<OPObject*> ObjectsToAdvect;                                     ///< Stores pointers to objects which have advection capability
    std::vector<OPObject*> ObjectsToRead;                                       ///< Stores pointers to objects which have raw data reading capability

    std::vector<std::string> PhaseNames;                                        ///< Names of thermodynamic phases
    std::vector<std::string> ElementNames;                                      ///< Names of chemical elements
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases
    std::vector<double> PhaseEquilibriumDensities;                              ///< Equilibrium densities for all phases

    //std::vector<ActivationModes> PhaseActivationModes;                          ///< Phase activation modes. The phase can be enabled or disabled
    //Table<bool> PhaseInteractions;                                              ///< Indicates which phase pairs are allowed to interact
    std::string OutPutDir;                                                         ///< Directory name for the VTK files
    std::string VTKDir;                                                         ///< Directory name for the VTK files
    std::string RawDataDir;                                                     ///< Directory name for the raw data files
    std::string InputRawDataDir;                                                ///< Directory name for the raw data files
    std::string TextDir;                                                        ///< Directory name for the text files

    Settings& operator= (const Settings& rhs);                                  ///< Assignment operator
#ifndef WIN32
    MetaData Meta;
#endif
    std::vector<GridParameters> GridHistory;                                    ///< Stores grid parameters evolution history over the simulation time. Syntax: {{TimeStep1, Nx1, Ny1, Nz1}, ... ,{TimeStepN, NxN, NyN, NzN}}

 protected:
 private:
};

}// namespace openphase

#endif
