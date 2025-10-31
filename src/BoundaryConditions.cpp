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

 *   File created :   2014
 *   Main contributors :   Oleg Shchyglo; Raphael Schiedung
 *
 */

#include "BoundaryConditions.h"
#include "Settings.h"

using  namespace std;
namespace openphase
{

BoundaryConditionTypes BoundaryConditions::TranslateBoundaryConditions(string Key)
{
    BoundaryConditionTypes myBC;
    bool check = false;
    if(Key == "PERIODIC")
    {
        myBC = BoundaryConditionTypes::Periodic;
        check = true;
    }
    if(Key == "NOFLUX")
    {
        myBC = BoundaryConditionTypes::NoFlux;
        check = true;
    }
    if(Key == "FREE")
    {
        myBC = BoundaryConditionTypes::Free;
        check = true;
    }
    if(Key == "FIXED")
    {
        myBC = BoundaryConditionTypes::Fixed;
        check = true;
    }
    if(Key == "MIRROR")
    {
        myBC = BoundaryConditionTypes::Mirror;
        check = true;
    }

    if (!check)
    {
        ConsoleOutput::WriteExit(string("Illegal Boundary Conditions Input. ") +
                                 string("Legal boundary conditions are ") +
                                 string("Periodic, NoFlux, Free, Fixed or Mirror."),
                                 thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }
    return myBC;
}

void BoundaryConditions::Initialize(Settings& Settings, std::string ObjectNameSuffix)
{
    thisclassname = "BoundaryConditions";
    thisobjectname = thisclassname + ObjectNameSuffix;

    BC0X = BoundaryConditionTypes::Periodic;
    BCNX = BoundaryConditionTypes::Periodic;
    BC0Y = BoundaryConditionTypes::Periodic;
    BCNY = BoundaryConditionTypes::Periodic;
    BC0Z = BoundaryConditionTypes::Periodic;
    BCNZ = BoundaryConditionTypes::Periodic;

#ifdef MPI_PARALLEL
    MPIperiodicX = false;
    MPIperiodicY = false;
    MPIperiodicZ = false;
#endif

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void BoundaryConditions::ReadInput(const string InputFileName)
{
    ConsoleOutput::WriteLineInsert("BoundaryConditions input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    std::string filetype = FileInterface::getFileExtension(InputFileName); 
	if (filetype == "opi")
	{
		fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

		if (!inp)
		{
		    ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
		    OP_Exit(EXIT_FAILURE);
		};

		std::stringstream data;
		data << inp.rdbuf();
		ReadInput(data);

		inp.close();
	}
	else
	if (filetype == "json")
	{
		ReadJSON(InputFileName);
	}
	else
	{
		std::cerr << "Filetype " << filetype << " not recognized. Filetype must be opi or json." << std::endl;
		OP_Exit(1);
	}
}

void BoundaryConditions::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    string BC0Xstring = FileInterface::ReadParameterK(inp, moduleLocation, string("BC0X"));
    string BCNXstring = FileInterface::ReadParameterK(inp, moduleLocation, string("BCNX"));
    string BC0Ystring = FileInterface::ReadParameterK(inp, moduleLocation, string("BC0Y"));
    string BCNYstring = FileInterface::ReadParameterK(inp, moduleLocation, string("BCNY"));
    string BC0Zstring = FileInterface::ReadParameterK(inp, moduleLocation, string("BC0Z"));
    string BCNZstring = FileInterface::ReadParameterK(inp, moduleLocation, string("BCNZ"));

    BC0X = TranslateBoundaryConditions(BC0Xstring);
    BCNX = TranslateBoundaryConditions(BCNXstring);
    BC0Y = TranslateBoundaryConditions(BC0Ystring);
    BCNY = TranslateBoundaryConditions(BCNYstring);
    BC0Z = TranslateBoundaryConditions(BC0Zstring);
    BCNZ = TranslateBoundaryConditions(BCNZstring);

#ifdef MPI_PARALLEL
    if(MPI_3D_DECOMPOSITION)
    {
        Setup_MPIX();
        Setup_MPIY();
        Setup_MPIZ();
    }
    else
    {
        Setup_MPI();
    }
#endif
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void BoundaryConditions::ReadJSON(const string InputFileName)
{
    std::ifstream f(InputFileName);
	json data = json::parse(f);
	if (data.contains(thisclassname))
	{
		json bc = data[thisclassname];

		string BC0Xstring = FileInterface::ReadParameter<std::string>(bc, {"BC0X"});
		string BCNXstring = FileInterface::ReadParameter<std::string>(bc, {"BCNX"});
		string BC0Ystring = FileInterface::ReadParameter<std::string>(bc, {"BC0Y"});
		string BCNYstring = FileInterface::ReadParameter<std::string>(bc, {"BCNY"});
		string BC0Zstring = FileInterface::ReadParameter<std::string>(bc, {"BC0Z"});
		string BCNZstring = FileInterface::ReadParameter<std::string>(bc, {"BCNZ"});
		
		std::transform(BC0Xstring.begin(), BC0Xstring.end(), BC0Xstring.begin(), ::toupper);
		std::transform(BCNXstring.begin(), BCNXstring.end(), BCNXstring.begin(), ::toupper);
		std::transform(BC0Ystring.begin(), BC0Ystring.end(), BC0Ystring.begin(), ::toupper);
		std::transform(BCNYstring.begin(), BCNYstring.end(), BCNYstring.begin(), ::toupper);
		std::transform(BC0Zstring.begin(), BC0Zstring.end(), BC0Zstring.begin(), ::toupper);
		std::transform(BCNZstring.begin(), BCNZstring.end(), BCNZstring.begin(), ::toupper);

		BC0X = TranslateBoundaryConditions(BC0Xstring);
		BCNX = TranslateBoundaryConditions(BCNXstring);
		BC0Y = TranslateBoundaryConditions(BC0Ystring);
		BCNY = TranslateBoundaryConditions(BCNYstring);
		BC0Z = TranslateBoundaryConditions(BC0Zstring);
		BCNZ = TranslateBoundaryConditions(BCNZstring);
	}
#ifdef MPI_PARALLEL
    if(MPI_3D_DECOMPOSITION)
    {
        Setup_MPIX();
        Setup_MPIY();
        Setup_MPIZ();
    }
    else
    {
        Setup_MPI();
    }
#endif
    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

#ifdef MPI_PARALLEL
bool BoundaryConditions::Setup_MPI()
{
    if(BC0X == BoundaryConditionTypes::Periodic or
       BCNX == BoundaryConditionTypes::Periodic)
    {
        MPIperiodicX = true;
        BC0X = BoundaryConditionTypes::MPIcomm;
        BCNX = BoundaryConditionTypes::MPIcomm;
        return true;
    }
    if (MPI_RANK > 0) BC0X = BoundaryConditionTypes::MPIcomm;
    if (MPI_RANK < MPI_SIZE - 1) BCNX = BoundaryConditionTypes::MPIcomm;
    MPIperiodicX = false;
    return false;
}

bool BoundaryConditions::Setup_MPIX()
{
    if(MPI_CART_SIZE[0] > 1)
    {
        if(BC0X == BoundaryConditionTypes::Periodic or
           BCNX == BoundaryConditionTypes::Periodic)
        {
            MPIperiodicX = true;
            BC0X = BoundaryConditionTypes::MPIcomm;
            BCNX = BoundaryConditionTypes::MPIcomm;
            return true;
        }

        if (MPI_CART_RANK[0] > 0) BC0X = BoundaryConditionTypes::MPIcomm;
        if (MPI_CART_RANK[0] < MPI_CART_SIZE[0] - 1) BCNX = BoundaryConditionTypes::MPIcomm;
        MPIperiodicX = false;
    }
    return false;
}

bool BoundaryConditions::Setup_MPIY()
{
    if(MPI_CART_SIZE[1] > 1)
    {
        if(BC0Y == BoundaryConditionTypes::Periodic or
           BCNY == BoundaryConditionTypes::Periodic)
        {
            MPIperiodicY = true;
            BC0Y = BoundaryConditionTypes::MPIcomm;
            BCNY = BoundaryConditionTypes::MPIcomm;
            return true;
        }
    
        if (MPI_CART_RANK[1] > 0) BC0Y = BoundaryConditionTypes::MPIcomm;
        if (MPI_CART_RANK[1] < MPI_CART_SIZE[1] - 1) BCNY = BoundaryConditionTypes::MPIcomm;
        MPIperiodicY = false;
    }
    return false;
}

bool BoundaryConditions::Setup_MPIZ()
{
    if(MPI_CART_SIZE[2] > 1)
    {
        if(BC0Z == BoundaryConditionTypes::Periodic or
        BCNZ == BoundaryConditionTypes::Periodic)
        {
            MPIperiodicZ = true;
            BC0Z = BoundaryConditionTypes::MPIcomm;
            BCNZ = BoundaryConditionTypes::MPIcomm;
            return true;
        }
     
        if (MPI_CART_RANK[2] > 0) BC0Z = BoundaryConditionTypes::MPIcomm;
        if (MPI_CART_RANK[2] < MPI_CART_SIZE[2] - 1) BCNZ = BoundaryConditionTypes::MPIcomm;
        MPIperiodicZ = false;
    }
    return false;
}
#endif

}// namespace openphase

