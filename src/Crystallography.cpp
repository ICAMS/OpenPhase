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

 *   File created :   2015
 *   Main contributors :   Philipp Engels; Hesham Salama
 *
 */

#include "Crystallography.h"
#include "Settings.h"

namespace openphase
{
using namespace std;

void Crystallography::ReadInput(string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Crystal Symmetry input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inpF(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inpF)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    stringstream inp;
    inp << inpF.rdbuf();
    inpF.close();

    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    string Symmetry = FileInterface::ReadParameterK(inp, moduleLocation, "CrystalSymmetry", false, "CUBIC");

    // Symmetry matrices for crystal lattice. Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambridge University Press 1998
    if (Symmetry == "CUBIC")
    {
        // Load cubic parameters (class 432)
        CrystalSym =  LatticeSystems::Cubic;

        nsym = 24;
        Theta_min =  0 * (Pi / 180);
        Theta_max = 45 * (Pi / 180);
        Phi_min =    0 * (Pi);
        Phi_max = acos(sqrt(1.0 / (2.0 + (pow(tan(Theta_max),2.0)))));

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set(1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(0, 0, 1, 1, 0, 0, 0, 1, 0);
        CrystalSymmetries[2].set(0, 1, 0, 0, 0, 1, 1, 0, 0);
        CrystalSymmetries[3].set(0,-1, 0, 0, 0, 1,-1, 0, 0);
        CrystalSymmetries[4].set(0,-1, 0, 0, 0,-1, 1, 0, 0);
        CrystalSymmetries[5].set(0, 1, 0, 0, 0,-1,-1, 0, 0);

        CrystalSymmetries[6].set( 0, 0,-1, 1, 0, 0, 0,-1, 0);
        CrystalSymmetries[7].set( 0, 0,-1,-1, 0, 0, 0, 1, 0);
        CrystalSymmetries[8].set( 0, 0, 1,-1, 0, 0, 0,-1, 0);
        CrystalSymmetries[9].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[10].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[11].set( 1, 0, 0, 0,-1, 0, 0, 0,-1);

        CrystalSymmetries[12].set( 0, 0,-1, 0,-1, 0,-1, 0, 0);
        CrystalSymmetries[13].set( 0, 0, 1, 0,-1, 0, 1, 0, 0);
        CrystalSymmetries[14].set( 0, 0, 1, 0, 1, 0,-1, 0, 0);
        CrystalSymmetries[15].set( 0, 0,-1, 0, 1, 0, 1, 0, 0);
        CrystalSymmetries[16].set(-1, 0, 0, 0, 0,-1, 0,-1, 0);
        CrystalSymmetries[17].set( 1, 0, 0, 0, 0,-1, 0, 1, 0);

        CrystalSymmetries[18].set( 1, 0, 0, 0, 0, 1, 0,-1, 0);
        CrystalSymmetries[19].set(-1, 0, 0, 0, 0, 1, 0, 1, 0);
        CrystalSymmetries[20].set( 0,-1, 0,-1, 0, 0, 0, 0,-1);
        CrystalSymmetries[21].set( 0, 1, 0,-1, 0, 0, 0, 0,-1);
        CrystalSymmetries[22].set( 0, 1, 0, 1, 0, 0, 0, 0,-1);
        CrystalSymmetries[23].set( 0,-1, 0, 1, 0, 0, 0, 0, 1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Cubic");
    }
    else if (Symmetry == "HEXAGONAL")
    {
        //  Load hexagonal parameters (class 622)
        CrystalSym = LatticeSystems::Hexagonal;

        nsym = 12;
        Theta_min =  0 * (Pi / 180);
        Theta_max = 30 * (Pi / 180);
        Phi_min =    0 * (Pi / 180);
        Phi_max =         Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-0.5, a, 0, -a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[2].set(-0.5,-a, 0,  a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[3].set( 0.5, a, 0, -a,  0.5, 0, 0, 0, 1);
        CrystalSymmetries[4].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0.5,-a, 0,  a,  0.5, 0, 0, 0, 1);

        CrystalSymmetries[6].set(-0.5, -a, 0, -a, 0.5, 0, 0, 0, -1);
        CrystalSymmetries[7].set( 1, 0, 0, 0, -1, 0, 0, 0, -1);
        CrystalSymmetries[8].set(-0.5,  a, 0,  a, 0.5, 0, 0, 0, -1);
        CrystalSymmetries[9].set( 0.5,  a, 0,  a,-0.5, 0, 0, 0, -1);
        CrystalSymmetries[10].set(-1, 0, 0, 0, 1, 0, 0, 0, -1);
        CrystalSymmetries[11].set(0.5, -a, 0, -a,-0.5, 0, 0, 0, -1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Hexagonal");
    }
    else if (Symmetry == "TETRAGONAL")
    {
        //  Load tetragonal parameters (class 422)
        CrystalSym = LatticeSystems::Tetragonal;

        nsym = 8;
        Theta_min =  0 * (Pi / 180);
        Theta_max = 45 * (Pi / 180);
        Phi_min =    0 * (Pi / 180);
        Phi_max =         Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[2].set( 1, 0, 0, 0,-1, 0, 0, 0,-1);
        CrystalSymmetries[3].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[4].set( 0, 1, 0,-1, 0, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0,-1, 0, 1, 0, 0, 0, 0, 1);
        CrystalSymmetries[6].set( 0, 1, 0, 1, 0, 0, 0, 0,-1);
        CrystalSymmetries[7].set( 0,-1, 0,-1, 0, 0, 0, 0,-1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Tetragonal");
    }
    else if (Symmetry == "TRIGONAL" or Symmetry == "RHOMBOHEDRAL")
    {
        //  Load trigonal parameters (class 32)
        CrystalSym = LatticeSystems::Trigonal;

        nsym = 6;
        Theta_min =  0 * (Pi / 180);
        Theta_max = 60 * (Pi / 180);
        Phi_min =    0 * (Pi / 180);
        Phi_max =         Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-0.5, a, 0,-a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[2].set(-0.5,-a, 0, a, -0.5, 0, 0, 0, 1);
        CrystalSymmetries[3].set( 0.5, a, 0, a, -0.5, 0, 0, 0,-1);
        CrystalSymmetries[4].set(-1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[5].set( 0.5,-a, 0,-a, -0.5, 0, 0, 0,-1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Trigonal");
    }
    else if (Symmetry == "ORTHORHOMBIC")
    {
        //  Load orthorhombic parameters (class 22)
        CrystalSym = LatticeSystems::Orthorhombic;

        nsym = 4;
        Theta_min =  0 * (Pi / 180);
        Theta_max = 90 * (Pi / 180);
        Phi_min =    0 * (Pi / 180);
        Phi_max =         Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
        CrystalSymmetries[2].set( 1, 0, 0, 0,-1, 0, 0, 0, 1);
        CrystalSymmetries[3].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Orthorhombic");
    }
    else if (Symmetry == "MONOCLINIC")
    {
        //  Load monoclinic parameters (class 2)
        CrystalSym = LatticeSystems::Monoclinic;

        nsym = 2;
        Theta_min =   0 * (Pi / 180);
        Theta_max = 180 * (Pi / 180);
        Phi_min =     0 * (Pi / 180);
        Phi_max =          Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set( 1, 0, 0, 0, 1, 0, 0, 0, 1);
        CrystalSymmetries[1].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Monoclinic");
    }
    else if (Symmetry == "TRICLINIC")
    {
        //  Load triclinic parameters (class 1)
        CrystalSym = LatticeSystems::Triclinic;

        nsym = 1;
        Theta_min =   0 * (Pi / 180);
        Theta_max = 360 * (Pi / 180);
        Phi_min =     0 * (Pi / 180);
        Phi_max =          Pi / 2;

        CrystalSymmetries.Allocate(nsym);

        CrystalSymmetries[0].set(1, 0, 0, 0, 1, 0, 0, 0, 1);

        ConsoleOutput::WriteSimple("CrystalSymmetry::Triclinic");
    }
    else
    {
        nsym = 0;
        ConsoleOutput::WriteWarning("No or wrong CrystalSymmetry specified!\nThe default \"Cubic\" model is used!", thisclassname, "ReadInput()");
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Crystallography::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Crystallography";
    thisobjectname = thisclassname + ObjectNameSuffix;

    // Symmetry matrices for cubic lattice. Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambridge University Press 1998

    SymmetriesCubic.Allocate(24);

    SymmetriesCubic[0].set(1, 0, 0, 0, 1, 0, 0, 0, 1);
    SymmetriesCubic[1].set(0, 0, 1, 1, 0, 0, 0, 1, 0);
    SymmetriesCubic[2].set(0, 1, 0, 0, 0, 1, 1, 0, 0);
    SymmetriesCubic[3].set(0,-1, 0, 0, 0, 1,-1, 0, 0);
    SymmetriesCubic[4].set(0,-1, 0, 0, 0,-1, 1, 0, 0);
    SymmetriesCubic[5].set(0, 1, 0, 0, 0,-1,-1, 0, 0);

    SymmetriesCubic[6].set( 0, 0,-1, 1, 0, 0, 0,-1, 0);
    SymmetriesCubic[7].set( 0, 0,-1,-1, 0, 0, 0, 1, 0);
    SymmetriesCubic[8].set( 0, 0, 1,-1, 0, 0, 0,-1, 0);
    SymmetriesCubic[9].set(-1, 0, 0, 0, 1, 0, 0, 0,-1);
    SymmetriesCubic[10].set(-1, 0, 0, 0,-1, 0, 0, 0, 1);
    SymmetriesCubic[11].set( 1, 0, 0, 0,-1, 0, 0, 0,-1);

    SymmetriesCubic[12].set( 0, 0,-1, 0,-1, 0,-1, 0, 0);
    SymmetriesCubic[13].set( 0, 0, 1, 0,-1, 0, 1, 0, 0);
    SymmetriesCubic[14].set( 0, 0, 1, 0, 1, 0,-1, 0, 0);
    SymmetriesCubic[15].set( 0, 0,-1, 0, 1, 0, 1, 0, 0);
    SymmetriesCubic[16].set(-1, 0, 0, 0, 0,-1, 0,-1, 0);
    SymmetriesCubic[17].set( 1, 0, 0, 0, 0,-1, 0, 1, 0);

    SymmetriesCubic[18].set( 1, 0, 0, 0, 0, 1, 0,-1, 0);
    SymmetriesCubic[19].set(-1, 0, 0, 0, 0, 1, 0, 1, 0);
    SymmetriesCubic[20].set( 0,-1, 0,-1, 0, 0, 0, 0,-1);
    SymmetriesCubic[21].set( 0, 1, 0,-1, 0, 0, 0, 0,-1);
    SymmetriesCubic[22].set( 0, 1, 0, 1, 0, 0, 0, 0,-1);
    SymmetriesCubic[23].set( 0,-1, 0, 1, 0, 0, 0, 0, 1);

    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}
}// namespace openphase
