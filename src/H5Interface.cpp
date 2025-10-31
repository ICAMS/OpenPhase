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

 *   File created :  
 *   Main contributors :   Marvin Tegeler
 *
 */

#include "H5Interface.h"
#include "Settings.h"
#ifdef H5OP
#include "../HighFive/include/highfive/H5Easy.hpp"
#include "../tinyxml2/tinyxml2.h"
#endif

#include <fstream>
#include <string>
#include <vector>
namespace openphase
{
#ifdef H5OP
using namespace tinyxml2;
#endif
#ifndef XMLCheckResult
    #define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { printf("Error: %i\n", a_eResult); return a_eResult; }
#endif

void H5Interface::OpenFile(const std::string InputFileName, const std::string OutputFileName)
{
    H5InputFileName = InputFileName;
    H5OutputFileName = OutputFileName;

    std::ifstream in(H5InputFileName);
    std::ifstream out(H5OutputFileName);
    if(in.good() && !out.good())
    {
        std::string command = "cp " + H5InputFileName + " " + H5OutputFileName;
        //system(command.c_str());
    }
    in.close();
    out.close();
}

void H5Interface::WriteSimulationSettings(const std::string InputFileName)
{
    #ifdef H5OP
    H5Easy::File file(H5OutputFileName, H5Easy::File::OpenOrCreate);
    {
        std::fstream inp(InputFileName.c_str(), std::ios::in | std::ios_base::binary);
        if (!inp)
        {
            Info::WriteExit("File " + InputFileName + " could not be opened", "H5", "WriteSimulationSettings()");
            OP_Exit(EXIT_FAILURE);
        };
        std::vector<std::string> inputlines;
        std::string str;
        while(std::getline(inp, str)){
            std::cout << str << std::endl;
            inputlines.push_back(str);
        }
        if (!file.exist("/SimulationSettings")) {
            // Create the HDF5 group path:
            file.createGroup("/SimulationSettings");
        }
        //H5Easy::dump(file, "/ProjectInput", inputlines);
        H5Easy::DataSet ds = H5Easy::dump(file, "/SimulationSettings/ProjectInput", inputlines, H5Easy::DumpMode::Overwrite);

        size_t el = ds.getElementCount();

        // Error checking:
        if(el == 0) {
            ConsoleOutput::WriteWarning("Zero elements were written to the HDF5 file", thisclassname, "WriteSimulationSettings()");
        }
    }
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteSimulationSettings()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}
void H5Interface::WriteOPID(const std::string InputFileName)
{
    #ifdef H5OP
    H5Easy::File file(H5OutputFileName, H5Easy::File::OpenOrCreate);
    {
        std::fstream inp(InputFileName.c_str(), std::ios::in | std::ios_base::binary);
        if (!inp)
        {
            Info::WriteExit("File " + InputFileName + " could not be opened", "H5", "WriteOPID()");
            OP_Exit(EXIT_FAILURE);
        };
        std::vector<std::string> inputlines;
        std::string str;
        while(std::getline(inp, str)){
            std::cout << str << std::endl;
            inputlines.push_back(str);
        }
        if (!file.exist("/SimulationSettings")) {
            // Create the HDF5 group path:
            file.createGroup("/SimulationSettings");
        }
        //H5Easy::dump(file, "/ProjectInput", inputlines);
        H5Easy::DataSet ds = H5Easy::dump(file, "/SimulationSettings/OPID", inputlines, H5Easy::DumpMode::Overwrite);

        size_t el = ds.getElementCount();

        // Error checking:
        if(el == 0) {
            ConsoleOutput::WriteWarning("Zero elements were written to the HDF5 file", thisclassname, "WriteOPID()");
        }
    }
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteOPID()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

void H5Interface::getProjectInput(std::stringstream& data)
{
    #ifdef H5OP
    H5Easy::File file(H5InputFileName, H5Easy::File::OpenOrCreate);
    if (!file.exist("/SimulationSettings")) {
        Info::WriteExit("/SimulationSettings not found.", "H5", "getProjectInput()");
        OP_Exit(EXIT_FAILURE);
    }
    if (!file.exist("/SimulationSettings/ProjectInput")) {
        Info::WriteExit("/SimulationSettings/ProjectInput not found.", "H5", "getProjectInput()");
        OP_Exit(EXIT_FAILURE);
    }
    std::vector<std::string> inputlines;
    inputlines = H5Easy::load<std::vector<std::string> >(file, "/SimulationSettings/ProjectInput");
    data.clear();
    for (int i = 0; i < inputlines.size(); ++i)
    {
        data << inputlines[i] << std::endl;
    }
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "getProjectInput()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

void H5Interface::getOPID(std::stringstream& data)
{
    #ifdef H5OP
    H5Easy::File file(H5InputFileName, H5Easy::File::OpenOrCreate);
    if (!file.exist("/SimulationSettings")) {
        Info::WriteExit("/SimulationSettings not found.", "H5", "getOPID()");
        OP_Exit(EXIT_FAILURE);
    }
    if (!file.exist("/SimulationSettings/OPID")) {
        Info::WriteExit("/SimulationSettings/OPID not found.", "H5", "getOPID()");
        OP_Exit(EXIT_FAILURE);
    }
    std::vector<std::string> inputlines;
    inputlines = H5Easy::load<std::vector<std::string> >(file, "/SimulationSettings/OPID");
    data.clear();
    for (int i = 0; i < inputlines.size(); ++i)
    {
        data << inputlines[i] << std::endl;
    }
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "getOPID()");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

void H5Interface::WriteVisualization(
        int tStep,
        const Settings& locSettings,
        std::vector<Field_t> ListOfFields,
        const int resolution)
    {
        #ifdef H5OP
        const long int Nx = resolution*locSettings.Nx;
        const long int Ny = resolution*locSettings.Ny;
        const long int Nz = resolution*locSettings.Nz;

        std::stringstream xdmffilename;
        xdmffilename << H5OutputFileName << ".xdmf";

        std::ifstream f(xdmffilename.str().c_str());
        if(!f.good())
        {
            XMLDocument xmlDoc;

            xmlDoc.InsertEndChild( xmlDoc.NewDeclaration() );

            XMLUnknown * pdoctype = xmlDoc.NewUnknown("DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []");
            xmlDoc.InsertEndChild(pdoctype);

            XMLElement * pXdmf = xmlDoc.NewElement("Xdmf");
            pXdmf->SetAttribute("Version", "2.0");
            xmlDoc.InsertEndChild(pXdmf);
            XMLElement * pDomain = xmlDoc.NewElement("Domain");
            pXdmf->InsertEndChild(pDomain);

            XMLElement * pTempGrid = xmlDoc.NewElement("Grid");
            pTempGrid->SetAttribute("Name", "CellTime");
            pTempGrid->SetAttribute("GridType", "Collection");
            pTempGrid->SetAttribute("CollectionType", "Temporal");
            pDomain->InsertEndChild(pTempGrid);

            xmlDoc.SaveFile(xdmffilename.str().c_str());
        }

        XMLDocument xmlDoc;

        XMLError eResult = xmlDoc.LoadFile(xdmffilename.str().c_str());

        XMLNode * pXdmf = xmlDoc.LastChild();

        XMLNode * pDomain = pXdmf->LastChild();

        XMLNode * pTempGrid = pDomain->LastChild();

        XMLElement * pStrucGrid = xmlDoc.NewElement("Grid");
        pStrucGrid->SetAttribute("Name", "Structured Grid");
        pStrucGrid->SetAttribute("GridType", "Uniform");
        pTempGrid->InsertEndChild(pStrucGrid);

        XMLElement * pTime = xmlDoc.NewElement("Time");
        pTime->SetAttribute("Value", tStep);
        pStrucGrid->InsertEndChild(pTime);

        XMLElement * pTop = xmlDoc.NewElement("Topology");
        pTop->SetAttribute("TopologyType", "3DCoRectMesh");
        std::stringstream dims;
        dims << Nx << " " << Ny << " " << Nz;
        pTop->SetAttribute("Dimensions", dims.str().c_str());
        pStrucGrid->InsertEndChild(pTop);

        XMLElement * pGeo = xmlDoc.NewElement("Geometry");
        pGeo->SetAttribute("GeometryType", "ORIGIN_DXDYDZ");
        pStrucGrid->InsertEndChild(pGeo);

        XMLElement * pOrigin = xmlDoc.NewElement("DataItem");
        pOrigin->SetAttribute("Name", "Origin");
        pOrigin->SetAttribute("Dimensions", 3);
        pOrigin->SetAttribute("NumberType", "Double");
        pOrigin->SetAttribute("Precision", 4);
        pOrigin->SetAttribute("Format", "XML");
        pOrigin->SetText("0 0 0");
        pGeo->InsertEndChild(pOrigin);

        XMLElement * pSpacing = xmlDoc.NewElement("DataItem");
        pSpacing->SetAttribute("Name", "Spacing");
        pSpacing->SetAttribute("Dimensions", 3);
        pSpacing->SetAttribute("NumberType", "Double");
        pSpacing->SetAttribute("Precision", 4);
        pSpacing->SetAttribute("Format", "XML");
        pSpacing->SetText("0.1 0.1 0.1");
        pGeo->InsertEndChild(pSpacing);

        std::vector<XMLElement*> Attribute(ListOfFields.size());
        std::vector<XMLElement*> Data(ListOfFields.size());
        int i = 0;
        for (auto Field : ListOfFields)
        {
            Attribute[i] = xmlDoc.NewElement("Attribute");
            Attribute[i]->SetAttribute("Name", Field.Name.c_str());

            if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                Attribute[i]->SetAttribute("AttributeType", "Vector");
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                Attribute[i]->SetAttribute("AttributeType", "Vector");
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                Attribute[i]->SetAttribute("AttributeType", "Tensor6");
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                Attribute[i]->SetAttribute("AttributeType", "Tensor6");
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                Attribute[i]->SetAttribute("AttributeType", "Matrix");
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                Attribute[i]->SetAttribute("AttributeType", "Matrix");
            }
            else
            {
                Attribute[i]->SetAttribute("AttributeType", "Scalar");
            }
            Attribute[i]->SetAttribute("Center", "Node");
            pStrucGrid->InsertEndChild(Attribute[i]);

            Data[i] = xmlDoc.NewElement("DataItem");

            std::stringstream ndims;
            ndims << dims.str();

            if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                ndims << " 3";
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                ndims << " 6";
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                ndims << " 6";
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                ndims << " 6";
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                ndims << " 3 3";
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                ndims << " 6 6";
            }
            Data[i]->SetAttribute("Dimensions", ndims.str().c_str());
            Data[i]->SetAttribute("NumberType", "Double");
            Data[i]->SetAttribute("Precision", 16);
            Data[i]->SetAttribute("Format", "HDF");
            std::stringstream path;
            path << H5OutputFileName << ":/Visualization/" << Field.Name << "/" << std::to_string(tStep);
            Data[i]->SetText(path.str().c_str());

            Attribute[i]->InsertEndChild(Data[i]);
            ++i;
        }

        xmlDoc.SaveFile(xdmffilename.str().c_str());

        WritePointData( tStep,ListOfFields,
             Nx,  Ny, Nz);

        #else
        ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteVisualization()");
        OP_Exit(EXIT_H5_ERROR);
        #endif
    }
}
