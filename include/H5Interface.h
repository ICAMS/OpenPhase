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

#ifndef H5Interface_H
#define H5Interface_H

#ifdef H5OP
#include "../HighFive/include/highfive/H5Easy.hpp"
#endif
#include <string>
#include "Settings.h"

namespace openphase
{

class OP_EXPORTS H5Interface
{
public:
    static constexpr auto thisclassname = "H5Interface";                        ///< Object's implementation class name

    struct Field_t                                                              ///< Short hand for a field's name and function which will be written to file
    {
        std::string Name;
        std::function<std::any(const int, const int, const int)> Function;

        Field_t(const std::string& inp_Name,
                const std::function<std::any(const int, const int, const int)>& inp_Function) :
            Name(inp_Name), Function(inp_Function) {};
    };
    void OpenFile(const std::string InputFileName, const std::string OutputFileName);
    void WriteSimulationSettings(const std::string InputFileName);
    void WriteOPID(const std::string InputFileName);
    void getProjectInput(std::stringstream& data);
    void getOPID(std::stringstream& data);
    void WriteVisualization(        int tStep,
        const Settings& locSettings,
        std::vector<Field_t> ListOfFields,
        const int resolution);

    //template <class T>
    void WriteCheckPoint(int tStep, std::string name, std::vector<double>& data)
    {
        #ifdef H5OP
        H5Easy::File file(H5OutputFileName, H5Easy::File::OpenOrCreate);
        if (!file.exist("/CheckPoints")) {
            // Create the HDF5 group path:
            file.createGroup("/CheckPoints");
        }
        std::stringstream check1;
        check1 << "/CheckPoints/" << name;
        if (!file.exist(check1.str().c_str())) {
            file.createGroup(check1.str().c_str());
        }
        std::stringstream check2;
        check2 << "/CheckPoints/" << name << "/" << tStep;
        H5Easy::DataSet ds = H5Easy::dump(file, check2.str(), data, H5Easy::DumpMode::Overwrite);

        size_t el = ds.getElementCount();
        std::cout << check2.str() << " written " << el << std::endl;
        // Error checking:
        if(el == 0) {
            ConsoleOutput::WriteWarning("Zero elements were written to the HDF5 file",thisclassname,"WriteCheckPoint()");
        }
        #else
        ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"",thisclassname,"WriteCheckPoint()");
        OP_Exit(EXIT_H5_ERROR);
        #endif
    }

    //template <class T>
    void ReadCheckPoint(int tStep, std::string name, std::vector<double>& data)
    {
        #ifdef H5OP
        H5Easy::File file(H5InputFileName, H5Easy::File::OpenOrCreate);
        if (!file.exist("/CheckPoints")) {
            Info::WriteExit("/CheckPoints not found.", "H5", "ReadCheckPoint()");
            OP_Exit(EXIT_FAILURE);
        }
        std::stringstream check1;
        check1 << "/CheckPoints/" << name;
        if (!file.exist(check1.str().c_str())) {
            Info::WriteExit(check1.str()+" not found.", "H5", "ReadCheckPoint()");
            OP_Exit(EXIT_FAILURE);
        }
        std::stringstream check2;
        check2 << "/CheckPoints/" << name << "/" << tStep;
        if (!file.exist(check2.str().c_str())) {
            Info::WriteExit(check2.str()+" not found.", "H5", "ReadCheckPoint()");
            OP_Exit(EXIT_FAILURE);
        }
        data = H5Easy::load<std::vector<double> >(file, check2.str());
        #else
        ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"",thisclassname,"ReadCheckPoint()");
        OP_Exit(EXIT_H5_ERROR);
        #endif
    }

    template <class Funktion_t>
    void ForEach (int tStep, std::string Name,
            const long int Nx, const  long int Ny, const long int Nz,
            Funktion_t Function)
    {
        #ifdef H5OP
        H5Easy::File file(H5OutputFileName, H5Easy::File::OpenOrCreate);
        std::vector<float> vdata;
        int it = 0;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            for (int n = 0; n < Function(i,j,k).size(); ++n)
            vdata.push_back(Function(i,j,k)[n]);
        }
        if (!file.exist("/Visualization")) {
            // Create the HDF5 group path:
            file.createGroup("/Visualization");
        }
        if (!file.exist("/Visualization/"+Name)) {
            // Create the HDF5 group path:
            file.createGroup("/Visualization/"+Name);
        }
        //H5Easy::dump(file, "/PhaseField", data);
        H5Easy::DataSet ds = H5Easy::dump(file,"/Visualization/"+Name+"/"+std::to_string(tStep), vdata, H5Easy::DumpMode::Overwrite);

        size_t el = ds.getElementCount();

        // Error checking:
        if(el == 0) {
            ConsoleOutput::WriteWarning("Zero elements were written to the HDF5 file",thisclassname,"ForEach()");
        }
        #else
        ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"",thisclassname,"ForEach()");
        OP_Exit(EXIT_H5_ERROR);
        #endif
    }

    template <class container_t>
    void WritePointData(int tStep,
            container_t ListOfFields,
            const long int Nx, const  long int Ny, const long int Nz,
            const int precision = 16)
    {
        #ifdef H5OP
        for (auto Field : ListOfFields)
        {
            if (Field.Function(0,0,0).type() == typeid(int))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<double> data;
                data.push_back(std::any_cast<double>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(size_t))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<double> data;
                data.push_back(std::any_cast<size_t>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(double))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<double> data;
                data.push_back(std::any_cast<double>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector3>(Field.Function(i,j,k)).writeBinary();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector6>(Field.Function(i,j,k)).writeBinary();});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStrain>(Field.Function(i,j,k)).writeBinary();});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStress>(Field.Function(i,j,k)).writeBinary();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix3x3>(Field.Function(i,j,k)).writeBinary();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                ForEach(tStep,Field.Name,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix6x6>(Field.Function(i,j,k)).writeBinary();});
            }
        }
        #else
        ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"",thisclassname,"WritePointData()");
        OP_Exit(EXIT_H5_ERROR);
        #endif
    }

std::string H5InputFileName;
std::string H5OutputFileName;
};

}

#endif
