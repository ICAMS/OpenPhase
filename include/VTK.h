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
 *   Main contributors :   Philipp Engels, Raphael Schiedung, Muhammad Adil Ali,
 *                         Marvin Tegeler, Oleg Shchyglo
 *
 */

#ifndef VTK_H
#define VTK_H

#include "Includes.h"
#include "Settings.h"
#include "ElasticProperties.h"
#include "Tools.h"
#include "../external/WinBase64/base64.h"
#include "../external/miniz.h"

namespace openphase
{

enum class VTKDataTypes                                                         ///< VTK data types
{
    PDScalars,                                                                  ///< Scalar
    PDVectors,                                                                  ///< Vector
    PDTensors                                                                   ///< Matrix/Tensor
};

typedef struct
{
    uint32_t nblocks;
    uint32_t blocksize;
    uint32_t lastblocksize;
    uint32_t compressedsize;
} header_singleblock;

struct OP_EXPORTS VTK                                                           ///< static class to write VTK data in xml format
{

    struct Field_t                                                              ///< Short hand for a field's name and function which will be written to file
    {
        std::string Name;
        std::function<std::any(const int, const int, const int)> Function;

        Field_t(const std::string& inp_Name,
                const std::function<std::any(const int, const int, const int)>& inp_Function) :
            Name(inp_Name), Function(inp_Function) {};
    };

    template <typename T, class function_t>
    static void ForEach (T& buffer,
            const long int Nx, const  long int Ny, const long int Nz,
            function_t Function)
    {
        for(long int k = 0; k < Nz; ++k)
        for(long int j = 0; j < Ny; ++j)
        for(long int i = 0; i < Nx; ++i)
        {
            buffer << Function(i,j,k) << "\n";
        }
        buffer << "</DataArray>\n";
    }

    template <typename T, class function_t>
    static void ForEachCompressed (T& buffer,
            const long int Nx, const  long int Ny, const long int Nz,
            function_t Function)
    {
        std::vector<float> vdata;
        //int it = 0;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            for (size_t n = 0; n < Function(i,j,k).size(); ++n)
            vdata.push_back(Function(i,j,k)[n]);
        }
        uint32_t datasize = vdata.size()*sizeof(float);
        size_t cdatasize;
        float* data = new float[vdata.size()];
        memcpy(data,vdata.data(),datasize);
        char* cpoints = VTK::compress_data(&cdatasize,(const char*)data,datasize);
        header_singleblock b64header = { 1, datasize, datasize, static_cast<unsigned int>(cdatasize)};
        char* b64string = VTK::encode_b64(cpoints, cdatasize);
        char* b64lenstring = VTK::encode_b64((const char*)&b64header, sizeof(header_singleblock));
        buffer << b64lenstring << b64string << std::endl;

        delete[] (data);
        free (cpoints);

        free (b64string);
        free (b64lenstring);

        buffer << "</DataArray>" << std::endl;
    }

    template <typename T>
    static void WriteFieldHeader(
            T& buffer,
            std::string Name,
            std::string Type,
            size_t NComponents = 1,
            std::string format = "ascii")
    {
        buffer << "<DataArray type = \"" << Type
               << "\" Name = \"" << Name
               << "\" NumberOfComponents=\"" << NComponents
               << "\" format=\""<< format << "\">\n";
    }

    template <typename T, class container_t>
    static void WritePointData(
            T& buffer,
            container_t ListOfFields,
            const long int Nx, const  long int Ny, const long int Nz,
            const int precision = 16)
    {
        // Use type with highest dimensions to write the data to file
        bool check = true;
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dMatrix3x3) or
            Field.Function(0,0,0).type()==typeid(dMatrix6x6))
        {
            buffer << "<PointData Tensors= \"TensorData\">\n";
            check = false;
            break;
        }
        if(check)
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dVector3) or
            Field.Function(0,0,0).type()==typeid(dVector6) or
            Field.Function(0,0,0).type()==typeid(vStress) or
            Field.Function(0,0,0).type()==typeid(vStrain))
        {
            buffer << "<PointData Vectors= \"VectorData\">\n";
            check = false;
            break;
        }
        if(check)
        {
            buffer << "<PointData Scalars= \"ScalarData\">\n";
        }
        for (auto Field : ListOfFields)
        {
            if (Field.Function(0,0,0).type() == typeid(int))
            {
                buffer << std::fixed;
                buffer << std::setprecision(0);
                WriteFieldHeader(buffer, Field.Name, "Int32");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<int>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(size_t))
            {
                buffer << std::fixed;
                buffer << std::setprecision(0);
                WriteFieldHeader(buffer, Field.Name, "UInt64");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<size_t>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(double))
            {
                //out << std::scientific; //NOTE: this results in unnecessarily large files
                buffer << std::defaultfloat;
                buffer << std::setprecision(precision);
                WriteFieldHeader(buffer, Field.Name, "Float64");
                ForEach(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<double>(Field.Function(i,j,k));});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 3);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector3>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dVector6>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStrain>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 6);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<vStress>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 9);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix3x3>(Field.Function(i,j,k)).write(precision);});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float64", 36);
                ForEach(buffer,Nx,Ny,Nz,[&Field, precision](int i,int j,int k){return std::any_cast<dMatrix6x6>(Field.Function(i,j,k)).write(precision);});
            }
        }
    }

    template <typename T, class container_t>
    static void WritePointDataCompressed(
            T& buffer,
            container_t ListOfFields,
            const long int Nx, const  long int Ny, const long int Nz,
            const int precision = 16)
    {
        // Use type with highest dimensions to write the data to file
        bool check = true;
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dMatrix3x3) or
            Field.Function(0,0,0).type()==typeid(dMatrix6x6))
        {
            buffer << "<PointData Tensors= \"TensorData\">\n";
            check = false;
            break;
        }
        if(check)
        for (auto Field : ListOfFields)
        if (Field.Function(0,0,0).type()==typeid(dVector3) or
            Field.Function(0,0,0).type()==typeid(dVector6) or
            Field.Function(0,0,0).type()==typeid(vStress) or
            Field.Function(0,0,0).type()==typeid(vStrain))
        {
            buffer << "<PointData Vectors= \"VectorData\">\n";
            check = false;
            break;
        }
        if(check)
        {
            buffer << "<PointData Scalars= \"ScalarData\">\n";
        }

        for (auto Field : ListOfFields)
        {
            if (Field.Function(0,0,0).type() == typeid(int))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 1, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<float> data;
                data.push_back(std::any_cast<int>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(size_t))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 1, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<float> data;
                data.push_back(std::any_cast<size_t>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(double))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 1, "binary");
                buffer << std::setprecision(precision) << std::scientific;
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){
                std::vector<float> data;
                data.push_back(std::any_cast<double>(Field.Function(i,j,k)));
                return data; });
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 3, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<dVector3>(Field.Function(i,j,k)).writeCompressed();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dVector6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 6, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<dVector6>(Field.Function(i,j,k)).writeCompressed();});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStrain))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 6, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<vStrain>(Field.Function(i,j,k)).writeCompressed();});
            }
            else if (Field.Function(0,0,0).type() == typeid(vStress))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 6, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<vStress>(Field.Function(i,j,k)).writeCompressed();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix3x3))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 9, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<dMatrix3x3>(Field.Function(i,j,k)).writeCompressed();});
            }
            else if (Field.Function(0,0,0).type() == typeid(dMatrix6x6))
            {
                WriteFieldHeader(buffer, Field.Name, "Float32", 36, "binary");
                ForEachCompressed(buffer,Nx,Ny,Nz,[&Field](int i,int j,int k){return std::any_cast<dMatrix6x6>(Field.Function(i,j,k)).writeCompressed();});
            }
        }
    }

    /// Writes a list of fields which returns a double into a VTK file
    static void Write(
            const std::string Filename,
            const Settings& locSettings,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);

    static void WriteCompressed(
            const std::string Filename,
            const Settings& locSettings,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);

    static void WriteDistorted(
            const std::string Filename,
            const Settings& locSettings,
            const ElasticProperties& EP,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);
            
     static void WriteDistortedCompressed(
            const std::string Filename,
            const Settings& locSettings,
            const ElasticProperties& EP,
            std::vector<Field_t> ListOfFields,
            const int precision = 16,
            const int resolution = 1);
        
    /// Short for a list of undefined data
    //typedef std::vector<std::function<void(std::stringstream& buffer)>> ListOfData_t;

    ///// Writes VTK data for file WriteVTKData function is provided
    //static void WriteData(
    //        const std::string Filename,
    //        const Settings& locSettings,
    //        ListOfData_t ListOfData);

    template <typename T>
    static void WriteHeader(T& buffer, const int Nx, const int Ny, const int Nz, int resolution = 1)
    {
        buffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
        buffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        buffer << "<StructuredGrid WholeExtent=\""
                 << 0 << " " << resolution * Nx - 1 << " "
                 << 0 << " " << resolution * Ny - 1 << " "
                 << 0 << " " << resolution * Nz - 1 << "\">\n";
    }

    template <typename T>
    static void WriteHeaderCompressed(T& buffer, const int Nx, const int Ny, const int Nz, int resolution = 1)
    {
        buffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << std::endl;
        buffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
        buffer << "<StructuredGrid WholeExtent=\""
                 << 0 << " " << resolution * Nx-1 << " "
                 << 0 << " " << resolution * Ny-1 << " "
                 << 0 << " " << resolution * Nz-1 << "\"> " << std::endl;
    }

    template <typename T>
    static void WriteBeginPointData(T& buffer, const std::vector<VTKDataTypes> PointDataTypes)
    {
        // Method requires vector that contains the types of the given point data.
        // Each type has to be given once.

        // Example vector initialization for an output with scalar and tensor data:
        // ' std::vector<int> DataTypes {PDScalars, PDTensors}; '

        buffer << "<PointData ";

        for(auto it = PointDataTypes.cbegin(); it != PointDataTypes.cend(); ++it)
        {
            switch(*it)
            {
                case VTKDataTypes::PDScalars:
                {
                    buffer << " Scalars= \"ScalarData\"";
                    break;
                }
                case VTKDataTypes::PDVectors:
                {
                    buffer << " Vectors= \"VectorData\"";
                    break;
                }
                case VTKDataTypes::PDTensors:
                {
                    buffer << " Tensors= \"TensorData\"";
                    break;
                }
                default:
                {
                    break;
                }
            }
        }
        buffer << ">\n";
    }

    template <typename T>
    static void WriteEndPointData(T& buffer)
    {
        buffer << "</PointData>\n";
    }

    template <typename T>
    static void WriteCoordinates(T& buffer, const Settings& locSettings, int resolution = 1)
    {
        const long int Nx      = get_Nx      (resolution, locSettings);
        const long int Ny      = get_Ny      (resolution, locSettings);
        const long int Nz      = get_Nz      (resolution, locSettings);
        const long int TotalNx = get_TotalNx (resolution, locSettings);
        const long int TotalNy = get_TotalNy (resolution, locSettings);
        const long int TotalNz = get_TotalNz (resolution, locSettings);
        const long int OffsetX = get_OffsetX (resolution, locSettings);
        const long int OffsetY = get_OffsetY (resolution, locSettings);
        const long int OffsetZ = get_OffsetZ (resolution, locSettings);

        buffer << "<Points>\n";
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        buffer << std::fixed;
        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        double loc_offsetX = 0.0;
        double loc_offsetY = 0.0;
        double loc_offsetZ = 0.0;

        if(resolution == 2)
        {
            a = locSettings.Grid.dNx/*(TotalNx-1)*/ ? 0.5*(TotalNx-1)/TotalNx : 0;
            b = locSettings.Grid.dNy/*(TotalNy-1)*/ ? 0.5*(TotalNy-1)/TotalNy : 0;
            c = locSettings.Grid.dNz/*(TotalNz-1)*/ ? 0.5*(TotalNz-1)/TotalNz : 0;

//            a = (TotalNx-1) ? 0.5 : 0;
//            b = (TotalNy-1) ? 0.5 : 0;
//            c = (TotalNz-1) ? 0.5 : 0;
//
//            loc_offsetX = (TotalNx-1) ? 0.25 : 0;
//            loc_offsetY = (TotalNy-1) ? 0.25 : 0;
//            loc_offsetZ = (TotalNz-1) ? 0.25 : 0;

            buffer << std::setprecision(2);
        }
        else
        {
            buffer << std::setprecision(0);
        }
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            buffer << (i + OffsetX)*a - loc_offsetX << " " << (j + OffsetY)*b - loc_offsetY << " " << (k + OffsetZ)*c - loc_offsetZ << "\n";
        }
        buffer << "</DataArray>\n";
        buffer << "</Points>\n";
    }

    template <typename T>
    static void WriteCoordinatesCompressed(T& buffer, const Settings& locSettings, int resolution = 1)
    {
    	const long int Nx      = get_Nx      (resolution, locSettings);
        const long int Ny      = get_Ny      (resolution, locSettings);
        const long int Nz      = get_Nz      (resolution, locSettings);
        const long int TotalNx = get_TotalNx (resolution, locSettings);
        const long int TotalNy = get_TotalNy (resolution, locSettings);
        const long int TotalNz = get_TotalNz (resolution, locSettings);
        const long int OffsetX = get_OffsetX (resolution, locSettings);
        const long int OffsetY = get_OffsetY (resolution, locSettings);
        const long int OffsetZ = get_OffsetZ (resolution, locSettings);
    
        buffer << "<Points>" << std::endl;
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" header_type=\"UInt32\" format=\"binary\">" << std::endl;

        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        double loc_offsetX = 0.0;
        double loc_offsetY = 0.0;
        double loc_offsetZ = 0.0;

        if(resolution == 2)
        {
            a = locSettings.Grid.dNx/*(TotalNx-1)*/ ? 0.5*(TotalNx-1)/TotalNx : 0;
            b = locSettings.Grid.dNy/*(TotalNy-1)*/ ? 0.5*(TotalNy-1)/TotalNy : 0;
            c = locSettings.Grid.dNz/*(TotalNz-1)*/ ? 0.5*(TotalNz-1)/TotalNz : 0;

//            a = (TotalNx-1) ? 0.5 : 0;
//            b = (TotalNy-1) ? 0.5 : 0;
//            c = (TotalNz-1) ? 0.5 : 0;
//
//            loc_offsetX = (TotalNx-1) ? 0.25 : 0;
//            loc_offsetY = (TotalNy-1) ? 0.25 : 0;
//            loc_offsetZ = (TotalNz-1) ? 0.25 : 0;
        }
        else
        {

        }

        int size = Nx*Ny*Nz*3;
        double* points = new double[size];
        int it = 0;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            points[it] = (i + OffsetX)*a - loc_offsetX; it++;
            points[it] = (j + OffsetY)*b - loc_offsetY; it++;
            points[it] = (k + OffsetZ)*c - loc_offsetZ; it++;
        }
        uint32_t datasize = size*sizeof(double);
        size_t cdatasize;
        char* cpoints = compress_data(&cdatasize,(const char*)points,datasize);
        header_singleblock b64header = { 1, datasize, datasize, static_cast<unsigned int>(cdatasize)};
        char* b64string = encode_b64(cpoints, cdatasize);
        char* b64lenstring = encode_b64((const char*)&b64header, sizeof(header_singleblock));
        buffer << b64lenstring << b64string << std::endl;
        buffer << "</DataArray>" << std::endl;
        buffer << "</Points>" << std::endl;
        delete[] (points);
        free(cpoints);
        free(b64string);
        free(b64lenstring);
    }

    template <typename T>
    static void WriteCoordinatesDistorted(T& buffer, const  ElasticProperties& EP, const Settings& locSettings, int resolution = 1)
    {
        dMatrix3x3 locDefGrad = EP.AverageDeformationGradient();
        //dMatrix3x3 locDefGrad = Tools::AlignBaseAxes(EP.AverageDeformationGradient);

        const long int Nx      = get_Nx      (resolution, locSettings);
        const long int Ny      = get_Ny      (resolution, locSettings);
        const long int Nz      = get_Nz      (resolution, locSettings);
        const long int TotalNx = get_TotalNx (resolution, locSettings);
        const long int TotalNy = get_TotalNy (resolution, locSettings);
        const long int TotalNz = get_TotalNz (resolution, locSettings);
        const long int OffsetX = get_OffsetX (resolution, locSettings);
        const long int OffsetY = get_OffsetY (resolution, locSettings);
        const long int OffsetZ = get_OffsetZ (resolution, locSettings);

        buffer << "<Points>\n";
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        buffer << std::fixed;

        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        if(resolution == 2)
        {
            a = locSettings.Grid.dNx/*(TotalNx-1)*/ ? 0.5*(TotalNx-1)/TotalNx : 0;
            b = locSettings.Grid.dNy/*(TotalNy-1)*/ ? 0.5*(TotalNy-1)/TotalNy : 0;
            c = locSettings.Grid.dNz/*(TotalNz-1)*/ ? 0.5*(TotalNz-1)/TotalNz : 0;
        }

        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            double x = i * a;
            double y = j * b;
            double z = k * c;
            dVector3 coordinates {(-0.5*TotalNx + x + a*OffsetX),
                                  (-0.5*TotalNy + y + b*OffsetY),
                                  (-0.5*TotalNz + z + c*OffsetZ)};
            //int l = i;//int(round(x));
            //int m = j;//int(round(y));
            //int n = k;//int(round(z));
            coordinates = locDefGrad*coordinates + EP.Displacements.at(x,y,z);
            buffer << std::setprecision(4);
            buffer << 0.5*TotalNx + coordinates[0] << " "
                   << 0.5*TotalNy + coordinates[1] << " "
                   << 0.5*TotalNz + coordinates[2] << "\n";
        }
        buffer << "</DataArray>\n";
        buffer << "</Points>\n";
    }
    
    template <typename T>
    static void WriteCoordinatesDistortedCompressed(T& buffer, const  ElasticProperties& EP, const Settings& locSettings, int resolution = 1)
    {
        dMatrix3x3 locDefGrad = EP.AverageDeformationGradient();
        //dMatrix3x3 locDefGrad = Tools::AlignBaseAxes(EP.AverageDeformationGradient);

        const long int Nx      = get_Nx      (resolution, locSettings);
        const long int Ny      = get_Ny      (resolution, locSettings);
        const long int Nz      = get_Nz      (resolution, locSettings);
        const long int TotalNx = get_TotalNx (resolution, locSettings);
        const long int TotalNy = get_TotalNy (resolution, locSettings);
        const long int TotalNz = get_TotalNz (resolution, locSettings);
        const long int OffsetX = get_OffsetX (resolution, locSettings);
        const long int OffsetY = get_OffsetY (resolution, locSettings);
        const long int OffsetZ = get_OffsetZ (resolution, locSettings);

        buffer << "<Points>\n";
        buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" header_type=\"UInt32\" format=\"binary\">" << std::endl;
        buffer << std::fixed;

        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        if(resolution == 2)
        {
            a = locSettings.Grid.dNx/*(TotalNx-1)*/ ? 0.5*(TotalNx-1)/TotalNx : 0;
            b = locSettings.Grid.dNy/*(TotalNy-1)*/ ? 0.5*(TotalNy-1)/TotalNy : 0;
            c = locSettings.Grid.dNz/*(TotalNz-1)*/ ? 0.5*(TotalNz-1)/TotalNz : 0;
        }
		int size = Nx*Ny*Nz*3;	
       	double* points = new double[size];
        int it = 0;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            double x = i * a;
            double y = j * b;
            double z = k * c;
            dVector3 coordinates {(-0.5*TotalNx + x + a*OffsetX),
                                  (-0.5*TotalNy + y + b*OffsetY),
                                  (-0.5*TotalNz + z + c*OffsetZ)};
            coordinates = locDefGrad*coordinates + EP.Displacements.at(x,y,z);
            points[it] = 0.5*TotalNx + coordinates[0]; it++;
            points[it] = 0.5*TotalNy + coordinates[1]; it++;
            points[it] = 0.5*TotalNz + coordinates[2]; it++;
        }
        uint32_t datasize = size*sizeof(double);
        size_t cdatasize;
        char* cpoints = compress_data(&cdatasize,(const char*)points,datasize);
        header_singleblock b64header = { 1, datasize, datasize, static_cast<unsigned int>(cdatasize)};
        char* b64string = encode_b64(cpoints, cdatasize);
        char* b64lenstring = encode_b64((const char*)&b64header, sizeof(header_singleblock));
        buffer << b64lenstring << b64string << std::endl;
        buffer << "</DataArray>\n";
        buffer << "</Points>\n";
    }

    template <typename T>
    static void WriteToFile(T& buffer, const std::string Filename)
    {
        buffer << "</StructuredGrid>\n";
        buffer << "</VTKFile>\n";

        std::ofstream vtk_file(Filename.c_str());
        vtk_file << buffer.rdbuf();
        vtk_file.close();
    }
    template <typename T>
    static void CloseFile(T& buffer)
    {
        buffer << "</StructuredGrid>\n";
        buffer << "</VTKFile>\n";
        buffer.close();
    }

    static char* encode_b64(const char* data, size_t size);
    static char* compress_data(size_t* csize, const char* data, size_t size);

private:

    static long int get_Nx (const int resolution, const Settings& locSettings) 
    {
        if (resolution == 1) return locSettings.Grid.Nx;
        else return (locSettings.Grid.dNx + 1)*locSettings.Grid.Nx;
    }
    static long int get_Ny (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.Ny;
        else return (locSettings.Grid.dNy + 1)*locSettings.Grid.Ny;
    }
    static long int get_Nz (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.Nz;
        else return (locSettings.Grid.dNz + 1)*locSettings.Grid.Nz;
    }
    static long int get_TotalNx (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.TotalNx;
        else return (locSettings.Grid.dNx + 1)*locSettings.Grid.TotalNx;
    }
    static long int get_TotalNy (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.TotalNy;
        else return (locSettings.Grid.dNy + 1)*locSettings.Grid.TotalNy;
    }
    static long int get_TotalNz (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.TotalNz;
        else return (locSettings.Grid.dNz + 1)*locSettings.Grid.TotalNz;
    }
    static long int get_OffsetX (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.OffsetX;
        else return (locSettings.Grid.dNx + 1)*locSettings.Grid.OffsetX;
    }
    static long int get_OffsetY (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.OffsetY;
        else return (locSettings.Grid.dNy + 1)*locSettings.Grid.OffsetY;
    }
    static long int get_OffsetZ (const int resolution, const Settings& locSettings)
    {
        if (resolution == 1) return locSettings.Grid.OffsetZ;
        else return (locSettings.Grid.dNz + 1)*locSettings.Grid.OffsetZ;
    }
};
}// namespace openphase
#endif
