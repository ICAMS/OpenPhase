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

#include "VTK.h"
#include "Settings.h"
#include "../external/WinBase64/base64.h"
#include "../external/miniz.h"
namespace openphase
{

using namespace std;

vector<string> voigtcon {"xx", "yy", "zz", "yz", "xz", "xy"};
vector<string> matrixcon  {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};

// Returns newly allocated memory filled with b64 encoded data. Remember to free this returned memory!
char* VTK::encode_b64(const char* data, size_t size)
{
    char* b64string;
//#ifdef _WIN32
    string result = base64_encode((unsigned char*)data, static_cast<unsigned int>(size));
    b64string = (char*)malloc(sizeof(char) * (result.size()+1));
    strcpy(b64string, result.c_str());
//#else
//  size_t expected_len = ceil(size*1.34)+3; // no idea why +3...
//  b64string = (char*)malloc(expected_len);
//  size_t stringlen;
//  base64_encode(data,size,b64string,&stringlen,0);
//  //printf("Allocated: %zu, got %zu for %zu\n",expected_len,stringlen,size);
//// remember we have to add the terminating \0
//  if (stringlen >= expected_len) printf("WARNING BUFFER TOO SMALL (USED %zu/%zu)\n", stringlen + 1, expected_len);
//  b64string[stringlen] = '\0';
//#endif
    return b64string;
}

// Returns memory and compresses into that. Remember to free this memory!
char* VTK::compress_data(size_t* csize, const char* data, size_t size)
{
    size_t comp_size = compressBound(size);
    //printf("Maximum Compressed Size: %zu (%zu)\n",comp_size,size);
    char* comp = (char*)malloc(comp_size);
    if(!comp) printf("OUT OF MEMORY!\n");
    int stat = compress2((unsigned char*)comp, (unsigned long*)&comp_size,(const unsigned char*)data,size,10);
    //printf("Compressed to: %zu (%zu)\n",comp_size,size);
    if(stat!=Z_OK) printf("COMPRESS FAILED!\n");
    *csize = comp_size;
    return comp;
}

void VTK::Write(
    const std::string Filename,
    const Settings& locSettings,
    std::vector<Field_t> ListOfFields,
    const int precision,
    const int resolution)
{
}

void VTK::WriteCompressed(
        const std::string Filename,
        const Settings& locSettings,
        std::vector<Field_t> ListOfFields,
        const int precision, const int resolution)
{
    const long int Nx = get_Nx(resolution, locSettings);
    const long int Ny = get_Ny(resolution, locSettings);
    const long int Nz = get_Nz(resolution, locSettings);

#ifndef MPI_PARALLEL
    ofstream vtk_file(Filename.c_str());
    VTK::WriteHeaderCompressed(vtk_file, Nx, Ny, Nz);
    {
        VTK::WritePointDataCompressed(vtk_file,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(vtk_file);
    VTK::WriteCoordinatesCompressed(vtk_file, locSettings, resolution);
    VTK::CloseFile(vtk_file);
#else
    std::stringstream buffer;
    std::stringstream hbuffer;
    std::stringstream tbuffer;

    const long int TotalNx = get_TotalNx(resolution, locSettings);
    const long int TotalNy = get_TotalNy(resolution, locSettings);
    const long int TotalNz = get_TotalNz(resolution, locSettings);
    const long int OffsetX = get_OffsetX(resolution, locSettings);
    const long int OffsetY = get_OffsetY(resolution, locSettings);
    const long int OffsetZ = get_OffsetZ(resolution, locSettings);

    buffer << "<Piece Extent=\""
           << OffsetX << " " << Nx + OffsetX-1 << " "
           << OffsetY << " " << Ny + OffsetY-1 << " "
           << OffsetZ << " " << Nz + OffsetZ-1 << "\">\n";
    {
        WritePointData(buffer,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, locSettings, resolution);
    buffer << "</Piece>\n";

    hbuffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    hbuffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    hbuffer << "<StructuredGrid WholeExtent=\""
            << 0 << " " << TotalNx-1 << " "
            << 0 << " " << TotalNy-1 << " "
            << 0 << " " << TotalNz-1 << "\"> \n";
    tbuffer << "</StructuredGrid> \n";
    tbuffer << "</VTKFile> \n";

    op_mpi_write_vtk(Filename, buffer, hbuffer, tbuffer);
 #endif
}

void VTK::WriteDistorted(
    const std::string Filename,
    const Settings& locSettings,
    const ElasticProperties& EP,
    std::vector<Field_t> ListOfFields,
    const int precision,
    const int resolution)
{
}

void VTK::WriteDistortedCompressed(
    const std::string Filename,
    const Settings& locSettings,
    const ElasticProperties& EP,
    std::vector<Field_t> ListOfFields,
    const int precision,
    const int resolution)
{
    const long int Nx = get_Nx(resolution, locSettings);
    const long int Ny = get_Ny(resolution, locSettings);
    const long int Nz = get_Nz(resolution, locSettings);

#ifndef MPI_PARALLEL
    ofstream vtk_file(Filename.c_str());
    VTK::WriteHeaderCompressed(vtk_file, Nx, Ny, Nz);
    {
        WritePointDataCompressed(vtk_file,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(vtk_file);
    VTK::WriteCoordinatesDistortedCompressed(vtk_file,EP,locSettings,resolution);
    VTK::CloseFile(vtk_file);
#else
    std::stringstream buffer;
    std::stringstream hbuffer;
    std::stringstream tbuffer;
    
    const long int TotalNx = get_TotalNx(resolution, locSettings);
    const long int TotalNy = get_TotalNy(resolution, locSettings);
    const long int TotalNz = get_TotalNz(resolution, locSettings);
    const long int OffsetX = get_OffsetX(resolution, locSettings);
    const long int OffsetY = get_OffsetY(resolution, locSettings);
    const long int OffsetZ = get_OffsetZ(resolution, locSettings);

    buffer << "<Piece Extent=\""
           << OffsetX << " " << Nx-1 + OffsetX << " "
           << OffsetY << " " << Ny-1 + OffsetY << " "
           << OffsetZ << " " << Nz-1 + OffsetZ << "\">\n";
    {
        WritePointData(buffer,ListOfFields,Nx,Ny,Nz,precision);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinatesDistorted(buffer,EP,locSettings,resolution);
    buffer << "</Piece>\n";

    hbuffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    hbuffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    hbuffer << "<StructuredGrid WholeExtent=\""
            << 0 << " " << TotalNx-1 << " "
            << 0 << " " << TotalNy-1 << " "
            << 0 << " " << TotalNz-1 << "\"> \n";
    tbuffer << "</StructuredGrid> \n";
    tbuffer << "</VTKFile> \n";

    op_mpi_write_vtk(Filename, buffer, hbuffer, tbuffer);
#endif
}

}// namespace openphase
