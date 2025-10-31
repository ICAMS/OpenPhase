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

 *   File created :   2018
 *   Main contributors :   Oleg Shchyglo; Hesham Salama
 *
 */

#include "Tools/MicrostructureAnalysis.h"
#include "PhaseField.h"
#include "SymmetryVariants.h"
#include "Tools.h"

namespace openphase
{
using namespace std;

void MicrostructureAnalysis::GrainsSurfaceArea(const PhaseField& Phase)
{
    int Nthreads = 1;

    #ifdef _OPENMP
    Nthreads = omp_get_max_threads();
    #endif
    const size_t size = Phase.FieldsProperties.size();
    vector<vector<double>> ThreadsSurfaceVolume(Nthreads);
    for(int t = 0; t < Nthreads; t++)
    {
        ThreadsSurfaceVolume[t].resize(size, 0.0);
    }
    // Calculate surface volumes in each OpenMP chunk
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        int thread = 0;

        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif

        for(auto it = Phase.Fields(i,j,k).cbegin();
                 it != Phase.Fields(i,j,k).cend(); ++it)
        if(it->value != 0.0 and it->value != 1.0)
        {
            ThreadsSurfaceVolume[thread][it->index] += 1;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    vector<double> SurfaceVolume(size);

    // Add surface volumes from different OpenMP chunks
    for(size_t idx = 0; idx < size; idx++)
    {
        for(int t = 0; t < Nthreads; t++)
        {
            SurfaceVolume[idx] += ThreadsSurfaceVolume[t][idx];
        }
    }

    // Update SurfaceVolume across MPI domains
#ifdef MPI_PARALLEL
    size_t loc_size = size;
    size_t max_size = size;

    OP_MPI_Allreduce(&loc_size, &max_size, 1, OP_MPI_UNSIGNED_LONG_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
    if (max_size > loc_size)
    {
        SurfaceVolume.resize(max_size);
    }

    for(size_t idx = 0; idx < SurfaceVolume.size(); idx++)
    {
        double loc_surface_volume = SurfaceVolume[idx];
        OP_MPI_Allreduce(&loc_surface_volume, &(SurfaceVolume[idx]), 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
    }
#endif

    // Calculate surface area from surface volume
    for(size_t idx = 0; idx < SurfaceVolume.size(); idx++)
    {
        SurfaceVolume[idx] /= Phase.Grid.iWidth;
    }
}

double MicrostructureAnalysis::GrainBoundaryStatistics(PhaseField& Phase, std::vector<dVector3> Facets, double DegreeTolerance)
{
    int totalInterface = 0;
    int totalFacet = 0;
    for(size_t n = 0; n < Facets.size(); n++)
    {
        Facets[n].normalize();
    }

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if(Phase.Fields(i,j,k).interface())
        {
            NodeAB<dVector3,dVector3> locNormals = Phase.Normals(i,j,k);
            if(locNormals.size() == 1)
            for(auto alpha = locNormals.begin();
                     alpha != locNormals.end(); alpha++)
            {
                totalInterface++;

                for(size_t n = 0; n < Facets.size(); n++)
                {
                    dVector3 abNorm;
                    abNorm = alpha->value1;
                    abNorm.normalize();
                    double cosTheta = Facets[n]*abNorm * 0.9999999;
                    double angle = acos(cosTheta) * 180.0/Pi;
                    if(fabs(angle) < DegreeTolerance and abNorm.abs() > 1.0 - 0.0000001 and abNorm.abs() < 1.0 + 0.0000001)
                    {
                        totalFacet++;
                    }
                }
            }
        }
    }
    STORAGE_LOOP_END
    if (totalInterface > 0)
    {
        cout << " (" << totalFacet << "/" << totalInterface << ") = ";
        return (double)totalFacet/(double)totalInterface;
    }
    else
    {
        return 0;
    }
}

void MicrostructureAnalysis::WriteEBSDDataQuaternions(const PhaseField& Phase, const int tStep, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        Quaternion  tempQuat;

        size_t pIndex = 0;
        // selecting the quaternion of the majority phase field.
        double value = 0.0;
        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        if(alpha->value > value)
        {
            tempQuat = Phase.FieldsProperties[alpha->index].Orientation;
            value = alpha->value;
            pIndex = Phase.FieldsProperties[alpha->index].Phase;
        }

        tempQuat.normalize();
        outbuffer << index       << "\t"
                  << pIndex + 1  << "\t"
                  << i*scale     << "\t"
                  << j*scale     << "\t"
                  << k*scale     << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = FileInterface::MakeFileName(DefaultRawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const int tStep,
                            const char Axis, const int Position, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t";
    if(Axis != 'X')
    {
        outbuffer << "X\t";
    }
    if(Axis != 'Y')
    {
        outbuffer << "Y\t";
    }
    if(Axis != 'Z')
    {
        outbuffer << "Z\t";
    }

    outbuffer << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if((Axis == 'X' and i == Position) or
           (Axis == 'Y' and j == Position) or
           (Axis == 'Z' and k == Position))
        {
            Quaternion  tempQuat;
            size_t pIndex = 0;

            // selecting the quaternion of the majority phase field.
            double value = 0.0;
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            if(alpha->value > value)
            {
                pIndex = Phase.FieldsProperties[alpha->index].Phase;
                tempQuat = Phase.FieldsProperties[alpha->index].Orientation;
                value = alpha->value;
            }

            tempQuat.normalize();

            outbuffer << index       << "\t"
                      << pIndex + 1  << "\t";
            if(Axis != 'X')
            {
                outbuffer << i*scale << "\t";
            }
            if(Axis != 'Y')
            {
                outbuffer << j*scale << "\t";
            }
            if(Axis != 'Z')
            {
                outbuffer << k*scale << "\t";
            }
            outbuffer << tempQuat[0] << "\t"
                      << tempQuat[1] << "\t"
                      << tempQuat[2] << "\t"
                      << tempQuat[3] << endl;
            index++;
        }
    }
    STORAGE_LOOP_END
    stringstream sliceInd;
    sliceInd << Axis << "-" << Position << "_";
    string FileName = FileInterface::MakeFileName(DefaultRawDataDir, "EBSD_" + sliceInd.str(), tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternions(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        Quaternion  tempQuat;
        int pIndex = 0;

        // selecting the quaternion of the majority phase field.
        double value = 0.0;
        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha != Phase.Fields(i, j, k).cend(); ++alpha)
        if(alpha->value > value)
        {
            pIndex = Phase.FieldsProperties[alpha->index].Phase;
            int variant = Phase.FieldsProperties[alpha->index].Variant;
            Quaternion locQ;
            dMatrix3x3 locRM = SV(pIndex,variant);
            locQ.set(locRM);
            tempQuat = Phase.FieldsProperties[alpha->index].Orientation + locQ;
            value = alpha->value;
        }

        tempQuat.normalize();

        outbuffer << index       << "\t"
                  << pIndex + 1  << "\t"
                  << i*scale     << "\t"
                  << j*scale     << "\t"
                  << k*scale     << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END

    string FileName = FileInterface::MakeFileName(DefaultRawDataDir, "EBSD_", tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

void MicrostructureAnalysis::WriteEBSDDataQuaternionsSlice(const PhaseField& Phase, const SymmetryVariants& SV, const int tStep,
                            const char Axis, const int Position, double scale)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t";
    if(Axis != 'X')
    {
        outbuffer << "X\t";
    }
    if(Axis != 'Y')
    {
        outbuffer << "Y\t";
    }
    if(Axis != 'Z')
    {
        outbuffer << "Z\t";
    }

    outbuffer << "Quat_real\t"
              << "Quat_i\t"
              << "Quat_j\t"
              << "Quat_k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if((Axis == 'X' and i == Position) or
           (Axis == 'Y' and j == Position) or
           (Axis == 'Z' and k == Position))
        {
            Quaternion  tempQuat;
            size_t pIndex = 0;

            // selecting the quaternion of the majority phase field.
            double value = 0.0;
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            if(alpha->value > value)
            {
                pIndex = Phase.FieldsProperties[alpha->index].Phase;
                size_t variant = Phase.FieldsProperties[alpha->index].Variant;
                Quaternion locQ;
                dMatrix3x3 locRM = SV(pIndex,variant);
                locQ.set(locRM);
                tempQuat = Phase.FieldsProperties[alpha->index].Orientation + locQ;
                value = alpha->value;
            }

            tempQuat.normalize();

            outbuffer << index       << "\t"
                      << pIndex + 1  << "\t";
            if(Axis != 'X')
            {
                outbuffer << i*scale << "\t";
            }
            if(Axis != 'Y')
            {
                outbuffer << j*scale << "\t";
            }
            if(Axis != 'Z')
            {
                outbuffer << k*scale << "\t";
            }
            outbuffer << tempQuat[0] << "\t"
                      << tempQuat[1] << "\t"
                      << tempQuat[2] << "\t"
                      << tempQuat[3] << endl;
            index++;
        }
    }
    STORAGE_LOOP_END
    stringstream sliceInd;
    sliceInd << Axis << "-" << Position << "_";
    string FileName = FileInterface::MakeFileName(DefaultRawDataDir, "EBSD_" + sliceInd.str(), tStep, ".dat");

    ofstream ebsd_file(FileName.c_str());
    ebsd_file << outbuffer.rdbuf();
    ebsd_file.close();
}

// Reference sample direction (sd) input as an integer between 1 and 3. Options are [100], [010], and [001], respectively.
void MicrostructureAnalysis::WriteEBSDVTK(PhaseField& Phase, Crystallography& CR, const Settings& locSettings, EulerConvention locConvention, const int sd, const int tStep)
{
    Storage3D<dVector3, 0> tempRGB;
    tempRGB.Allocate(Phase.Grid, 1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        EulerAngles tempAng;
        Quaternion tempQuat;

        size_t locPF = 0;
        double locVal = 0.0;
        for (auto beta = Phase.Fields(i,j,k).cbegin();
                beta != Phase.Fields(i,j,k).cend();  ++beta)
        {
            if(beta->value > locVal)
            {
                locVal = beta->value;
                locPF = beta->index;
            }
        }
        tempQuat =  Phase.FieldsProperties[locPF].Orientation;
        tempAng.set(tempQuat, locConvention, false);

        tempRGB(i,j,k) = IPFColor(sd, tempAng, CR);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    WriteVTKRGB(tempRGB, locSettings, "sRGB"+std::to_string(sd), tStep);
}

dVector3 MicrostructureAnalysis::IPFColor(size_t sd, EulerAngles& tempEuler, Crystallography& CR)
{
    size_t index = 0;
    double Red = 0.0;
    double Green = 0.0;
    double Blue = 0.0;
    double MaxRGB = 0.0;

    double Theta = 0.0;
    double Phi = 0.0;
    double Phi_max2 = 0.0;

    dVector3 ref_dir;
    dVector3 hkl;
    dMatrix3x3 R;
    dMatrix3x3 tempOM;

    dVector3 RGBint;
    dVector3 RGB;

    double phi1 = tempEuler.Q[0];
    double PHI = tempEuler.Q[1];
    double phi2 = tempEuler.Q[2];

    // Assign reference sample direction
    switch (sd)
    {
        case 1: // 100
        {
            ref_dir = {1,0,0};
            break;
        }

        case 2: // 010
        {
            ref_dir = {0,1,0};
            break;
        }

        case 3: // 001
        {
            ref_dir = {0,0,1};
            break;
        }
    };

    // Start of main routine //
    // Assign black RGB values for bad data points (nsym = 0)
    if (CR.nsym == 0)
    {
        RGB.set_to_zero();
    }

    // Assign black RGB value for Euler angles outside of allowable range
    else if (phi1 > 2.0*Pi || PHI > Pi || phi2 > 2.0*Pi)
    {
        RGB.set_to_zero();
    }

    //  Routine for valid set of Euler angles
    else
    {
        // Construct 3X3 orientation matrix from Euler Angles
        tempOM(0,0) = std::cos(phi1) * std::cos(phi2) - std::sin(phi1) * std::cos(PHI) * std::sin(phi2);
        tempOM(0,1) = std::sin(phi1) * std::cos(phi2) + std::cos(phi1) * std::cos(PHI) * std::sin(phi2);
        tempOM(0,2) = std::sin(phi2) * std::sin(PHI);
        tempOM(1,0) = -std::cos(phi1) * std::sin(phi2) - std::sin(phi1) * std::cos(PHI) * std::cos(phi2);
        tempOM(1,1) = -std::sin(phi1) * std::sin(phi2) + std::cos(phi1) * std::cos(PHI) * std::cos(phi2);
        tempOM(1,2) = std::cos(phi2) * std::sin(PHI);
        tempOM(2,0) = std::sin(phi1) * std::sin(PHI);
        tempOM(2,1) = -std::cos(phi1) * std::sin(PHI);
        tempOM(2,2) = std::cos(PHI);

        //Sorting Euler angles into standard stereographic triangle (SST)
        index = 0;
        while(index < CR.nsym)
        {
            // Form orientation matrix
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    R(i,j) = 0.0;
                    for (size_t k = 0; k < 3; ++k)
                    {
                        R(i,j) += CR.CrystalSymmetries[index](i,k) * tempOM(k,j);
                    }
                }
            }

            // Multiple orientation matrix by reference sample direction
            for (size_t i = 0; i < 3; ++i)
            {
                hkl[i] = 0;
                for (size_t j = 0; j < 3; ++j)
                {
                    hkl[i] += R(i,j) * ref_dir[j];
                }
            }

            // Convert to spherical coordinates (ignore "r" variable since r=1)
            Theta = abs(atan2(hkl[1], hkl[0]));
            Phi = acos(abs(hkl[2]));

            // Continue if Theta and Phi values are within the SST
            if (Theta >= CR.Theta_min && Theta < CR.Theta_max && Phi >= CR.Phi_min && Phi < CR.Phi_max)
            {
                break;
            }

            // Increment to next symmetry operator if not in SST
            else
            {
                index++;
            }
        }

        //  Adjust maximum Phi value to ensure it falls within the SST (cubic materials only)
        if(CR.nsym == 24)
        {
            Phi_max2 = acos(sqrt(1.0 / (2.0 + (pow(tan(Theta),2.0)))));
        }
        else
        {
            Phi_max2 = Pi / 2.0;
        }

        // Calculate the RGB color values and make adjustments to maximize colorspace
        Red = abs(1.0 - (Phi / Phi_max2));
        Blue = abs((Theta - CR.Theta_min) / (CR.Theta_max - CR.Theta_min));
        Green = 1.0 - Blue;

        Blue *= (Phi / Phi_max2);
        Green *= (Phi / Phi_max2);

        // Check for negative RGB values before taking square root
        if(Red < 0 || Green < 0 || Blue < 0)
        {
            string msg = "RGB component values must be positive!";
            ConsoleOutput::WriteWarning(msg,"Tools()","Euler2rgb");
        }

        RGB[0] = sqrt(Red);
        RGB[1] = sqrt(Green);
        RGB[2] = sqrt(Blue);

        // Find maximum value of red, green, or blue
        MaxRGB = max({RGB[0], RGB[1], RGB[2]});

        // Normalize position of SST center point
        RGB /= MaxRGB;
    }
    return RGB;
}

void MicrostructureAnalysis::WriteSymmetryVariantsStatistics(size_t pIndex, const PhaseField& Phase, const SymmetryVariants& SV, const int tStep)
{
    string FileName = DefaultTextDir + "SymmetryVariantsStatistics.dat";

    if(tStep == 0)
    {
        stringstream header;
        header << std::left << std::setw(12) << "tStep";
        for(size_t n = 0; n < SV.Nvariants(pIndex); n++)
        {
            stringstream variantN;
            variantN << "V" << n;
            header << std::right << std::setw(6) << variantN.str();
        }
        header << std::right << std::setw(8) << "Total" << endl;

        fstream out_file(FileName.c_str(), ios::out);
        out_file << header.rdbuf();
        out_file.close();
    }

    vector<int> Npfs(SV.Nvariants(pIndex),0);
    vector<double> Volumes(SV.Nvariants(pIndex),0.0);

    for(size_t n = 0; n < Phase.FieldsProperties.size(); n++)
    if(Phase.FieldsProperties[n].Exist and Phase.FieldsProperties[n].Phase == pIndex)
    {
        size_t variant = Phase.FieldsProperties[n].Variant;
        Npfs[variant] += 1.0;
        Volumes[variant] += Phase.FieldsProperties[n].Volume;
    }
    stringstream outbuffer;

    outbuffer << std::left << std::setw(12) << tStep;
    int counter = 0;
    for(size_t n = 0; n < Npfs.size(); n++)
    {
        counter += Npfs[n];
        outbuffer << std::right << std::setw(6) << Npfs[n] /*<< " " << Volumes[n]*/;
    }
    outbuffer << std::right << std::setw(8) << counter << endl;

    fstream out_file(FileName.c_str(),ios::app);
    out_file << outbuffer.rdbuf();
    out_file.close();
}

void MicrostructureAnalysis::WriteGrainsStatistics(const PhaseField& Phase, const int tStep)
{
/// ++++++++++++++++++++++++++++++++++++++++++++
    size_t nPFs = Phase.FieldsProperties.size();

    map<size_t, size_t> AllGrains;
    map <size_t,size_t> pairs;

    pairs.clear();
    AllGrains.clear();

    for (int i = 1; i < Phase.Grid.Nx+1; ++i)
    for (int j = 1; j < Phase.Grid.Ny+1; ++j)
    for (int k = 1; k < Phase.Grid.Nz+1; ++k)
    {
      if (Phase.Fields(i, j, k).interface() && Phase.Fields(i, j, k).size() == 2)
      {
        int idx1 = Phase.Fields(i, j, k).cbegin()->index;
        int idx2 = (Phase.Fields(i, j, k).cbegin()+1)->index;
        pairs[nPFs*idx1 + idx2] = 1;
        pairs[nPFs*idx2 + idx1] = 1;
      }
      for(auto n = Phase.Fields(i, j, k).cbegin(); n < Phase.Fields(i, j, k).cend();++n)
      {
          AllGrains[n->index] += n->value;
      }
    }

    double AveSize = double(Phase.Grid.LocalNumberOfCells())/AllGrains.size();
/// ++++++++++++++++++++++++++++++++++++++++++++ ///

    ofstream Fl1;
    ofstream Fl2;
    ofstream Fl3;
    if(tStep)
    {
        Fl1.open( DefaultTextDir + "SizeAveInfo.dat", ios::app);
        Fl2.open( DefaultTextDir + "SizeDetails.dat", ios::app);
        Fl3.open( DefaultTextDir + "NeighboInfo.dat", ios::app);
    }
    else
    {
        Fl1.open( DefaultTextDir + "SizeAveInfo.dat", ios::out);
        Fl2.open( DefaultTextDir + "SizeDetails.dat", ios::out);
        Fl3.open( DefaultTextDir + "NeighboInfo.dat", ios::out);
    }
    Fl1.precision(10);
    Fl2.precision(10);
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl1 << tStep << " " << AllGrains.size() << " " << AveSize << " " << AveSize*Phase.Grid.CellVolume() << endl;
    Fl1.close();
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl2 << tStep << " " << Phase.FieldsProperties.size() << " ";
    for (size_t i = 0; i < nPFs; i++)
    {
        Fl2 << Phase.FieldsProperties[i].Volume << " ";
    }
    Fl2 << endl;
    Fl2.close();
    /// ++++++++++++++++++++++++++++++++++++++++++++
    Fl3 << tStep << " " << Phase.FieldsProperties.size() << " ";
    vector<size_t> sumpairs(nPFs);
    sumpairs.assign(nPFs, 0);
    for(size_t i = 0; i < nPFs; ++i)
    {
        for(size_t j = 0; j < nPFs; ++j)
        {
            if ((i != j) && (pairs.find(i*nPFs+j) != pairs.end())) sumpairs[i] += pairs.find(i*nPFs+j)->second;
        }
        Fl3 << sumpairs[i] << " ";
    }
    Fl3 << endl;
    Fl3.close();
}

dVector3 MicrostructureAnalysis::FindValuePosition(PhaseField& Phi, size_t index, double value, dVector3 start_position, dVector3 direction, double tolerance)
{
    direction.normalize();
    dVector3 current_position = start_position;
    double current_value = Phi.Fields.at(start_position[0],start_position[1],start_position[2]).get_value(index);
    double step = 1.0;

    double local_difference = 1.0;

    while(current_position[0] >= 0 and current_position[0] < Phi.Grid.Nx and
          current_position[1] >= 0 and current_position[1] < Phi.Grid.Ny and
          current_position[2] >= 0 and current_position[2] < Phi.Grid.Nz and
          local_difference > tolerance)
    {
        double local_value = current_value;
        current_position += direction*step;
        current_value = Phi.Fields.at(current_position[0],current_position[1],current_position[2]).get_value(index);
        if((current_value < value and local_value > value) or
           (current_value > value and local_value < value))
        {
            current_position -= direction*step;
            step *= 0.5;
            local_difference = fabs(local_value - current_value);
            current_value = local_value;
        }
    }
    return current_position;
};

}// namespace openphase
