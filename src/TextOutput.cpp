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
 *   Main contributors :   Matthias Stratmann; Johannes Goerler
 *
 */

#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "ElasticProperties.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Temperature.h"
#include "Nucleation.h"
#include "TextOutput.h"
#include "Composition.h"

namespace openphase
{
using namespace std;

void TextOutput::WriteValue(double Value, string filename, double time)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            ofstream file(filename, ios::out);
            file << "time"
                 << separator << "Value\n";
            file.close();
        }

        ofstream file(filename, ios::app);

        file << time << separator << Value << "\n";
        file.close();
    }
}

void TextOutput::WriteMultipleValues(std::vector<string> Names, std::vector<double> value, string filename, double time)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            ofstream file(filename, ios::out);
            file << "time";
            for (auto& N : Names)
            {
                file << separator << N;
            }
            file << "\n";
            file.close();
        }

        ofstream file(filename, ios::app);
        file << time << scientific;
        for (auto& V : value)
        {
          file << " , " << V;
        }
        file << "\n";
        file.close();
    }
}

void TextOutput::WritedVectorNValues(std::vector<string> Names, dVectorN value, string filename, double time)
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
#endif
    {
        if (!std::filesystem::exists(filename))
        {
            ofstream file(filename, ios::out);
            file << "time";
            for (auto& N : Names)
            {
                file << separator << N;
            }
            file << "\n";
            file.close();
        }

        ofstream file(filename, ios::app);
        file << time << scientific;
        for (size_t i = 0; i < value.size(); i++ )
        {
          file << " , " << value[i];
        }
        file << "\n";
        file.close();
    }
}

void TextOutput::AveragePlasticStrain(ElasticProperties& EP, string filename,
                                double timeOrStrain, size_t precision)
{
    /** This function will create tabulated data on the average strain in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    vStrain avgStrain = EP.AveragePlasticStrain();

    static vStrain OldStrain;
    static double OldTimeStep = 0.0;
    vStrain Rate;
    if (timeOrStrain)
    {
        Rate = (avgStrain - OldStrain)/(timeOrStrain - OldTimeStep);
    }
    OldStrain = avgStrain;
    OldTimeStep = timeOrStrain;

    double averagePEEQ = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.DeformationGradientsTotal,0,reduction(+:averagePEEQ))
    {
        averagePEEQ += EP.PlasticStrains(i,j,k).PEEQ();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    double  tempaverage = averagePEEQ ;
    OP_MPI_Reduce(&tempaverage,&averagePEEQ ,1, OP_MPI_DOUBLE, OP_MPI_SUM,0, OP_MPI_COMM_WORLD);
    if(MPI_RANK == 0)
#endif
    {
        averagePEEQ /= double(EP.Grid.TotalNumberOfCells());
    if (!std::filesystem::exists(filename))
    {
        ofstream file(filename, ios::out);
        file << "time"
             << separator << "Epsilon_0"
             << separator << "Epsilon_1"
             << separator << "Epsilon_2"
             << separator << "Epsilon_3"
             << separator << "Epsilon_4"
             << separator << "Epsilon_5"
             << separator << "PEEQ"
             << separator << "Rate_0"
             << separator << "Rate_1"
             << separator << "Rate_2"
             << separator << "Rate_3"
             << separator << "Rate_4"
             << separator << "Rate_5\n";
        file.close();
    }

    ofstream file(filename, ios::app);
    file.precision(precision);
    file << timeOrStrain << scientific
         << separator << avgStrain[0]
         << separator << avgStrain[1]
         << separator << avgStrain[2]
         << separator << avgStrain[3]
         << separator << avgStrain[4]
         << separator << avgStrain[5]
         << separator << averagePEEQ
         << separator << Rate[0]
         << separator << Rate[1]
         << separator << Rate[2]
         << separator << Rate[3]
         << separator << Rate[4]
         << separator << Rate[5] << "\n";
    file.close();
    }
}
void TextOutput::maxElasticRotation(ElasticProperties& EP, PhaseField& Phase,Orientations& OR, int tStep, string Filename)
{
    double maxRotation = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields,0, reduction(max:maxRotation))
    {
        double locangle = 0;
        dVector3 locaxis;
        Tools::getAxisAngle(EP.DeformationGradientsTotal(i, j, k), locaxis, locangle);
        locangle *= 180.0 / Pi;
        maxRotation = std::max(maxRotation, locangle);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#ifdef MPI_PARALLEL
    double tempmaxRotation = maxRotation;
    OP_MPI_Allreduce(&tempmaxRotation,&maxRotation, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
#endif
    WriteValue(maxRotation, Filename, tStep);
}

void TextOutput::LineConcentration(Composition& Cx, PhaseField& Phi,
                                   string filename, double timestep,
                                   string type, string axis,
                                   int x, int y, int z)
{
    /** This function will create a separate file each time it is called, with
    a different file name according to the current time step. In this file, the
    tabulated composition data will be written down along a single line along a
    given axis with type = "X", "Y" or "Z". The line can be positioned with a
    point on it given with "x", "y" and "z". The output can be weight fraction
    with type = "WF", weight percent "WP", mole fraction "MF" and mole percent
    with "MP".*/

    int mode = -1;
    int dir = -1;
    int tablength = 16;
    std::stringstream ss;

    if(type == "WF")
    {
        mode = 0;
    }
    else if(type == "WP")
    {
        mode = 1;
    }
    else if(type == "MF")
    {
        mode = 2;
    }
    else if(type == "MP")
    {
        mode = 3;
    }
    else
    {
        cout << "Undefined type " << type
             << " in TextOutput::LineConcentration()! Chose \"WF\", \"WP\", "
             << "\"MF\" or \"MP\". Now exiting!\n";
        OP_Exit(EXIT_FAILURE);
    }

    if(axis == "X")
    {
        dir = 0;
    }
    else if(axis == "Y")
    {
        dir = 1;
    }
    else if(axis == "Z")
    {
        dir = 2;
    }
    else
    {
        cout << "Undefined axis direction " << axis
             << " in TextOutput::LineConcentration()! Chose \"X\", \"Y\" or"
             << "\"Z\". Now exiting!\n";
        OP_Exit(EXIT_FAILURE);
    }

    ss << fixed << std::setw(9) << std::setfill('0') << int(timestep);
    std::string s = ss.str();
    string name = filename + s + ".opd";

    ofstream file(name, ios::out);                                              //Create file and write header (or overwrite)

    file << std::setw(5) << std::setfill(' ');
    switch (dir)
    {
        case 0:
        {
            file << "x";
            break;
        }
        case 1:
        {
            file << "y";
            break;
        }
        case 2:
        {
            file << "z";
            break;
        }
    }

    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
    {
        string temp;
        string name(Cx.Component[comp].Name);

        if(mode == 0)
        {
            temp = "wf." + name;
        }
        else if(mode == 1)
        {
            temp = "wp." + name;
        }
        else if(mode == 2)
        {
            temp = "mf." + name;
        }
        else if(mode == 3)
        {
            temp = "mp." + name;
        }

        file << std::setw(tablength) << std::setfill(' ') << temp;
    }

    file << "\n";

    switch(dir)                                                                 // Write data to the file
    {
        case 0:
        {
            for(int x = 0; x < Cx.Grid.Nx; x++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << x;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z,{comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z,{comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << "\n";
            }
            break;
        }
        case 1:
        {
            for(int y = 0; y < Cx.Grid.Ny; y++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << y;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z,{comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z,{comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << "\n";
            }
            break;
        }
        case 2:
        {
            for(int z = 0; z < Cx.Grid.Nz; z++)
            {
                file << fixed << std::setw(5) << std::setfill(' ') << z;
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double val = 0.0;
                    switch(mode)
                    {
                        case 0:
                        {//WF
                            val = Cx.WeightFractionsTotal(x,y,z)({comp});
                            break;
                        }
                        case 1:
                        {//WP
                            val = Cx.WeightFractionsTotal(x,y,z)({comp})*100.0;
                            break;
                        }
                        case 2:
                        {//MF
                            val = Cx.MoleFractionsTotal(x,y,z,{comp});
                            break;
                        }
                        case 3:
                        {//MP
                            val = Cx.MoleFractionsTotal(x,y,z,{comp})*100.0;
                            break;
                        }
                    }
                    file << scientific << std::setw(tablength) << std::setfill(' ')
                         << val;
                }
                file << "\n";
            }
            break;
        }
    }
    file.close();                                                               //Close file
}

NodeAB<double,double> TextOutput::IntegrateDG(PhaseField& Phase, DrivingForce& dG, std::string filename, double RealTime)
{
    NodeAB<double,double> dg;
    NodeAB<double,double> num_point;
    STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0)
    {
        if (Phase.Fields(i,j,k).interface())
        {
            for(auto it = dG.Force(i,j,k).begin();
                     it != dG.Force(i,j,k).end(); ++it)
            {
                size_t indexA = it->indexA;
                size_t indexB = it->indexB;

                size_t pIndexA = Phase.FieldsProperties[indexA].Phase;
                size_t pIndexB = Phase.FieldsProperties[indexB].Phase;
                double phase = 1.; //sqrt(Phase.Fields(i,j,k).get_value(indexA)*Phase.Fields(i,j,k).get_value(indexB));
                dg.add_asym1(pIndexA, pIndexB, phase*it->raw);
                num_point.add_sym1(pIndexA, pIndexB,phase);
            }
        }
    }
    STORAGE_LOOP_END
    #ifdef MPI_PARALLEL
    for (size_t i = 0; i < Phase.Nphases; ++i)
    for (size_t j = 0; j < i; ++j)
    {
        double value = dg.get_asym1(i,j);
        double rvalue = 0;
        OP_MPI_Allreduce(&value,&rvalue, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        dg.set_asym1(i,j,rvalue);
        double numpoints = num_point.get_sym1(i,j);
        double rnumpoints = 0;
        OP_MPI_Allreduce(&numpoints,&rnumpoints, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        num_point.set_asym1(i,j,rnumpoints);
    }
    #endif
    for (size_t i = 0; i < Phase.Nphases; ++i)
    for (size_t j = 0; j < i; ++j)
    {
        double value = dg.get_asym1(i,j);
        double numpoints = std::abs(num_point.get_sym1(i,j));
        if (numpoints > 0.)
            dg.set_asym1(i,j,value/numpoints);
        else
            dg.set_asym1(i,j,0.);
    }
    #ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    #endif
    {
        if (!std::filesystem::exists(filename) or RealTime == 0.)
        {
            ofstream file(filename, ios::out | ios::trunc);
            file << "Time";
            for (size_t i = 0; i < Phase.Nphases; ++i)
            for (size_t j = i+1; j < Phase.Nphases; ++j)
            {
                file << separator << "dG_" << i << "_" << j;
            }
            file << "\n";
            file.close();
        }

        ofstream file(filename, ios::app);
        file << RealTime;
        for (size_t i = 0; i < Phase.Nphases; ++i)
        for (size_t j = i+1; j < Phase.Nphases; ++j)
        {
            file << separator << dg.get_asym1(i,j);
        }
        file << "\n";
        file.close();
    }
    return dg;
}

NodeAB<double,double> TextOutput::IntegrateDG(PhaseField& Phase, DrivingForce& dG, std::string filename, NodeAB<double,double> dgold, double RealTime)
{
    NodeAB<double,double> dg;
    NodeAB<double,double> num_point;
    STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0)
    {
        if (Phase.Fields(i,j,k).interface())
        {
            for(auto it = dG.Force(i,j,k).begin();
                     it != dG.Force(i,j,k).end(); ++it)
            {
                size_t indexA = it->indexA;
                size_t indexB = it->indexB;

                size_t pIndexA = Phase.FieldsProperties[indexA].Phase;
                size_t pIndexB = Phase.FieldsProperties[indexB].Phase;
                double phase = 1.; //sqrt(Phase.Fields(i,j,k).get_value(indexA)*Phase.Fields(i,j,k).get_value(indexB));
                dg.add_asym1(pIndexA, pIndexB, phase*it->raw);
                num_point.add_sym1(pIndexA, pIndexB,phase);
            }
        }
    }
    STORAGE_LOOP_END
    #ifdef MPI_PARALLEL
    for (size_t i = 0; i < Phase.Nphases; ++i)
    for (size_t j = 0; j < i; ++j)
    {
        double value = dg.get_asym1(i,j);
        double rvalue = 0;
        OP_MPI_Allreduce(&value,&rvalue, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        dg.set_asym1(i,j,rvalue);
        double numpoints = num_point.get_sym1(i,j);
        double rnumpoints = 0;
        OP_MPI_Allreduce(&numpoints,&rnumpoints, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
        num_point.set_asym1(i,j,rnumpoints);
    }
    #endif
    for (size_t i = 0; i < Phase.Nphases; ++i)
    for (size_t j = 0; j < i; ++j)
    {
        double value = dg.get_asym1(i,j);
        double numpoints = std::abs(num_point.get_sym1(i,j));
        if (numpoints > 0.)
            dg.set_asym1(i,j,value/numpoints);
        else
            dg.set_asym1(i,j,0.);
    }
    #ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    #endif
    {
        if (!std::filesystem::exists(filename) or RealTime == 0.)
        {
            ofstream file(filename, ios::out);
            file << "Time";
            for (size_t i = 0; i < Phase.Nphases; ++i)
            for (size_t j = i+1; j < Phase.Nphases; ++j)
            {
                file << separator << "dG_" << i << "_" << j;
            }
            file << "\n";
            file.close();
        }

        ofstream file(filename, ios::app);
        file << RealTime;
        for (size_t i = 0; i < Phase.Nphases; ++i)
        for (size_t j = i+1; j < Phase.Nphases; ++j)
        {
            file << separator << dg.get_asym1(i,j)-dgold.get_asym1(i,j);
        }
        file << "\n";
        file.close();
    }
    return dg;
}
}// namespace openphase
