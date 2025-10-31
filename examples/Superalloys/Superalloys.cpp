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

 *   File created :   2017
 *   Main contributors : Muhammad Adil Ali; Johannes Görler
 *
 */

#include "Mechanics.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "Initializations.h"
#include "InterfaceProperties.h"
#include "TextOutput.h"
#include "BoundaryConditions.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "Composition.h"
#include "RunTimeControl.h"
#include "Tools/TimeInfo.h"
#include "Tools.h"
#include "VTK.h"

using namespace std;
using namespace openphase;

void WriteEnergy(double EE, double IE, RunTimeControl& RTC, string OutFile);

int main(int argc, char *argv[])
{
    string InputFile = "ProjectInput.opi";
    Settings                            OPSettings(InputFile);
    RunTimeControl                      RTC(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    ElasticProperties                   EP(OPSettings);
    ElasticitySolverSpectral            ES(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    ConsoleOutput::WriteLineInsert("Initialization stage done!");

    if(RTC.Restart)
    {
        ConsoleOutput::WriteLineInsert("Restart data being read! ");
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC,RTC.tStart);
        EP.Read(OPSettings, BC, RTC.tStart);
        EP.SetEffectiveProperties(Phi);
        ES.Solve(EP, BC, RTC.dt);
        ConsoleOutput::WriteLineInsert("Done");
    }
    else
    {
        Initializations::Single(Phi, 0, BC);
//        Initializations::Sphere(Phi, 1, 5, (OPSettings.Nx+1)/2, (OPSettings.Ny+1)/2, (OPSettings.Nz+1)/2, BC);
//        int particleIdx1 = Initializations::Sphere(Phi, 1, 5, OPSettings.Nx/2, OPSettings.Ny/2, OPSettings.Nz/2, BC);
//        int particleIdx2 = Initializations::Sphere(Phi, 1, 5, 0, 0, 0, BC);

        ConsoleOutput::WriteLineInsert("Planting particles on a regular grid with random offset! ");
        Initializations::QuasiRandomNuclei(Phi, 1, {8, 8, 8}, {16, 16, 16}, {4, 4, 4}, 0.05);
        ConsoleOutput::WriteLineInsert("Done!");
        Tx.SetInitial(BC);
        Cx.SetInitialMoleFractions(Phi);
        EP.SetEffectiveProperties(Phi);
    }

    ConsoleOutput::WriteLineInsert("Entering the Time Loop!!!");
    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        dG.Clear();
        IP.Set(Phi, BC);
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        EP.SetEffectiveProperties(Phi);
        ES.Solve(EP, BC, RTC.dt);
        EP.CalculateDrivingForce(Phi, dG);
        dG.Average(Phi, BC);

        if (RTC.WriteVTK())
        {
            dG.WriteVTK(OPSettings, RTC.tStep, 1, 0);
        }
        dG.MergePhaseFieldIncrements(Phi, IP);
        Phi.NormalizeIncrements(BC, RTC.dt);
        Phi.MergeIncrements(BC, RTC.dt);
        DF.SolveDiffusion(Phi, Cx, Tx, BC, RTC.dt);
        Tx.Set(BC,Phi, RTC.SimulationTime, RTC.dt);

        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings, RTC.tStep);
            Cx.Write(OPSettings, RTC.tStep);
            Tx.Write(OPSettings, RTC.tStep);
            EP.Write(OPSettings, RTC.tStep);
        }

        if (RTC.WriteVTK())
        {
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            EP.WriteStressesVTK(OPSettings, RTC.tStep);
            EP.WriteElasticStrainsVTK(OPSettings, RTC.tStep);
            Phi.WriteAverageVolume(RTC.tStep, 1);
            Cx.WriteStatistics(OPSettings, RTC.tStep, RTC.SimulationTime);
        }
        if (RTC.WriteToScreen())
        {
            double I_En = DO.AverageEnergyDensity(Phi, IP);
            double E_En = EP.AverageEnergyDensity(Phi);
            double T_En = E_En + I_En;
            ConsoleOutput::WriteTimeStep(RTC);
            ConsoleOutput::WriteStandard("Elastic energy density",E_En);
            ConsoleOutput::WriteStandard("Interface energy density",I_En);
            ConsoleOutput::WriteLine("-");
            ConsoleOutput::WriteStandard("Total energy density",T_En);
            ConsoleOutput::WriteLine("=");
            Phi.PrintVolumeFractions();
            dG.PrintDiagnostics();
            WriteEnergy(E_En, I_En, RTC, "Energy");
            ConsoleOutput::WriteLine();
        }
    }
   simulation_end();
   return 0;
}

void WriteEnergy(double EE, double IE, RunTimeControl& RTC, string OutFile)
{
    ofstream output_file;
    stringstream FileName;
    FileName<<OutFile.c_str()<<".txt";
    if (RTC.tStep == RTC.tStart)
    {
        remove(FileName.str().c_str());
        output_file.open(FileName.str().c_str(), ios::out);
        output_file << setw(12) << left << "time"
                    << setw(20) << left << "Elastic Energy"
                    << setw(20) << left << "Interface Energy"
                    << setw(20) << left << "Total Energy";
        output_file<< endl;
        output_file.close();
    }
    output_file.open(FileName.str().c_str(), ios::app);
    output_file << setw(12) << left << RTC.SimulationTime
                << setw(20) << left << EE
                << setw(20) << left << IE
                << setw(20) << left << IE + EE;
    output_file<< endl;
    output_file.close();
}
