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

 *   File created :   2010
 *   Main contributors :   Oleg Shchyglo; Dmitry Medvedev
 *
 */

#include "Settings.h"
#include "RunTimeControl.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Composition.h"
#include "Temperature.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Tools/TimeInfo.h"
#include "Tools/MicrostructureAnalysis.h"
#include "EquilibriumPartitionDiffusionBinary.h"

using namespace std;
using namespace openphase;

/*********** <<< The Main >>> ***********/
int main(int argc, char *argv[])
{
#ifdef MPI_PARALLEL
    int provided = 0;
    OP_MPI_Init_thread(&argc, &argv, OP_MPI_THREAD_FUNNELED, &provided);
    OP_MPI_Comm_rank(OP_MPI_COMM_WORLD, &MPI_RANK);
    OP_MPI_Comm_size(OP_MPI_COMM_WORLD, &MPI_SIZE);
    {
#endif

    Settings                            OPSettings;
    OPSettings.ReadInput();

    RunTimeControl                      RTC(OPSettings);
    PhaseField                          Phi(OPSettings);
    DoubleObstacle                      DO(OPSettings);
    InterfaceProperties                 IP(OPSettings);
    EquilibriumPartitionDiffusionBinary DF(OPSettings);
    Composition                         Cx(OPSettings);
    Temperature                         Tx(OPSettings);
    DrivingForce                        dG(OPSettings);
    BoundaryConditions                  BC(OPSettings);
    BoundaryConditions                  BC_C(OPSettings);
    TimeInfo                            Timer(OPSettings, "Execution Time Statistics");

    BC_C.BCNZ = BoundaryConditionTypes::Fixed;

    ofstream fit_out;
    ofstream tip_velocity;

    if(RTC.Restart)
    {
        cout << "Restart data being read!" << endl;
        cout << "Restart time step: " << RTC.tStart << endl;
        Phi.Read(OPSettings, BC, RTC.tStart);
        Cx.Read(OPSettings, BC, RTC.tStart);
        Tx.Read(OPSettings, BC, RTC.tStart);
        cout << "Done reading restart parameters!" << endl;
    }
    else
    {
        ignore_result(system("rm -rf DendriteProperties"));
        ignore_result(system("mkdir DendriteProperties"));

        tip_velocity.open("tStep-tip_velocity-tip_temperature.dat", ios::out);
        tip_velocity << "tStep \t tip_velocity \t tip_temperature" << endl;
        tip_velocity.close();

        Initializations::Single(Phi, 0, BC);

        Phi.PlantGrainNucleus(1, 0, 0, 0);
        Cx.SetInitialMoleFractions(Phi);
        Tx.SetInitial(BC);
        IP.Set(Phi, BC);
    }
    cout << "Initialization stage done!" << endl;

    int x1 = 0;
    int y1 = 0;
    int z1 = (OPSettings.Grid.Nz)/3 + Phi.Grid.iWidth;

    // Evaluating the dendrite tip and trunk radius
    double old_tip_pos = 0.0;
    double old_time = 0.0;
    // End of evaluating the dendrite tip and trunk radius

    cout << "Diffusion max time step  : " << DF.ReportMaximumTimeStep(Tx) << endl;
    cout << "Phase-field max time step: " << IP.ReportMaximumTimeStep() << endl;
    //RTC.dt = DF.ReportMaximumTimeStep(Tx);

    cout << "Simulation time step     : " << RTC.dt << endl;

    cout << "Entering the Time Loop!!!" << endl;

    for(RTC.tStep = RTC.tStart; RTC.tStep <= RTC.nSteps; RTC.IncrementTimeStep())
    {
        Timer.SetStart();
        dG.Clear();
        Timer.SetTimeStamp("Clear driving forces");
        IP.Set(Phi, BC);
        Timer.SetTimeStamp("Set IPs");
        //DF.CalculateInterfaceMobility(Phi, Tx, BC, IP);
        DO.CalculatePhaseFieldIncrements(Phi, IP, dG);
        Timer.SetTimeStamp("Get PsiDot");
        DF.CalculateDrivingForce(Phi, Cx, Tx, dG);
        Timer.SetTimeStamp("Chemical Driving Force");
        dG.Average(Phi, BC);
        Timer.SetTimeStamp("Driving Force Average");
        dG.MergePhaseFieldIncrements(Phi, IP);
        Timer.SetTimeStamp("Driving Force merge to Psi");
        Phi.NormalizeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Psi Normalize");
        DF.SolveDiffusion(Phi, Cx, Tx, BC_C, RTC.dt);
        Timer.SetTimeStamp("Solve diffusion");
        Tx.Set(BC, Phi, RTC.SimulationTime, RTC.dt);
        Timer.SetTimeStamp("Set temperature");
        Phi.MergeIncrements(BC, RTC.dt);
        Timer.SetTimeStamp("Merge Phase Fields");

        if (Phi.Fields(0, 0, OPSettings.Grid.Nz/3).get_value(0) <= 0.4)
        {
            Phi.MoveFrame(0,0,1, BC);
            Tx.MoveFrame(0,0,1, BC);
            Cx.MoveFrame(0,0,1, BC_C);
            old_tip_pos -= 1;
        }
        Timer.SetTimeStamp("Moving Frame");

        // Evaluating the dendrite tip and trunk radius
        {
            double dx        = OPSettings.Grid.dx;
            double delta     = 1.0e-4;
            int    tip_pos   = old_tip_pos;
            double tip_vel   = 0.0;
            bool   tip_found = false;
            for(int z = (OPSettings.Grid.Nz)/4; z < OPSettings.Grid.Nz - 1; z++)
            if(Phi.Fields(0, 0, z).get_value(0) > 0.5 - delta and
               Phi.Fields(0, 0, z).get_value(0) < 0.5 + delta and
               z != old_tip_pos)
            {
                tip_pos = z;
                tip_vel = (tip_pos - old_tip_pos)*dx/(RTC.SimulationTime - old_time);
                old_tip_pos = z;
                old_time = RTC.SimulationTime;
                tip_found = true;
                tip_velocity.open("tStep-tip_velocity-tip_temperature.dat", ios::app);
                tip_velocity << RTC.tStep << ", " << tip_vel << ", " << Tx.at(0,0,z) << endl;
                tip_velocity.close();
                break;
            }
            if (tip_found)
            {
                z1 = tip_pos;
                vector<double> parabola (tip_pos, 0.0);
                vector<double>   radius (tip_pos, 0.0);
                vector<double>   radius2(tip_pos, 0.0);
                vector<double> tip_position(tip_pos, 0.0);
                for(int z = tip_pos-1; z >= 0; z--)
                {
                    parabola[z] = MicrostructureAnalysis::FindValuePosition(Phi,1,0.5,(dVector3){0,0,(double)z},(dVector3){1,0,0})[0];
                }
                for(int z = tip_pos-1; z > 0; z--)
                {
                    double xx1 = -parabola[z] - 1;
                    double xx2 = 0.0;
                    double xx3 = parabola[z];
                    double zz1 = z;
                    double zz2 = tip_pos;
                    double zz3 = z;
                    radius[z] = dx/(2.0*((-(xx2*zz1) + xx3*zz1 + xx1*zz2 - xx3*zz2 - xx1*zz3 + xx2*zz3)/
                                        ((xx1 - xx2)*(xx2 - xx3)*(xx1 - xx3))));
                    xx1 = parabola[z-1];
                    xx3 = parabola[z];
                    zz1 = z-1;
                    zz3 = z;
                    radius2[z] = dx/(2.0*(zz3 - zz1)/(pow(xx1,2) - pow(xx3,2)));
                    tip_position[z] = -((pow(xx3,2)*zz1 - pow(xx1,2)*zz3)/(pow(xx1,2) - pow(xx3,2)));
                }
                fstream parabola_out(string("DendriteProperties/DendriteContour_") + to_string(RTC.tStep)+ string(".txt"), ios::out);
                parabola_out << tip_pos << " " << 0 << endl;
                for(int z = tip_pos-1; z >= 0; z--)
                {
                    parabola_out << z << " " << parabola[z] << endl;
                }
                parabola_out.close();
                fstream curvature_radius_out(string("DendriteProperties/DendriteTipRadius_") + to_string(RTC.tStep)+ string(".txt"), ios::out);
                curvature_radius_out << "TipPosition" << " " << "TipRadius"<< " " << "TipPosition2" << " "<< "TipRadius2" << " " << "dTipRadius_dz" << " " << "d2TipRadius_dz2" << " " << "d2TipRadius_dz2-5point" << endl;
                curvature_radius_out << tip_pos << " " << 0 << " " << tip_pos << " " << 0 << " " << 0 << " " << 0 << endl;
                for(int z = tip_pos-1; z > 1; z--)
                {
                    curvature_radius_out << z << " " << radius[z]<< " " << tip_position[z] <<  " " << radius2[z] << " " << radius[z-1] - radius[z+1] << " " << radius[z-1]+radius[z+1]-2.0*radius[z] << endl;
                }
                curvature_radius_out.close();
            }
        }
        // End of evaluating the dendrite tip and trunk curvature

        //  Output to VTK file
        if (RTC.WriteVTK())
        {
            // Write data in VTK format
            Phi.WriteVTK(OPSettings, RTC.tStep);
            Cx.WriteVTK(OPSettings, RTC.tStep);
            Tx.WriteVTK(OPSettings, RTC.tStep);
            IP.WriteVTK(OPSettings, RTC.tStep);
            dG.WriteVTK(OPSettings, Phi, RTC.tStep);
            Cx.WriteStatistics(OPSettings,RTC.tStep, RTC.SimulationTime);
        }
        Timer.SetTimeStamp("Write VTK files");

        // Write raw data
        if (RTC.WriteRawData())
        {
            Phi.Write(OPSettings,RTC.tStep);
            Cx.Write(OPSettings,RTC.tStep);
            Tx.Write(OPSettings,RTC.tStep);
        }
        Timer.SetTimeStamp("Write raw data");

        //  Output to screen
        if(RTC.WriteToScreen())
        {
            std::string message  = ConsoleOutput::GetStandard("Interface energy density", DO.AverageEnergyDensity(Phi, IP));
                        message += ConsoleOutput::GetStandard("Interface energy", DO.Energy(Phi, IP));
            ConsoleOutput::WriteTimeStep(RTC, message);
            // Statistics
            Phi.PrintPointStatistics(x1,y1,z1);
            Cx.PrintPointStatistics(x1,y1,z1);
            Tx.PrintPointStatistics(x1,y1,z1);
            dG.PrintDiagnostics();
            Phi.PrintPFVolumes();
            Timer.PrintWallClockSummary();
        }
    } //end time loop
#ifdef MPI_PARALLEL
    }
    OP_MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
