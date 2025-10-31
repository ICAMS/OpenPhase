#include "BoundaryConditions.h"
#include "ConsoleOutput.h"
#include "DrivingForce.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"
#include "Initializations.h"
#include "PhaseField.h"
#include "RunTimeControl.h"
#include "Settings.h"

using namespace std;
using namespace openphase;

int main(int argc, char *argv[])
{
    Settings                    OPSettings;
    OPSettings.ReadInput();

    RunTimeControl              RTC(OPSettings);
    BoundaryConditions          BC(OPSettings);
    PhaseField                  Phi(OPSettings);
    ElasticProperties           EP(OPSettings);
    ElasticitySolverSpectral    ES(OPSettings);

    DrivingForce                dG(OPSettings);

    double centerx = (OPSettings.Grid.dNx > 0) ? (OPSettings.Grid.Nx+1)/2.0 : 0;
    double centery = (OPSettings.Grid.dNy > 0) ? (OPSettings.Grid.Ny+1)/2.0 : 0;
    double centerz = (OPSettings.Grid.dNz > 0) ? (OPSettings.Grid.Nz+1)/2.0 : 0;

    int iRadius = 10;

    Initializations::Single(Phi, 0, BC);
    Initializations::Sphere(Phi, 1, iRadius, centerx, centery, centerz, BC);
    ConsoleOutput::WriteStandard("Initial Microstructure", "Set");
    EP.SetEffectiveProperties(Phi);
    ConsoleOutput::WriteStandard("Grains Properties", "Set");

    ConsoleOutput::WriteStandard("Calculation", "Started");

    int niterations = ES.Solve(EP, BC, RTC.dt);

    //EP.CalculateDrivingForce(Phi,dG);
    //dG.Average(Phi,OPSettings);
    //dG.WriteVTK(OPSettings,0,0,1);

    ConsoleOutput::WriteStandard("Solver iterations", niterations);
    ConsoleOutput::WriteStandard("Elastic energy", EP.Energy(Phi));

    Phi.WriteDistortedVTK(OPSettings,EP,0,true);
    EP.WriteStressesVTK(OPSettings,0);
    EP.WriteTotalStrainsVTK(OPSettings,0);

    fstream out("LineStresses.dat", ios::out);
    string dlm = ",";

    out << std::fixed;
    out << "X" << dlm << "Sigma_11" << dlm << "Sigma_33" << dlm << "Esh_Sigma_11" << dlm << "Esh_Sigma_33" << endl;

    double C11 = EP.PhaseElasticConstants[0](0,0);
    double C12 = EP.PhaseElasticConstants[0](0,1);
    double nu = C12/(C12 + C11);

    vStrain locEigenStrain = EP.EigenStrains(centerx,centery,centerz);

    double epsilon = locEigenStrain[0];
    double sigma_0 = 2.0/3.0*epsilon*(C11 + 2.0*C12)*(1.0 - 2.0*nu)/(1.0 - nu);
    for (int i = centerx; i < centerx + iRadius; i++)
    {
        /* First Piola-Kirchhoff stress:*/
        //dMatrix3x3 locStress = EP.DeformationGradientsTotal(i, centery, centerz) * EP.Stresses(i, centery, centerz).tensor();
        /* Second Piola-Kirchhoff stress:*/
        dMatrix3x3 locStress = EP.Stresses(i, centery, centerz).tensor();

        double sigma_11 = -sigma_0;
        double sigma_33 = -sigma_0;
        out << double(i - centerx)/(OPSettings.Grid.Nx-1)   << dlm
            << locStress(0,0)                          << dlm
            << locStress(2,2)                          << dlm
            << sigma_11                                << dlm
            << sigma_33                                << endl;
    }
    for (int i = centerx + iRadius; i < OPSettings.Grid.Nx; i++)
    {
        /* First Piola-Kirchhoff stress:*/
        //dMatrix3x3 locStress = EP.DeformationGradientsTotal(i, centery, centerz) * EP.Stresses(i, centery, centerz).tensor();
        /* Second Piola-Kirchhoff stress:*/
        dMatrix3x3 locStress = EP.Stresses(i, centery, centerz).tensor();

        double sigma_11 = -sigma_0*pow(iRadius/double(i-centerx), 3);
        double sigma_33 = 0.5*sigma_0*pow(iRadius/double(i-centerx), 3);

        out << double((i - centerx))/(OPSettings.Grid.Nx-1) << dlm
            << locStress(0,0)                          << dlm
            << locStress(2,2)                          << dlm
            << sigma_11                                << dlm
            << sigma_33                                << endl;
    }
    out.close();

    ConsoleOutput::WriteLine();
    simulation_end();
    return 0;
}
