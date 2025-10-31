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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
 *
 */

#ifndef SOLIDBODY_H_INCLUDED
#define SOLIDBODY_H_INCLUDED

#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Velocities.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <array>
#include <random>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <tuple>

#include "Initializations.h"

#include "FlowMixture.h"
#include "EnergyTransport.h"

using namespace std;

namespace openphase {
class EnergyTransport;
class FlowMixture;
 }
namespace openphase  {
struct Circle 
{
	double x, z, r;
};

struct Sphere 
{
    double x, y, z, r;
};

class SolidBody
{
    public:
	SolidBody(const std::string InputFileName = DefaultInputFileName)                                  ///< Constructor
    {
		ReadInput(InputFileName);
    }
	void ReadInput(string InputFile);
	void DistributeRandomSolidBodies(PhaseField& Phase, Settings& locSettings, string dir);

	bool isOverlapping(const Circle& c1, const Circle& c2);
	vector<Circle> generateNonOverlappingCircles(double width, double height,int N,double minRadius, double maxRadius, double dx); 
	void saveCirclesToFile(const std::vector<Circle>& circles, const std::string& filename, double dx);
	vector<Circle> loadCirclesFromFile(const std::string& filename);

	void saveSpheresToFile(const std::vector<Sphere>& spheres, const std::string& filename, double dx);
	vector<Sphere> loadSpheresFromFile(const std::string& filename);
	vector<Sphere> generateNonOverlappingSpheres(double Nx, double Ny, double Nz, int N, double minRadius, double maxRadius, double dx);

	// Find 3 nearest neighbors (indices) for point i
	std::array<int, 3> find_3_nearest(const std::vector<dVector3>& pts, size_t idx, BoundaryConditions& BC);
	double quad_area(const dVector3& a, const dVector3& b, const dVector3& c, const dVector3& d);
	dVector3 quad_centroid(const dVector3& a, const dVector3& b, const dVector3& c, const dVector3& d);
	bool is_close(const dVector3& a, const dVector3& b, double tol=0.4) ;
	double dist3(const dVector3& a, const dVector3& b);
	double triangle_area(const dVector3& a, const dVector3& b, const dVector3& c);
	std::vector<double> surface_elements_3d(const std::vector<dVector3>& points, size_t kneigh = 6);
	double curve_element_2d(const vector<dVector3>& points, size_t i, const BoundaryConditions& BC);

	double BilinearInterPolation(double xp, double zp, double Tx1z1, double Tx1z2, double Tx2z1, double Tx2z2);
	double Lagrange_polynomial4 (double T1, double T2, double T3,  double T4, double T5, double f1, double f2, double f3,  double f4);
	double TrilinearInterPolation(double ic, double jc, double kc,  
                                  double q000, double q100, double q010, double q110, 
                                  double q001, double q101, double q011, double q111);

	void CalculateLocalNusseltNumber(PhaseField& Phase, EnergyTransport& ET, FlowMixture& FM, FlowSolverLBM& FL, double Length,
                                        double Tin, double dn, double Time, string dir, string fname, BoundaryConditions& BC);
	void CalculateDragCoeff(PhaseField& Phase, FlowSolverLBM& FL, FlowMixture& FM, double A_proj, double Rho0, double U0, double dn, double Time, 
	        size_t axis, string dir, string fname, BoundaryConditions& BC);
	bool isOverlapping(const Sphere& s1, const Sphere& s2); 
	double CalculateInterpolatedTemperature(EnergyTransport& ET, double ic, double jc, double kc);
	dVector3 CalculateInterpolatedVelocity(FlowMixture& FM,double ic, double jc, double kc);
	double CalculateInterpolatedPressure(FlowSolverLBM& FL, double ic, double jc, double kc);
	double CalculateInterpolatedViscosity(FlowSolverLBM& FL, double ic, double jc, double kc);

    GridParameters Grid;                                                        ///< Simulation grid parameters
    size_t nParticles;
	int FirstRow;
	int nRows;
	int nCols;
	double Porosity;
	double PartDiameter;

	double nDist;
	
	bool   Do_randDist;
	double X0DistZone;
	double XNDistZone;
	double Y0DistZone;
	double YNDistZone;
	double Z0DistZone;
	double ZNDistZone;
	double MinRadious;
	double MaxRadious;
	double Clearance;

	vector<Circle>  rand_Circles;
	vector<Sphere>  rand_Spheres;

	template <typename T>
	void writeData(const std::string& directory, const std::string& definition, const T& parameter)
	{
    	namespace fs = std::filesystem;
		
    	// Create directory if it doesn't exist
    	if (!fs::exists(directory)) 
    	{
    	    if (!fs::create_directories(directory)) {
    	        std::cerr << "Failed to create directory: " << directory << std::endl;
    	        return;
    	    }
    	}
	
    	std::string filePath = directory + "/simData.txt";
    	std::ofstream outFile(filePath, std::ios::app); // append mode
	
    	if (!outFile) 
    	{
    	    std::cerr << "Failed to open file for writing: " << filePath << std::endl;
    	    return;
    	}
	
    	constexpr int columnWidth = 50; // Adjust as needed
    	outFile << std::left << std::setw(columnWidth) << (definition + ":") << parameter << std::endl;
    	outFile.close();
	}

};
}
#endif
