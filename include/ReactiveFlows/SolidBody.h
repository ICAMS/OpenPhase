/*
 *  This file is part of the OpenPhase (R) software library.
 *
 *  Copyright (c) 2009-2026 Ruhr-Universitaet Bochum,
 *                Universitaetsstrasse 150, D-44801 Bochum, Germany
 *            AND 2018-2026 OpenPhase Solutions GmbH,
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
 *
 *  File created :   2025
 *  Main contributors :   Oleg Shchyglo; Reza Namdar
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

#include "MixtureFlow.h"
#include "Energy.h"

using namespace std;

namespace openphase {
class Energy;
class MixtureFlow;
 }
namespace openphase  {

struct Body
{
    double x, y, z, r;
};

class OP_EXPORTS SolidBody
{
    public:
	SolidBody(Settings& locSettings, const std::string InputFileName = DefaultInputFileName)          ///< Constructor
    {
		ReadInput(InputFileName);
		Initialize(locSettings);
    }
	void ReadInput(string InputFile);
	void Initialize(Settings& locSettings);

	bool isOverlapping(const Body& a, const Body& b, int nDim, int DistPlaneAxis) const;
	std::vector<Body> generateNonOverlappingBodies(int nDim,double Lx, double Ly, double Lz,      // physical box lengths
	                                               int N, double minRadius, double maxRadius,
                                                       double dx, int DistPlaneAxis, double DistPlaneW0);
	void saveBodiesToFile(const std::vector<Body>& bodies, const std::string& filename);
	std::vector<Body> loadBodiesFromFile(const std::string& filename);
	void DistributeSolidBodies(int DistPlaneAxis, double DistPlaneW0, const std::string& dir);

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
	double TrilinearInterPolation(double ic, double jc, double kc,  
                                  double q000, double q100, double q010, double q110, 
                                  double q001, double q101, double q011, double q111);

	void CalculateLocalNusseltNumber(PhaseField& Phase, Energy& EN, MixtureFlow& MF, FlowSolverLBM& FL, double Length,
                                        double Tin, double dn, double Time, string dir, string fname, BoundaryConditions& BC);
	//void CalculateDragCoeff(PhaseField& Phase, FlowSolverLBM& FL, MixtureFlow& MF, double A_proj, double Rho0, double U0, double dn, double Time,
	 //       size_t axis, string dir, string fname, BoundaryConditions& BC);
	void CalculateDragCoeff(PhaseField& Phase, MixtureFlow& MF, FlowSolverLBM& FL, double rhoRef, double Uref, double pRef,
                            const dVector3& dragDir, double dn, double Time, std::string dir, std::string fname,  BoundaryConditions& BC);

	double CalculateInterpolatedTemperature(Energy& EN, double ic, double jc, double kc);
	dVector3 CalculateInterpolatedVelocity(MixtureFlow& MF,double ic, double jc, double kc);
	double CalculateInterpolatedPressure(FlowSolverLBM& FL, double ic, double jc, double kc);
	double CalculateInterpolatedViscosity(FlowSolverLBM& FL, double ic, double jc, double kc);
    void CalculateDistanceField(PhaseField& Phase);
	void WriteVTKDistanceField(Settings& locSettings, const int tStep, const int precision = 16);
	void UpdateFields();
	void SetBoundaryConditions(BoundaryConditions& BC);
	dVector3 DistanceToInterface(const dVector3& X1, const dVector3& X2);
	void AdvectPhaseField(PhaseField &Phase, FlowSolverLBM &FL, MixtureFlow &MF, BoundaryConditions &BC, double dt);
	void CalculateDistancePeriodic(dVector3 A, dVector3 B, dVector3 &dist, double Nx, double Ny, double Nz);
	void WriteVTKSolidVelocity(Settings& locSettings, const int tStep, const int precision = 16);

    GridParameters Grid;                                                        ///< Simulation grid parameters
    size_t nParticles;
	int FirstRow;
	int nRows;
	int nCols;
	double Porosity;
	double PartDiameter;

	double nDist;
	double nFallBack;
	
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


    Storage3D<double,0> DistanceField;                                                   ///< Distance Field
    Storage3D<dVector3,0> NormField;                                                   ///< Normal vector Field
	string IntProf;

	vector<Body>  rand_Bodies;

		// --- scalar masked trilinear ---
	template <class ValueAt , class keepCorner>
	double MaskedTrilinearScalar(double x, double y, double z,
	                             ValueAt&& value_at,
	                             keepCorner&& keep_corner, const double &fallback)
	{
	    constexpr double kTiny = 1e-30;

	    const int i0 = (int)std::floor(x);
	    const int j0 = (int)std::floor(y);
	    const int k0 = (int)std::floor(z);
	    const int i1 = i0 + Grid.dNx;
	    const int j1 = j0 + Grid.dNy;
	    const int k1 = k0 + Grid.dNz;

	    const double fx = x - i0;
	    const double fy = y - j0;
	    const double fz = z - k0;

	    struct Corner { int i,j,k; double w; };
	    const Corner C[8] = {
	        {i0,j0,k0,(1-fx)*(1-fy)*(1-fz)},
	        {i1,j0,k0,( fx)*(1-fy)*(1-fz)},
	        {i0,j1,k0,(1-fx)*( fy)*(1-fz)},
	        {i1,j1,k0,( fx)*( fy)*(1-fz)},
	        {i0,j0,k1,(1-fx)*(1-fy)*( fz)},
	        {i1,j0,k1,( fx)*(1-fy)*( fz)},
	        {i0,j1,k1,(1-fx)*( fy)*( fz)},
	        {i1,j1,k1,( fx)*( fy)*( fz)}
	    };

	    double num = 0.0, den = 0.0;
	    for (const auto& c : C)
	    {
	        if (!keep_corner(c.i, c.j, c.k)) continue;
	        num += c.w * value_at(c.i,c.j,c.k);
	        den += c.w;
	    }

	    if (den > kTiny) return num / den;
		return fallback;
	}

	// --- vector3 masked trilinear ---
	template <class ValueAt , class keepCorner>
	dVector3 MaskedTrilinearVector3(double x, double y, double z,
	                                       ValueAt&& value_at,
	                                       keepCorner&& keep_corner, const dVector3 &fallback)
	{
	    constexpr double kTiny = 1e-30;

	    const int i0 = (int)std::floor(x);
	    const int j0 = (int)std::floor(y);
	    const int k0 = (int)std::floor(z);
	    const int i1 = i0 + Grid.dNx;
	    const int j1 = j0 + Grid.dNy;
	    const int k1 = k0 + Grid.dNz;

	    const double fx = x - i0;
	    const double fy = y - j0;
	    const double fz = z - k0;

	    struct Corner { int i,j,k; double w; };
	    const Corner C[8] = {
	        {i0,j0,k0,(1-fx)*(1-fy)*(1-fz)},
	        {i1,j0,k0,( fx)*(1-fy)*(1-fz)},
	        {i0,j1,k0,(1-fx)*( fy)*(1-fz)},
	        {i1,j1,k0,( fx)*( fy)*(1-fz)},
	        {i0,j0,k1,(1-fx)*(1-fy)*( fz)},
	        {i1,j0,k1,( fx)*(1-fy)*( fz)},
	        {i0,j1,k1,(1-fx)*( fy)*( fz)},
	        {i1,j1,k1,( fx)*( fy)*( fz)}
	    };

	    dVector3 num;
	    double den = 0.0;

	    for (const auto& c : C)
	    {
	        if (!keep_corner(c.i, c.j, c.k)) continue;
	        const dVector3 v = value_at(c.i,c.j,c.k);
	        for (size_t s=0; s<3; ++s) num[s] += c.w * v[s];
	        den += c.w;
	    }

	    if (den > kTiny)
	    {
	        const double invDen = 1.0 / den;
	        for (size_t s=0; s<3; ++s) num[s] *= invDen;
	        return num;
	    }
		return fallback;
	}
	// --- tensor masked trilinear (Tensor<double,1> of size nSpecies) ---
	template <class ValueAt , class keepCorner>
	Tensor<double,1> MaskedTrilinearTensor(double x, double y, double z,
	                                       ValueAt&& value_at,
	                                       keepCorner&& keep_corner,
	                                       size_t nSpecies, const Tensor<double,1>& fallback)
	{
	    constexpr double kTiny = 1e-30;

	    const int i0 = (int)std::floor(x);
	    const int j0 = (int)std::floor(y);
	    const int k0 = (int)std::floor(z);
	    const int i1 = i0 + Grid.dNx;
	    const int j1 = j0 + Grid.dNy;
	    const int k1 = k0 + Grid.dNz;

	    const double fx = x - i0;
	    const double fy = y - j0;
	    const double fz = z - k0;

	    struct Corner { int i,j,k; double w; };
	    const Corner C[8] = {
	        {i0,j0,k0,(1-fx)*(1-fy)*(1-fz)},
	        {i1,j0,k0,( fx)*(1-fy)*(1-fz)},
	        {i0,j1,k0,(1-fx)*( fy)*(1-fz)},
	        {i1,j1,k0,( fx)*( fy)*(1-fz)},
	        {i0,j0,k1,(1-fx)*(1-fy)*( fz)},
	        {i1,j0,k1,( fx)*(1-fy)*( fz)},
	        {i0,j1,k1,(1-fx)*( fy)*( fz)},
	        {i1,j1,k1,( fx)*( fy)*( fz)}
	    };

	    Tensor<double,1> num({nSpecies});
	    for (size_t s=0; s<nSpecies; ++s) num[s] = 0.0;
	    double den = 0.0;

	    for (const auto& c : C)
	    {
	        if (!keep_corner(c.i, c.j, c.k)) continue;
	        const Tensor<double,1> v = value_at(c.i,c.j,c.k);
	        for (size_t s=0; s<nSpecies; ++s) num[s] += c.w * v[s];
	        den += c.w;
	    }

	    if (den > kTiny)
	    {
	        const double invDen = 1.0 / den;
	        for (size_t s=0; s<nSpecies; ++s) num[s] *= invDen;
	        return num;
	    }
		return fallback;
	}

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
