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
 *   Main contributors :   Raphael Schiedung
 *
 */

#ifndef NOISE_H
#define NOISE_H

#include "Includes.h"
#include "fftw3.h"

namespace openphase
{
class DrivingForce;
class PhaseField;
class Settings;
class Temperature;

class Noise : public OPObject                                                   ///< Calculates a noise, which can be added to the DrivingForce
{
 public: 
    ~Noise(void);                                                               ///< Destructor
    Noise(void){};                                                              ///< Simple constructor, but Initialize and ReadInput have to be called additionally
    Noise(Settings& locSettings,
          const std::string InputFileName = DefaultInputFileName)               ///< Constructor includes the call of Initialize and ReadInput
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    };
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Allocates memory, initializes the settings
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input from input file
    void UpdateRaw(int tStep);                                                  ///< Updates Raw Driving force
    void AddToDrivingForce(int tStep, Temperature& Temp, DrivingForce& DF);     ///< Adds Fields to driving force
    void AddToDrivingForce(int tStep, DrivingForce& DF);                        ///< Adds Fields to driving force
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input from input file
    void WriteVTK(const Settings& locSettings,
                  const int tStep,
                  const int precision = 16);                                    ///< Writes the noise acting on the driving force from phase field with indexB onto phase field with indexA in VTK format

    double CutOffWaveLength;                                                    ///< Cut off wave length of noise
    double FixAmpl;                                                             ///< Cut off wave length of noise
    double gamma;                                                               ///< Macroscopic response of the system
    int    BCells;                                                              ///< Number of boundary cells around the computation domain
    int    Nphases;                                                             ///< Number of phase-fields

    GridParameters Grid;                                                        ///< Simulation grid parameters

    int    RandomSeed;                                                          ///< Defines the random number table (used to make results reproducible)
    int    TimeSteps;                                                           ///< Number of iteration after which a new noise will be scrabbled

    //const double kBoltzmann;                                                  ///< Physical constant

 protected:

    void Generate(void);                                                        ///< Generates a new noise field (it will be interpolated between new and old noise field)

    Storage3D< double, 0 >           Raw;                                       ///< Real raw noise, which will be merge into the driving force
    fftw_complex*                    fftw_In;                                   ///< FFTW input pointer
    fftw_complex*                    fftw_Out;                                  ///< FFTW out pointer
    fftw_plan                        FFTBackward;                               ///< Plan for the execution of FFTW
    std::complex<double>*            RandomFourier;                             ///< Pointer to storage of Fourier-space random number distribution
    std::complex<double>*            RandomReal;                                ///< Pointer to storage of real space random number distribution
    std::default_random_engine       generator;                                 ///< Random number generator
    std::normal_distribution<double> distribution;                              ///< Random number distribution
 private:
};
}
#endif
