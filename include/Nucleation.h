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
 *   Main contributors :   Oleg Shchyglo; Alexander Monas; Marvin Tegeler;
 *                         Matthias Stratmann
 *
 */

#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "Includes.h"
#include "H5Interface.h"

namespace openphase
{
class PhaseField;
class Settings;
class InterfaceProperties;
class DrivingForce;
class BoundaryConditions;
class Orientations;
class Temperature;
class SymmetryVariants;
class H5Interface;
class ThermodynamicProperties;

enum class NucleiGenerationModes
{
    Static,                                                                     ///< Nucleation seeds are generated once and should be regenerated if new seed positions are required
    Dynamic                                                                     ///< Nucleation seeds are generated on the fly for each nucleation event
};

enum class NucleiOrientationModes                                               ///< Possible orientation modes for new grains
{
    Reference,                                                                  ///< Use frame of reference orientation
    Parent,                                                                     ///< Use parent grain orientation
    Random,                                                                     ///< Generate random orientation
};

enum class NucleiLocationModes                                                  ///< Possible location modes for new grains
{
    Bulk,                                                                       ///< In the bulk of the parent phase
    BulkAndGrainBoundaries,                                                     ///< In the bulk and grain boundaries of the parent phase
    BulkAndInterfaces,                                                          ///< In the bulk of the parent phase and interfaces containing parent phase
    GrainBoundaries,                                                            ///< In the entire region of the grain boundaries of the parent phase
    PhaseBoundaries,                                                            ///< In the boundaries of the parent phase with any other phase
    Junctions,                                                                  ///< In the triple and higher order junctions containing 100% matrix phase
    Interfaces,                                                                 ///< In all interfaces, e.g. grain and phase boundaries involving matrix phase
    XBottom,                                                                    ///< At the bottom of x-axis
    XTop,                                                                       ///< At the top of x-axis
    YBottom,                                                                    ///< At the bottom of y-axis
    YTop,                                                                       ///< At the top of y-axis
    ZBottom,                                                                    ///< At the bottom of z-axis
    ZTop                                                                        ///< At the top of z-axis
};

enum class NucleiSizeDistributions                                              ///< Available seed size distributions
{
    Normal,                                                                     ///< Normal (Gaussian) distribution
    Cauchy,                                                                     ///< Cauchy distribution
    Uniform,                                                                    ///< Uniform distribution
    FixedRadius,                                                                ///< Fixed size seeds
    None                                                                        ///< Seed size is set to half of the interface width
};

enum class NucleiVariantsModes                                                  ///< Variants selection modes
{
    Random,                                                                     ///< Random variant selection
    LowestEnergy                                                                ///< Highest driving force variant selection
};

class Nucleus                                                                   ///< Nucleus properties
{
 public:
    size_t index;                                                               ///< Phase-field index of the nucleus
    size_t phase;                                                               ///< Thermodynamic phase of the nucleus
    size_t parentPhase;                                                         ///< Parent Thermodynamic phase for the nucleus
    size_t variant;                                                             ///< Symmetry (or other) variant
    size_t time_stamp;                                                          ///< Time step of the nucleation event

    double radius;                                                              ///< Effective seed particle radius
    double dGmin;                                                               ///< Minimum driving force for nucleation
    double dGnuc;                                                               ///< Driving force of the nucleus

    Quaternion orientation;                                                     ///< Nucleus orientation
    iVector3 position;                                                          ///< Nucleus position

    bool planted;                                                               ///< True if nucleus is planted successfully (and growing)

    Nucleus() :
        index(0),
        phase(0),
        parentPhase(0),
        variant(0),
        time_stamp(0),
        radius(0.0),
        dGmin(0.0),
        dGnuc(0.0),
        planted(false)
    {
       orientation.set_to_zero();
       position.set_to_zero();
    }
    void read(std::fstream& inp)
    {
        inp >> index;
        inp >> phase;
        inp >> parentPhase;
        inp >> variant;

        inp >> time_stamp;

        inp >> radius;
        inp >> dGmin;
        inp >> dGnuc;

        orientation.read(inp);
        position.read(inp);

        inp >> planted;
    }
    void write(std::fstream& out) const
    {
        out << index << std::endl;
        out << phase << std::endl;
        out << parentPhase << std::endl;
        out << variant << std::endl;
        out << time_stamp << std::endl;

        out << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);

        out << radius << std::endl;
        out << dGmin << std::endl;
        out << dGnuc << std::endl;

        orientation.write(out);
        out << std::defaultfloat;
        position.write(out);

        out << planted << std::endl;
    }
};

struct NucleationParameters                                                     ///< Pairwise nucleation settings ("phase alpha" in "phase beta")
{
    bool   Allowed;                                                             ///< True if nucleation is allowed
    bool   Generated;                                                           ///< True if nucleation seeds were generated
    bool   TrueRadius;                                                          ///< If true use nucleus true radius instead of interface width
    bool   RelativeDensity;                                                     ///< Relative nucleation density if "true" or absolute otherwise
    size_t NucleiPhase;                                                         ///< Index of nucleating phase
    size_t MatrixPhase;                                                         ///< Index of matrix phase
    double Tmin;                                                                ///< Upper temperature limit for nucleation
    double Tmax;                                                                ///< Lower temperature limit for nucleation
    double Density;                                                             ///< Average nucleation density
    size_t Nsites;                                                              ///< Number of requested nucleation sites
    size_t Nseeds;                                                              ///< Number of generated nucleation seeds
    size_t Nvariants;                                                           ///< Number of permitted crystallographic variants per nucleation site
    double ShieldingRadius;                                                     ///< Shielding radius for the nuclei to avoid numerical artifacts
    double SeedRadius;                                                          ///< Radius for fixed size seeds
    double SeedRadiusMIN;                                                       ///< Minimum seed radius
    double SeedRadiusMAX;                                                       ///< Maximum seed radius

    NucleiLocationModes LocationMode;                                           ///< Nuclei location mode
    NucleiOrientationModes OrientationMode;                                     ///< Nuclei orientation mode
    NucleiSizeDistributions Distribution;                                       ///< Nuclei size distribution mode
    NucleiVariantsModes VariantsMode;                                           ///< Nuclei variant selection mode

    std::mt19937_64 SizeGenerator;                                              ///< Nuclei size generator
    std::mt19937_64 VariantsGenerator;                                          ///< Variants generator

    std::mt19937_64 PositionGeneratorX;                                         ///< Nuclei x-position generator
    std::mt19937_64 PositionGeneratorY;                                         ///< Nuclei y-position generator
    std::mt19937_64 PositionGeneratorZ;                                         ///< Nuclei z-position generator

    std::mt19937_64 OrientationGenerator1;                                      ///< Nuclei orientation generator Q1
    std::mt19937_64 OrientationGenerator2;                                      ///< Nuclei orientation generator Q2
    std::mt19937_64 OrientationGenerator3;                                      ///< Nuclei orientation generator Q3

    std::uniform_int_distribution <int> VariantSelector;                        ///< Nuclei symmetry variants distribution

    std::uniform_int_distribution <int> PositionDistributionX;                  ///< Nuclei x-position distribution
    std::uniform_int_distribution <int> PositionDistributionY;                  ///< Nuclei y-position distribution
    std::uniform_int_distribution <int> PositionDistributionZ;                  ///< Nuclei z-position distribution

    std::uniform_real_distribution <double> OrientationDistributionQ;           ///< Nuclei orientation distribution in 3D (quaternion)
    std::uniform_real_distribution <double> OrientationDistributionA;           ///< Nuclei orientation distribution in 2D (Euler angle)

    std::normal_distribution <double> SizeDistributionNormal;                   ///< Normal nuclei size distribution
    std::cauchy_distribution <double> SizeDistributionCauchy;                   ///< Cauchy nuclei size distribution
    std::uniform_real_distribution <double> SizeDistributionUniform;            ///< Uniform nuclei size distribution

    std::vector <Nucleus> GeneratedNuclei;                                      ///< Generated nuclei storage
    std::vector <Nucleus> PlantedNuclei;                                        ///< Planted nuclei storage

    NucleationParameters():
        Allowed(false),
        Generated(false),
        TrueRadius(false),
        RelativeDensity(false),
        NucleiPhase(0),
        MatrixPhase(0),
        Tmin(0.0),
        Tmax(0.0),
        Density(0.0),
        Nsites(0ul),
        Nseeds(0ul),
        Nvariants(0ul),
        ShieldingRadius(0.0),
        SeedRadius(1.0),
        SeedRadiusMIN(1.0),
        SeedRadiusMAX(1.0),
        LocationMode(NucleiLocationModes::Bulk),
        OrientationMode(NucleiOrientationModes::Reference),
        Distribution(NucleiSizeDistributions::Normal),
        VariantsMode(NucleiVariantsModes::Random){};

    double     SetSeedRadius();
    iVector3   SetSeedPosition();
    Quaternion SetSeedOrientation();
    bool       CheckLocation(const PhaseField& Phase,
                             const iVector3 loc_position) const;                ///< Returns "true" if seed position fulfills nucleation constraints, "false" otherwise
    bool       PositionNotShielded(const iVector3 position) const;              ///< Returns "true" if seed position is outside of the "shielding" radius, "false" otherwise.
};

class OP_EXPORTS Nucleation : public OPObject                                   ///< Handles the nucleation of new phases and grains
{
  public:

    Nucleation(){};
    Nucleation(Settings& locSettings, const std::string InputFileName = DefaultInputFileName);
    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input values from file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input values from file
    void MoveFrame(const int dX, const int dY, const int dZ,
                   const BoundaryConditions& BC) override;                      ///< Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondingly.

    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC) override; ///< Changes the mesh size while keeping the data.

    void GenerateNucleationSites(PhaseField& Phase, Temperature& Tx);           ///< Generates nuclei
    void PlantNuclei(PhaseField& Phi, int tstep);                               ///< Plants generated nuclei according to their nucleation parameters
    void CheckNuclei(PhaseField& Phi, InterfaceProperties& IP, DrivingForce& dG, int tstep); ///< Checks planted nuclei, removes unstable ones //TODO update description method is const

    void WriteStatistics(const Settings& locSettings, const int tStep) const;   ///< Writes nucleation statistics to a file
    void PrintStatistics(void);                                                 ///< Prints nucleation statistics to screen.

    bool Read(const Settings& locSettings,
              const BoundaryConditions& BC, const int tStep) override;          ///< Read stored nucleation information from a file
    bool Write(const Settings& locSettings, const int tStep) const override;    ///< Write nucleation information to a file

    void ReadH5(H5Interface& H5, const int tStep);                              ///< Read stored nucleation information from a file
    void WriteH5(H5Interface& H5, const int tStep);                             ///< Write nucleation information to a file
    
    int NucleateEvery;                                                          ///< How frequently nucleation attempts should be performed

  private:
    size_t Nphases;                                                             ///< Number of thermodynamic phases
    std::vector<size_t> Nvariants;                                              ///< Number of crystallographic (symmetry/translation/...) variants
    std::vector<AggregateStates> PhaseAggregateStates;                          ///< Aggregate states of all phases

    GridParameters Grid;                                                        ///< Simulation grid parameters

    int NumberOfAttempts;                                                       ///< Maximum number of attempts to generate each nucleation site
    bool NucleiPlanted;                                                         ///< True if PlantNuclei() has been called. It is reset to false in CheckNuclei()

    Matrix< NucleationParameters > Parameters;                                  ///< Stores nucleation parameters for each pair of phases (A in B)

    void CalculateNumberOfSeeds(PhaseField& Phase, Temperature& Tx);            ///< Generates random nucleation sites and gives them weights according to the chosen distribution parameters
    void GenerateSeeds(PhaseField& Phase, Temperature& Tx);                     ///< Assigns seeds properties to generated nuclei
    void Clear();                                                               ///< Clears the seeds storage
    void ClearOutOfRange(const Temperature& Tx);                                ///< Clears the seeds storage if phase pair goes out of the nucleation range
    void SeedRandomGenerators(int RandomNumberSeedInput);                       ///< Generates new seeds for random number generators
};

}// namespace openphase

#endif // NUCLEATION_H
