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

#include "Nucleation.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperature.h"
#include "DrivingForce.h"
#include "InterfaceProperties.h"
#include "Orientations.h"
#include "ConsoleOutput.h"
#include "SymmetryVariants.h"

namespace openphase
{
using namespace std;

double NucleationParameters::SetSeedRadius()
{
    double radius = 0.0;

    switch(Distribution)
    {
        case NucleiSizeDistributions::Normal:
        {
            radius = 0.5*SizeDistributionNormal(SizeGenerator);
            break;
        }
        case NucleiSizeDistributions::Cauchy:
        {
            radius = 0.5*SizeDistributionCauchy(SizeGenerator);
            break;
        }
        case NucleiSizeDistributions::Uniform:
        {
            radius = 0.5*SizeDistributionUniform(SizeGenerator);
            break;
        }
        case NucleiSizeDistributions::FixedRadius:
        case NucleiSizeDistributions::None:
        {
            radius = SeedRadius;
            break;
        }
        default:
        {
            radius = 0.0;
            break;
        }
    }

    if (radius == 0.0)
    {
        ConsoleOutput::WriteWarning("Seed radius has been set to zero. Check input Parameters!", "NucleationParametersEXP", "SetSeedRadius");
    }
    return radius;
}
iVector3 NucleationParameters::SetSeedPosition()
{
    iVector3 position{-1,-1,-1};

    if(LocationMode == NucleiLocationModes::XBottom)
    {
        position[0] = 0;
    }
    else if(LocationMode == NucleiLocationModes::XTop)
    {
        position[0] = PositionDistributionX.b() - 1;
    }
    else
    {
        position[0] = PositionDistributionX(PositionGeneratorX);
    }

    if(LocationMode == NucleiLocationModes::YBottom)
    {
        position[1] = 0;
    }
    else if(LocationMode == NucleiLocationModes::YTop)
    {
        position[1] = PositionDistributionY.b() - 1;
    }
    else
    {
        position[1] = PositionDistributionY(PositionGeneratorY);
    }

    if(LocationMode == NucleiLocationModes::ZBottom)
    {
        position[2] = 0;
    }
    else if(LocationMode == NucleiLocationModes::ZTop)
    {
        position[2] = PositionDistributionZ.b() - 1;
    }
    else
    {
        position[2] = PositionDistributionZ(PositionGeneratorZ);
    }
    return position;
}

Quaternion NucleationParameters::SetSeedOrientation()
{
    Quaternion loc_orientation;

    int active_dimensions = (PositionDistributionX.b() != 0)
                          + (PositionDistributionY.b() != 0)
                          + (PositionDistributionZ.b() != 0);

    if(OrientationMode == NucleiOrientationModes::Random)
    {
        switch(active_dimensions)
        {
            case 1: // No rotation in 1D
            {
                loc_orientation.set(1.0, 0.0, 0.0, 0.0);
                break;
            }
            case 2: // Only in-plain rotation in 2D
            {
                double a1 = 0.0;
                double a2 = 0.0;
                double a3 = 0.0;

                if(PositionDistributionX.b() == 0) a1 = OrientationDistributionA(OrientationGenerator1);
                if(PositionDistributionY.b() == 0) a2 = OrientationDistributionA(OrientationGenerator2);
                if(PositionDistributionZ.b() == 0) a3 = OrientationDistributionA(OrientationGenerator3);

                EulerAngles ph1({a1,a2,a3},XYZ);

                loc_orientation = ph1.getQuaternion().normalized();
                break;
            }
            case 3: // Full rotation freedom in 3D
            {
                double u1 = OrientationDistributionQ(OrientationGenerator1);
                double u2 = OrientationDistributionQ(OrientationGenerator2);
                double u3 = OrientationDistributionQ(OrientationGenerator3);

                Quaternion locQuaternion;
                locQuaternion.set(sqrt(1.0-u1)*sin(2.0*Pi*u2),
                                  sqrt(1.0-u1)*cos(2.0*Pi*u2),
                                  sqrt(u1)*sin(2.0*Pi*u3),
                                  sqrt(u1)*cos(2.0*Pi*u3));

                loc_orientation = locQuaternion.normalized();
                break;
            }
            default: // No rotation in zero dimensions.
            {
                loc_orientation.set(1.0, 0.0, 0.0, 0.0);
                break;
            }
        }
    }
    return loc_orientation;
}

bool NucleationParameters::CheckLocation(const PhaseField& Phase,
                                            const iVector3 loc_position) const
{
    bool enable_seed = false;
    if (Phase.Grid.PositionInLocalBounds(loc_position) and
        PositionNotShielded(loc_position))
    {
        int i = loc_position.get_x() - Phase.Grid.OffsetX;
        int j = loc_position.get_y() - Phase.Grid.OffsetY;
        int k = loc_position.get_z() - Phase.Grid.OffsetZ;

        switch (LocationMode)
        {
            case NucleiLocationModes::XBottom:
            case NucleiLocationModes::XTop:
            case NucleiLocationModes::YBottom:
            case NucleiLocationModes::YTop:
            case NucleiLocationModes::ZBottom:
            case NucleiLocationModes::ZTop:
            case NucleiLocationModes::BulkAndInterfaces:
            {
                if(Phase.Fractions(i,j,k,{MatrixPhase}) > 0.25)
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::Bulk:
            {
                if( !Phase.Fields(i,j,k).interface() and
                    (Phase.Fractions(i,j,k,{MatrixPhase}) == 1.0))
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::BulkAndGrainBoundaries:
            {
                if(Phase.Fractions(i,j,k,{MatrixPhase}) == 1.0)
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::GrainBoundaries:
            {
                if( Phase.Fields(i,j,k).interface() and
                   (Phase.Fractions(i,j,k,{MatrixPhase}) == 1.0))
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::PhaseBoundaries:
            {
                if( Phase.Fields(i,j,k).interface() and
                   (Phase.Fractions(i,j,k,{MatrixPhase}) < 1.0) and
                   (Phase.Fractions(i,j,k,{MatrixPhase}) > 0.25))
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::Junctions:
            {
                if((Phase.Fractions(i,j,k,{MatrixPhase}) == 1.0) and
                   (Phase.Fields(i,j,k).size() > 2))
                {
                    enable_seed = true;
                }
                break;
            }
            case NucleiLocationModes::Interfaces:
            {
                if( Phase.Fields(i,j,k).interface() and
                   (Phase.Fractions(i,j,k,{MatrixPhase}) > 0.25))
                {
                    enable_seed = true;
                }
                break;
            }
        }
    }
    return enable_seed;
}

bool NucleationParameters::PositionNotShielded(const iVector3 position) const
{
    bool my_return = true;
    double ShieldingRadiusSquare = ShieldingRadius*ShieldingRadius;

    long int TotalNx = PositionDistributionX.b();
    long int TotalNy = PositionDistributionY.b();
    long int TotalNz = PositionDistributionZ.b();

    for (auto it = GeneratedNuclei.begin(); it != GeneratedNuclei.end(); ++it)
    {
        int x = position.get_x();
        int y = position.get_y();
        int z = position.get_z();
        const double xdis = std::min(std::fabs(x - it->position[0]), std::min( std::fabs(x - it->position[0] + TotalNx), std::fabs(x - it->position[0] - TotalNx)));
        const double ydis = std::min(std::fabs(y - it->position[1]), std::min( std::fabs(y - it->position[1] + TotalNy), std::fabs(y - it->position[1] - TotalNy)));
        const double zdis = std::min(std::fabs(z - it->position[2]), std::min( std::fabs(z - it->position[2] + TotalNz), std::fabs(z - it->position[2] - TotalNz)));
        const double distanceSquare = xdis*xdis + ydis*ydis + zdis*zdis;
        if(distanceSquare < ShieldingRadiusSquare)
        {
            my_return = false;
        }
    }
    return my_return;
}

//=============================================================================

Nucleation::Nucleation(Settings& locSettings, const std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Nucleation::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "Nucleation";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Grid = locSettings.Grid;

    Nphases = locSettings.Nphases;
    Nvariants = locSettings.Nvariants;

    PhaseAggregateStates = locSettings.PhaseAggregateStates;

    for(size_t n = 0; n < Nvariants.size(); n++)
    {
        Nvariants[n] = max(Nvariants[n], (size_t)1);
    }

    Parameters.Allocate(Nphases, Nphases);

    NucleateEvery = 1;
    NumberOfAttempts = 100;
    NucleiPlanted = false;

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    {
        Parameters(n, m).Allowed = false;
        Parameters(n, m).Generated = false;
        Parameters(n, m).NucleiPhase = n;
        Parameters(n, m).MatrixPhase = m;
        Parameters(n, m).Nsites  = 0;
        Parameters(n, m).Nseeds  = 0;
        Parameters(n, m).Density = 0;
        Parameters(n, m).PositionDistributionX.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNx - 1));
        Parameters(n, m).PositionDistributionY.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNy - 1));
        Parameters(n, m).PositionDistributionZ.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNz - 1));

        Parameters(n, m).OrientationDistributionQ.param(uniform_real_distribution<double>::param_type(0, 1));
        Parameters(n, m).OrientationDistributionA.param(uniform_real_distribution<double>::param_type(0.0, 2.0*Pi));
    }

    locSettings.AddForRemeshing(*this);
    locSettings.AddForReading(*this);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void Nucleation::ReadInput(const std::string InputFileName)
{
    ConsoleOutput::WriteLineInsert("Nucleation input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in | ios_base::binary);

    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    };
    std::stringstream data;
    data << inp.rdbuf();
    ReadInput(data);

    inp.close();
}

void Nucleation::ReadInput(std::stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    if(moduleLocation == -1)
    {
        ConsoleOutput::WriteWarning("Input for module \"" + thisclassname + "\" is not found in the input file! Nucleation is disabled!",thisclassname, "ReadInput()");
    }
    else
    {
        const long int RandomNumberSeedInput = FileInterface::ReadParameterI(inp, moduleLocation, "RandomNumberSeed", false, 1);
        SeedRandomGenerators(RandomNumberSeedInput);

        NucleateEvery = FileInterface::ReadParameterI(inp, moduleLocation, "NucleateEvery", false, 1);
        NumberOfAttempts = FileInterface::ReadParameterI(inp, moduleLocation, "NumberOfAttempts", false, 100);

        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = 0; m < Nphases; m++)
        {
            stringstream converter;
            converter << n << "_" << m;
            string counter = converter.str();

            Parameters(n, m).Allowed = FileInterface::ReadParameterB(inp, moduleLocation, string("Allowed_") + counter, false, "NO");

            if(Parameters(n, m).Allowed)
            {
                string location_input = FileInterface::ReadParameterK(inp, moduleLocation, string("Location_") + counter, false, "BULK");

                bool valid_location_input = false;

                // Setting nucleation sites location
                if(location_input == "BULK")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::Bulk;
                    valid_location_input = true;
                }
                if(location_input == "BULKANDGRAINBOUNDARIES")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::BulkAndGrainBoundaries;
                    valid_location_input = true;
                }
                if(location_input == "BULKANDINTERFACES")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::BulkAndInterfaces;
                    valid_location_input = true;
                }
                if(location_input == "GB" or location_input == "GRAINBOUNDARIES")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::GrainBoundaries;
                    valid_location_input = true;
                }
                if(location_input == "PHASEBOUNDARIES")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::PhaseBoundaries;
                    valid_location_input = true;
                }
                if(location_input == "INTERFACE" or location_input == "INTERFACES")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::Interfaces;
                    valid_location_input = true;
                }
                if(location_input == "JUNCTIONS")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::Junctions;
                    valid_location_input = true;
                }
                if(location_input == "XBOTTOM")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::XBottom;
                    valid_location_input = true;
                }
                if(location_input == "XTOP")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::ZBottom;
                    valid_location_input = true;
                }
                if(location_input == "YBOTTOM")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::YBottom;
                    valid_location_input = true;
                }
                if(location_input == "YTOP")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::YTop;
                    valid_location_input = true;
                }
                if(location_input == "ZBOTTOM")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::ZBottom;
                    valid_location_input = true;
                }
                if(location_input == "ZTOP")
                {
                    Parameters(n, m).LocationMode = NucleiLocationModes::ZTop;
                    valid_location_input = true;
                }
                if(valid_location_input)
                {
                    // Reading seeds number or density
                    if (FileInterface::ReadParameterS(inp, moduleLocation, string("IMode_") + counter, false, "SITES") == "SITES")
                    {
                        Parameters(n, m).Nsites  = FileInterface::ReadParameterI(inp, moduleLocation, string("Nsites_") + counter, false, 0);
                    }
                    if(Parameters(n, m).Nsites == 0)
                    {
                        Parameters(n, m).Density = FileInterface::ReadParameterD(inp, moduleLocation, string("Density_") + counter, true, 0.0);

                        Parameters(n, m).RelativeDensity = FileInterface::ReadParameterB(inp, moduleLocation, {string("RelativeDensity_") + counter, string("RelDensity_") + counter}, false, true);

                        if(Parameters(n, m).Density == 0.0)
                        {
                            ConsoleOutput::WriteWarning("The nucleation sites Nsites_" + counter +
                                               " and the nucleation density Density_" + counter +
                                               " for phase pair " + counter + " are zero!\n" +
                                               "No nuclei will be generated for this phase pair!",
                                               thisclassname, "ReadInput()");
                        }
                    }

                    // Reading nucleation temperature range
                    Parameters(n, m).Tmin = FileInterface::ReadParameterD(inp, moduleLocation, string("Tmin_") + counter);
                    Parameters(n, m).Tmax = FileInterface::ReadParameterD(inp, moduleLocation, string("Tmax_") + counter);

                    Parameters(n, m).ShieldingRadius = FileInterface::ReadParameterD(inp, moduleLocation, string("ShieldingRadius_") + counter, false, 5.0);

                    // Reading seeds size distribution
                    string distribution_input = "NONE";
                    if(PhaseAggregateStates[n] == AggregateStates::Solid)
                    {
                        distribution_input = FileInterface::ReadParameterK(inp, moduleLocation, string("Distribution_") + counter, false, "NONE");
                    }
                    bool distribution_set = false;
                    if(distribution_input == "NONE")
                    {
                        Parameters(n, m).Distribution = NucleiSizeDistributions::None;
                        Parameters(n, m).SeedRadius   = 0.5*Grid.iWidth*Grid.dx;
                        Parameters(n, m).SeedRadiusMIN = 0;
                        Parameters(n, m).SeedRadiusMAX = Grid.iWidth*Grid.dx;
                        distribution_set = true;
                    }
                    if(distribution_input == "FIXEDSIZE" or distribution_input == "FIXEDRADIUS")
                    {
                        Parameters(n, m).Distribution = NucleiSizeDistributions::FixedRadius;
                        Parameters(n, m).SeedRadius   = FileInterface::ReadParameterD(inp, moduleLocation, string("SeedRadius_") + counter, true, 0.0);
                        Parameters(n, m).SeedRadiusMIN = 0;
                        Parameters(n, m).SeedRadiusMAX = (Grid.Nx+Grid.Ny+Grid.Nz)*Grid.dx;
                        distribution_set = true;
                    }
                    if(distribution_input == "NORMAL")
                    {
                        Parameters(n, m).Distribution = NucleiSizeDistributions::Normal;
                        double mu    = FileInterface::ReadParameterD(inp, moduleLocation, string("Center_") + counter, true, 0.0);
                        double sigma = FileInterface::ReadParameterD(inp, moduleLocation, string("Deviation_") + counter, true, 0.0);
                        Parameters(n, m).SeedRadiusMIN = std::max(0.0, mu - 3.0*sigma);
                        Parameters(n, m).SeedRadiusMAX = mu + 3.0*sigma;
                        Parameters(n, m).SizeDistributionNormal.param(normal_distribution<double>::param_type(mu, sigma));
                        distribution_set = true;
                    }
                    if(distribution_input == "CAUCHY")
                    {
                        Parameters(n, m).Distribution = NucleiSizeDistributions::Cauchy;
                        double mu    = FileInterface::ReadParameterD(inp, moduleLocation, string("Center_") + counter, true, 0.0);
                        double sigma = FileInterface::ReadParameterD(inp, moduleLocation, string("HalfWidth_") + counter, true, 0.0);
                        Parameters(n, m).SeedRadiusMIN = std::max(0.0, mu - 10.0*sigma);
                        Parameters(n, m).SeedRadiusMAX = mu + 10.0*sigma;
                        Parameters(n, m).SizeDistributionCauchy.param(cauchy_distribution<double>::param_type(mu, sigma));
                        distribution_set = true;
                    }
                    if(distribution_input == "UNIFORM")
                    {
                        Parameters(n, m).Distribution = NucleiSizeDistributions::Uniform;
                        Parameters(n, m).SeedRadiusMIN = FileInterface::ReadParameterD(inp, moduleLocation,string("SeedRadiusMIN_") + counter, true, 0.0);
                        Parameters(n, m).SeedRadiusMAX = FileInterface::ReadParameterD(inp, moduleLocation,string("SeedRadiusMAX_") + counter, true, 0.0);
                        Parameters(n, m).SizeDistributionUniform.param(uniform_real_distribution<double>::param_type(Parameters(n, m).SeedRadiusMIN, Parameters(n, m).SeedRadiusMAX));
                        distribution_set = true;

                        if (Parameters(n, m).SeedRadiusMIN >= Parameters(n, m).SeedRadiusMAX)
                        {
                            ConsoleOutput::WriteExit("Wrong input! Minimum seed radius must be smaller than maximum nucleation radius", thisclassname, "ReadInput()");
                            OP_Exit(EXIT_FAILURE);
                        }
                        if (Parameters(n, m).SeedRadiusMIN < 0.0 or Parameters(n, m).SeedRadiusMAX < 0.0)
                        {
                            ConsoleOutput::WriteExit("Wrong input! Minimum and Maximum seed radius must be positive", thisclassname, "ReadInput()");
                            OP_Exit(EXIT_FAILURE);
                        }
                    }

                    Parameters(n, m).TrueRadius = FileInterface::ReadParameterB(inp, moduleLocation,string("TrueRadius_") + counter, false, false);

                    if(!distribution_set)
                    {
                        //Wrong seed size distribution input
                        string message  = "Wrong or no input for the seeds size distribution!";
                        ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
                        OP_Exit(EXIT_FAILURE);
                    }

                    // Reading Nuclei Orientation Mode
                    string orientations_input = "REFERENCE";
                    if(PhaseAggregateStates[n] == AggregateStates::Solid)
                    {
                        orientations_input = FileInterface::ReadParameterK(inp,moduleLocation, {string("OrientationMode_") + counter, string("Orientation_") + counter});
                    }

                    bool orientations_set = false;
                    if(orientations_input == "RANDOM")
                    {
                        Parameters(n, m).OrientationMode = NucleiOrientationModes::Random;
                        orientations_set = true;
                    }
                    if(orientations_input == "PARENT")
                    {
                        Parameters(n, m).OrientationMode = NucleiOrientationModes::Parent;
                        orientations_set = true;
                    }
                    if(orientations_input == "REFERENCE")
                    {
                        Parameters(n, m).OrientationMode = NucleiOrientationModes::Reference;
                        orientations_set = true;
                    }
                    if(!orientations_set)
                    {
                        //Wrong seed orientation input
                        string message  = "Wrong or no input for the seeds orientation!";
                        ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
                        OP_Exit(EXIT_FAILURE);
                    }

                    if(Nvariants[n] > 1)
                    {
                        Parameters(n, m).Nvariants = FileInterface::ReadParameterI(inp, moduleLocation, string("Variants_") + counter, false, 1);

                        string variants_mode = FileInterface::ReadParameterK(inp,moduleLocation, string("VariantsMode_")+counter, false, "RANDOM");
                        bool variants_mode_set = false;
                        if(variants_mode == "RANDOM")
                        {
                            Parameters(n, m).VariantsMode = NucleiVariantsModes::Random;
                            Parameters(n, m).VariantSelector.param(uniform_int_distribution<int>::param_type(0, Nvariants[n] - 1));
                            variants_mode_set = true;
                        }
                        if(variants_mode == "LOWESTENERGY")
                        {
                            Parameters(n, m).VariantsMode = NucleiVariantsModes::LowestEnergy;
                            variants_mode_set = true;
                        }
                        if(!variants_mode_set)
                        {
                            //Wrong variants mode input
                            string message  = "Wrong or no input for the variants mode!";
                            ConsoleOutput::WriteExit(message, thisclassname, "ReadInput()");
                            OP_Exit(EXIT_FAILURE);
                        }
                    }
                }
                else
                {
                    ConsoleOutput::WriteExit("Nucleation mode for phase pair " + counter +
                                    " could not be read", thisclassname, "ReadInput()");
                    OP_Exit(EXIT_NUCMODE_ERROR);
                }
            }
        }
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void Nucleation::Remesh(int newNx, int newNy, int newNz,
                        const BoundaryConditions& BC)
{
    Grid.SetDimensions(newNx, newNy, newNz);

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed)
    {
        Parameters(n, m).PositionDistributionX.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNx - 1));
        Parameters(n, m).PositionDistributionY.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNy - 1));
        Parameters(n, m).PositionDistributionZ.param(uniform_int_distribution<int>::param_type(0, Grid.TotalNz - 1));

        if(Parameters(n, m).Generated)
        for(auto it  = Parameters(n, m).GeneratedNuclei.begin();
                 it != Parameters(n, m).GeneratedNuclei.end(); )
        {
            if(!Grid.PositionInTotalBounds(it->position))
            {
                it = Parameters(n, m).GeneratedNuclei.erase(it);
            }
            else
            {
                it++;
            }
        }
    }

    ConsoleOutput::WriteStandard(thisclassname, "Remeshed");
}

void Nucleation::SeedRandomGenerators(int RandomNumberSeedInput)
{
    const size_t RandomNumberSeed = (RandomNumberSeedInput < 0)
               ? static_cast<size_t>(std::chrono::system_clock::now().time_since_epoch().count())
               : static_cast<size_t>(RandomNumberSeedInput);

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    //if(Parameters(n, m).Allowed)
    {
        size_t seed_multiplier = n + m*Nphases + 1;

        Parameters(n, m).SizeGenerator.seed(34784u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).VariantsGenerator.seed(5467u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).PositionGeneratorX.seed(253u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).PositionGeneratorY.seed(4958u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).PositionGeneratorZ.seed(54861u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).OrientationGenerator1.seed(45u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).OrientationGenerator2.seed(697u*seed_multiplier*RandomNumberSeed);
        Parameters(n, m).OrientationGenerator3.seed(2597u*seed_multiplier*RandomNumberSeed);
    }
}

void Nucleation::Clear()
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Generated)
    {
        Parameters(n, m).GeneratedNuclei.clear();
        Parameters(n, m).PlantedNuclei.clear();
    }
}

void Nucleation::ClearOutOfRange(const Temperature& Tx)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Generated and
       !(Parameters(n, m).Tmin <= Tx.Tmax and
         Parameters(n, m).Tmax >= Tx.Tmin))
    {
        Parameters(n, m).GeneratedNuclei.clear();
        //Parameters(n, m).NucleatedParticles.clear();
        Parameters(n, m).Generated = false;
    }
}

void Nucleation::MoveFrame(const int dX, const int dY, const int dZ,
                           const BoundaryConditions& BC)
{
    iVector3 frame_shift({dX,dY,dZ});

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Generated)
    {
        for (auto it  = Parameters(n, m).GeneratedNuclei.begin();
                  it != Parameters(n, m).GeneratedNuclei.end(); )
        {
            it->position += frame_shift;
            if(!(Grid.PositionInTotalBounds(it->position)))
            {
                it = Parameters(n, m).GeneratedNuclei.erase(it);
            }
            else
            {
                it++;
            }
        }
        for (auto it  = Parameters(n, m).PlantedNuclei.begin();
                  it != Parameters(n, m).PlantedNuclei.end(); it++)
        {
            it->position += frame_shift;
        }
    }
}

void Nucleation::GenerateNucleationSites(PhaseField& Phase, Temperature& Tx)
{
    CalculateNumberOfSeeds(Phase, Tx);
    GenerateSeeds(Phase, Tx);
    ClearOutOfRange(Tx);
}

void Nucleation::CalculateNumberOfSeeds(PhaseField& Phase, Temperature& Tx)
{
    double TotalVolume = Grid.TotalNumberOfCells();
    double RealUnitsVolume = TotalVolume * pow(Grid.dx, Grid.Active());

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Tmin <= Tx.Tmax and
       Parameters(n, m).Tmax >= Tx.Tmin)
    {
        double PhaseFractions_m = 0.0;

        if(Parameters(n, m).RelativeDensity)
        {
            for(size_t idx = 0; idx < Phase.FieldsProperties.size(); idx++)
            if(Phase.FieldsProperties[idx].Phase == m)
            {
                PhaseFractions_m += Phase.FieldsProperties[idx].Volume;
            }
            PhaseFractions_m /= TotalVolume;
        }
        else
        {
            PhaseFractions_m = 1.0;
        }

        if(Parameters(n, m).Density != 0.0)
        {
            double dNsites = Parameters(n, m).Density*RealUnitsVolume*PhaseFractions_m;
            if (dNsites > double(std::numeric_limits<size_t>::max()))
            {
                ConsoleOutput::WriteExit("Too large particles density!", thisclassname, "CalculateNumberOfSeeds");
                OP_Exit(EXIT_FAILURE);
            }
            assert(dNsites >= 0.0);
            size_t locNsites = std::floor(dNsites);

            if(locNsites > Parameters(n, m).Nsites)
            {
                Parameters(n, m).Nsites = locNsites;
                Parameters(n, m).Generated = false;
            }
        }

//        for (auto it  = Parameters(n,m).GeneratedParticles.begin();
//                  it != Parameters(n,m).GeneratedParticles.end(); )
//        {
//            if (it->planted)
//            {
//                //Parameters(n,m).NucleatedParticles.push_back(*it);
//                it = Parameters(n,m).GeneratedParticles.erase(it);
//            }
//            else
//            {
//                it++;
//            }
//        }
//
        Parameters(n,m).Nseeds = Parameters(n,m).GeneratedNuclei.size();// + Parameters(n,m).NucleatedParticles.size();
    }
}

void Nucleation::GenerateSeeds(PhaseField& Phase, Temperature& Tx)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       !Parameters(n, m).Generated and
       Parameters(n, m).Tmin <= Tx.Tmax and
       Parameters(n, m).Tmax >= Tx.Tmin)
    {
        Parameters(n, m).Generated = true;

        int attempts = 0;

        bool PrintMsg2 = false;
        while (Parameters(n, m).Nseeds < Parameters(n, m).Nsites and
               attempts < NumberOfAttempts)
        {
            bool enable_seed = false;

            int increment_part = 0;
            int increment_attempts = 0;

            double loc_radius = Parameters(n,m).SetSeedRadius();
            if(loc_radius >= Parameters(n,m).SeedRadiusMIN and
               loc_radius <= Parameters(n,m).SeedRadiusMAX)
            {
                iVector3 loc_position = Parameters(n,m).SetSeedPosition();

                enable_seed = Parameters(n,m).CheckLocation(Phase, loc_position);

    #ifdef MPI_PARALLEL
                int tmp_enable_seed = enable_seed;
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &tmp_enable_seed, 1, OP_MPI_INT, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                enable_seed = tmp_enable_seed;
    #endif
                if(enable_seed)
                {
                    Nucleus locSeed;
                    locSeed.position    = loc_position;
                    locSeed.radius      = loc_radius;
                    locSeed.orientation = Parameters(n, m).SetSeedOrientation();

                    std::stringstream message;
                    message << "Nucleation: Generated seed particle " << Parameters(n, m).Nseeds << " at ["
                            << loc_position.get_x() << ", "
                            << loc_position.get_y() << ", "
                            << loc_position.get_z() << "] and effective radius of "
                            << loc_radius;
                    ConsoleOutput::WriteSimple(message.str());

                    Parameters(n, m).GeneratedNuclei.push_back(locSeed);
                    increment_part = 1;
                    attempts = 0;
                    PrintMsg2 = true;
                }
                else // if no seed location found
                {
                    increment_attempts = 1;
                }
            }
            else // if radius outside [R_min,R_max]
            {
                increment_attempts = 1;
            }
            Parameters(n, m).Nseeds += increment_part;
            attempts += increment_attempts;
        } // end while loop
        if (PrintMsg2)
        {
            std::stringstream message2;
            message2 << "Nucleation: Generated " << Parameters(n, m).Nseeds
                     << " nucleation sites (of " << Parameters(n, m).Nsites
                     << ") for phase " << n
                     << " in phase " << m << ".";
            ConsoleOutput::Write(message2.str());
        }
        if ((Parameters(n, m).Nseeds == 0) and (Parameters(n, m).Nsites != 0))
        {
            std::stringstream message;
            message << "No seeds has been generated for " << Parameters(n, m).Nsites
                    << " nucleation sites for " << Phase.PhaseNames[n] << " in " << Phase.PhaseNames[m] << "\n";
            ConsoleOutput::WriteWarning(message.str(), thisclassname, "GenerateSeed");
        }
    }
}

void Nucleation::WriteStatistics(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif

    string FileName = locSettings.TextDir + dirSeparator + "NucleationStatistics.dat";

    fstream NucFile;

    NucFile.open(FileName, ios::out);

    NucFile.setf(ios::left);
    NucFile << setw(10) << "TimeStep";
    NucFile.setf(ios::right);
    NucFile << setw(10) << "PFindex"
            << setw(10) << "Nphase"
            << setw(10) << "Mphase"
            << setw(10) << "Variant"
            << setw(10) << "X"
            << setw(10) << "Y"
            << setw(10) << "Z"
            << setw(10) << "Q1"
            << setw(10) << "Q2"
            << setw(10) << "Q3"
            << setw(10) << "Q4"
            << setw(14) << "dGnuc"
            << setw(14) << "dGmin" << endl;

    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    //if(Parameters(n,m).Allowed)
    {
        for(auto it  = Parameters(n,m).PlantedNuclei.begin();
                 it != Parameters(n,m).PlantedNuclei.end(); it++)
        if(it->planted)
        {
            NucFile.setf(ios::left);
            NucFile << setw(10)  << setprecision(6) << it->time_stamp;
            NucFile.setf(ios::right);
            NucFile << setw(10)  << setprecision(6) << it->index
                    << setw(10)  << setprecision(6) << it->phase
                    << setw(10)  << setprecision(6) << it->parentPhase
                    << setw(10)  << setprecision(6) << it->variant
                    << setw(10)  << setprecision(6) << it->position[0]
                    << setw(10)  << setprecision(6) << it->position[1]
                    << setw(10)  << setprecision(6) << it->position[2]
                    << setw(10)  << setprecision(6) << it->orientation[0]
                    << setw(10)  << setprecision(6) << it->orientation[1]
                    << setw(10)  << setprecision(6) << it->orientation[2]
                    << setw(10)  << setprecision(6) << it->orientation[3]
                    << setw(14)  << setprecision(6) << it->dGnuc
                    << setw(14)  << setprecision(6) << it->dGmin << endl;
        }
    }
    NucFile.close();
#ifdef MPI_PARALLEL
    }
#endif
}

void Nucleation::PrintStatistics(void)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    if(Parameters(n, m).Allowed and
       Parameters(n, m).Generated)
    {
        if(Parameters(n, m).Nsites == 0)
        {
            std::stringstream message;
            message << "Nucleation: Too low particles density! No nucleation "
                    << "sites were generated for phase " << n << " in phase "
                    << m << ".";
            ConsoleOutput::Write(message.str());
        }
    }
}

void Nucleation::PlantNuclei(PhaseField& Phase, int tStep)
{
    for(size_t n = 0; n < Nphases; n++)
    for(size_t m = 0; m < Nphases; m++)
    for(auto ind  = Parameters(n, m).GeneratedNuclei.begin();
             ind != Parameters(n, m).GeneratedNuclei.end(); ind++)
    {
        if(!ind->planted)
        {
            // detect parent grain by max phase-field value
            size_t ParentGrainIndex = 0;
            if(Grid.PositionInLocalBounds(ind->position))
            {
                double locMax = 0.0;
                iVector3 locPosition = Grid.ConvertToLocal(ind->position);

                for(auto alpha  = Phase.Fields(locPosition[0], locPosition[1], locPosition[2]).begin();
                         alpha != Phase.Fields(locPosition[0], locPosition[1], locPosition[2]).end(); ++alpha)
                {
                    if(Phase.FieldsProperties[alpha->index].Phase == m and alpha->value > locMax)
                    {
                        locMax = alpha->value;
                        ParentGrainIndex = alpha->index;
                    }
                }
            }
#ifdef MPI_PARALLEL
            unsigned long locParentGrainIndex = ParentGrainIndex;
            OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locParentGrainIndex, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
            ParentGrainIndex = locParentGrainIndex;
#endif
            // select symmetry variants
            if(Nvariants[n] > 1)
            {
                vector<bool> locVariants(Nvariants[n], false);

                switch(Parameters(n,m).VariantsMode)
                {
                    case NucleiVariantsModes::Random:
                    {
                        size_t locNvariants = 0;
                        while(locNvariants < Parameters(n,m).Nvariants)
                        {
                            int locVariant = Parameters(n, m).VariantSelector(Parameters(n, m).VariantsGenerator);

                            if(!locVariants[locVariant])
                            {
                                locVariants[locVariant] = true;
                                locNvariants++;

                                size_t locIndex = Phase.PlantGrainNucleus(n, ind->position[0], ind->position[1], ind->position[2]);

                                if(Parameters(n, m).TrueRadius)
                                {
                                    Phase.FieldsProperties[locIndex].RefVolume = Phase.CalculateReferenceVolume(ind->radius/Grid.dx);
                                }

                                if(Parameters(n, m).OrientationMode == NucleiOrientationModes::Parent)
                                {
                                    Phase.FieldsProperties[locIndex].Orientation =
                                           Phase.FieldsProperties[ParentGrainIndex].Orientation;

                                    ind->orientation = Phase.FieldsProperties[ParentGrainIndex].Orientation;
                                }
                                else
                                {
                                    Phase.FieldsProperties[locIndex].Orientation = ind->orientation;
                                }

                                Phase.FieldsProperties[locIndex].Variant = locVariant;
                                Phase.FieldsProperties[locIndex].Parent  = ParentGrainIndex;
                                Phase.FieldsProperties[locIndex].RefVolume /= Parameters(n,m).Nvariants;
                            }
                        }
                        break;
                    }
                    case NucleiVariantsModes::LowestEnergy:
                    {
                        for(size_t locVariant = 0; locVariant < Nvariants[n]; locVariant++)
                        {
                            size_t locIndex = Phase.PlantGrainNucleus(n, ind->position[0], ind->position[1], ind->position[2]);

                            if(Parameters(n, m).TrueRadius)
                            {
                                Phase.FieldsProperties[locIndex].RefVolume = Phase.CalculateReferenceVolume(ind->radius/Grid.dx);
                            }

                            if(Parameters(n, m).OrientationMode == NucleiOrientationModes::Parent)
                            {
                                Phase.FieldsProperties[locIndex].Orientation =
                                       Phase.FieldsProperties[ParentGrainIndex].Orientation;

                                ind->orientation = Phase.FieldsProperties[ParentGrainIndex].Orientation;
                            }
                            else
                            {
                                Phase.FieldsProperties[locIndex].Orientation = ind->orientation;
                            }

                            Phase.FieldsProperties[locIndex].Variant = locVariant;
                            Phase.FieldsProperties[locIndex].Parent  = ParentGrainIndex;
                            Phase.FieldsProperties[locIndex].RefVolume /= Parameters(n,m).Nvariants;
                        }
                        break;
                    }
                }
            }
            else
            {
                size_t locIndex = Phase.PlantGrainNucleus(n, ind->position[0], ind->position[1], ind->position[2]);

                if(Parameters(n, m).TrueRadius)
                {
                    Phase.FieldsProperties[locIndex].RefVolume = Phase.CalculateReferenceVolume(ind->radius/Grid.dx);
                }

                if(Parameters(n, m).OrientationMode == NucleiOrientationModes::Parent)
                {
                    Phase.FieldsProperties[locIndex].Orientation =
                           Phase.FieldsProperties[ParentGrainIndex].Orientation;

                    ind->orientation = Phase.FieldsProperties[ParentGrainIndex].Orientation;
                }
                else
                {
                    Phase.FieldsProperties[locIndex].Orientation = ind->orientation;
                }

                Phase.FieldsProperties[locIndex].Parent = ParentGrainIndex;
            }
            //ind->planted = true;
            ind->time_stamp = tStep;

            NucleiPlanted = true;
        }
    }
}

void Nucleation::CheckNuclei(PhaseField& Phase, InterfaceProperties& IP, DrivingForce& dG, int tStep)
{
    if(NucleiPlanted)
    {
        NucleiPlanted = false;

        for(size_t n = 0; n < Nphases; n++)
        for(size_t m = 0; m < Nphases; m++)
        for(auto ind  = Parameters(n, m).GeneratedNuclei.begin();
                 ind != Parameters(n, m).GeneratedNuclei.end(); ind++)
        {
            if(!ind->planted)
            {
                /** Minimum driving force barrier the nucleus has to overcome.*/
                double dGmin = 0.0;
                if(ind->radius > 0.0)
                {
                    dGmin = 2.0*IP.InterfaceEnergy(n, m).MaxEnergy/ind->radius;
                }

                size_t locNvariants = 0;
                NodeA<double> dGloc;
                if(Grid.PositionInLocalBounds(ind->position))
                {
                    iVector3 locPosition = Grid.ConvertToLocal(ind->position);

                    for(auto it  = dG.Force(locPosition[0], locPosition[1], locPosition[2]).begin();
                             it != dG.Force(locPosition[0], locPosition[1], locPosition[2]).end(); )
                    {
                        bool remove = false;
                        if(Phase.FieldsProperties[it->indexA].Phase == n and
                           Phase.FieldsProperties[it->indexA].Stage == GrainStages::Seed)
                        {
                            if(it->average < dGmin)
                            {
                                remove = true;
                            }
                            else
                            {
                                dGloc.add_value(it->indexA, it->average);
                                locNvariants ++;
                            }
                        };

                        if(Phase.FieldsProperties[it->indexB].Phase == n and
                           Phase.FieldsProperties[it->indexB].Stage == GrainStages::Seed)
                        {
                            if(it->average > -dGmin)
                            {
                                remove = true;
                            }
                            else
                            {
                                dGloc.add_value(it->indexB, -it->average);
                                locNvariants ++;
                            }
                        }
                        if(Phase.FieldsProperties[it->indexA].Phase == n and
                          (Phase.FieldsProperties[it->indexA].Stage == GrainStages::Seed or
                           Phase.FieldsProperties[it->indexB].Stage == GrainStages::Seed) and
                           Phase.FieldsProperties[it->indexB].Phase == n)
                        {
                            remove = true;
                        };

                        if(remove)
                        {
                            it = dG.Force(locPosition[0], locPosition[1], locPosition[2]).erase(it);
                        }
                        else
                        {
                            it++;
                        }
                    }
                }

    #ifdef MPI_PARALLEL
                unsigned long tmp_locNvariants = locNvariants;
                OP_MPI_Allreduce(OP_MPI_IN_PLACE, &tmp_locNvariants, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                locNvariants = tmp_locNvariants;
    #endif
                if(!locNvariants)
                {
                    ind->planted = false;
                }
                else
                {
                    ind->planted = true;

                    vector<size_t> PFidx(locNvariants,0);// Phase-field indices
                    vector<size_t> Vidx(locNvariants,0); // Variants indices
                    vector<double> dGmax(locNvariants,0.0);// Driving forces at the nucleation event

                    if(Grid.PositionInLocalBounds(ind->position))
                    {
                        switch(Parameters(n,m).VariantsMode)
                        {
                            case NucleiVariantsModes::Random:
                            {
                                /* Nothing to do here.*/
                                break;
                            }
                            case NucleiVariantsModes::LowestEnergy:
                            {
                                for(auto it1  = dGloc.begin(); it1 != dGloc.end(); it1++)
                                for(auto it2  = it1 + 1; it2 != dGloc.end(); it2++)
                                if(it1->value < it2->value)
                                {
                                    double tmp_value = it1->value;
                                    size_t tmp_idx = it1->index;

                                    it1->value = it2->value;
                                    it1->index = it2->index;

                                    it2->value = tmp_value;
                                    it2->index = tmp_idx;
                                }

                                while(dGloc.size() > Parameters(n,m).Nvariants)
                                {
                                    dGloc.erase(dGloc.end()-1);
                                }

                                break;
                            }
                        }
                        size_t idx = 0;
                        for(auto it = dGloc.begin(); it != dGloc.end(); it++)
                        {
                            PFidx[idx] = it->index;
                            Vidx[idx] = Phase.FieldsProperties[it->index].Variant;
                            dGmax[idx] = it->value;
                            idx++;
                        }
                    }
                    locNvariants = dGloc.size();
    #ifdef MPI_PARALLEL
                    unsigned long tmpNvariants = locNvariants;
                    OP_MPI_Allreduce(OP_MPI_IN_PLACE, &tmpNvariants, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                    locNvariants = tmpNvariants;

                    for(size_t it = 0; it < locNvariants; it++)
                    {
                        unsigned long locPF = PFidx[it];
                        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locPF, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                        PFidx[it] = locPF;

                        unsigned long locVariant = Vidx[it];
                        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locVariant, 1, OP_MPI_UNSIGNED_LONG, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                        Vidx[it] = locVariant;

                        double locdGmax = dGmax[it];
                        OP_MPI_Allreduce(OP_MPI_IN_PLACE, &locdGmax, 1, OP_MPI_DOUBLE, OP_MPI_MAX, OP_MPI_COMM_WORLD);
                        dGmax[it] = locdGmax;
                    }
    #endif
                    for(size_t it = 0; it < locNvariants; it++)
                    {
                        Nucleus locNucleus = *ind;

                        locNucleus.phase = n;
                        locNucleus.parentPhase = m;
                        locNucleus.variant = Vidx[it];
                        locNucleus.index = PFidx[it];
                        locNucleus.dGmin = dGmin;
                        locNucleus.dGnuc = dGmax[it];
                        locNucleus.time_stamp = tStep;

                        Parameters(n,m).PlantedNuclei.push_back(locNucleus);
                    }
                }
            }
        }
    }
}

bool Nucleation::Write(const Settings& locSettings, const int tStep) const
{
#ifdef MPI_PARALLEL
    if(MPI_RANK == 0)
    {
#endif
    string FileName = FileInterface::MakeFileName(locSettings.RawDataDir, "Nucleation_",
                                                             tStep, ".dat");
    fstream out(FileName.c_str(), ios::out);

    if (!out)
    {
        ConsoleOutput::WriteWarning("File/" + FileName + "could not be created",
                                                     thisclassname, "Write(");
        return false;
    }

    for(size_t n = 0; n < Nphases; ++n)
    {
        for(size_t m = 0; m < Nphases; ++m)
        {
            out << Parameters(n,m).GeneratedNuclei.size() << " " ;
            out << Parameters(n,m).PlantedNuclei.size() << " " ;
        }
        out << endl;
    }
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < Parameters(n,m).GeneratedNuclei.size(); ++i)
    {
        Parameters(n,m).GeneratedNuclei[i].write(out);
    }
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < Parameters(n,m).PlantedNuclei.size(); ++i)
    {
        Parameters(n,m).PlantedNuclei[i].write(out);
    }
    out.close();
#ifdef MPI_PARALLEL
    }
#endif
    return true;
}

bool Nucleation::Read(const Settings& locSettings, const BoundaryConditions& BC, const int tStep)
{
    string FileName = FileInterface::MakeFileName(locSettings.InputRawDataDir, "Nucleation_",
                                                             tStep, ".dat");
    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        ConsoleOutput::WriteWarning("File/" + FileName + "could not be opened",
                                                    thisclassname,"Read()");
        return false;
    };

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    {
        int tmpGen;
        inp >> tmpGen;

        if(tmpGen)
        {
            Parameters(n,m).GeneratedNuclei.resize(tmpGen);
            Parameters(n,m).Nsites = tmpGen;
        }

        int tmpNuc;
        inp >> tmpNuc;

        if(tmpNuc)
        {
            Parameters(n,m).PlantedNuclei.resize(tmpNuc);
        }

        if(tmpGen != 0 or tmpNuc != 0)
        {
            Parameters(n, m).Generated = true;
        }
    }

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < Parameters(n,m).GeneratedNuclei.size(); ++i)
    {
        Parameters(n,m).GeneratedNuclei[i].read(inp);
    }

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < Parameters(n,m).PlantedNuclei.size(); ++i)
    {
        Parameters(n,m).PlantedNuclei[i].read(inp);
    }
    inp.close();
    ConsoleOutput::WriteStandard(thisclassname, "Binary input loaded");
    return true;
}

void Nucleation::WriteH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;

    for(size_t n = 0; n < Nphases; ++n)
    {
        for(size_t m = 0; m < Nphases; ++m)
        {
            dbuffer.push_back(GeneratedNuclei(n,m).size());
        }
    }
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedNuclei(n,m).size(); ++i)
    {
        dbuffer.push_back(GeneratedNuclei(n,m)[i].x);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].y);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].z);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].Q[0]);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].Q[1]);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].Q[2]);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].Q[3]);
        dbuffer.push_back(GeneratedNuclei(n,m)[i].radius);
    }
    H5.WriteCheckPoint(tStep, "Nucleation", dbuffer);
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "WriteH5");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}

void Nucleation::ReadH5(H5Interface& H5, const int tStep)
{
    #ifdef H5OP
    std::vector<double> dbuffer;
    H5.ReadCheckPoint(tStep, "Nucleation", dbuffer);
    int i = 0;
    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    {
        int tmp;
        tmp = dbuffer[i]; ++i;
        if(tmp)
        {
            GeneratedNuclei(n,m).resize(tmp);
        }
    }

    for(size_t n = 0; n < Nphases; ++n)
    for(size_t m = 0; m < Nphases; ++m)
    for(size_t i = 0; i < GeneratedNuclei(n,m).size(); ++i)
    {
        GeneratedNuclei(n,m)[i].x = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].y = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].z = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].Q[0] = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].Q[1] = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].Q[2] = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].Q[3] = dbuffer[i]; ++i;
        GeneratedNuclei(n,m)[i].radius= dbuffer[i]; ++i;
    }
    #else
    ConsoleOutput::WriteExit("OpenPhase is not compiled with HDF5 support, use: make Settings=\"H5\"", thisclassname, "ReadH5");
    OP_Exit(EXIT_H5_ERROR);
    #endif
}
} //namespace openphase
