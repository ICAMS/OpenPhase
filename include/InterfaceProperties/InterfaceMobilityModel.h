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

 *   File created :   2019
 *   Main contributors :   Oleg Shchyglo; Alexander Monas; Efim Borukhovich
 *
 */

#ifndef INTERFACEMOBILITYMODEL_H
#define INTERFACEMOBILITYMODEL_H

#include "Includes.h"
#include "Tools.h"

namespace openphase
{
enum class InterfaceMobilityModels
{
    Ext,                                                                        ///< Interface mobility should be calculated externally
    Iso,                                                                        ///< Isotropic interface mobility
    Cubic,                                                                      ///< Cubic anisotropy for solid phases during solidification
    HexBoettger,                                                                ///< Hexagonal anisotropy for solid phases during solidification (after Boettger et al.)
    HexSun,                                                                     ///< Hexagonal anisotropy for solid phases during solidification (after Sun et al.)
    HexYang,                                                                    ///< Hexagonal anisotropy for solid phases during solidification (after Yang et al.)
    Faceted                                                                     ///< Interface energy with strong anisotropy for faceted crystals
};

class InterfaceMobilityModel                                                    ///< Interface mobility models implementation class
{
 public:
    InterfaceMobilityModels Model;

    InterfaceMobilityModel()
    {
        Model = InterfaceMobilityModels::Iso;
        Mobility = 0.0;
        MaxMobility = 0.0;
        ActivationEnergy = 0.0;
        Epsilon1 = 0.0;
        Epsilon2 = 0.0;
        Epsilon3 = 0.0;
        Epsilon4 = 0.0;
        Nfacets  = 0;
    };
    void ReadInputIso(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the isotropic model
    {
        Mobility = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        ActivationEnergy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::Iso;
    };
    void ReadInputCubic(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the cubic model
    {
        Mobility = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM") + counter);
        ActivationEnergy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::Cubic;
    };
    void ReadInputHexBoettger(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Boettger et al
    {
        Mobility = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM") + counter);
        ActivationEnergy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::HexBoettger;
    };
    void ReadInputHexSun(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Sun et al
    {
        Mobility = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM20") + counter, false, -0.026);
        Epsilon2 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM40") + counter, false,  0.0);
        Epsilon3 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM60") + counter, false,  0.0);
        Epsilon4 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM66") + counter, false,  0.003);

        ActivationEnergy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexSun;
    };
    void ReadInputHexYang(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the hexagonal model by Yang et al
    {
        Mobility = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Mu") + counter);
        ActivationEnergy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Q") + counter, false, 0.0);

        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM1") + counter, false,-0.02);
        Epsilon2 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM2") + counter, false, 0.15);
        Epsilon3 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonM3") + counter, false, 0.15);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexYang;
    };
    void ReadJSONIso(json& j, size_t alpha, size_t beta)///< Read parameters of the isotropic model
    {
        Mobility = FileInterface::ReadParameter<double>(j, {"Mu", alpha, beta});
        ActivationEnergy = FileInterface::ReadParameter<double>(j, {"Q", alpha, beta}, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::Iso;
    };
    void ReadJSONCubic(json& j, size_t alpha, size_t beta)///< Read parameters of the cubic model
    {
        Mobility = FileInterface::ReadParameter<double>(j, {"Mu", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonM", alpha, beta});
        ActivationEnergy = FileInterface::ReadParameter<double>(j, {"Q", alpha, beta}, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::Cubic;
    };
    void ReadJSONHexBoettger(json& j, size_t alpha, size_t beta)///< Read parameters of the hexagonal model by Boettger et al
    {
        Mobility = FileInterface::ReadParameter<double>(j, {"Mu", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonM", alpha, beta});
        ActivationEnergy = FileInterface::ReadParameter<double>(j, {"Q", alpha, beta}, 0.0);

        MaxMobility = Mobility;//*std::min(1.0, (1.0 - Epsilon1));
        Model = InterfaceMobilityModels::HexBoettger;
    };
    void ReadJSONHexSun(json& j, size_t alpha, size_t beta)///< Read parameters of the hexagonal model by Sun et al
    {
        Mobility = FileInterface::ReadParameter<double>(j, {"Mu", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonM20", alpha, beta}, -0.026);
        Epsilon2 = FileInterface::ReadParameter<double>(j, {"EpsilonM40", alpha, beta}, 0.0);
        Epsilon3 = FileInterface::ReadParameter<double>(j, {"EpsilonM60", alpha, beta}, 0.0);
        Epsilon4 = FileInterface::ReadParameter<double>(j, {"EpsilonM66", alpha, beta}, 0.003);

        ActivationEnergy = FileInterface::ReadParameter<double>(j, {"Q", alpha, beta}, 0.0);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexSun;
    };
    void ReadJSONHexYang(json& j, size_t alpha, size_t beta)///< Read parameters of the hexagonal model by Yang et al
    {
        Mobility = FileInterface::ReadParameter<double>(j, {"Mu", alpha, beta});
        ActivationEnergy = FileInterface::ReadParameter<double>(j, {"Q", alpha, beta}, 0.0);

        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonM1", alpha, beta},-0.02);
        Epsilon2 = FileInterface::ReadParameter<double>(j, {"EpsilonM2", alpha, beta}, 0.15);
        Epsilon3 = FileInterface::ReadParameter<double>(j, {"EpsilonM3", alpha, beta}, 0.15);

        MaxMobility = Mobility;
        Model = InterfaceMobilityModels::HexYang;
    };
    void ReadInputFaceted(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the Faceted model
    {
        bool FacetFamilies = true;
        Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("NfacetFamilies") + counter, false, 0);
        if(Nfacets == 0)
        {
            Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("Nfacets") + counter);
            FacetFamilies = false;
        }

        FacetVectors.resize(Nfacets);
        FacetMobility.resize(Nfacets,0.0);
        FacetEpsilon.resize(Nfacets,0.0);

        for(size_t n = 0 ; n < Nfacets; n++)
        {
            std::stringstream converter;
            converter << "_" << n;
            std::string counter1 = converter.str();

            iVector3 locFacet = FileInterface::ReadParameterV3I(inp, moduleLocation, std::string("Facet") + counter + counter1, true);

            if(FacetFamilies)
            {
                FacetVectors[n] = Tools::findIndexPermutations(locFacet);
            }
            else
            {
                FacetVectors[n].push_back(locFacet);
            }

            FacetMobility[n]  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("MuF") + counter + counter1);
            FacetEpsilon[n] = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonMF") + counter + counter1);
        }

        for(size_t n = 0; n < Nfacets; n++)
        {
            if (FacetMobility[n] > MaxMobility)
            {
                MaxMobility = FacetMobility[n];
            }
        }
        Model = InterfaceMobilityModels::Faceted;
    };

    double Calculate(dVector3& locNormal)                                       ///< Calculate anisotropic interface mobility for cubic symmetry grains
    {
        double locMobility = 0.0;
        switch(Model)
        {
            case InterfaceMobilityModels::Ext:
            {
                break;
            }
            case InterfaceMobilityModels::Iso:
            {
                locMobility = MaxMobility;
                break;
            }
            case InterfaceMobilityModels::Cubic:
            {
                locMobility = Mobility*(1.0 - Epsilon1 * (1.5 - 2.5*(pow(locNormal[0], 4) +
                                                                     pow(locNormal[1], 4) +
                                                                     pow(locNormal[2], 4))));
                break;
            }
            case InterfaceMobilityModels::HexBoettger:
            {
                locMobility = Mobility*(1.0 + Epsilon1 *
                            (      std::pow(locNormal[0], 6) -
                                   std::pow(locNormal[1], 6) -
                            15.0 * std::pow(locNormal[0], 4) * locNormal[1]*locNormal[1] +
                            15.0 * std::pow(locNormal[1], 4) * locNormal[0]*locNormal[0] +
                            (5.0 * std::pow(locNormal[2], 4) -
                             5.0 * std::pow(locNormal[2], 2) +
                                   std::pow(locNormal[2], 6))));
                break;
            }
            case InterfaceMobilityModels::HexSun:
            {
                locMobility = Mobility*(1.0
                        -Epsilon1*sqrt(5.0/16.0/Pi)*(3.0*locNormal[2]*locNormal[2] - 1.0)
                        -Epsilon2*3.0/16.0/sqrt(Pi)*(35.0*pow(locNormal[2],4) -
                                                     30.0*locNormal[2]*locNormal[2] + 3.0)
                        -Epsilon3*sqrt(13.0/Pi)/32.0*(231.0*pow(locNormal[2],6) -
                                                      315.0*pow(locNormal[2],4) +
                                                      105.0*locNormal[2]*locNormal[2] - 5.0)
                        -Epsilon4*sqrt(6006.0/Pi)/64.0*(pow(locNormal[0],6) -
                                                        15.0*pow(locNormal[0],4)*locNormal[1]*locNormal[1] +
                                                        15.0*locNormal[0]*locNormal[0]*pow(locNormal[1],4) -
                                                        pow(locNormal[1],6)));
                break;
            }
            case InterfaceMobilityModels::HexYang:
            {
                locMobility = Mobility*(1.0
                            - Epsilon1 * pow((3.0*locNormal[2]*locNormal[2] - 1.0), 2)
                            - Epsilon2 * pow((locNormal[0]*locNormal[0]*locNormal[0] - 3.0*locNormal[0]*locNormal[1]*locNormal[1]),2)
                            * pow((9.0*locNormal[2]*locNormal[2] - 1.0 + Epsilon3),2));
                break;
            }
            case InterfaceMobilityModels::Faceted:
            {
                double Inclination = Pi;
                double MobilityA = 0.0;
                double EpsilonMA = 0.0;

                // Calculating inclination angle between facet normal and interface normal
                for(size_t n = 0; n < Nfacets; n++)
                for(size_t j = 0; j < FacetVectors[n].size(); j++)
                {
                    dVector3 locFacet = FacetVectors[n][j].normalized();
                    double dotVec = acos(locNormal * locFacet);

                    if(fabs(dotVec) <= fabs(Inclination))
                    {
                        Inclination = dotVec;
                        MobilityA = FacetMobility[n];
                        EpsilonMA = FacetEpsilon[n];
                    }
                }

                locMobility = MobilityA * (EpsilonMA + (1.0 - EpsilonMA) * fabs(tan(Inclination)) * tanh(pow(fabs(tan(Inclination)),-1)));
                break;
            }
        }
        return locMobility;
    };

    size_t Nfacets;                                                             ///< Number of facet types, e.g. <100>, <111>, <110>, etc...
    std::vector<std::vector<iVector3>> FacetVectors;                            ///< Set of planes for different facet types
    std::vector<double> FacetMobility;                                          ///< Interface mobility for a facet family
    std::vector<double> FacetEpsilon;                                           ///< Interface mobility anisotropy for a facet family

    double Mobility;                                                            ///< Interface mobility
    double MaxMobility;                                                         ///< Maximum interface mobility for a phase pair
    double Epsilon1;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon2;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon3;                                                            ///< Interface mobility anisotropy parameter
    double Epsilon4;                                                            ///< Interface mobility anisotropy parameter
    double ActivationEnergy;                                                    ///< Interface mobility activation energy
};
}// namespace openphase
#endif
