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

#ifndef INTERFACEENERGYMODEL_H
#define INTERFACEENERGYMODEL_H

#include "Includes.h"
#include "Tools.h"

namespace openphase
{

enum class InterfaceEnergyModels
{
    Ext,                                                                        ///< Interface energy should be calculated externally
    Iso,                                                                        ///< Isotropic interface energy
    Cubic,                                                                      ///< Cubic anisotropy for solid phases during solidification
    CubicFull,
    HexBoettger,                                                                ///< Hexagonal anisotropy for solid phases during solidification (after Boettger et al.)
    HexSun,                                                                     ///< Hexagonal anisotropy for solid phases during solidification (after Sun et al.)
    HexYang,                                                                    ///< Hexagonal anisotropy for solid phases during solidification (after Yang et al.)
    Faceted,                                                                    ///< Interface energy with strong anisotropy for faceted crystals
    FacetedEXP,                                                                 ///< Interface energy with strong anisotropy for faceted crystals (experimental)
    FacetedFull,
    Disorientation                                                              ///< Interface energy with strong anisotropy considering grains disorientation
};

enum class InterfaceEnergyModelType                                             ///< Interface energy model types
{
    Stiffness,                                                                  ///< Stiffness model
    Energy                                                                      ///< Energy model with the driving force contribution
};

class InterfaceEnergyModel                                                      ///< Interface energy models implementation class
{
 public:
    InterfaceEnergyModels Model;
    InterfaceEnergyModelType Type;

    InterfaceEnergyModel()
    {
        Model = InterfaceEnergyModels::Iso;
        Type  = InterfaceEnergyModelType::Stiffness;

        Energy = 0.0;
        MinEnergy = 0.0;
        MaxEnergy = 0.0;
        Sigma0   = 0.0;
        Epsilon1 = 0.0;
        Epsilon2 = 0.0;
        Epsilon3 = 0.0;
        Epsilon4 = 0.0;
        Nfacets  = 0;
    }
    void ReadInputIso(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the isotropic model
    {
        Energy = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::Iso;
        Type  = InterfaceEnergyModelType::Energy;
    }
    void ReadInputCubic(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the cubic model
    {
        Energy   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE") + counter);

        MaxEnergy = Energy*(1.0 + Epsilon1*1.41);
        Model = InterfaceEnergyModels::Cubic;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadInputCubicFull(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the cubic model
    {
        Sigma0    = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma0") + counter);
        Epsilon1  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE") + counter);
        Energy    = Sigma0*(1.0 + Epsilon1/3.0); //Minimal interface energy?TODO
        MaxEnergy = Sigma0*(1.0 + Epsilon1);

        Model = InterfaceEnergyModels::CubicFull;
        Type  = InterfaceEnergyModelType::Energy;
    }
    void ReadInputHexBoettger(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Boettger et al
    {
        Energy   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE") + counter);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexBoettger;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadInputHexSun(std::stringstream& inp, int moduleLocation, std::string counter)///< Read parameters of the hexagonal model by Sun et al
    {
        Energy   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE20") + counter, false, -0.026);
        Epsilon2 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE40") + counter, false,  0.0);
        Epsilon3 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE60") + counter, false,  0.0);
        Epsilon4 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE66") + counter, false,  0.003);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexSun;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadInputHexYang(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the hexagonal model by Yang et al
    {
        Energy   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        Epsilon1 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE1") + counter, false,-0.02);
        Epsilon2 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE2") + counter, false, 0.15);
        Epsilon3 = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonE3") + counter, false, 0.15);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexYang;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadJSONIso(json& j, size_t alpha, size_t beta)///< Read parameters of the isotropic model
    {
        Energy = FileInterface::ReadParameter<double>(j, {"Sigma", alpha, beta});
        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::Iso;
        Type  = InterfaceEnergyModelType::Energy;
    }
    void ReadJSONCubic(json& j, size_t alpha, size_t beta)///< Read parameters of the cubic model
    {
        Energy   = FileInterface::ReadParameter<double>(j, {"Sigma", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonE", alpha, beta});

        MaxEnergy = Energy*(1.0 + Epsilon1*1.41);
        Model = InterfaceEnergyModels::Cubic;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadJSONCubicFull(json& j, size_t alpha, size_t beta)///< Read parameters of the cubic model
    {
        Sigma0    = FileInterface::ReadParameter<double>(j, {"Sigma0", alpha, beta});
        Epsilon1  = FileInterface::ReadParameter<double>(j, {"EpsilonE", alpha, beta});
        Energy    = Sigma0*(1.0 + Epsilon1/3.0); //Minimal interface energy?TODO
        MaxEnergy = Sigma0*(1.0 + Epsilon1);

        Model = InterfaceEnergyModels::CubicFull;
        Type  = InterfaceEnergyModelType::Energy;
    }
    void ReadJSONHexBoettger(json& j, size_t alpha, size_t beta)///< Read parameters of the hexagonal model by Boettger et al
    {
        Energy   = FileInterface::ReadParameter<double>(j, {"Sigma", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonE", alpha, beta});

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexBoettger;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadJSONHexSun(json& j, size_t alpha, size_t beta)///< Read parameters of the hexagonal model by Sun et al
    {
        Energy   = FileInterface::ReadParameter<double>(j, {"Sigma", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonE20", alpha, beta}, -0.026);
        Epsilon2 = FileInterface::ReadParameter<double>(j, {"EpsilonE40", alpha, beta}, 0.0);
        Epsilon3 = FileInterface::ReadParameter<double>(j, {"EpsilonE60", alpha, beta}, 0.0);
        Epsilon4 = FileInterface::ReadParameter<double>(j, {"EpsilonE66", alpha, beta}, 0.003);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexSun;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    void ReadJSONHexYang(json& j, size_t alpha, size_t beta) ///< Read parameters of the hexagonal model by Yang et al
    {
        Energy   = FileInterface::ReadParameter<double>(j, {"Sigma", alpha, beta});
        Epsilon1 = FileInterface::ReadParameter<double>(j, {"EpsilonE1", alpha, beta},-0.02);
        Epsilon2 = FileInterface::ReadParameter<double>(j, {"EpsilonE2", alpha, beta}, 0.15);
        Epsilon3 = FileInterface::ReadParameter<double>(j, {"EpsilonE3", alpha, beta}, 0.15);

        MaxEnergy = Energy;
        Model = InterfaceEnergyModels::HexYang;
        Type  = InterfaceEnergyModelType::Stiffness;
    }
    
    void ReadInputFaceted(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the Faceted model
    {
        bool FacetFamilies = true;
        Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("NfacetFamilies") + counter, false, 0);
        if(Nfacets == 0)
        {
            Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("Nfacets") + counter);

            if(Nfacets == 0)
            {
                ConsoleOutput::WriteExit("Zero Nfacets specified!\nThe input should contain non-zero number of facets!", "InterfaceEnergyModel", "ReadInputFaceted()");
                OP_Exit(EXIT_FAILURE);
            }
            FacetFamilies = false;
        }

        FacetVectors.resize(Nfacets);
        FacetEnergy.resize(Nfacets,0.0);
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

            FacetEnergy[n]  = FileInterface::ReadParameterD(inp, moduleLocation, std::string("SigmaF") + counter + counter1);
            FacetEpsilon[n] = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonEF") + counter + counter1);
        }

        for(size_t n = 0; n < Nfacets; n++)
        {
            Energy += FacetEnergy[n]/Nfacets;

            if(FacetEnergy[n] > MaxEnergy)
            {
                MaxEnergy = FacetEnergy[n];
            }
        }
        Model = InterfaceEnergyModels::Faceted;
    };
    void ReadInputFacetedEXP(std::stringstream& inp, int moduleLocation, std::string counter) ///< Read parameters of the Faceted model
    {
        bool FacetFamilies = true;
        Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("NfacetFamilies") + counter, false, 0);
        if(Nfacets == 0)
        {
            Nfacets = FileInterface::ReadParameterI(inp, moduleLocation, std::string("Nfacets") + counter);

            if(Nfacets == 0)
            {
                ConsoleOutput::WriteExit("Zero Nfacets specified!\nThe input should contain non-zero number of facets!", "InterfaceEnergyModel", "ReadInputFaceted()");
                OP_Exit(EXIT_FAILURE);
            }
            FacetFamilies = false;
        }

        Energy   = FileInterface::ReadParameterD(inp, moduleLocation, std::string("Sigma") + counter);
        MaxEnergy = Energy;

        FacetVectors.resize(Nfacets);
        FacetEpsilon.resize(Nfacets,0.0);
        FacetPower.resize(Nfacets,2.0);

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
            FacetEpsilon[n] = FileInterface::ReadParameterD(inp, moduleLocation, std::string("EpsilonEF") + counter + counter1);
            FacetPower[n] = FileInterface::ReadParameterD(inp, moduleLocation, std::string("PowerEF") + counter + counter1);
            if(FacetPower[n] < 2)
            {
                ConsoleOutput::WriteExit("Facet power parameter must be >= 2", "InterfaceEnergyModel", "ReadInputFacetedEXP()");
                OP_Exit(EXIT_FAILURE);
            }
        }
        Model = InterfaceEnergyModels::FacetedEXP;
    };
    void ReadInputFacetedFull(std::stringstream& inp, int moduleLocation, std::string counter)
    {
        ReadInputFaceted(inp,moduleLocation,counter);
        Model = InterfaceEnergyModels::FacetedFull;
        Type  = InterfaceEnergyModelType::Energy;
    };
    void ReadInputDisorientation(std::stringstream& inp, int moduleLocation, std::string counter)
    {
        Model = InterfaceEnergyModels::Disorientation;
    };

    double CalculateDisorientation(Quaternion& OrientationA, Quaternion& OrientationB)
    {
        double constant = 0.85;                                                 // Shear Modulus*Burgers vector/4*Pi*(1-v)
        double misorientation = Tools::getDisorientationCubic(OrientationA, OrientationB);
        double locEnergy = 0.0;

        if (misorientation <= 0.349066)
        {
            locEnergy = misorientation*constant*(0.8 - log(misorientation));
        }
        else
        {
            locEnergy = 0.349066*constant*(0.8 - log(0.349066));
        }
        return locEnergy;
    }

    double Calculate(dVector3& locNormal)                                       ///< Calculate anisotropic interface energy
    {
        double locEnergy = 0.0;
        switch(Model)
        {
            case InterfaceEnergyModels::Ext:
            {
                break;
            }
            case InterfaceEnergyModels::Iso:
            {
                locEnergy = MaxEnergy;
                break;
            }
            case InterfaceEnergyModels::Cubic:                                  //
            {
                locEnergy = Energy*(1.0 + Epsilon1 * (1.5 - 2.5*(pow(locNormal[0], 4) +
                                                                 pow(locNormal[1], 4) +
                                                                 pow(locNormal[2], 4))));
                break;
            }
            case InterfaceEnergyModels::CubicFull:
            {
                const double nx4 = std::pow(locNormal[0],4);
                const double ny4 = std::pow(locNormal[1],4);
                const double nz4 = std::pow(locNormal[2],4);
                locEnergy = Sigma0*(1.0 + Epsilon1*(nx4+ny4+nz4));
                break;
            }
            case InterfaceEnergyModels::HexBoettger:                            // Bottger, B.; Eiken, J.; Steinbach, I. Phase field simulation of equiaxed solidification in technical alloys. Acta Mater. 2006, 54, 2697–2704. [http://doi.org/10.1016/j.actamat.2006.02.008]
            {
                locEnergy = Energy*(1.0 - Epsilon1 *
                            (      std::pow(locNormal[0], 6) -
                                   std::pow(locNormal[1], 6) -
                            15.0 * std::pow(locNormal[0], 4) * locNormal[1]*locNormal[1] +
                            15.0 * std::pow(locNormal[1], 4) * locNormal[0]*locNormal[0] +
                            (5.0 * std::pow(locNormal[2], 4) -
                             5.0 * std::pow(locNormal[2], 2) +
                                   std::pow(locNormal[2], 6))));
                break;
            }
            case InterfaceEnergyModels::HexSun:                                 // Sun, D.; Mendelev, M.; Becker, C.; Kudin, K.; Haxhimali, T.; Asta, M.; Hoyt, J.; Karma, A.; Srolovitz, D. Crystal-melt interfacial free energies in hcp metals: A molecular dynamics study of Mg. Phys. Rev. B 2006, 73, 024116. [http://doi.org/10.1103/PhysRevB.73.024116]
            {
                locEnergy = Energy*(1.0
                          + Epsilon1*sqrt(5.0/16.0/Pi)*(3.0*locNormal[2]*locNormal[2] - 1.0)
                          + Epsilon2*3.0/16.0/sqrt(Pi)*(35.0*pow(locNormal[2],4) -
                                                        30.0*locNormal[2]*locNormal[2] + 3.0)
                          + Epsilon3*sqrt(13.0/Pi)/32.0*(231.0*pow(locNormal[2],6) -
                                                         315.0*pow(locNormal[2],4) +
                                                         105.0*locNormal[2]*locNormal[2] - 5.0)
                          + Epsilon4*sqrt(6006.0/Pi)/64.0*(pow(locNormal[0],6) -
                                                           15.0*pow(locNormal[0],4)*locNormal[1]*locNormal[1] +
                                                           15.0*locNormal[0]*locNormal[0]*pow(locNormal[1],4) -
                                                           pow(locNormal[1],6)));
                break;
            }
            case InterfaceEnergyModels::HexYang:                                // Yang, M.; Xiong, S.; Guo, Z. Characterisation of the 3-D dendrite morphology of magnesium alloys using synchrotron X-ray tomography and 3-D phase-field modelling. Acta Mater. 2015, 92, 8–17. [http://doi.org/10.1016/j.actamat.2015.03.044]
            {
                locEnergy = Energy*(1.0
                          + Epsilon1 * pow((3.0*locNormal[2]*locNormal[2] - 1.0), 2)
                          + Epsilon2 * pow((locNormal[0]*locNormal[0]*locNormal[0] - 3.0*locNormal[0]*locNormal[1]*locNormal[1]),2)
                          * pow((9.0*locNormal[2]*locNormal[2] - 1.0 + Epsilon3),2));
                break;
            }
            case InterfaceEnergyModels::FacetedEXP:                                // Inspired by M. Salvalaglio et al. Cryst. Growth Des. 2015, 15, 2787−2794 DOI: 10.1021/acs.cgd.5b00165
            {
                locEnergy = Energy;
                for(size_t n = 0; n < Nfacets; n++)
                for(size_t j = 0; j < FacetVectors[n].size() ; j++)
                {
                    dVector3 locFacet = FacetVectors[n][j].normalized();
                    double cos_inc  = locNormal * locFacet;
                    if(cos_inc > 0.0)
                    {
                        // part of the energy function
                        locEnergy -= Energy*FacetEpsilon[n]*pow(cos_inc,FacetPower[n]);
                        // second order derivative
                        locEnergy += Energy*FacetEpsilon[n]*FacetPower[n]*pow(cos_inc,FacetPower[n]);
                        //// !!!The contribution below is zero near the bottom of the cusp and can be ignored as irrelevant!!!
                        //locEnergy -= FacetEnergy[n]*FacetEpsilon[n]*FacetPower[n]*(FacetPower[n]-1)*pow(cos_inc,FacetPower[n]-2)*(1.0 - cos_inc*cos_inc);
                    }
                }
                if(locEnergy < 0.0)
                {
                    ConsoleOutput::WriteWarning("Negative interface stiffness encountered for InterfaceEnergyModels::FacetedEXP!", "InterfaceEnergyModels","Calculate()");
                }

                break;
            }
            case InterfaceEnergyModels::Faceted:                                // H. F. M. A. Salama, J. Kundin, O. Shchyglo, V. Mohles, K. Marquardt, I. Steinbach. Role of inclination dependence of grain boundary energy on the microstructure evolution during grain growth, Acta Materialia, 188, 641-651, (2020) [http://dx.doi.org/10.1016/j.actamat.2020.02.043]
            {
                double CosInc2  = 0.0;
                double EnergyA  = 0.0;
                double EpsilonA = 0.0;

                for(size_t n = 0; n < Nfacets; n++)
                for(size_t j = 0; j < FacetVectors[n].size() ; j++)
                {
                    dVector3 locFacet = FacetVectors[n][j].normalized();
                    double cos_inc  = locNormal * locFacet;
                    double cos_inc2 = cos_inc*cos_inc;

                    if(cos_inc2 > CosInc2)
                    {
                        CosInc2  = cos_inc2;
                        EnergyA  = FacetEnergy[n];
                        EpsilonA = FacetEpsilon[n];
                    }
                }
                double SinInc2 = 1.0 - CosInc2;
                locEnergy = EnergyA * EpsilonA * EpsilonA * pow(sqrt(SinInc2 + EpsilonA*EpsilonA*CosInc2),-3.0);
                break;
            }
            case InterfaceEnergyModels::FacetedFull:
            {
                double Inclination = Pi;
                double EnergyA = 0.0;
                double EpsilonA = 0.0;

                // Calculating inclination angle between facet normal and interface normal
                for(size_t n = 0; n < Nfacets; n++)
                for(size_t j = 0; j < FacetVectors[n].size() ; j++)
                {
                    const dVector3 locFacet = FacetVectors[n][j].normalized();
                    const double dotVec = std::acos(locNormal * locFacet);

                    if(std::fabs(dotVec) <= std::fabs(Inclination))
                    {
                        Inclination = dotVec;
                        EnergyA     = FacetEnergy[n];
                        EpsilonA    = FacetEpsilon[n];
                    }
                }

                const double sin_inc = std::sin(Inclination);
                const double cos_inc = std::cos(Inclination);
                locEnergy = EnergyA * std::sqrt(sin_inc*sin_inc + EpsilonA*EpsilonA*cos_inc*cos_inc);

                break;
            }
            case InterfaceEnergyModels::Disorientation:
            {
                // The disorientation model is implemented in its own method
                break;
            }
        }
        assert(locEnergy > 0.0 && "Negative interface energy reduce anisotropy");
        return locEnergy;
    };

    dVector3 Derivative(const dVector3& locNormal) const
    {
        switch(Model)
        {
            case InterfaceEnergyModels::CubicFull:
            {
                const double nx3 = std::pow(locNormal[0],3);
                const double ny3 = std::pow(locNormal[1],3);
                const double nz3 = std::pow(locNormal[2],3);
                return {Sigma0*Epsilon1*4*nx3,Sigma0*Epsilon1*4*ny3,Sigma0*Epsilon1*4*nz3};
            }
            case InterfaceEnergyModels::FacetedFull:
            {
                double Inclination = Pi;
                double EnergyA = 0.0;
                double EpsilonA = 0.0;
                dVector3 Facet = {0.0,0.0,0.0};

                // Calculating inclination angle between facet normal and interface normal
                for(size_t n = 0; n < Nfacets; n++)
                for(size_t j = 0; j < FacetVectors[n].size() ; j++)
                {
                    const dVector3 locFacet = FacetVectors[n][j].normalized();
                    const double dotVec = std::acos(locNormal * locFacet);

                    if(std::fabs(dotVec) <= std::fabs(Inclination))
                    {
                        Inclination = dotVec;
                        EnergyA     = FacetEnergy[n];
                        EpsilonA    = FacetEpsilon[n];
                        Facet       = locFacet;
                    }
                }

                const double sin_inc = std::sin(Inclination);
                const double cos_inc = std::cos(Inclination);
                const double dEnergy_dInclination = EnergyA * (1.0-EpsilonA*EpsilonA)*sin_inc*cos_inc/std::sqrt(sin_inc*sin_inc + EpsilonA*EpsilonA*cos_inc*cos_inc);
                const dVector3 dInclination_dn = Facet*(-1.0)/std::sqrt(1-(locNormal*Facet)*(locNormal*Facet));

                return dInclination_dn*dEnergy_dInclination;
            }
            default:
                return {0.0,0.0,0.0};
        }
    };
    size_t Nfacets;                                                             ///< Number of facet types, e.g. <100>, <111>, <110>, etc...
    std::vector<std::vector<iVector3>> FacetVectors;                            ///< Set of planes for different facet types
    std::vector<double> FacetEnergy;                                            ///< Interface energy for a facet family
    std::vector<double> FacetEpsilon;                                           ///< Interface energy anisotropy parameter for a facet family
    std::vector<double> FacetPower;                                             ///< Experimental interface energy anisotropy parameter

    double Energy;                                                              ///< Interface energy
    double MaxEnergy;                                                           ///< Maximum interface energy for a phase pair
    double MinEnergy;                                                           ///< Minimum interface energy for a phase pair
    double Sigma0;                                                              ///< Interface energy parameter
    double Epsilon1;                                                            ///< Interface energy anisotropy parameter
    double Epsilon2;                                                            ///< Interface energy anisotropy parameter
    double Epsilon3;                                                            ///< Interface energy anisotropy parameter
    double Epsilon4;                                                            ///< Interface energy anisotropy parameter
};
}// namespace openphase
#endif
