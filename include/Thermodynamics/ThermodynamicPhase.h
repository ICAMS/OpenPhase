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

 *   File created :   2015
 *   Main contributors :   Matthias Stratmann; Oleg Shchyglo
 *
 */

#ifndef THERMODYNAMICPHASE_H
#define THERMODYNAMICPHASE_H

#include "Includes.h"
#include "Thermodynamics/SublatticeModel.h"

namespace openphase
{
class BoundaryConditions;
class ElasticProperties;
class PhaseField;
class SublatticeModel;
class Settings;
class Composition;
enum class DiffusionModels;

/***/
class OP_EXPORTS ThermodynamicPhase                                                        ///  Stores the properties of the thermodynamic phase
{
 public:
    //Constructors and structures
	~ThermodynamicPhase(void);                                                  ///< Destructor
    ThermodynamicPhase(void);                                                   ///< Constructor
    ThermodynamicPhase(const ThermodynamicPhase& RHS)                           ///< Copy constructor
    {
        /**This constructor allows allocating and resizing of a vector of this
        phase, used in ChemicalProperties as std::vector<ThermodynamicPhase>*/
        Cmax = RHS.Cmax;
        Cmin = RHS.Cmin;
        Name = RHS.Name;
        Tmax = RHS.Tmax;
        Tmin = RHS.Tmin;
        Index = RHS.Index;
        Active = RHS.Active;
        Initial = RHS.Initial;
        ConsIndex = RHS.ConsIndex;
        Sublattice = RHS.Sublattice;
        TotInitial = RHS.TotInitial;
        MolarVolume = RHS.MolarVolume;
        Number = RHS.Number;
        State = RHS.State;
		TQName = RHS.TQName;
    }
    ThermodynamicPhase& operator=(const ThermodynamicPhase& RHS)                ///< Assignment operator
    {
        /**This operator allows allocating and resizing of a vector of this
        phase, used in ChemicalProperties as std::vector<ThermodynamicPhase>*/
        if(this != &RHS)
        {
            Cmax = RHS.Cmax;
            Cmin = RHS.Cmin;
            Name = RHS.Name;
            Tmax = RHS.Tmax;
            Tmin = RHS.Tmin;
            Index = RHS.Index;
            Active = RHS.Active;
            Initial.DeAllocate();
            Initial = RHS.Initial;
            ConsIndex.DeAllocate();
            ConsIndex = RHS.ConsIndex;
            Sublattice = RHS.Sublattice;
            TotInitial = RHS.TotInitial;
            MolarVolume.DeAllocate();
            MolarVolume = RHS.MolarVolume;
            Number = RHS.Number;
            State = RHS.State;
            TQName = RHS.TQName;
        }
        return *this;
    }
    struct CompositionLimits
    {
        double Min;
        double Max;
    };

    double Nmoles(Tensor<double,1> Cx);
    //void Initialize(PhaseInput& phinput);                                       ///< Initialize this class

    int Index;                                                                  ///< Phase index in e.g. TQ list of phases.
    size_t Number;                                                              ///< Phase index as a position in the phase-vector
    size_t ReferenceElement;
    bool Active;                                                                ///< Indicates if current phase is considered in the current simulation
    double Tmin;                                                                ///< Indicates min temperature at which the phase is stable
    double Tmax;                                                                ///< Indicates max temperature at which the phase is stable
    AggregateStates State;
    Tensor<int,   1> ConsIndex;                                                 ///< Indicator where a sublattice starts and where it ends
    Tensor<double,1> MolarVolume;                                               ///< Molar volume for each phase //TODO: Storage3D
    Tensor<double,1> Cmin;                                                      ///< Minimum composition value for each phase
    Tensor<double,1> Cmax;                                                      ///< Maximum composition value for each phase
    Tensor<double,1> Initial;                                                   ///< Initial composition of components in all phases
    Tensor<double,1> Initial2;                                                  ///< Initial composition of components in all phases
    std::string Name;                                                           ///< Phase name
    std::string TQName;                                                         ///< Phase name alternative (only used for data output)
    std::vector<double> TotInitial;                                             ///< Initial amount of each component
    std::vector<Element>         Component;                                     ///< List of all chemical components present in the phase
    std::vector<SublatticeModel> Sublattice;                                    ///< Sublattices composing the phase
    std::string thisclassname;                                                  ///< Name of this class
    //TODO: create constant values for and make them private:
    double Nsites;                                                              ///< Returns the sum of site coefficients of all sublattices
    double Nsubsites;                                                           ///< Returns the sum of site coefficients of all substitutional sublattices
    size_t Ncons;                                                               ///< Number of constituents in the present phase
    size_t Nsubs;                                                               ///< Number of sublattices in the present phase
    size_t Ncomp;                                                               ///< Number of components
    double nUfractions;                                                         ///< Amount, to which all u-fractions sum up to in the given phase.
    std::vector<size_t> SublIdx;                                                ///< A vector which contains the first constituent in each sublattice
    std::vector<int> ConsIdx;                                                   ///< A vector that holds the Index of each Constituent
    std::vector<double> NsitesWithI;                                            ///< Summation of all sites of those sublattices, that include a certain element!
    std::vector<std::string> ConNames;                                          ///< A vector of strings with all the Names of each constituent
    std::vector<CompositionLimits> SiteFractions;
    bool isInterstitial(size_t comp);                                           ///<
    bool isEndmember(int CIdxI);
    bool ispairwithVA(int CIdxI);
    bool isStoichiometric(int CIdxI);
    bool get_isStoichiometric(void);
    bool StoichiometricFlag;
    size_t Idx2Cons(size_t sub, int Idx);
    std::vector<std::vector<double> > dMAdYi;
    bool AnalyticChemicalPotentials;
    /*SublatticeModel& Add_Sublattice(unsigned int n_sites)
    {
        SublatticeModel locSublattice(n_sites);
        unsigned int old_size = Constituent.size();
        Constituent.push_back(locElement);
        return Constituent[old_size];
    }*/

 protected:
 private:
    double getNsites(void);                                                     ///< Returns the sum of site coefficients of all sublattices
    size_t getNcons(void);                                                      ///< Number of constituents in the present phase
    size_t getNsubs(void);                                                      ///< Number of sublattices in the present phase
    size_t getNcomp(void);                                                      ///< Number of components
    std::vector<size_t> getSublIdx(void);                                       ///< A vector which contains the first constituent in each sublattice
    std::vector<int> getConsIdx(void);                                          ///< A vector that holds the Index of each Constituent
    std::vector<double> getNsitesWithI(void);                                   ///< Summation of all sites of those sublattices, that include a certain element!
    std::vector<std::string> getConNames(void);                                 ///< A vector of strings with all the Names of each constituent
    std::vector<std::vector<double> > getdMAdYi(void);
    bool getAnalyticChemicalPotentials(void);
    double getSubstitutionalSites(void);
};

}
#endif//THERMODYNAMICPHASE_H
