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

#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include "Thermodynamics/ThermodynamicPhase.h"
#include "Settings.h"
#include "Composition.h"

namespace openphase
{
using namespace std;

ThermodynamicPhase::~ThermodynamicPhase(void)
{
    /**Destructor*/
}

ThermodynamicPhase::ThermodynamicPhase(void)
{
    /**Constructor*/
}
/*
void ThermodynamicPhase::Initialize(PhaseInput& phinput)
{
    Component = phinput.Components;
    Number = phinput.Index;
    vector<string> RefElements = phinput.RefElements;

    //This function is called from the ChemicalProperties::Initialize() and
    //will set the dimension of each storage and fill all parameters, that can
    //be taken from the read in input.

    thisclassname = "ThermodynamicPhase";
    Ncomp = getNcomp();
    Nsubs = getNsubs();
    int ReferenceElementIndex = -1;

    if(RefElements.size() > Number)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            if(Component[comp].Name == RefElements[Number])
            {
                Component[comp].Major = true;
                ReferenceElement = comp;
                ReferenceElementIndex = Component[comp].Index;
            }
        }
    }

    if(!(ConsIndex.IsAllocated()))
    {
        ConsIndex.Allocate({Nsubs+1});
    }
    else
    {
        ConsIndex.Reallocate({Nsubs+1});
    }
    ConsIndex({0}) = 0;
    for(size_t s = 0; s < Nsubs; s++)
    {
        Sublattice[s].Initialize();
        ConsIndex({s+1}) = ConsIndex({s}) + Sublattice[s].Ncons;

        //Check if Reference-Element is present and store position of said
        //   element, otherwise take first element as reference element of this
        //   sublattice! 

        size_t locReferenceElement = 0;
        for(size_t i = 0; i < Sublattice[s].Ncons; i++)
        if(Sublattice[s].Constituent[i].Index == ReferenceElementIndex)
        {
            locReferenceElement = i;
        }
        Sublattice[s].Reference = locReferenceElement;
    }
    Ncons = getNcons();
    Nsites = getNsites();
    Nsubsites = getSubstitutionalSites();
    if (Nsubsites == 0) Nsubsites = 1;
    nUfractions = Nsites/Nsubsites;
    SublIdx = getSublIdx();
    ConsIdx = getConsIdx();
    ConNames = getConNames();
    NsitesWithI = getNsitesWithI();
    dMAdYi = getdMAdYi();
    AnalyticChemicalPotentials = getAnalyticChemicalPotentials();

    if(!(Cmin.IsAllocated()))
    {
        Cmin.Allocate({Ncons});
        Cmax.Allocate({Ncons});
        MolarVolume.Allocate({Ncons});
        Initial.Allocate({Ncons});
        Initial2.Allocate({Ncons});
    }
    else
    {
        Cmin.Reallocate({Ncons});
        Cmax.Reallocate({Ncons});
        MolarVolume.Reallocate({Ncons});
        Initial.Reallocate({Ncons});
        Initial2.Reallocate({Ncons});
    }
    //Composition.Allocate(Nx, Ny, Nz, {Ncons}, boundarysize);
    //MoleFractions.Allocate(Nx, Ny, Nz, {Ncomp}, boundarysize);

    //SublDot.resize(Nsubs);
    //for(int sub = 0; sub < Nsubs; sub++)
    //{
    //    int subNcons = Sublattice[sub].Ncons;
    //    //SublDot[sub].Allocate(Nx,Ny,Nz,{subNcons,subNcons,3,3,3},boundarysize);
    //}
    SiteFractions.resize(Ncons);
    for(size_t con = 0; con < Ncons; con++)
    {
        SiteFractions[con].Min = 0.0;
        SiteFractions[con].Max = 1.0;
    }

    if(phinput.diffMod == DiffusionModels::EquilibriumPartitioning)
    {
        for(size_t comp = 0; comp < Ncomp; comp++)
        {
            Component[comp].isStoichiometric = false;
        }
    }
}*/

double ThermodynamicPhase::getSubstitutionalSites()
{
    double result = 0.0;
    for(size_t sub = 0; sub < Nsubs; sub++)
    {
        bool substitutional = true;
        for(size_t cons = 0; cons < Sublattice[sub].Ncons; cons++)
        if(Sublattice[sub].Constituent[cons].Index == -1)
        substitutional = false;

        if(substitutional)
        result += Sublattice[sub].Site;
    }

    return result;
}

size_t ThermodynamicPhase::getNcons(void)
{
    /**This function will return the total number of constituents for this
    phase.*/
    size_t tmp = 0;
    for(size_t s = 0; s < Nsubs; s++)
    {
        tmp += Sublattice[s].Ncons;
    }
    return tmp;
}

size_t ThermodynamicPhase::getNsubs(void)
{
    /**This function will return the number of sublattices for this phase.*/
    return Sublattice.size();
}

size_t ThermodynamicPhase::getNcomp(void)
{
    /** This function returns the number of components (!not in the present
    phase, but in the total system) as there is no access to
    ChemicalProperties.Ncomp from within this class.*/
    return Component.size();
}

double ThermodynamicPhase::getNsites(void)
{
    /** This function returns the sum of site coefficients of all
    sublattices. For the phase (A,B)_1 (A,B,C)_3 (C)_5 this function will
    return 1 + 3 + 5 = 9.*/
    double locNsites = 0.0;
    for(size_t n = 0; n < Nsubs; n++)
    {
        locNsites += Sublattice[n].Site;
    }
    return locNsites;
}

double ThermodynamicPhase::Nmoles(const Tensor<double,1> Cx)
{
    /** This function returns the number of moles in the present phase (not
    weighted with the phase field). Does not return the overall number of
    moles over all phases in the point [x,y,z].*/

    size_t counter = 0;
    double tempNmoles = 0.0;
    for(size_t n = 0; n < Nsubs; n++)
    {
        double tempSite = Sublattice[n].Site;
        double temp = 0.0;
        for(size_t i = 0; i < Sublattice[n].Ncons; i++)
        {
            if(Sublattice[n].Constituent[i].Index >= 0)
            {
                temp += Cx({counter});
            }
            counter++;
        }
        tempNmoles += temp*tempSite;
    }
    return tempNmoles;
}

std::vector<std::string> ThermodynamicPhase::getConNames(void)
{
    /** This vector of strings holds all the Names of each constituent in
    this phase. Can be used for printouts. */
    std::vector<std::string> tempNames;
    for(size_t sub = 0; sub < Nsubs; sub++)
    for(size_t con = 0; con < Sublattice[sub].Ncons; con++)
    {
        tempNames.push_back(Sublattice[sub].Constituent[con].Name);
    }
    return tempNames;
}

std::vector<size_t> ThermodynamicPhase::getSublIdx(void)
{
    /** For the example sublattice (A,B,C)(D,E)(F)(G,H), SublIdx will hold
    the values [0,3,5,6,8]. This helps when looping over all constituents of
    a phase, but where the sublattice has to be identified. Constituent 0-2,
    3-4, 5, 6-7 all share the same sublattice. Use as:
            for(int sub = 0; sub < Nsubs; sub++)
            for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)*/
    std::vector<size_t> tempSublIdx(Nsubs+1);
    tempSublIdx[0] = 0;
    for(size_t sub = 0; sub < Nsubs; sub++)
    {
        tempSublIdx[sub+1] = tempSublIdx[sub]+Sublattice[sub].Ncons;
    }
    return tempSublIdx;
}

std::vector<int> ThermodynamicPhase::getConsIdx(void)
{
    /**This function returns a vector of the size of all Constituents of the
    phase, that holds the Index of each Constituent to identify each
    Constituent as a Component. If you want to know the Index of one
    Constituent 'con', you have to use ConsIdx[con].*/
    std::vector<int> tempConsIdx(Ncons,-1);
    for(size_t com = 0; com < Ncomp; com++)
    {
        size_t counter = 0;
        for(size_t sub = 0; sub < Nsubs; sub++)
        for(size_t con = 0; con < Sublattice[sub].Ncons; con++)
        {
            if(Component[com].Index == Sublattice[sub].Constituent[con].Index)
            {
                tempConsIdx[counter] = Component[com].Index;
            }
            counter++;
        }
    }
    return tempConsIdx;
}

size_t ThermodynamicPhase::Idx2Cons(size_t sub, int Idx)
{
    /**This function returns the constituent number of the given component index
    Idx*/
    int result = -1;
    for(size_t con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    if(ConsIdx[con] == Idx)
    result = con;
    if(result < 0)
    {
        cout << "Error in Idx2Cons("<<sub<<","<<Idx<<")" << endl;
        OP_Exit(EXIT_FAILURE);
    }
    return result;
}

std::vector<double> ThermodynamicPhase::getNsitesWithI(void)
{
    /**For the example phase (A,B)_1 (A)_3 (A,B)_5 this function will
    hold the following values: NsitesWithI[A] = 9, NsitesWithI[B] = 6,
    NsitesWithI[C] = 0. For the total summation of all Sites, see Nsites.*/
    std::vector<double> temp(Ncomp+1, 0.0);
    std::vector<int> locConsIdx = ConsIdx;
    for(size_t sub = 0; sub < Nsubs; sub++)
    for(size_t con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    {
        if(locConsIdx[con] >= 0)
        {
            temp[locConsIdx[con]]  += Sublattice[sub].Site;
        }
        else
        {
            temp[Ncomp]  += Sublattice[sub].Site;
        }
    }
    return temp;
}

bool ThermodynamicPhase::isInterstitial(size_t comp)
{
    /**This boolean value will identify if the component is an interstitial
    component (true, if present on sublattice with vacancies) or an
    substitutional component (false, if not present on any sublattice with
    vacancies).*/
    bool temp = false;
    for(size_t sub = 0; sub < Nsubs; sub++)
    if(Sublattice[sub].hasVacancies
       and Sublattice[sub].isElementPresent(Component[comp].Index))
    {
        temp = true;
    }
    return temp;
}

bool ThermodynamicPhase::isEndmember(int CIdxI)
{
    /**This boolean value will identify if the component I is an endmember in
    this phase. That means if component I exists on every sublattice, which
    doesn't have vacancies. CIdxI has to be the index of the component!!*/

    bool Iispresent = false;
    bool Iisendmember = false;

    for(size_t sub = 0; sub < Nsubs; sub++)
    if(Sublattice[sub].isElementPresent(CIdxI))
    Iispresent = true;

    if(Iispresent)
    {
        Iisendmember = true;
        for(size_t sub = 0; sub < Nsubs; sub++)
        if((not Sublattice[sub].isElementPresent(CIdxI))
        and(not Sublattice[sub].hasVacancies)
        /*and(Sublattice[sub].Ncons > 1)*/)
        {// no I nor Va on sublattice -> I not an endmember (if Ncons > 1)
            Iisendmember = false;
        }
    }
    return Iisendmember;
}

bool ThermodynamicPhase::ispairwithVA(int CIdxI)
{
    /**This boolean value will identify if the component I shares its sublattice
    with vacancy VA on at least one sublattice. CIdxI has to be the index of the
    component!!*/

    bool IshareswithVA = false;

    for(size_t sub = 0; sub < Nsubs; sub++)
    if((Sublattice[sub].isElementPresent(CIdxI)) and Sublattice[sub].isElementPresent(-1))
    IshareswithVA = true;

    return IshareswithVA;
}

bool ThermodynamicPhase::isStoichiometric(int CIdxI)
{
    /**This boolean value will identify if the component I is a stoichiometric
    component, as it shares its sublattice without any other component.*/

    bool aloneinsubl = true;
    bool nointerstitials = true;

    bool isStoich = false;

    for(size_t sub = 0; sub < Nsubs; sub++)
    {
        if((Sublattice[sub].isElementPresent(CIdxI)) and Sublattice[sub].Ncons > 1)
        {
            aloneinsubl = false;
        }

        if(Sublattice[sub].hasVacancies and Sublattice[sub].Ncons > 1)
        {
            nointerstitials = false;
        }
    }

    if(aloneinsubl and nointerstitials)
    isStoich = true;

    return isStoich;
}

bool ThermodynamicPhase::get_isStoichiometric(void)
{
    /***/

    bool isStoich = false;

    for(size_t n = 0; n < Ncons; n++)
    if(ConsIdx[n] > -1 and isStoichiometric(ConsIdx[n]))
    {
        isStoich = true;
    }

    return isStoich;
}

vector<vector<double> > ThermodynamicPhase::getdMAdYi(void)
{
    /**This matrix is populated with values for dM_A/dY_i, which are necessary
    for calculation of chemical potentials*/

    vector<vector<double> > result;

    for(size_t comp = 0; comp < Ncomp; comp++)
    {
        vector<double> tempcons(Ncons, 0.0);

        size_t counter = 0;
        for(size_t sub = 0; sub < Nsubs; sub++)
        for(size_t con = 0; con < Sublattice[sub].Ncons; con++)
        {
            if(Sublattice[sub].Constituent[con].Index == Component[comp].Index)
            {
                tempcons[counter] = Sublattice[sub].Site;
            }
            counter++;
        }

        result.push_back(tempcons);
    }

    return result;
}

bool ThermodynamicPhase::getAnalyticChemicalPotentials(void)
{
    /**This boolean determines, whether chemical potentials can be calculated
    analytically or only by numeric calculation*/

    bool result = true;

    for(size_t A = 0; A < Ncomp; A++)
    {
        int CIdxI = Component[A].Index;
        if(!(isEndmember(CIdxI) or ispairwithVA(CIdxI)))
        result = false;

    }

    //check if all components are endmembers or share a sublattice with VA

    return result;
}

}
