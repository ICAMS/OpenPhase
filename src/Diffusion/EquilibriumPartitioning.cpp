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
 *  Main contributors :   Oleg Shchyglo; Marvin Tegeler; Matthias Stratmann
 *
 */

#include "Diffusion/EquilibriumPartitioning.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "PhaseField.h"
#include "Composition.h"
#include "Temperature.h"
#include "Thermodynamics/ThermodynamicPropertiesEQP.h"
#include "BoundaryConditions.h"
#include "NumericalMethods/SystemOfLinearEquationsSolvers.h"

namespace openphase
{

using namespace std;

struct Info_MinimumRelaxation
{
    int timesrelaxed;
    double Relaxation;
    iVector3 Position;
    std::vector<std::vector<std::vector<double> > > EquilibriumComp;
    std::vector<std::vector<std::vector<double> > > EquilibriumSlope;
    std::vector<double> TotalComposition;
    std::vector<double> Fractions;
    void reset()
    {
        timesrelaxed = 0;
        Relaxation = 1;
        Position.set_to_zero();
        EquilibriumComp.clear();
        EquilibriumSlope.clear();
        TotalComposition.clear();
        Fractions.clear();
    }
};

struct Info_MaximumDeviation
{
    double Deviation;
    iVector3 Position;
    int comp;
    std::vector<std::vector<double> > EquilibriumComp;
    std::vector<std::vector<double> > EquilibriumSlope;
    double TotalComposition;
    std::vector<double> Fractions;
    void reset()
    {
        Deviation = 0;
        Position.set_to_zero();
        comp = -1;
        EquilibriumComp.clear();
        EquilibriumSlope.clear();
        TotalComposition = 0;
        Fractions.clear();
    }
};

enum class PartitioningModels
{
    Eiken,
    Optimization,
    Steinbach,
    TwoStep
};

class OP_EXPORTS EquilibriumPartitioning_Impl : public OPObject                             ///< Composition partitioning solver
{
 public:

    EquilibriumPartitioning_Impl();
    ~EquilibriumPartitioning_Impl();
    EquilibriumPartitioning_Impl(Settings& locSettings,
                     const std::string InputFileName = DefaultInputFileName)    ///< Initializes storages, sets internal variables.
    {
        Initialize(locSettings);
        ReadInput(InputFileName);
    }

    void Initialize(Settings& locSettings, std::string ObjectNameSuffix = "") override; ///< Initializes storages, sets internal variables.
    void ReadInput(const std::string InputFileName) override;                   ///< Reads input parameters from the input file
    void ReadInput(std::stringstream& inp) override;                            ///< Reads input parameters from the input stream

    void CalculatePhaseConcentrations(PhaseField& Phi,
                                      Composition& Cx,
                                      Temperature& Tx,
                                      ThermodynamicPropertiesEQP& TP,
                                      BoundaryConditions& BC);                  ///< Calculates phase compositions for each present phase from the total composition

    size_t Nphases;                                                             ///< Number of thermodynamic phases
    size_t Ncomp;                                                               ///< Number of chemical components
    PartitioningModels Model;
    Table<bool> Partition;                                                      ///< Switches partitioning between pairs of phases on and off

 protected:
 private:

    void WriteStatistics();                                                     ///< Writes partitioning statistics in debug mode

    void CalculatePhaseConcentrationsTwoStep(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP);
    void CalculatePhaseConcentrationsOptimization(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP);
    void CalculatePhaseConcentrationsEiken(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP);
    void CalculatePhaseConcentrationsSteinbach(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP);

    Info_MinimumRelaxation MinRelax;
    Info_MaximumDeviation  MaxDev;
};

EquilibriumPartitioning::EquilibriumPartitioning() : impl_(new EquilibriumPartitioning_Impl) {}

EquilibriumPartitioning::~EquilibriumPartitioning() = default;

EquilibriumPartitioning::EquilibriumPartitioning(Settings& locSettings,
                 const std::string InputFileName)
                : impl_(new EquilibriumPartitioning_Impl)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void EquilibriumPartitioning::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    impl_->Initialize(locSettings, ObjectNameSuffix);
    thisclassname = impl_->thisclassname;
    thisobjectname = impl_->thisobjectname;
    initialized = impl_->initialized;
}

void EquilibriumPartitioning::ReadInput(const std::string InputFileName)
{
    impl_->ReadInput(InputFileName);
}

void EquilibriumPartitioning::ReadInput(std::stringstream& inp)
{
    impl_->ReadInput(inp);
}

void EquilibriumPartitioning::CalculatePhaseConcentrations(PhaseField& Phi,
                                        Composition& Cx,
                                        Temperature& Tx,
                                        ThermodynamicPropertiesEQP& TP,
                                        BoundaryConditions& BC)
{
    impl_->CalculatePhaseConcentrations(Phi,Cx,Tx,TP,BC);
}

EquilibriumPartitioning_Impl::EquilibriumPartitioning_Impl()
{
    Model = PartitioningModels::Optimization;
    Nphases = 0;
    Ncomp   = 0;
    MinRelax.reset();
    MaxDev.reset();
};

EquilibriumPartitioning_Impl::~EquilibriumPartitioning_Impl() = default;

void EquilibriumPartitioning_Impl::Initialize(Settings& locSettings, std::string ObjectNameSuffix)
{
    thisclassname = "EquilibriumPartitioning";
    thisobjectname = thisclassname + ObjectNameSuffix;

    Nphases = locSettings.Nphases;
    Ncomp   = locSettings.Ncomp;

    Partition.Allocate(Nphases, Nphases);

    initialized = true;
    ConsoleOutput::WriteStandard(thisclassname, "Initialized");
}

void EquilibriumPartitioning_Impl::ReadInput(const string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        ConsoleOutput::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        OP_Exit(EXIT_FAILURE);
    };

    std::stringstream data;
    data << inp.rdbuf();
    inp.close();

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteLineInsert(thisclassname+" input");
    ConsoleOutput::WriteStandard("Source", InputFileName);

    ReadInput(data);

    ConsoleOutput::WriteLine();
}

void EquilibriumPartitioning_Impl::ReadInput(stringstream& inp)
{
    int moduleLocation = FileInterface::FindModuleLocation(inp, thisclassname);

    // Reading partitioning switches
    for(size_t n =   0; n < Nphases-1; n++)
    for(size_t m = n+1; m < Nphases  ; m++)
    //if(n != m)
    {
        string counter = to_string(n) + "_" + to_string(m);
        Partition(n,m) = FileInterface::ReadParameterB(inp, moduleLocation, string("Partition_") + counter, false, true);
        Partition(m,n) = Partition(n,m);
    }

    string tmpModel = FileInterface::ReadParameterK(inp, moduleLocation, "PartitioningModel", false, "OPTIMIZATION");
    if (tmpModel == "OPTIMIZATION")
    {
        Model = PartitioningModels::Optimization;
    }
    else if (tmpModel == "EIKEN")
    {
        Model = PartitioningModels::Eiken;
    }
    else if (tmpModel == "STEINBACH")
    {
        Model = PartitioningModels::Steinbach;
    }
    else if (tmpModel == "TWOSTEP")
    {
        Model = PartitioningModels::TwoStep;
    }
    else
    {
        ConsoleOutput::WriteWarning("Wrong partitioning model specified!", thisclassname, "ReadInput()");
        OP_Exit(EXIT_FAILURE);
    }

    ConsoleOutput::WriteLine();
    ConsoleOutput::WriteBlankLine();
}

void EquilibriumPartitioning_Impl::CalculatePhaseConcentrations(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP,
                                            BoundaryConditions& BC)
{
#ifdef DEBUG
    MinRelax.EquilibriumComp.resize(Cx.Ncomp);
    MinRelax.EquilibriumSlope.resize(Cx.Ncomp);
    MinRelax.TotalComposition.resize(Cx.Ncomp);
    MinRelax.Fractions.resize(Cx.Nphases);
    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
    {
        MinRelax.EquilibriumComp[comp].resize(Cx.Nphases);
        MinRelax.EquilibriumSlope[comp].resize(Cx.Nphases); 
        for(size_t m = 0; m < Cx.Nphases; m++)
        {
            MinRelax.EquilibriumComp[comp][m].resize(Cx.Nphases);
            MinRelax.EquilibriumSlope[comp][m].resize(Cx.Nphases); 
        }
    }
    MaxDev.EquilibriumComp.resize(Cx.Nphases);
    MaxDev.EquilibriumSlope.resize(Cx.Nphases);
    MaxDev.Fractions.resize(Cx.Nphases);
    for(size_t m = 0; m < Cx.Nphases; m++)
    {
        MaxDev.EquilibriumComp[m].resize(Cx.Nphases);
        MaxDev.EquilibriumSlope[m].resize(Cx.Nphases); 
    }
#endif
    /** This function calculates the partitioning from the total composition
    to the separate phase compositions.*/
        
    switch(Model)
    {
        case PartitioningModels::Optimization:
        {
            CalculatePhaseConcentrationsOptimization(Phi,Cx,Tx,TP);
            break;
        }
        case PartitioningModels::Steinbach:
        {
            CalculatePhaseConcentrationsSteinbach(Phi,Cx,Tx,TP);
            break;
        }
        case PartitioningModels::Eiken:
        {
            CalculatePhaseConcentrationsEiken(Phi,Cx,Tx,TP);
            break;
        }
        case PartitioningModels::TwoStep:
        {
            CalculatePhaseConcentrationsTwoStep(Phi,Cx,Tx,TP);
            break;
        }
        default: std::cerr << "No Valid Partitioning Scheme selected." << std::endl;
    }

    Cx.SetBoundaryConditions(BC);
#ifdef DEBUG
   // WriteStatistics();
#endif
}

void EquilibriumPartitioning_Impl::WriteStatistics()
{
    std::cout << "EquilibriumPartitioning Statistics" << std::endl;
    std::cout << "Minimum Relaxation = " << MinRelax.Relaxation << std::endl;
    std::cout << "At position " << MinRelax.Position.print() << std::endl;
    for(size_t comp = 0; comp < MinRelax.EquilibriumComp.size(); comp++)
    {
        for(size_t m = 0; m < MinRelax.EquilibriumComp[comp].size(); m++)
        {
            for(size_t n = 0; n < MinRelax.EquilibriumComp[comp][m].size(); n++)
            {
                std::cout << "EquilibriumComp[" << comp << "][" << m << "][" << n << "] = " << MinRelax.EquilibriumComp[comp][m][n]<< std::endl;
            }
        }  
    }
    for(size_t comp = 0; comp < MinRelax.EquilibriumSlope.size(); comp++)
    {
        for(size_t m = 0; m < MinRelax.EquilibriumSlope[comp].size(); m++)
        {
            for(size_t n = 0; n < MinRelax.EquilibriumSlope[comp][m].size(); n++)
            {
                std::cout << "EquilibriumSlope[" << comp << "][" << m << "][" << n << "] = " << MinRelax.EquilibriumSlope[comp][m][n] << std::endl;
            }
        }  
    }
    for(size_t comp = 0; comp < MinRelax.TotalComposition.size(); comp++)
    {
        std::cout << "TotalComposition[" << comp << "] = " << MinRelax.TotalComposition[comp] << std::endl;
    }
    for(size_t m = 0; m < MinRelax.Fractions.size(); m++)
    {
        std::cout << "Fractions[" << m << "] = " << MinRelax.Fractions[m] << std::endl;
    }  
    std::cout << "Maximum Deviation = " << MaxDev.Deviation << std::endl;
    std::cout << "At position " << MaxDev.Position.print()  << std::endl;
    std::cout << "Component " << MaxDev.comp << std::endl;
    for(size_t m = 0; m < MaxDev.EquilibriumComp.size(); m++)
    {
        for(size_t n = 0; n < MaxDev.EquilibriumComp[m].size(); n++)
        {
            std::cout << "EquilibriumComp[" << m << "][" << n << "] = " << MaxDev.EquilibriumComp[m][n]<< std::endl;
        }
    }  
    for(size_t m = 0; m < MaxDev.EquilibriumSlope.size(); m++)
    {
        for(size_t n = 0; n < MaxDev.EquilibriumSlope[m].size(); n++)
        {
            std::cout << "EquilibriumSlope[" << m << "][" << n << "] = " << MaxDev.EquilibriumSlope[m][n] << std::endl;
        }
    }  
    std::cout << "TotalComposition = " << MaxDev.TotalComposition << std::endl;
    for(size_t m = 0; m < MaxDev.Fractions.size(); m++)
    {
        std::cout << "Fractions[" << m << "] = " << MaxDev.Fractions[m] << std::endl;
    }  
    MinRelax.reset();
    MaxDev.reset();
}

void EquilibriumPartitioning_Impl::CalculatePhaseConcentrationsTwoStep(PhaseField& Phi,
                                          Composition& Cx,
                                          Temperature& Tx,
                                          ThermodynamicPropertiesEQP& TP)
{
    /** This function calculates the partitioning from the total composition
    to the separate phase compositions in each point.*/

    Tensor<double,2> locPhaseMoleFractions({Cx.Nphases,Cx.Ncomp});
    locPhaseMoleFractions.set_to_zero();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,firstprivate(locPhaseMoleFractions))
    {
        Cx.MoleFractions(i,j,k).set_to_zero();

        if(Phi.Fields(i,j,k).interface())
        {
            const Tensor<EquilibriumData,2>& locExtrapolationData = TP.ExtrapolationData(i,j,k);

            /* In case of an interface, first check for multi-phase region.
            Checking appearance of every phase not only by its fraction
            but also checking if there is a seed in the grid point to allow
            partitioning for nuclei*/

            bool single_phase = true;

            for(auto alpha  = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t m = Phi.FieldsProperties[alpha->index].Phase;
                if(Phi.Fractions(i,j,k,{m}) != 1.0)
                {
                    single_phase = false;
                    break;
                }
            }

            if(single_phase)
            {
                /* If only a single thermodynamic phase is present, don't partition
                and set the phase composition to the total composition. */

                for(size_t m = 0; m < Cx.Nphases; m++)
                if(Phi.Fractions(i,j,k,{m}) != 0.0)
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    locPhaseMoleFractions({m,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                }
            }
            else
            {
                /* If more than one phase is present, partitioning has to be
                calculated for each phase pair and for each component
                separately.*/

                double relaxation = 1.0;
                bool outoflimits = true;

                do
                {
                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    {
                        // Set equilibrium concentrations

                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            locPhaseMoleFractions({m,comp}) = 0.0;

                            double denominator = 0.0;
                            for(size_t n = 0; n < Cx.Nphases; n++)
                            if (n != m)
                            {
                                if (locExtrapolationData({m,n}).is_set and
                                    Partition(m,n))
                                {
                                    locPhaseMoleFractions({m,comp}) += Phi.Fractions(i,j,k,{n})*(locExtrapolationData({m,n}).eq_composition_extra(comp,Tx(i,j,k))*relaxation + Cx.MoleFractionsTotal(i,j,k,{comp})*(1.0-relaxation));
                                    denominator += Phi.Fractions(i,j,k,{n});
                                }
                                else
                                {
                                    locPhaseMoleFractions({m,comp}) += Phi.Fractions(i,j,k,{n})*Cx.MoleFractionsTotal(i,j,k,{comp});
                                    denominator += Phi.Fractions(i,j,k,{n});
                                }
                            }
                            if (denominator != 0.0)
                            {
                                locPhaseMoleFractions({m,comp}) /= denominator;
                            }
                            else
                            {
                                locPhaseMoleFractions({m,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                            }
                        }
                        // Calculate composition delta
                        double Delta_c_tot = Cx.MoleFractionsTotal(i,j,k,{comp});
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            Delta_c_tot -= Phi.Fractions(i,j,k,{m})*locPhaseMoleFractions({m,comp});
                        }

                        // Calculate phase concentrations
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(Phi.Fractions(i,j,k,{m}) != 0.0 and
                          !Cx.Phase[m].Component[comp].isStoichiometric)
                        {
                            double numerator   = Phi.Fractions(i,j,k,{m});
                            double denominator = Phi.Fractions(i,j,k,{m});

                            for(size_t n = 0; n < Cx.Nphases; n++)
                            if(m != n and
                               locExtrapolationData({n,m}).is_set and
                               Partition(n,m) and
                               !Cx.Phase[n].Component[comp].isStoichiometric)
                            {
                                numerator   += Phi.Fractions(i,j,k,{n});
                                denominator += Phi.Fractions(i,j,k,{n})*(relaxation*locExtrapolationData({m,n}).part_coefficient({comp}) + (1.0 - relaxation)*1.0);
                            }
                            if(denominator != 0.0)
                            {
                                numerator /= denominator;
                            }
                            else
                            {
                                numerator = 0.0;
                            }
                            locPhaseMoleFractions({m,comp}) += Delta_c_tot*numerator;
                        }

                        // Calculate concentration deviation
                        double deviation = Cx.MoleFractionsTotal(i,j,k,{comp});
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            deviation -= locPhaseMoleFractions({m,comp})*Phi.Fractions(i,j,k,{m});
                        }

                        // Check and adjust mass conservation
#ifdef DEBUG
                        if (std::fabs(deviation) > MaxDev.Deviation)
                        {
                            #pragma omp critical
                            {
                                MaxDev.Deviation = std::fabs(deviation);
                                MaxDev.Position = {(int)i,(int)j,(int)k};
                                MaxDev.comp = comp;
                                for(size_t m = 0; m < Cx.Nphases; m++)
                                {
                                    for(size_t n = 0; n < Cx.Nphases; n++)
                                    if(m != n)
                                    {
                                        MaxDev.EquilibriumComp[m][n]  = locExtrapolationData({m,n}).eq_composition_extra(comp,Tx(i,j,k));
                                        MaxDev.EquilibriumSlope[m][n] = locExtrapolationData({m,n}).eq_slope({comp});
                                    }
                                    MaxDev.Fractions[m] = Phi.Fractions(i,j,k,{m});
                                }
                                MaxDev.TotalComposition = Cx.MoleFractionsTotal(i,j,k,{comp});
                            }
                        }
#endif
                        if (fabs(deviation) > DBL_EPSILON)
                        {
                            double deltaCstoich = 0.0;
                            for(size_t m = 0; m < Cx.Nphases; m++)
                            if(Cx.Phase[m].Component[comp].isStoichiometric)
                            {
                                deltaCstoich += Cx.MoleFractions(i,j,k,{m,comp})*Phi.Fractions(i,j,k,{m});
                            }
                            double c_total_sans_stoich = Cx.MoleFractionsTotal(i,j,k,{comp}) - deltaCstoich;
                            double correction_factor = c_total_sans_stoich - deviation;
                            for(size_t m = 0; m < Cx.Nphases; m++)
                            if(Phi.Fractions(i,j,k,{m}) != 0.0 and
                               !Cx.Phase[m].Component[comp].isStoichiometric and
                               fabs(correction_factor) > DBL_EPSILON)
                            {
                                locPhaseMoleFractions({m,comp}) *= c_total_sans_stoich/correction_factor;
                            }
                        }
                    }

                    // Check composition limits
                    outoflimits = false;
                    for(size_t m = 0; m < Cx.Nphases; m++)
                    if(Phi.Fractions(i,j,k,{m}) != 0.0)
                    {
                        for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                        if (locPhaseMoleFractions({m,comp}) < Cx.Phase[m].Component[comp].Min or
                            locPhaseMoleFractions({m,comp}) > Cx.Phase[m].Component[comp].Max)
                        {
                            outoflimits = true;
                            break;
                        }
                    }
                    if (outoflimits and relaxation > DBL_EPSILON)
                    {
                        // Repeat partitioning with reduced partitioning coefficients
                        relaxation *= 0.9;
#ifdef DEBUG
                        if (relaxation < MinRelax.Relaxation)
                        {
                            #pragma omp critical
                            {
                                MinRelax.Relaxation = relaxation;
                                MinRelax.Position = {(int)i,(int)j,(int)k};
                                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                                {
                                    for(size_t m = 0; m < Cx.Nphases; m++)
                                    {
                                        for(size_t n = 0; n < Cx.Nphases; n++)
                                        if(m != n)
                                        {
                                            MinRelax.EquilibriumComp[comp][m][n]  = locExtrapolationData({m,n}).eq_composition_extra(comp,Tx(i,j,k));
                                            MinRelax.EquilibriumSlope[comp][m][n] = locExtrapolationData({m,n}).eq_slope({comp});
                                        }
                                    }
                                    MinRelax.TotalComposition[comp] = Cx.MoleFractionsTotal(i,j,k,{comp});
                                }
                                for(size_t m = 0; m < Cx.Nphases; m++)
                                {
                                    MinRelax.Fractions[m] = Phi.Fractions(i,j,k,{m});
                                }
                            }
                        }
#endif
                    }
                    else
                    {
                        // Limit composition to prescribed limits and calculate reference element concentration
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(Phi.Fractions(i,j,k,{m}) != 0.0)
                        {
                            double unity = 1.0;
                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(!Cx.Phase[m].Component[comp].isReference)
                            {
                                if(locPhaseMoleFractions({m,comp}) < Cx.Phase[m].Component[comp].Min)
                                {
                                    locPhaseMoleFractions({m,comp}) = Cx.Phase[m].Component[comp].Min;
                                }

                                if(locPhaseMoleFractions({m,comp}) > Cx.Phase[m].Component[comp].Max)
                                {
                                    locPhaseMoleFractions({m,comp}) = Cx.Phase[m].Component[comp].Max;
                                }
                                unity -= locPhaseMoleFractions({m,comp});
                            }

                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(Cx.Phase[m].Component[comp].isReference)
                            {
                                locPhaseMoleFractions({m,comp})
                                        = max(min(unity,Cx.Phase[m].Component[comp].Max)
                                                       ,Cx.Phase[m].Component[comp].Min);
                            }
                        }
                    }
                } while (outoflimits and relaxation > DBL_EPSILON);
            }
            Cx.MoleFractions(i,j,k) = locPhaseMoleFractions;
        }
        else
        {
            /* If only one phase is present, the phase composition will be set to
            the total composition. All other phase compositions were set to zero
            before.*/

            size_t pIndex
                = Phi.FieldsProperties[Phi.Fields(i,j,k).front().index].Phase;

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                Cx.MoleFractions(i,j,k,{pIndex, comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitioning_Impl::CalculatePhaseConcentrationsOptimization(PhaseField& Phi,
                                        Composition& Cx,
                                        Temperature& Tx,
                                        ThermodynamicPropertiesEQP& TP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        Cx.MoleFractions(i,j,k).set_to_zero();

        if(Phi.Fields(i,j,k).interface())
        {
            /* In case of an interface, first check for dual-phase re
            Next: Checking appearance of every phase not by its fraction but by
            its entry in Phi.Fields (to allow partitioning for nuclei)*/

            double locTx = Tx(i,j,k);

            Tensor<bool,1> ispresent({Cx.Nphases});

            ispresent.set_to_value(false);

            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t m = Phi.FieldsProperties[alpha->index].Phase;
                ispresent({m}) = true;
            }

            int num_phases = 0;
            for(size_t m = 0; m < Cx.Nphases; m++)
            if(ispresent({m}))
            {
                num_phases++;
            }
            double relaxation = 1.0;
            bool outoflimits = true;
            /* If at least two phases are present, partitioning has to be
            calculated for each phase pair and for each component
            separately.*/

            Tensor<double,2> Phase({Cx.Nphases,Cx.Ncomp});
            Phase.set_to_zero();

            if(num_phases == 1)
            {
                /* If only a single thermodynamic phase is present, don't partition,
                but set the total composition to the phase composition. */

                for(size_t m = 0; m < Cx.Nphases; m++)
                if(ispresent({m}))
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    Phase({m,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                }
            }
            else
            {
                do
                {
                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    {
                        double maxSlope = 1.0;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        for(size_t n = 0; n < Cx.Nphases; n++)
                        if(m != n)
                        {
                            maxSlope = std::max(maxSlope, std::fabs(TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp})));
                        }
                        double r = 1e8*maxSlope*relaxation;
                        if ( maxSlope == 0.0) std::cout << "No Slopes available!" << std::endl;
                        Matrix<double> A;
                        std::vector<double> B(Cx.Nphases,0.0);
                        std::vector<double> X(Cx.Nphases,0.0);
                        A.Allocate(Cx.Nphases, Cx.Nphases);
                        A.set_to_value(0.0);
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(!Cx.Phase[m].Component[comp].isStoichiometric and Phi.Fractions(i,j,k,{m}) > 0.0)
                        {
                            A(m,m) = r*Phi.Fractions(i,j,k,{m})*Phi.Fractions(i,j,k,{m});
                            B[m] = r*Cx.MoleFractionsTotal(i,j,k,{comp})*Phi.Fractions(i,j,k,{m});
                            for(size_t n = 0; n < Cx.Nphases; n++)
                            if(ispresent({n}) and Phi.Fractions(i,j,k,{n}) > 0.0)
                            if(m != n)
                            {
                                if (TP.ExtrapolationData(i,j,k)({m,n}).is_set)
                                {
                                    A(m,m) += std::fabs(TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp}));
                                }
                                A(m,n) += r*Phi.Fractions(i,j,k,{m})*Phi.Fractions(i,j,k,{n});
                                if (TP.ExtrapolationData(i,j,k)({m,n}).is_set)
                                {
                                    B[m] += TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx)*std::fabs(TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp}));
                                }
                            }
                        }
                        else
                        {
                            A(m,m) = 1.0;
                            for(size_t n = 0; n < Cx.Nphases; n++)
                            if(m != n)
                            {
                                if (TP.ExtrapolationData(i,j,k)({m,n}).is_set)
                                {
                                    B[m] = TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                    break;
                                }
                                else
                                {
                                    B[m] = Cx.MoleFractionsTotal(i,j,k,{comp});
                                }
                            }
                        }
                        /*std::cout << "A = " << std::endl;
                        std::cout << A.print() << std::endl;
                        std::cout << "B = " << std::endl;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            std::cout << B[m] << std::endl;
                        }*/
                        SystemOfLinearEquationsSolvers::Gauss(A,X,B);
                        /*std::cout << "X = " << std::endl;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            std::cout << X[m] << std::endl;
                        }*/
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            Phase({m,comp}) = X[m];
                            //std::cout << "Phase({"<<m<<","<<comp<<"}) = " << Phase({m,comp}) << std::endl;
                        }
                        double total = 0.0;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            total += Phase({m,comp})*Phi.Fractions(i,j,k,{m});
                        }
                        double deviation = Cx.MoleFractionsTotal(i,j,k,{comp}) - total;

                        total = 0.0;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}) and Phi.Fractions(i,j,k,{m}) > 0.0)
                        {
                            if (Cx.MoleFractionsTotal(i,j,k,{comp}) - deviation != 0.0)
                            {
                                Phase({m,comp}) *= Cx.MoleFractionsTotal(i,j,k,{comp})/(Cx.MoleFractionsTotal(i,j,k,{comp})-deviation);
                            }
                            total += Phase({m,comp})*Phi.Fractions(i,j,k,{m});
                        }
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}))
                        {
                            double unity = 1.0;
                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(!Cx.Phase[m].Component[comp].isReference)
                            {
                                if(Phase({m,comp}) < Cx.Phase[m].Component[comp].Min)
                                {
                                    Phase({m,comp}) = Cx.Phase[m].Component[comp].Min;
                                }

                                if(Phase({m,comp}) > Cx.Phase[m].Component[comp].Max)
                                {
                                    Phase({m,comp}) = Cx.Phase[m].Component[comp].Max;
                                }

                                unity -= Phase({m,comp});
                            }

                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(Cx.Phase[m].Component[comp].isReference)
                            {
                                Phase({m,comp})
                                        = max(min(unity,Cx.Phase[m].Component[comp].Max)
                                                       ,Cx.Phase[m].Component[comp].Min);
                            }
                        }

#ifdef DEBUG
                        if (std::fabs(deviation) > MaxDev.Deviation)
                        {
                            #pragma omp critical
                            {
                                MaxDev.Deviation = std::fabs(deviation);
                                MaxDev.Position = {(int)i,(int)j,(int)k};
                                MaxDev.comp = comp;
                                for(size_t m = 0; m < Cx.Nphases; m++)
                                {
                                    for(size_t n = 0; n < Cx.Nphases; n++)
                                    if(m != n)
                                    {
                                        MaxDev.EquilibriumComp[m][n]  = TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                        MaxDev.EquilibriumSlope[m][n] = TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                                    }
                                    MaxDev.Fractions[m] = Phi.Fractions(i,j,k,{m});
                                }
                                MaxDev.TotalComposition = Cx.MoleFractionsTotal(i,j,k,{comp});
                            }
                        }
#endif
                    }
                    outoflimits = false;
                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    {
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}))
                        {
                            if (Phase({m,comp}) < 0.0 or Phase({m,comp}) > 1.0)
                            {
                                outoflimits = true;
                            }
                        }
                    }
                    if (outoflimits)
                    {
                        relaxation *= 0.1;
                    }
                } while (outoflimits);
            }
            Cx.MoleFractions(i,j,k) = Phase;
        }
        else
        {
            /* If only one phase is present, the phase composition will be set to
            the total composition. All other phase compositions were set to zero
            before.*/

            size_t pIndex
                = Phi.FieldsProperties[Phi.Fields(i,j,k).front().index].Phase;

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                Cx.MoleFractions(i,j,k,{pIndex, comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void EquilibriumPartitioning_Impl::CalculatePhaseConcentrationsSteinbach(PhaseField& Phi,
                                            Composition& Cx,
                                            Temperature& Tx,
                                            ThermodynamicPropertiesEQP& TP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        Cx.MoleFractions(i,j,k).set_to_zero();

        if(Phi.Fields(i,j,k).interface())
        {
            /* In case of an interface, first check for dual-phase re
            Next: Checking appearance of every phase not by its fraction but by
            its entry in Phi.Fields (to allow partitioning for nuclei)*/

            double locTx = Tx(i,j,k);

            Tensor<bool,1> ispresent({Cx.Nphases});
            ispresent.set_to_value(false);

            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t m = Phi.FieldsProperties[alpha->index].Phase;
                ispresent({m}) = true;
            }

            int num_phases = 0;
            for(size_t m = 0; m < Cx.Nphases; m++)
            if(ispresent({m}))
            {
                num_phases++;
            }

            /* If at least two phases are present, partitioning has to be
            calculated for each phase pair and for each component
            separately.*/

            Tensor<double,2> Phase({Cx.Nphases,Cx.Ncomp});
            Phase.set_to_zero();

            if(num_phases == 1)
            {
                /* If only a single thermodynamic phase is present, don't partition,
                but set the total composition to the phase composition. */

                for(size_t m = 0; m < Cx.Nphases; m++)
                if(ispresent({m}))
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    Phase({m,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                }
            }
            else
            {
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    double maxSlope = 1.0;
                    for(size_t m = 0; m < Cx.Nphases; m++)
                    for(size_t n = 0; n < Cx.Nphases; n++)
                    if(m != n)
                    {
                        maxSlope = std::max(maxSlope, std::fabs(TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp})));
                    }
                    //double r = 1e8*maxSlope;
                    if ( maxSlope == 0.0) std::cout << "No Slopes available!" << std::endl;
                    Matrix<double> A;
                    std::vector<double> B(Cx.Nphases+1,0.0);
                    std::vector<double> X(Cx.Nphases+1,0.0);
                    A.Allocate(Cx.Nphases+1, Cx.Nphases+1);
                    A.set_to_value(0.);
                    maxSlope = 0;
                    for(size_t m = 0; m < Cx.Nphases; m++)
                    {
                        double slope = 1.0;
                        double ceq = 0.0;
                        double fractions = 0.0;
                        for(size_t n = 0; n < Cx.Nphases; n++)
                        if (TP.ExtrapolationData(i,j,k)({m,n}).is_set)
                        {
                            slope += Phi.Fractions(i,j,k,{n})*TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                            ceq += Phi.Fractions(i,j,k,{n})*TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                            fractions += Phi.Fractions(i,j,k,{n});
                        }
                        else
                        {
                            slope += Phi.Fractions(i,j,k,{n})*TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                            ceq += Phi.Fractions(i,j,k,{n})*Cx.MoleFractionsTotal(i,j,k,{comp});
                            fractions += Phi.Fractions(i,j,k,{n});
                        }
                        if (fractions > 0)
                        {
                            slope /= fractions;
                            ceq /= fractions;
                        }

                        A(m,m) = 1.0;
                        if (slope != 0)
                        {
                            A(m,Cx.Nphases) = 1.0/slope;
                        }
                        maxSlope = std::max(maxSlope, std::fabs(slope));
                        B[m] = ceq;
                        A(Cx.Nphases,m) = Phi.Fractions(i,j,k,{m});
                    }
                    B[Cx.Nphases] = Cx.MoleFractionsTotal(i,j,k,{comp});
                    /*std::cout << "A = " << std::endl;
                    std::cout << A.print() << std::endl;
                    std::cout << "B = " << std::endl;
                    for(size_t m = 0; m < Cx.Nphases+1; m++)
                    {
                        std::cout << B[m] << std::endl;
                    }*/
                    if (maxSlope > 0)
                    {
                        SystemOfLinearEquationsSolvers::Gauss(A,X,B);
                    }
                    else
                    {
                        X = B;
                    }
                    /*std::cout << "X = " << std::endl;
                    for(size_t m = 0; m < Cx.Nphases+1; m++)
                    {
                        std::cout << X[m] << std::endl;
                    }*/
                    for(size_t m = 0; m < Cx.Nphases; m++)
                    {
                        Phase({m,comp}) = X[m];
                        //std::cout << "Phase({"<<m<<","<<comp<<"}) = " << Phase({m,comp}) << std::endl;
                    }
#ifdef DEBUG
                    double total = 0.0;
                    for(size_t m = 0; m < Cx.Nphases; m++)
                    {
                        total += Phase({m,comp})*Phi.Fractions(i,j,k,{m});
                    }

                    double deviation = Cx.MoleFractionsTotal(i,j,k,{comp}) - total;
                    if (std::fabs(deviation) > MaxDev.Deviation)
                    {
                        #pragma omp critical
                        {
                            MaxDev.Deviation = std::fabs(deviation);
                            MaxDev.Position = {(int)i,(int)j,(int)k};
                            MaxDev.comp = comp;
                            for(size_t m = 0; m < Cx.Nphases; m++)
                            {
                                for(size_t n = 0; n < Cx.Nphases; n++)
                                {
                                    MaxDev.EquilibriumComp[m][n] = TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                    MaxDev.EquilibriumSlope[m][n] = TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                                }
                                MaxDev.Fractions[m] = Phi.Fractions(i,j,k,{m});
                            }
                            MaxDev.TotalComposition = Cx.MoleFractionsTotal(i,j,k,{comp});
                        }
                    }
#endif
                }
            }
            Cx.MoleFractions(i,j,k) = Phase;
        }
        else
        {
            /* If only one phase is present, the phase composition will be set to
            the total composition. All other phase compositions were set to zero
            before.*/

            size_t pIndex
                = Phi.FieldsProperties[Phi.Fields(i,j,k).front().index].Phase;

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                Cx.MoleFractions(i,j,k,{pIndex, comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void EquilibriumPartitioning_Impl::CalculatePhaseConcentrationsEiken(PhaseField& Phi,
                                        Composition& Cx,
                                        Temperature& Tx,
                                        ThermodynamicPropertiesEQP& TP)
{
    /** This function calculates the partitioning from the total composition
    to the separate phase compositions in each point.*/
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.MoleFractionsTotal,0,)
    {
        Cx.MoleFractions(i,j,k).set_to_zero();

        if(Phi.Fields(i,j,k).interface())
        {
            /* In case of an interface, first check for dual-phase re
            Next: Checking appearance of every phase not by its fraction but by
            its entry in Phi.Fields (to allow partitioning for nuclei)*/

            double locTx = Tx(i,j,k);

            Tensor<bool,1> ispresent({Cx.Nphases});
            ispresent.set_to_value(false);

            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                size_t m = Phi.FieldsProperties[alpha->index].Phase;
                ispresent({m}) = true;
            }

            int num_phases = 0;
            for(size_t m = 0; m < Cx.Nphases; m++)
            if(ispresent({m}))
            {
                num_phases++;
            }

            /* If at least two phases are present, partitioning has to be
            calculated for each phase pair and for each component
            separately.*/

            Tensor<double,2> Phase({Cx.Nphases,Cx.Ncomp});
            Phase.set_to_zero();
            double relaxation = 1.0;
            bool outoflimits = true;

            if(num_phases == 1)
            {
                /* If only a single thermodynamic phase is present, don't partition,
                but set the total composition to the phase composition. */

                for(size_t m = 0; m < Cx.Nphases; m++)
                if(ispresent({m}))
                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                {
                    Phase({m,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                }
            }
            else
            {
                do
                {
                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    {
                        Tensor<double,1> eqCx({Cx.Nphases});
                        eqCx.set_to_zero();

                        std::vector<double> c_eq(Cx.Nphases,0.0);

                        for(size_t m = 0; m < Cx.Nphases; m++)
                        {
                            double denominator = 0.0;
                            for(size_t n = 0; n < Cx.Nphases; n++)
                            if (n != m)
                            {
                                if (TP.ExtrapolationData(i,j,k)({m,n}).is_set)
                                {
                                    c_eq[m] += Phi.Fractions(i,j,k,{n})*(TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx)*relaxation + Cx.MoleFractionsTotal(i,j,k,{comp})*(1.-relaxation));
                                    denominator +=  Phi.Fractions(i,j,k,{n});
                                }
                                else
                                {
                                    c_eq[m] += Phi.Fractions(i,j,k,{n})*Cx.MoleFractionsTotal(i,j,k,{comp});
                                    denominator +=  Phi.Fractions(i,j,k,{n});
                                }
                            }
                            if (denominator > 0)
                            {
                                c_eq[m] /= denominator;
                            }
                            else
                            {
                                c_eq[m] = Cx.MoleFractionsTotal(i,j,k,{comp});
                            }
                        }

                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}))
                        {
                            if(!Cx.Phase[m].Component[comp].isStoichiometric)
                            {
                                double numerator = Cx.MoleFractionsTotal(i,j,k,{comp});                             ///< Sum of phase-pair equilibrium composition for alpha weighted with phase-fraction of beta phases
                                double denominator = 0.0;                           ///< Sum of phase-fractions of beta phases

                                for(size_t n = 0; n < Cx.Nphases; n++)
                                if(ispresent({n}))
                                {
                                    double ecm = 0.0;                                ///< Equilibrium composition for phase n. If no thermodynamic properties are available, set to Cx.Total
                                    double ecn = 0.0;
                                    if(TP.ExtrapolationData(i,j,k)({n,m}).is_set)
                                    {
                                        ecm = c_eq[m];
                                        ecn = c_eq[n];
                                    }
                                    else
                                    {
                                        ecm = Cx.MoleFractionsTotal(i,j,k,{comp});
                                        ecn = Cx.MoleFractionsTotal(i,j,k,{comp});
                                    }

                                    if(!Cx.Phase[n].Component[comp].isStoichiometric)
                                    {
                                        numerator -= Phi.Fractions(i,j,k,{n})*(ecn - ecm*TP.ExtrapolationData(i,j,k)({m,n}).part_coefficient({comp}));
                                        denominator += Phi.Fractions(i,j,k,{n})*TP.ExtrapolationData(i,j,k)({m,n}).part_coefficient({comp});
                                    }
                                    else
                                    {
                                        numerator -= Phi.Fractions(i,j,k,{n})*(ecn);
                                    }
                                }

                                if(denominator > DBL_EPSILON)
                                {
                                    eqCx({m}) = numerator/denominator;
                                }
                                else
                                {
                                    eqCx({m}) = Cx.MoleFractionsTotal(i,j,k,{comp});
                                }
                            }
                            else
                            {
                                double eqcomp = 0.0;
                                int count = 0;
                                for(size_t n = 0; n < Cx.Nphases; n++)
                                if(ispresent({n}) and n != m)
                                {
                                    if(TP.ExtrapolationData(i,j,k)({n,m}).is_set)
                                    {
                                        eqcomp += TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                    }
                                    else
                                    {
                                        eqcomp += Cx.MoleFractionsTotal(i,j,k,{comp});
                                    }
                                    ++count;
                                }
                                eqCx({m}) = eqcomp/count;
                            }
                        }
                        double total = 0.0;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}) and Phi.Fractions(i,j,k,{m}) > 0.0)
                        {
                            Phase({m,comp}) = eqCx({m});
                            total += eqCx({m})*Phi.Fractions(i,j,k,{m});
                        }
                        else
                        {
                            Phase({m,comp}) = c_eq[m];
                        } 
                        double deviation = Cx.MoleFractionsTotal(i,j,k,{comp}) - total;
#ifdef DEBUG
                        if (std::fabs(deviation) > MaxDev.Deviation)
                        {
                            #pragma omp critical
                            {
                                MaxDev.Deviation = std::fabs(deviation);
                                MaxDev.Position = {(int)i,(int)j,(int)k};
                                MaxDev.comp = comp;
                                for(size_t m = 0; m < Cx.Nphases; m++)
                                {
                                    for(size_t n = 0; n < Cx.Nphases; n++)
                                    {
                                        MaxDev.EquilibriumComp[m][n]  = TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                        MaxDev.EquilibriumSlope[m][n] = TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                                    }
                                    MaxDev.Fractions[m] = Phi.Fractions(i,j,k,{m});
                                }
                                MaxDev.TotalComposition = Cx.MoleFractionsTotal(i,j,k,{comp});
                            }
                        }
#endif
                        total = 0.0;
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}) and Phi.Fractions(i,j,k,{m}) > 0.0)
                        {
                            if (Cx.MoleFractionsTotal(i,j,k,{comp}) - deviation != 0.0)
                            {
                                Phase({m,comp}) *= Cx.MoleFractionsTotal(i,j,k,{comp})/(Cx.MoleFractionsTotal(i,j,k,{comp}) - deviation);
                            }
                            total += Phase({m,comp})*Phi.Fractions(i,j,k,{m});
                        }
                    }

                    outoflimits = false;
                    for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                    {
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}))
                        {
                            if (Phase({m,comp}) < 0.0 or Phase({m,comp}) > 1.0)
                            {
                                outoflimits = true;
                            }
                        }
                    }
                    if (outoflimits)
                    {
                        relaxation *= 0.9;
#ifdef DEBUG
                        if (relaxation < MinRelax.Relaxation)
                        {
                            #pragma omp critical
                            {
                                MinRelax.Relaxation = relaxation;
                                MinRelax.Position = {(int)i,(int)j,(int)k};
                                for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                                {
                                    for(size_t m = 0; m < Cx.Nphases; m++)
                                    {
                                        for(size_t n = 0; n < Cx.Nphases; n++)
                                        {
                                            MinRelax.EquilibriumComp[comp][m][n]  = TP.ExtrapolationData(i,j,k)({m,n}).eq_composition_extra(comp,locTx);
                                            MinRelax.EquilibriumSlope[comp][m][n] = TP.ExtrapolationData(i,j,k)({m,n}).eq_slope({comp});
                                        }
                                    }
                                    MinRelax.TotalComposition[comp] = Cx.MoleFractionsTotal(i,j,k,{comp});
                                }
                                for(size_t m = 0; m < Cx.Nphases; m++)
                                {
                                    MinRelax.Fractions[m] = Phi.Fractions(i,j,k,{m});
                                }  
                            }
                        }
#endif
                    }
                    else
                    {
                        for(size_t m = 0; m < Cx.Nphases; m++)
                        if(ispresent({m}))
                        {
                            double unity = 1.0;
                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(!Cx.Phase[m].Component[comp].isReference)
                            {
                                if(Phase({m,comp}) < Cx.Phase[m].Component[comp].Min)
                                {
                                    Phase({m,comp}) = Cx.Phase[m].Component[comp].Min;
                                }

                                if(Phase({m,comp}) > Cx.Phase[m].Component[comp].Max)
                                {
                                    Phase({m,comp}) = Cx.Phase[m].Component[comp].Max;
                                }

                                unity -= Phase({m,comp});
                            }

                            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
                            if(Cx.Phase[m].Component[comp].isReference)
                            {
                                Phase({m,comp})
                                        = max(min(unity,Cx.Phase[m].Component[comp].Max)
                                                       ,Cx.Phase[m].Component[comp].Min);
                            }
                        }
                    }
                } while (outoflimits);
                Cx.MoleFractions(i,j,k) = Phase;
            }
        }
        else
        {
            /* If only one phase is present, the phase composition will be set to
            the total composition. All other phase compositions were set to zero
            before.*/

            size_t pIndex
                = Phi.FieldsProperties[Phi.Fields(i,j,k).front().index].Phase;

            for(size_t comp = 0; comp < Cx.Ncomp; comp++)
            {
                Cx.MoleFractions(i,j,k,{pIndex,comp}) = Cx.MoleFractionsTotal(i,j,k,{comp});
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

}// namespace openphase
