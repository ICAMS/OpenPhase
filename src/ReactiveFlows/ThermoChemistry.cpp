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
 *  Main contributors :   Reza Namdar
 *
 */
 
#ifdef CANTERA

#include "ReactiveFlows/ThermoChemistry.h" 

using namespace std;
using namespace openphase;

void ThermoChemistry::ReadInput(string InputFile)
{
	ConsoleOutput::WriteBlankLine();
	ConsoleOutput::WriteLineInsert("ThermoChemistry");
	ConsoleOutput::WriteStandard("Source", InputFile);
	std::fstream inp(InputFile, std::ios::in);
	if (!inp)
	{
	    std::stringstream message;
	    message << "File \"" << InputFile << "\" could not be opened";
	    throw std::runtime_error(message.str());
	};
	std::stringstream inp_data;
	inp_data << inp.rdbuf();
	inp.close();
	int moduleLocation   = FileInterface::FindModuleLocation(inp_data, "ThermoChemistry");
	ReactionMechanism    = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("ReactionMechanism"),false,"bferch4.yaml");
	GasPhaseName   		 = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("GasPhaseName"),false,"CH4_BFER_mix");
	GasTransportName   	 = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("GasTransportName"),false,"mixture-averaged");
	InitialMixture	     = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("InitialMixture"),false,"O2:0.232, N2:0.768");
	STOICH   			 = FileInterface::ReadParameterB(inp_data, moduleLocation, std::string("STOICH"),false, false);
	EquiValRatio         = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("EquiValRatio"), false, 1.0);
	FuelName   		     = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("FuelName"),false,"CH4");
	Oxidizer		     = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("Oxidizer"),false,"O2:0.232, N2:0.768");
	nPhases        		 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("nPhases"), false, 1.0);
}

Tensor<double, 1> ThermoChemistry::GettingGasFuelMixture(Cantera::ThermoPhase& gas, const std::string FuelName, 
												const std::string Oxidizer, double ER)
{
    gas.setEquivalenceRatio(ER, FuelName, Oxidizer, Cantera::ThermoBasis::mass);
    size_t nSp = gas.nSpecies();
    Tensor<double, 1> FM({nSp});
    gas.getMassFractions(FM.data());
    return FM;
}

Tensor<double, 1> ThermoChemistry::GettingGasMixture(Cantera::ThermoPhase& gas, const double T0, const double P0)
{
	gas.setState_TPY(T0, P0, InitialMixture);
    size_t nSp = gas.nSpecies();
    Tensor<double, 1> FM({nSp});
    gas.getMassFractions(FM.data());
    return FM;
}

#ifdef ONNX
void ThermoChemistry::CalculateGasKinetics(Energy &EN, Species& SP, MixtureFlow &MF,
                                       FlowSolverLBM& FL, NeuralNetwork& NN)
{
    double R = 8.314462618;
	SP.CalculateCpSpAndhSp(EN, FL);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.T, EN.T.Bcells(),)
    {
        if(!FL.Obstacle(i,j,k))
        {
            double MMw = SP.CalculateMeanMolarMass(i, j, k, "new");

            double T_d = EN.T(i,j,k);
			MF.DensityMix(i,j,k) = EN.CalculateIdealGasDensity(FL.Pth, R/MMw, T_d);

            // Build NN inputs
            float T = static_cast<float>(T_d);

            std::vector<float> x_prop;
            x_prop.reserve(SP.nSpecies + 1);
            for (size_t ig = 0; ig < SP.nSpecies; ig++)
                x_prop.push_back(static_cast<float>(SP.MassFrac(i,j,k,{ig})));
            x_prop.push_back(T);

            std::vector<float> x_react;
            x_react.reserve(SP.nSpecies + 1);
			for (size_t ig = 0; ig < SP.nSpecies-1; ig++)
                x_react.push_back(static_cast<float>(SP.MassFrac(i,j,k,{ig})));
            x_react.push_back(T);
            x_react.push_back(1.0f / T);

            auto yprop  = NN.PredictProp(x_prop);
            auto yreact = NN.PredictReact(x_react); 

            // Use outputs (adjust indices to your model definition!)
			for (size_t ig = 0; ig < SP.nSpecies; ig++)
                SP.DSp(i,j,k,{ig}) = yprop[ig];

            EN.KMix(i,j,k) = yprop[6];
            MF.ViscosityMix(i,j,k)= yprop[7] / MF.DensityMix(i,j,k);

			for (size_t ig = 0; ig < SP.nSpecies-1; ig++)
                SP.WSp(i,j,k,{ig}) = yreact[ig];
			SP.WSp(i,j,k,{SP.nSpecies-1})=0.0;

            EN.HRR(i,j,k) = 0.0;
            double CpMix = 0.0;
			for (size_t ig = 0; ig < SP.nSpecies; ig++)
            {
                SP.WSp(i,j,k,{ig}) *= (1000.0 * SP.MwSp({ig}));
                EN.HRR(i,j,k) -= (SP.hSp(i,j,k,{ig}) * SP.WSp(i,j,k,{ig}));
                CpMix += SP.MassFrac(i,j,k,{ig}) * SP.CpSp(i,j,k,{ig});
            }
            EN.CpMix(i,j,k) = CpMix;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
#endif

void ThermoChemistry::CalculateGasKinetics(Energy &EN, Species& SP, MixtureFlow &MF, FlowSolverLBM& FL,
						Cantera::ThermoPhase &gas, Cantera::Transport &transp, Cantera::Kinetics &kin, BoundaryConditions& BC) 
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.T, 0,)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			gas.setState_TPY(EN.T(i,j,k), FL.Pth, SP.MassFrac(i,j,k).data());
			MF.DensityMix(i,j,k)=gas.density();
			MF.ViscosityMix(i,j,k)=transp.viscosity()/MF.DensityMix(i,j,k);
			EN.CpMix(i,j,k)=gas.cp_mass();
			EN.KMix(i,j,k)=transp.thermalConductivity();
			gas.getPartialMolarCp(SP.CpSp(i,j,k).data());
			transp.getMixDiffCoeffsMass(SP.DSp(i,j,k).data());
			kin.getNetProductionRates(SP.WSp(i,j,k).data()); //kmol/m^3/s
			gas.getPartialMolarEnthalpies(SP.hSp(i,j,k).data());
	        EN.HRR(i,j,k) = 0.0;
			for (size_t ig = 0; ig < SP.nSpecies; ig++)
			{
				SP.WSp(i,j,k,{ig})  *= SP.MwSp({ig})*1000.0;
            	SP.CpSp(i,j,k,{ig}) /= SP.MwSp({ig})*1000.0;
				SP.hSp(i,j,k,{ig})  /= SP.MwSp({ig})*1000.0;
            	EN.HRR(i,j,k)       -=  (SP.hSp(i,j,k,{ig})*SP.WSp(i,j,k)({ig}));
			}
    	}
	}
	OMP_PARALLEL_STORAGE_LOOP_END
		
	if (EN.Grid.dNx) BC.SetX(MF.DensityMix);
    if (EN.Grid.dNy) BC.SetY(MF.DensityMix);
    if (EN.Grid.dNz) BC.SetZ(MF.DensityMix);

	if (EN.Grid.dNx) BC.SetX(MF.ViscosityMix);
    if (EN.Grid.dNy) BC.SetY(MF.ViscosityMix);
    if (EN.Grid.dNz) BC.SetZ(MF.ViscosityMix);

	if (EN.Grid.dNx) BC.SetX(EN.CpMix);
    if (EN.Grid.dNy) BC.SetY(EN.CpMix);
    if (EN.Grid.dNz) BC.SetZ(EN.CpMix);

	if (EN.Grid.dNx) BC.SetX(EN.KMix);
    if (EN.Grid.dNy) BC.SetY(EN.KMix);
    if (EN.Grid.dNz) BC.SetZ(EN.KMix);

	if (EN.Grid.dNx) BC.SetX(SP.DSp);
    if (EN.Grid.dNy) BC.SetY(SP.DSp);
    if (EN.Grid.dNz) BC.SetZ(SP.DSp);

	if (EN.Grid.dNx) BC.SetX(SP.CpSp);
    if (EN.Grid.dNy) BC.SetY(SP.CpSp);
    if (EN.Grid.dNz) BC.SetZ(SP.CpSp);

	if (EN.Grid.dNx) BC.SetX(SP.hSp);
    if (EN.Grid.dNy) BC.SetY(SP.hSp);
    if (EN.Grid.dNz) BC.SetZ(SP.hSp);
}

void ThermoChemistry::CalculateGasKinetics(PhaseField& Phase, Energy &EN, Species& SP, MixtureFlow &MF, SolidBody& SB, FlowSolverLBM& FL,
						Cantera::ThermoPhase &gas, Cantera::Transport &transp, Cantera::Kinetics &kin, BoundaryConditions& BC) 
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EN.T, 0,)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			gas.setState_TPY(EN.T(i,j,k), FL.Pth, SP.MassFrac(i,j,k).data());
			MF.DensityMix(i,j,k)=gas.density();
			MF.ViscosityMix(i,j,k)=transp.viscosity()/MF.DensityMix(i,j,k);
			EN.CpMix(i,j,k)=gas.cp_mass();
			EN.KMix(i,j,k)=transp.thermalConductivity();
			gas.getPartialMolarCp(SP.CpSp(i,j,k).data());
			transp.getMixDiffCoeffsMass(SP.DSp(i,j,k).data());
			kin.getNetProductionRates(SP.WSp(i,j,k).data()); //kmol/m^3/s
			gas.getPartialMolarEnthalpies(SP.hSp(i,j,k).data());
	        EN.HRR(i,j,k) = 0.0;
			for (size_t ig = 0; ig < SP.nSpecies; ig++)
			{
				SP.WSp(i,j,k,{ig})  *= SP.MwSp({ig})*1000.0;
            	SP.CpSp(i,j,k,{ig}) /= SP.MwSp({ig})*1000.0;
				SP.hSp(i,j,k,{ig})  /= SP.MwSp({ig})*1000.0;
            	EN.HRR(i,j,k)             -=  (SP.hSp(i,j,k,{ig})*SP.WSp(i,j,k)({ig}));
			}
    	}
		else if(Phase.Fields(i,j,k).interface())
		{
            double Tint = EN.InterfaceTemperature(Phase,FL,SB,i,j,k);
            Tensor<double, 1> Yint = SP.InterfaceSpecies(FL,SB,i,j,k);
			gas.setState_TPY(Tint, FL.Pth, Yint.data());
			MF.DensityMix(i,j,k)=gas.density();
			MF.ViscosityMix(i,j,k)=transp.viscosity()/MF.DensityMix(i,j,k);
			EN.CpMix(i,j,k)=gas.cp_mass();
			EN.KMix(i,j,k)=transp.thermalConductivity();
			gas.getPartialMolarCp(SP.CpSp(i,j,k).data());
			transp.getMixDiffCoeffsMass(SP.DSp(i,j,k).data());
			gas.getPartialMolarEnthalpies(SP.hSp(i,j,k).data());
			for (size_t ig = 0; ig < SP.nSpecies; ig++)
			{
            	SP.CpSp(i,j,k,{ig}) /= SP.MwSp({ig})*1000.0;
				SP.hSp(i,j,k,{ig})  /= SP.MwSp({ig})*1000.0;
			}
		}
	}
	OMP_PARALLEL_STORAGE_LOOP_END

	if (Phase.Grid.dNx) BC.SetX(MF.DensityMix);
    if (Phase.Grid.dNy) BC.SetY(MF.DensityMix);
    if (Phase.Grid.dNz) BC.SetZ(MF.DensityMix);

	if (Phase.Grid.dNx) BC.SetX(MF.ViscosityMix);
    if (Phase.Grid.dNy) BC.SetY(MF.ViscosityMix);
    if (Phase.Grid.dNz) BC.SetZ(MF.ViscosityMix);

	if (Phase.Grid.dNx) BC.SetX(EN.CpMix);
    if (Phase.Grid.dNy) BC.SetY(EN.CpMix);
    if (Phase.Grid.dNz) BC.SetZ(EN.CpMix);

	if (Phase.Grid.dNx) BC.SetX(EN.KMix);
    if (Phase.Grid.dNy) BC.SetY(EN.KMix);
    if (Phase.Grid.dNz) BC.SetZ(EN.KMix);

	if (Phase.Grid.dNx) BC.SetX(SP.DSp);
    if (Phase.Grid.dNy) BC.SetY(SP.DSp);
    if (Phase.Grid.dNz) BC.SetZ(SP.DSp);

	if (Phase.Grid.dNx) BC.SetX(SP.CpSp);
    if (Phase.Grid.dNy) BC.SetY(SP.CpSp);
    if (Phase.Grid.dNz) BC.SetZ(SP.CpSp);

	if (Phase.Grid.dNx) BC.SetX(SP.hSp);
    if (Phase.Grid.dNy) BC.SetY(SP.hSp);
    if (Phase.Grid.dNz) BC.SetZ(SP.hSp);
}

vector<vector<double>> ThermoChemistry::GettingGasSpeciesPolyConstants(Cantera::ThermoPhase &gas)
{
	vector<vector<double>> speciesdata;
	for(size_t k=0; k < gas.nSpecies(); k++)
	{
		std::shared_ptr<Cantera::Species> sp = gas.species(k);
		std::shared_ptr<Cantera::SpeciesThermoInterpType> sptype = sp->thermo;
		size_t index;
		int type;
		double Tmin, Tmax, P_ref;
		vector<double> cantCoeffs(sptype->nCoeffs());
		sptype->reportParameters(index,type,Tmin,Tmax,P_ref,cantCoeffs.data());

		vector<double> locdata;
		nCoeffs = cantCoeffs.size();
		for(size_t ii=0;ii<sptype->nCoeffs();ii++)
		{
			locdata.push_back(cantCoeffs[ii]);
		}
		speciesdata.push_back(locdata);
	}
	std::cout<<"Number of poly nomial constants used to determine specific heat capacity and heat of formation: "<<nCoeffs<<endl;
	return speciesdata;
}

Tensor<string, 1> ThermoChemistry::GettingGasSpeciesNames(Cantera::ThermoPhase& gas)
{
    size_t nSp = gas.nSpecies();
    Tensor<string, 1> GasMixtureNames({nSp});
	for (size_t iSp=0;iSp<nSp;iSp++)
	{
		GasMixtureNames({iSp}) = gas.speciesName(iSp);
	}
	return GasMixtureNames;
}

Tensor<double, 1> ThermoChemistry::GettingAirMixture()
{
	Tensor<double, 1> AM({GasSpeciesNames.size()});
	for (size_t iSp=0;iSp<GasSpeciesNames.size();iSp++)
	{
		string Sp=GasSpeciesNames({iSp});
		if(Sp=="N2")
		{
			AM({iSp})=0.765;
		}
		else if(Sp=="O2")
		{
			AM({iSp})=0.235;
		}
		else
		{
			AM({iSp})=0.0;
		}
	}
	return AM;
}
Tensor<double, 1> ThermoChemistry::GettingBackFlowMixture()
{
	Tensor<double, 1> BFM({GasSpeciesNames.size()});
	for (size_t iSp=0;iSp<GasSpeciesNames.size();iSp++)
	{
		string Sp=GasSpeciesNames({iSp});
		if(Sp=="N2")
		{
			BFM({iSp})=0.765;
		}
		else if(Sp=="O2")
		{
			BFM({iSp})=0.235;
		}
		else
		{
			BFM({iSp})=0.0;
		}
	}
	return BFM;
}
Tensor<double, 1> ThermoChemistry::GettingBurntGasMixture(Cantera::ThermoPhase& gas, const double T0, const double P0)
{
	gas.setState_TP(T0, P0);
	gas.equilibrate("HP");
	Tensor<double, 1> BM({gas.nSpecies()});
	gas.getMassFractions(BM.data());
	BurntGasTemp = gas.temperature();
	return BM;
}
size_t ThermoChemistry::GettingGasFuelIndex(Cantera::ThermoPhase& gas, string FuelName)
{
	size_t FI=0;
	size_t nSp=gas.nSpecies();
	for (size_t ic=0; ic<nSp; ic++)
    {
		if(gas.speciesName(ic)==FuelName)
		{
			FI = ic;
		}
    }
	return FI;
}
Tensor<double, 1> ThermoChemistry::GettingGasMolecularweight(Cantera::ThermoPhase& gas)
{
    Tensor<double, 1> MW({gas.nSpecies()});
	gas.getMolecularWeights(MW.data());
	return MW;
}

void ThermoChemistry::WriteSpeciesNamesInSettingsInputFile(string SettingsInputFile) 
{
	double nb=1.0;
	#ifdef MPI_PARALLEL
    	if(MPI_RANK==0)
		{
			fstream WD;
    		std::string tempfilename=SettingsInputFile;
    		WD.open(tempfilename, ios_base::app);
    		string Element="$Comp_";

			WD<<"	"<<endl;
			for (size_t iSp=0;iSp<GasSpeciesNames.size();iSp++)
			{
				WD<<Element + to_string(iSp) + "	 Chemical Component " + to_string(iSp) + "	 : " + GasSpeciesNames({iSp})<<endl;
			}
		}
		OP_MPI_Allreduce(OP_MPI_IN_PLACE, &nb, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#else
    	fstream WD;
    	std::string tempfilename=SettingsInputFile;
    	WD.open(tempfilename, ios_base::app);
    	string Element="$Comp_";

		WD<<"	"<<endl;
		for (size_t iSp=0; iSp<GasSpeciesNames.size(); iSp++)
		{
			WD<<Element + to_string(iSp) + "	 Chemical Component " + to_string(iSp) + "	 : " + GasSpeciesNames({iSp})<<endl;
		}
	#endif
}

void ThermoChemistry::ExportingGasSpeciesData(Species& SP)
{
	for (size_t iSp=0;iSp<SP.nSpecies;iSp++)
	{
		for (size_t ic = 0; ic < SP.nCoeffs; ic++)
		{
			SP.PolyNomCoeffs({iSp,ic}) = Coeffs.at(iSp).at(ic);
		}
		SP.MwSp({iSp}) = MW({iSp})/1000.0; 
	}
}

void ThermoChemistry::InitializeGasPhase(Cantera::ThermoPhase& gas, const double T0, const double P0, string SettingsInputFile)
{
    GasSpeciesNames.Allocate({gas.nSpecies()});
    GasMixture.Allocate({gas.nSpecies()});
    AirMixture.Allocate({gas.nSpecies()});
    BurntGasMixture.Allocate({gas.nSpecies()});
    BackFlowMixture.Allocate({gas.nSpecies()});
    MW.Allocate({gas.nSpecies()});

	if(!STOICH) GasMixture = GettingGasMixture(gas, T0, P0);
	if(STOICH)  GasMixture = GettingGasFuelMixture(gas, FuelName, Oxidizer, EquiValRatio);
	GasSpeciesNames  		  = GettingGasSpeciesNames(gas);
	Coeffs		  	 		  = GettingGasSpeciesPolyConstants(gas);
	AirMixture	  	 		  = GettingAirMixture();
	BackFlowMixture	 		  = GettingBackFlowMixture();
	BurntGasMixture  		  = GettingBurntGasMixture(gas, T0, P0);
	MW			  	 		  = GettingGasMolecularweight(gas);
	GasFuelIndex 	 		  = GettingGasFuelIndex(gas, FuelName);
	WriteSpeciesNamesInSettingsInputFile(SettingsInputFile);
}

#endif
