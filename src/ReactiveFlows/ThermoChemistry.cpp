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

 *   File created :   2025
 *   Main contributors :   Oleg Shchyglo; Reza Namdar
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
	ReactionMechanism    = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("ReactionMechanism"),false,"CH4_BFER.yaml");
	PhaseName   		 = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("PhaseName"),false,"CH4_BFER_mix");
	TransportName   	 = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("TransportName"),false,"mixture-averaged");
	FuelName   		     = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("FuelName"),false,"CH4");
	Oxidizer   		     = FileInterface::ReadParameterS(inp_data, moduleLocation, std::string("Oxidizer"),false,"air");
	ER        			 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("EquivRatio"), false, 1.0);
	nPhases        		 = FileInterface::ReadParameterD(inp_data, moduleLocation, std::string("nPhases"), false, 1.0);
}

vector<double> ThermoChemistry::GettingFuelMixture(std::shared_ptr<Cantera::ThermoPhase> & gas, string FuelName, string Oxidizer, double ER)
{
	string oxidizer_cant;
	if(Oxidizer=="air")
	{
		oxidizer_cant = "O2:0.235, N2:0.765";
	}
	gas->setEquivalenceRatio(ER, FuelName, oxidizer_cant, Cantera::ThermoBasis::mass);
    int nSp = gas->nSpecies();
    double FuelMixture[nSp];
	gas->getMassFractions(FuelMixture);
	vector<double> FM;
	for (int iSp=0;iSp<nSp;iSp++)
	{
		FM.push_back(FuelMixture[iSp]);
	}
	return FM;
}

void ThermoChemistry::UpdateProperties(EnergyTransport& ET, SpeciesTransport& ST, FlowSolverLBM& FL,
				shared_ptr<Cantera::ThermoPhase> &gas, shared_ptr<Cantera::Transport> &transp, shared_ptr<Cantera::Kinetics> &kin) 
{
	//double R = 8.314462618; // J/(mol*K)
	//size_t nr = kin->nReactions();
	size_t MixtureComp=0;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ET.Tx, ET.Tx.Bcells(),)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			gas->setState_TPY(ET.Tx(i,j,k), FL.Pth, ST.MassFractions(i,j,k).data());
			FL.DensityWetting(i,j,k)({MixtureComp})=gas->density();
			FL.nut(i,j,k)({MixtureComp})=transp->viscosity()/FL.DensityWetting(i,j,k)({MixtureComp});
			ET.Cp_Mixture(i,j,k)=gas->cp_mass();
			ET.K_Mixture(i,j,k)=transp->thermalConductivity();
			if(ST.Species)
			{
				transp->getMixDiffCoeffsMass(ST.MassDiff_Species(i,j,k).data());
				kin->getNetProductionRates(ST.W_Species(i,j,k).data()); //kmol/m^3/s
				for (size_t ic = 0; ic < ST.nSpecies; ic++)
				{
            		ST.W_Species(i,j,k)({ic}) *= (1000.0*ST.MolecularWeight({ic}));
				}

	        	ST.HRR(i,j,k) = 0.0;
				for (size_t ic = 0; ic < ST.nSpecies; ic++)
				{
            		ST.HRR(i,j,k) -=  (ST.h_Species(i,j,k)({ic})*ST.W_Species(i,j,k)({ic}));
				}			 
			}
    	}
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

void ThermoChemistry::UpdateAllProperties(EnergyTransport& ET, SpeciesTransport& ST, FlowSolverLBM& FL,
				shared_ptr<Cantera::ThermoPhase> &gas, shared_ptr<Cantera::Transport> &transp, shared_ptr<Cantera::Kinetics> &kin) 
{
	//double R = 8.314462618; // J/(mol*K)
	//size_t nr = kin->nReactions();
	size_t MixtureComp=0;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,ET.Tx, ET.Tx.Bcells(),)
    {
		if(!FL.Obstacle(i,j,k))
    	{
			gas->setState_TPY(ET.Tx(i,j,k), FL.Pth, ST.MassFractions(i,j,k).data());
			FL.DensityWetting(i,j,k)({MixtureComp})=gas->density();
			FL.nut(i,j,k)({MixtureComp})=transp->viscosity()/FL.DensityWetting(i,j,k)({MixtureComp});
			ET.Cp_Mixture(i,j,k)=gas->cp_mass();
			ET.K_Mixture(i,j,k)=transp->thermalConductivity();
			if(ST.Species)
			{
				gas->getPartialMolarCp(ST.Cp_Species(i,j,k).data());
				transp->getMixDiffCoeffsMass(ST.MassDiff_Species(i,j,k).data());
				kin->getNetProductionRates(ST.W_Species(i,j,k).data()); //kmol/m^3/s

				std::vector<double> h(ST.nSpecies);
				gas->getPartialMolarEnthalpies(h.data());
				for (size_t ic = 0; ic < ST.nSpecies; ic++)
				{
            		ST.W_Species(i,j,k)({ic}) *= (1000.0*ST.MolecularWeight({ic}));
            		ST.Cp_Species(i,j,k)({ic}) /= (1000.0*ST.MolecularWeight({ic}));
					h[ic] /= (1000.0*ST.MolecularWeight({ic}));
				}

	        	ST.HRR(i,j,k) = 0.0;
				for (size_t ic = 0; ic < ST.nSpecies; ic++)
				{
            		ST.HRR(i,j,k) -=  (h[ic]*ST.W_Species(i,j,k)({ic}));
				}
							 
			}
    	}
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

vector<vector<double>> ThermoChemistry::GettingSpeciesPolyConstants(std::shared_ptr<Cantera::ThermoPhase> & gas)
{
	vector<vector<double>> speciesdata;
	for(size_t k=0; k < gas->nSpecies(); k++)
	{
		std::shared_ptr<Cantera::Species> sp = gas->species(k);
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

vector<string> ThermoChemistry::GettingSpeciesNames(std::shared_ptr<Cantera::ThermoPhase> & gas)
{
    int nSp = gas->nSpecies();
    vector<string> MixtureNames;
	for (int iSp=0;iSp<nSp;iSp++)
	{
		MixtureNames.push_back(gas->speciesName(iSp));
	}
	return MixtureNames;
}

vector<double> ThermoChemistry::GettingAirMixture()
{
	vector<double> AM;
	for (size_t iSp=0;iSp<SpeciesNames.size();iSp++)
	{
		string Sp=SpeciesNames.at(iSp);
		if(Sp=="N2")
		{
			AM.push_back(0.765);
		}
		else if(Sp=="O2")
		{
			AM.push_back(0.235);
		}
		else
		{
			AM.push_back(0.0);
		}
	}
	return AM;
}
vector<double> ThermoChemistry::GettingBackFlowMixture()
{
	vector<double> BFM;
	for (size_t iSp=0;iSp<SpeciesNames.size();iSp++)
	{
		string Sp=SpeciesNames.at(iSp);
		if(Sp=="N2")
		{
			BFM.push_back(0.765);
		}
		else if(Sp=="O2")
		{
			BFM.push_back(0.235);
		}
		else
		{
			BFM.push_back(0.0);
		}
	}
	return BFM;
}
vector<double> ThermoChemistry::GettingBurntMixture(std::shared_ptr<Cantera::ThermoPhase> & gas, double Temp, double Pressure)
{
	gas->setState_TP(Temp,Pressure);
	gas->equilibrate("HP");
	size_t nSp=gas->nSpecies();
    double BurntMixture[nSp];
	gas->getMassFractions(BurntMixture);
	BurntTemp = gas->temperature();
	vector<double> BM;
	for (size_t iSp=0;iSp<nSp;iSp++)
	{
		BM.push_back(BurntMixture[iSp]);
	}
	return BM;
}
size_t ThermoChemistry::GettingFuelIndex(std::shared_ptr<Cantera::ThermoPhase> & gas, string FuelName)
{
	size_t FI=0;
	size_t nSp=gas->nSpecies();
	for (size_t ic=0; ic<nSp; ic++)
    {
		if(gas->speciesName(ic)==FuelName)
		{
			FI = ic;
		}
    }
	return FI;
}
vector<double> ThermoChemistry::GettingMolecularweight(std::shared_ptr<Cantera::ThermoPhase> & gas)
{
	size_t nSp=gas->nSpecies();
    vector<double> MW(nSp);
	gas->getMolecularWeights(MW.data());
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
			for (size_t iSp=0;iSp<SpeciesNames.size();iSp++)
			{
				WD<<Element + to_string(iSp) + "	 Chemical Component " + to_string(iSp) + "	 : " + SpeciesNames.at(iSp)<<endl;
			}
		}
		OP_MPI_Allreduce(OP_MPI_IN_PLACE, &nb, 1, OP_MPI_DOUBLE, OP_MPI_SUM, OP_MPI_COMM_WORLD);
	#else
    	fstream WD;
    	std::string tempfilename=SettingsInputFile;
    	WD.open(tempfilename, ios_base::app);
    	string Element="$Comp_";

		WD<<"	"<<endl;
		for (size_t iSp=0;iSp<SpeciesNames.size();iSp++)
		{
			WD<<Element + to_string(iSp) + "	 Chemical Component " + to_string(iSp) + "	 : " + SpeciesNames.at(iSp)<<endl;
		}
	#endif

}

void ThermoChemistry::ExportingData(SpeciesTransport& ST)
{
	for (size_t iSp=0;iSp<ST.nSpecies;iSp++)
	{
		for (size_t ic = 0; ic < ST.nCoeffs; ic++)
		{
			ST.PolyNomCoeffs({iSp,ic}) = Coeffs.at(iSp).at(ic);
		}
		ST.MolecularWeight({iSp}) = MW.at(iSp)/1000.0; 
	}
}

void ThermoChemistry::Initialize(std::shared_ptr<Cantera::ThermoPhase> & gas, string SettingsInputFile)
{
	FuelMixture	  	=	  GettingFuelMixture(gas,FuelName, Oxidizer, ER);
	SpeciesNames  	=	  GettingSpeciesNames(gas);
	Coeffs		  	=	  GettingSpeciesPolyConstants(gas);
	AirMixture	  	=	  GettingAirMixture();
	BackFlowMixture	=	  GettingBackFlowMixture();
	BurntMixture  	=	  GettingBurntMixture(gas, 300.0,Cantera::OneAtm);
	MW			  	=	  GettingMolecularweight(gas);
	FuelIndex 	  	=  	  GettingFuelIndex(gas, FuelName);
	WriteSpeciesNamesInSettingsInputFile(SettingsInputFile);
}

#endif