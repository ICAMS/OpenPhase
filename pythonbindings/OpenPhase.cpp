#include <pybind11/pybind11.h>
#include "Settings.h"
#include "RunTimeControl.h"
#include "BoundaryConditions.h"
#include "InterfaceProperties.h"
#include "DoubleObstacle.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Initializations.h"
#include "ElasticProperties.h"
#include "ElasticitySolverSpectral.h"
#include "UserDrivingForce.h"
#include "Temperature.h"
#include "Nucleation.h"
#include "TextOutput.h"
#include "HeatDiffusion.h"
#include "HeatSources.h"
#include "Globals.h"
#include "Composition.h"
#include "Crystallography.h"
#include "InterfaceRegularization.h"
#include "AdvectionHR.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "EquilibriumPartitionDiffusionBinary.h"
#include "Tools/TimeInfo.h"

namespace py = pybind11;

PYBIND11_MODULE(OpenPhase, m) {
    py::class_<openphase::Settings>(m, "Settings")
        .def(py::init<const std::string &>())
        .def("ReadInput",  static_cast<void (openphase::Settings::*)(std::string)>(&openphase::Settings::ReadInput));
        
    py::class_<openphase::GridParameters>(m, "GridParameters")
        .def(py::init<const std::string &>())
        .def("ReadInput",  static_cast<void (openphase::GridParameters::*)(std::string)>(&openphase::GridParameters::ReadInput))
        .def_readwrite("Nx", &openphase::GridParameters::Nx)
        .def_readwrite("Ny", &openphase::GridParameters::Ny)
        .def_readwrite("Nz", &openphase::GridParameters::Nz);
        
    py::enum_<openphase::AggregateStates>(m, "AggregateStates")
        .value("Liquid", openphase::AggregateStates::Liquid)
        .value("Solid", openphase::AggregateStates::Solid)
        .value("Gas", openphase::AggregateStates::Gas)
        .export_values();
        
    py::class_<openphase::RunTimeControl>(m, "RunTimeControl")
        .def(py::init<>()) 
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::RunTimeControl::Initialize)  
        .def("ReadInput",  static_cast<void (openphase::RunTimeControl::*)(std::string)>(&openphase::RunTimeControl::ReadInput))
        .def("WriteVTK", &openphase::RunTimeControl::WriteVTK)
        .def("WriteRawData", &openphase::RunTimeControl::WriteRawData)
        .def("WriteToScreen", &openphase::RunTimeControl::WriteToScreen)
        .def("IncrementTimeStep", &openphase::RunTimeControl::IncrementTimeStep)
        .def_readwrite("dt", &openphase::RunTimeControl::dt)
        .def_readwrite("MaxTimeStep", &openphase::RunTimeControl::MaxTimeStep)
        .def_readwrite("TimeStep", &openphase::RunTimeControl::TimeStep)
        .def_readwrite("StartTimeStep", &openphase::RunTimeControl::StartTimeStep)
        .def_readwrite("RestartSwitch", &openphase::RunTimeControl::RestartSwitch);
        
    py::class_<openphase::BoundaryConditions>(m, "BoundaryConditions")
        .def(py::init<>())
        .def(py::init<openphase::Settings&, const std::string &>())      
        .def("Initialize", &openphase::BoundaryConditions::Initialize)
        .def("ReadInput",  static_cast<void (openphase::BoundaryConditions::*)(std::string)>(&openphase::BoundaryConditions::ReadInput));      
        
    py::class_<openphase::PhaseField>(m, "PhaseField")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::PhaseField::Initialize)
        .def("ReadInput",  static_cast<void (openphase::PhaseField::*)(std::string)>(&openphase::PhaseField::ReadInput))
        .def("MergeIncrements",
         &openphase::PhaseField::MergeIncrements,
         py::arg("BC"),
         py::arg("dt"),
         py::arg("finalize") = true,
         py::arg("clear") = true)
        .def("NormalizeIncrements",&openphase::PhaseField::NormalizeIncrements) 
        .def("PlantGrainNucleus", &openphase::PhaseField::PlantGrainNucleus)
        .def("WriteVTK",
         &openphase::PhaseField::WriteVTK,
         py::arg("locSettings"),
         py::arg("tStep"),
         py::arg("CurvatureOutput") = false,
         py::arg("precision") = 16)
    	.def("Read", static_cast<bool (openphase::PhaseField::*)(std::string)>(&openphase::PhaseField::Read))
		.def("Read", static_cast<bool (openphase::PhaseField::*)(const openphase::Settings&, const openphase::BoundaryConditions&, int)>(&openphase::PhaseField::Read))
        .def("Write", static_cast<bool (openphase::PhaseField::*)(const std::string&) const>(&openphase::PhaseField::Write))
        .def("Write", static_cast<bool (openphase::PhaseField::*)(const openphase::Settings&, const int) const>(&openphase::PhaseField::Write))
        .def("PrintVolumeFractions", &openphase::PhaseField::PrintVolumeFractions)
        .def("WriteDistortedVTK",&openphase::PhaseField::WriteDistortedVTK)
        .def("Advect", [](openphase::PhaseField& self,
				          openphase::AdvectionHR& Adv,
				          const openphase::Velocities& Vel,
				          const openphase::BoundaryConditions& BC,
				          double dt, double tStep,
				          bool finalize) {
			self.Advect(Adv, Vel, BC, dt, tStep, finalize);
		}, py::arg("Adv"), py::arg("Vel"), py::arg("BC"), 
		   py::arg("dt"), py::arg("tStep"), py::arg("finalize") = true)

		.def("AdvectWithLBM", [](openphase::PhaseField& self,
				                 openphase::AdvectionHR& Adv,
				                 const openphase::Velocities& Vel,
				                 const openphase::BoundaryConditions& BC,
				                 openphase::FlowSolverLBM& LBM,
				                 double dt, double tStep,
				                 bool finalize) {
			self.Advect(Adv, Vel, BC, LBM, dt, tStep, finalize);
		}, py::arg("Adv"), py::arg("Vel"), py::arg("BC"), 
		   py::arg("LBM"), py::arg("dt"), py::arg("tStep"), py::arg("finalize") = true)
		.def_readwrite("FieldsProperties", &openphase::PhaseField::FieldsProperties);

        
    py::class_<openphase::InterfaceProperties>(m, "InterfaceProperties")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::InterfaceProperties::Initialize)  
        .def("ReadInput",  static_cast<void (openphase::InterfaceProperties::*)(std::string)>(&openphase::InterfaceProperties::ReadInput))  
	    .def("Set", static_cast<void (openphase::InterfaceProperties::*)(const openphase::PhaseField&, const openphase::BoundaryConditions&)>(&openphase::InterfaceProperties::Set))
	    .def("Set", static_cast<void (openphase::InterfaceProperties::*)(const openphase::PhaseField&, const openphase::Temperature&, const openphase::BoundaryConditions&)>(&openphase::InterfaceProperties::Set))
	    .def("ReportMaximumTimeStep", &openphase::InterfaceProperties::ReportMaximumTimeStep);;	
	    
    py::class_<openphase::InterfaceRegularization>(m, "InterfaceRegularization")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::InterfaceRegularization::Initialize)  
        .def("ReadInput",  static_cast<void (openphase::InterfaceRegularization::*)(std::string)>(&openphase::InterfaceRegularization::ReadInput))
        .def("Average", &openphase::InterfaceRegularization::Average)
        .def("MergePhaseFieldIncrements", &openphase::InterfaceRegularization::MergePhaseFieldIncrements);
        
    py::class_<openphase::HeatSources>(m, "HeatSources")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::HeatSources::Initialize)  
        .def("ReadInput",  static_cast<void (openphase::HeatSources::*)(std::string)>(&openphase::HeatSources::ReadInput))
        .def("Apply", &openphase::HeatSources::Apply);
        
	py::class_<openphase::HeatDiffusion>(m, "HeatDiffusion")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::HeatDiffusion::Initialize)  
        .def("ReadInput",  static_cast<void (openphase::HeatDiffusion::*)(std::string)>(&openphase::HeatDiffusion::ReadInput))  
        .def("SetEffectiveProperties", static_cast<void (openphase::HeatDiffusion::*)(
			const openphase::PhaseField&,
			const openphase::Temperature&)>(
				&openphase::HeatDiffusion::SetEffectiveProperties)
		)
        .def("SolveImplicit", &openphase::HeatDiffusion::SolveImplicit);	       
			
	py::class_<openphase::DrivingForce>(m, "DrivingForce")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::DrivingForce::Initialize)
        .def("ReadInput",  static_cast<void (openphase::DrivingForce::*)(std::string)>(&openphase::DrivingForce::ReadInput))
        .def("Average", &openphase::DrivingForce::Average)
        .def("Clear", &openphase::DrivingForce::Clear)
        .def("MergePhaseFieldIncrements", &openphase::DrivingForce::MergePhaseFieldIncrements)
        .def("PrintDiagnostics", &openphase::DrivingForce::PrintDiagnostics)
        .def("WriteVTK",static_cast<void (openphase::DrivingForce::*)(const openphase::Settings&, const openphase::PhaseField&, const int, const int) const>(&openphase::DrivingForce::WriteVTK));      
        
    py::class_<openphase::DoubleObstacle>(m, "DoubleObstacle")
        .def(py::init<>())   
        .def(py::init<openphase::Settings&, const std::string &>())
        .def("Initialize", &openphase::DoubleObstacle::Initialize)  
	    .def("CalculatePhaseFieldIncrements",
         py::overload_cast<openphase::PhaseField&, openphase::InterfaceProperties&>(
             &openphase::DoubleObstacle::CalculatePhaseFieldIncrements))
        .def("CalculatePhaseFieldIncrements",
             py::overload_cast<openphase::PhaseField&, openphase::InterfaceProperties&, openphase::DrivingForce&>(
                 &openphase::DoubleObstacle::CalculatePhaseFieldIncrements))
        .def("CalculatePhaseFieldIncrements",
             py::overload_cast<openphase::PhaseField&, openphase::InterfaceProperties&, openphase::InterfaceRegularization&>(
                 &openphase::DoubleObstacle::CalculatePhaseFieldIncrements))
        .def("CalculatePhaseFieldIncrements",
             py::overload_cast<openphase::PhaseField&, openphase::InterfaceProperties&, openphase::DrivingForce&, openphase::InterfaceRegularization&>(
                 &openphase::DoubleObstacle::CalculatePhaseFieldIncrements))
	    .def("AverageEnergyDensity", &openphase::DoubleObstacle::AverageEnergyDensity);
	
    py::class_<openphase::Composition>(m, "Composition")
        .def(py::init<>())      
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::Composition::Initialize)
        .def("ReadInput",  static_cast<void (openphase::Composition::*)(std::string)>(&openphase::Composition::ReadInput))
        .def("SetBoundaryConditions",&openphase::Composition::SetBoundaryConditions)
        .def("CalculateMoleFractionsAverage",&openphase::Composition::CalculateMoleFractionsAverage)
        .def("SetInitialMoleFractions", &openphase::Composition::SetInitialMoleFractions)
        .def("Read", &openphase::Composition::Read)
        .def("Write", &openphase::Composition::Write)
        .def("WriteStatistics", &openphase::Composition::WriteStatistics)
        .def("WriteVTK",&openphase::Composition::WriteVTK)
		.def("Advect",&openphase::Composition::Advect);     
           
     py::class_<openphase::Temperature>(m, "Temperature")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::Temperature::Initialize)
        .def("ReadInput",  static_cast<void (openphase::Temperature::*)(std::string)>(&openphase::Temperature::ReadInput))
        .def("SetInitial", &openphase::Temperature::SetInitial)
        .def("WriteVTK",&openphase::Temperature::WriteVTK)
        .def("Write",&openphase::Temperature::Write)
        .def("Read",&openphase::Temperature::Read)
        .def("Advect",&openphase::Temperature::Advect)
        .def("PrintStatistics",&openphase::Temperature::PrintStatistics);
              
    py::class_<openphase::ElasticProperties>(m, "ElasticProperties")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::ElasticProperties::Initialize)
        .def("ReadInput",  static_cast<void (openphase::ElasticProperties::*)(std::string)>(&openphase::ElasticProperties::ReadInput))
        .def("SetEffectiveProperties", py::overload_cast<const openphase::PhaseField&>(&openphase::ElasticProperties::SetEffectiveProperties))
        .def("SetEffectiveProperties", py::overload_cast<const openphase::PhaseField&, const openphase::Composition&>(&openphase::ElasticProperties::SetEffectiveProperties))
        .def("SetEffectiveProperties", py::overload_cast<const openphase::PhaseField&, const openphase::Temperature&>(&openphase::ElasticProperties::SetEffectiveProperties))
        .def("SetEffectiveProperties", py::overload_cast<const openphase::PhaseField&, const openphase::Composition&, const openphase::Temperature&>(&openphase::ElasticProperties::SetEffectiveProperties))
        .def("SetEffectiveProperties", py::overload_cast<const openphase::PhaseField&, const openphase::InterfaceProperties&>(&openphase::ElasticProperties::SetEffectiveProperties))
        .def("CalculateDrivingForce",static_cast<void (openphase::ElasticProperties::*)(const openphase::PhaseField&, openphase::DrivingForce&)>(&openphase::ElasticProperties::CalculateDrivingForce))
        .def("WritePlasticStrainsVTK", &openphase::ElasticProperties::WritePlasticStrainsVTK)
        .def("WriteTotalRotationsVTK", &openphase::ElasticProperties::WriteTotalRotationsVTK);     
        
    py::class_<openphase::ElasticitySolverSpectral>(m, "ElasticitySolverSpectral")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", static_cast<void (openphase::ElasticitySolverSpectral::*)(openphase::Settings&, std::string)>(&openphase::ElasticitySolverSpectral::Initialize))
		.def("Initialize", static_cast<void (openphase::ElasticitySolverSpectral::*)(openphase::Settings&, const openphase::BoundaryConditions&, std::string)>(&openphase::ElasticitySolverSpectral::Initialize))
        .def("ReadInput",  static_cast<void (openphase::ElasticitySolverSpectral::*)(std::string)>(&openphase::ElasticitySolverSpectral::ReadInput))
        .def("Solve", static_cast<int (openphase::ElasticitySolverSpectral::*)(openphase::ElasticProperties&, openphase::BoundaryConditions&, double, std::function<bool()>)>(&openphase::ElasticitySolverSpectral::Solve));    
                
    py::class_<openphase::Nucleation>(m, "Nucleation")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::Nucleation::Initialize)
        .def("ReadInput",  static_cast<void (openphase::Nucleation::*)(std::string)>(&openphase::Nucleation::ReadInput))
        .def("GenerateNucleationSites", &openphase::Nucleation::GenerateNucleationSites)
        .def("PlantNuclei", &openphase::Nucleation::PlantNuclei)
        .def("CheckNuclei", &openphase::Nucleation::CheckNuclei)
        .def("Write", &openphase::Nucleation::Write)
        .def("Read", &openphase::Nucleation::Read); 
        
    py::class_<openphase::UserDrivingForce>(m, "UserDrivingForce")
        .def(py::init<>())    
        .def(py::init<openphase::Settings&, const std::string &>()) 
		.def("Initialize", &openphase::UserDrivingForce::Initialize)
        .def("ReadInput",  static_cast<void (openphase::UserDrivingForce::*)(std::string)>(&openphase::UserDrivingForce::ReadInput))
        .def("SetDrivingForce", static_cast<void (openphase::UserDrivingForce::*)(openphase::PhaseField&, openphase::DrivingForce&, openphase::Temperature&)>(&openphase::UserDrivingForce::SetDrivingForce));     
                  
    py::class_<openphase::Crystallography>(m, "Crystallography")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::Crystallography::Initialize)
        .def("ReadInput",  static_cast<void (openphase::Crystallography::*)(std::string)>(&openphase::Crystallography::ReadInput));  
        
    py::class_<openphase::AdvectionHR>(m, "AdvectionHR")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::AdvectionHR::Initialize)
        .def("ReadInput",  static_cast<void (openphase::AdvectionHR::*)(std::string)>(&openphase::AdvectionHR::ReadInput));  
        
        py::class_<openphase::FlowSolverLBM>(m, "FlowSolverLBM")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, double, const std::string &>())
		.def("Initialize", &openphase::FlowSolverLBM::Initialize)
        .def("ReadInput",  static_cast<void (openphase::FlowSolverLBM::*)(std::string)>(&openphase::FlowSolverLBM::ReadInput))
        .def("Solve", static_cast<void (openphase::FlowSolverLBM::*)(
			openphase::PhaseField&, openphase::Velocities&, const openphase::BoundaryConditions&)>(
			&openphase::FlowSolverLBM::Solve))
		.def("Solve", static_cast<void (openphase::FlowSolverLBM::*)(
			openphase::PhaseField&, const openphase::Composition&, openphase::Velocities&, const openphase::BoundaryConditions&)>(
			&openphase::FlowSolverLBM::Solve))
		.def("SetUniformVelocity", &openphase::FlowSolverLBM::SetUniformVelocity);   
    
    py::class_<openphase::Initializations>(m, "Initializations")
        .def_static("VoronoiTessellation",  &openphase::Initializations::VoronoiTessellation)   
        .def_static("Single",  &openphase::Initializations::Single);   
        
    py::class_<openphase::EquilibriumPartitionDiffusionBinary>(m, "EquilibriumPartitionDiffusionBinary")
        .def(py::init<>())     
        .def(py::init<openphase::Settings&, const std::string &>())
		.def("Initialize", &openphase::EquilibriumPartitionDiffusionBinary::Initialize)
        .def("ReadInput",  static_cast<void (openphase::EquilibriumPartitionDiffusionBinary::*)(std::string)>(&openphase::EquilibriumPartitionDiffusionBinary::ReadInput))
        .def("SolveDiffusion", &openphase::EquilibriumPartitionDiffusionBinary::SolveDiffusion)
        .def("CalculateDrivingForce", &openphase::EquilibriumPartitionDiffusionBinary::CalculateDrivingForce);  
        
    py::class_<openphase::TimeInfo>(m, "TimeInfo")
        .def(py::init<>())     
		.def(py::init<const openphase::Settings&, const std::string&, bool>(),
			py::arg("settings"), py::arg("name") = "Execution Time Statistics", py::arg("verbose") = false)
		.def("Initialize", &openphase::TimeInfo::Initialize,
			py::arg("settings"), py::arg("name") = "Execution Time Statistics", py::arg("verbose") = false)
		.def("SetStart", &openphase::TimeInfo::SetStart)
		.def("SetTimeStamp", &openphase::TimeInfo::SetTimeStamp)
		.def("SkipToHere", &openphase::TimeInfo::SkipToHere)
		.def("Reset", &openphase::TimeInfo::Reset)
		.def("PrintWallClockSummary", &openphase::TimeInfo::PrintWallClockSummary);
		   		
	py::class_<openphase::Grain>(m, "Grain")
		.def(py::init<>())
		.def_readwrite("Exist", &openphase::Grain::Exist)
		.def_readwrite("Mobile", &openphase::Grain::Mobile)
		.def_readwrite("Parent", &openphase::Grain::Parent)
		.def_readwrite("Phase", &openphase::Grain::Phase)
		.def_readwrite("Variant", &openphase::Grain::Variant)
		.def_readwrite("Stage", &openphase::Grain::Stage)
		.def_readwrite("State", &openphase::Grain::State)
		.def_readwrite("Density", &openphase::Grain::Density)
		.def_readwrite("RefVolume", &openphase::Grain::RefVolume)
		.def_readwrite("Volume", &openphase::Grain::Volume)
		.def_readwrite("MAXVolume", &openphase::Grain::MAXVolume)
		.def_readwrite("VolumeRatio", &openphase::Grain::VolumeRatio)
		.def_readwrite("Rcm", &openphase::Grain::Rcm)
		.def_readwrite("Vcm", &openphase::Grain::Vcm)
		.def_readwrite("Acm", &openphase::Grain::Acm)
		.def_readwrite("aVel", &openphase::Grain::aVel)
		.def_readwrite("aAcc", &openphase::Grain::aAcc)
		.def_readwrite("Force", &openphase::Grain::Force)
		.def_readwrite("Torque", &openphase::Grain::Torque)
		.def_readwrite("InertiaM", &openphase::Grain::InertiaM)
		.def_readwrite("Orientation", &openphase::Grain::Orientation)
		.def("Clear", &openphase::Grain::Clear);  // if defined

 	py::class_<openphase::GrainsProperties>(m, "GrainsProperties")
		.def(py::init<>())
		.def("size", &openphase::GrainsProperties::size)
		.def("Allocate", &openphase::GrainsProperties::Allocate)
		.def("Reallocate", &openphase::GrainsProperties::Reallocate)
		.def("__getitem__", [](openphase::GrainsProperties &gp, size_t i) -> openphase::Grain& {
		    if (i >= gp.size()) throw py::index_error();
		    return gp[i];
		}, py::return_value_policy::reference_internal)
		.def("__setitem__", [](openphase::GrainsProperties &gp, size_t i, const openphase::Grain& g) {
		    if (i >= gp.size()) throw py::index_error();
		    gp[i] = g;
		});
                            
}


