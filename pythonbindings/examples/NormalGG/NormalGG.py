import sys
import os
from pathlib import Path

lib_dir1 = Path(__file__).resolve().parents[3] / "lib" / "libOpenPhase.so"

link1 = Path("libOpenPhase.so")
if not link1.exists():
    os.symlink(lib_dir1, link1)
    print(f"Created symlink: {link1} → {lib_dir1}")
else:
    print(f"Symlink already exists: {link1}")
    
sys.path.append(str(Path(__file__).resolve().parents[2]))

from OpenPhase import (
    Settings,
    RunTimeControl,
    PhaseField,
    DoubleObstacle,
    InterfaceProperties,
    BoundaryConditions,
    Initializations,
    TimeInfo,
)
# Optional (only if available in bindings)
try:
    from OpenPhase import MicrostructureAnalysis, ConsoleOutput
except ImportError:
    MicrostructureAnalysis = None
    ConsoleOutput = None


def main():
    # --- Handle command line arguments ---
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        print("No Inputfile provided, trying to use default ProjectInput.opi")
        input_file = "ProjectInput.opi"

    # --- Initialize core simulation objects ---
    OPSettings = Settings(input_file)

    RTC = RunTimeControl(OPSettings, input_file)
    Phi = PhaseField(OPSettings, input_file)
    DO = DoubleObstacle(OPSettings, input_file)
    IP = InterfaceProperties(OPSettings, input_file)
    BC = BoundaryConditions(OPSettings, input_file)
    Timer = TimeInfo(OPSettings, "Execution Time Statistics")

    # --- Generate initial microstructure (Voronoi) ---
    number_of_grains = 200
    GrainsPhase = 0
    Initializations.VoronoiTessellation(Phi, BC, number_of_grains, GrainsPhase)

    print("Entering the Time Loop!!!")

    # --- Main time-stepping loop ---
    t = RTC.StartTimeStep
    while t <= RTC.MaxTimeStep:
        Timer.SetStart()

        IP.Set(Phi, BC)
        Timer.SetTimeStamp("IP.Set()")

        DO.CalculatePhaseFieldIncrements(Phi, IP)
        Timer.SetTimeStamp("CalculatePhaseFieldIncrements")

        Phi.NormalizeIncrements(BC, RTC.dt)
        Timer.SetTimeStamp("NormalizeIncrements")

        Phi.MergeIncrements(BC, RTC.dt)
        Timer.SetTimeStamp("MergeIncrements")

        # --- Output to VTK ---
        if RTC.WriteVTK():
            Phi.WriteVTK(OPSettings, RTC.TimeStep)
            if MicrostructureAnalysis:
                MicrostructureAnalysis.WriteGrainsStatistics(Phi, RTC.TimeStep)

        # --- Output raw data ---
        if RTC.WriteRawData():
            Phi.Write(OPSettings, RTC.TimeStep)

        # --- Output to screen ---
        if RTC.WriteToScreen():
            I_En = DO.AverageEnergyDensity(Phi, IP)
            if ConsoleOutput:
                msg = ConsoleOutput.GetStandard("Interface energy density", I_En)
                ConsoleOutput.WriteTimeStep(RTC, msg)
            else:
                print(f"[t={RTC.TimeStep}] Interface energy density: {I_En:.6e}")

            Timer.PrintWallClockSummary()

        # --- Increment timestep ---
        RTC.IncrementTimeStep()
        t = RTC.TimeStep

    print("Simulation finished successfully.")


if __name__ == "__main__":
    main()
