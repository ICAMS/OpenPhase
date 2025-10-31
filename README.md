# OpenPhase-Academic v2.0

<div align="center">
  <img src="documentation/UserManual/Figures/openphase-logo.png" alt="OpenPhase Academic Logo" width="300">
</div>

OpenPhase is an open-source software library for phase-field simulations of complex scientific problems. It is written in object-oriented C++ and features a modular structure where each module is dedicated to solving a single problem along the phase-field simulation strategy.

## Description

OpenPhase-Academic v2.0 provides a comprehensive framework for phase-field modeling, enabling researchers to simulate complex material science phenomena. The library is designed with a modular architecture, making it easy to combine different physical models and simulation components for your specific research needs.

### Simulated Phenomena

OpenPhase can simulate a wide range of physical phenomena, including:

- **Solidification processes**: Solidification in binary systems considering melt flow
- **Grain growth**: Structural transformation in solids including grain growth
- **Phase transformations**: Solid-solid phase transformations accompanied by large deformations and plasticity
- **Multiphase flows**: Phase transformations in fluid systems considering multiphase flow and wetting phenomena
- **Reactive flows**: Heat exchange in reactive gas flows through particle beds and solid combustion
- **Sintering processes**: Solid-state and liquid-phase sintering
- **Heat transfer**: Heat diffusion and heat sources with latent heat effects
- **Elastic deformation**: Small strain mechanics with spectral elasticity solvers
- **Chemical diffusion**: Binary diffusion
- **Magnetic and electric fields**: Magnetism and electrical potential calculations
- **Fracture mechanics**: Crack propagation and fracture field modeling

## Gallery of Simulations

Here are some visual examples of the complex phenomena that can be simulated using OpenPhase:

<div align="center">
  <img src="documentation/UserManual/Figures/oppenphase-examples.png" alt="OpenPhase Simulations Collage" width="900">
</div>

*A collage of various phase-field simulations performed with OpenPhase, showcasing diverse microstructural evolutions, fluid dynamics, and material processes.*

For more simulation videos and visualizations, visit the [Simulation Gallery](https://openphase.rub.de/gallery.html).


## Features

- **Modular Architecture**: Each module handles a specific aspect of phase-field simulations, making it easy to combine different physical models
- **Smart Microstructure Creation**: Easy creation of initial microstructures, realistic nucleation based on various properties
- **Extensive Physical Modules**: Diffusion, thermodynamics, mechanics, chemical, fluid-dynamics, magnetism, electrics, and more
- **All-in-One Design Workflow**: Modular structure which allows easy extensions of the library and simplifies the development of user programs
- **High Performance**: Parallel computing support via OpenMP and MPI hybrid parallelism
- **Comprehensive Examples**: Multiple example simulations demonstrating various capabilities
- **Benchmark Suite**: Validated test cases to verify installation and functionality
- **Cross-Platform**: Supports Linux and Windows (via WSL)
- **Python Bindings**: Python interface for easier scripting and integration
- **Documentation**: User manual and code documentation available

## Requirements

### Linux

**Prerequisites:**
- GNU GCC C++17 compliant compiler (GCC 9.0.0 or greater).
- FFTW3 library (www.fftw.org)

**Note:** The OpenPhase project is developed for Linux, but it is not restricted to use it in any other OS.

#### Installing G++ on Ubuntu

Step 1: Update the packages list:
```bash
sudo apt update
```

Step 2: Install C++ compiler by installation of the development package build-essential:
```bash
sudo apt install build-essential
```

Step 3: Check C++ compiler version:
```bash
g++ --version
```

#### Installing FFTW3 on Ubuntu

Installing fftw-dev package on Ubuntu:
```bash
sudo apt-get install fftw-dev libfftw3-3 libfftw3-dev
```

**Installing the MPI Library for computing Fast Fourier Transforms package on Ubuntu (MPI users):**
```bash
sudo apt-get install libfftw3-mpi-dev
```

### Windows (using WSL)

For Windows users, we recommend using Windows Subsystem for Linux (WSL). This allows you to use the same Linux-based installation and build process.

#### Installing WSL

1. **Open PowerShell as Administrator** and run:
   ```powershell
   wsl --install
   ```

2. **Restart your computer** when prompted.

3. **After restart**, WSL will finish installation. You'll be prompted to create a username and password for your Linux distribution.

4. **Update your Linux distribution**:
   ```bash
   sudo apt update && sudo apt upgrade -y
   ```

5. **Install required tools** in your WSL terminal:
   ```bash
   sudo apt update
   sudo apt install build-essential
   sudo apt-get install fftw-dev libfftw3-3 libfftw3-dev
   ```
   
   **For MPI support (optional):**
   ```bash
   sudo apt-get install libfftw3-mpi-dev
   ```

Once WSL is set up, you can follow the Linux installation instructions below.

## Installation

### Compiling the Library

To compile the OpenPhase library, run `make` in the distribution root directory:

```bash
make
```

This will compile the library and place it into the `./lib` directory. This will take a few minutes depending on your system.

**Note:** If you got a 'Compilation done' message on your terminal, now you are ready to run your first example.

### Compiling Examples

To compile example simulations:

```bash
make examples
```

Examples are located in the `./examples` directory. To run them, navigate into the specific example directory and follow the instructions in the README file.

### Compiling Benchmarks

To compile benchmarks to test the library:

```bash
make benchmarks
```

To run the benchmarks, enter the `./benchmarks` directory and follow the instructions in the README file.

### Changing the Compiler

To change the compiler, use:

```bash
make CXX=<preferred_compiler_name>
```

**Note:** Intel compilers are not supported. Use GNU GCC 9.0.0 or greater.

### Compilation Options

You can customize the build with additional options using:

```bash
make SETTINGS="<option1> <option2> ..."
```

Note that OpenMP parallelism is enabled by default.

Available options:
- `serial`: No OpenMP parallelization
- `mpi-parallel`: Compiles OpenPhase with MPI support
- `debug`: Compiles OpenPhase with debug information and enables additional output and checks via the DEBUG flag
- `silent`: Disables most of the console output from the library

Example:
```bash
make SETTINGS="debug serial"
```

## Project Structure

- `include/`: Header files containing the library API
- `src/`: Source code implementation
- `examples/`: Example simulations demonstrating various features
- `benchmarks/`: Test cases for validation
- `documentation/`: User manual and documentation
- `external/`: External dependencies

## Usage

Each simulation typically consists of:
1. A C++ source file implementing your simulation
2. An input file (`.opi` extension) with simulation parameters

The library assumes that each module that requires user input will read it from a dedicated input file. The input files (with the `.opi` file extension) are located in the corresponding example directories.

Example simulations are located in the `examples/` directory. Navigate to any example directory and follow the instructions in the README file.

To run an example simulation:
```bash
cd examples/NormalGG
./NormalGG
```

## Documentation

- **User Manual**: Available in the `documentation/UserManual/` directory
- **Code Documentation**: Generate using Doxygen by running `doxygen` in the project root. The documentation will be created in the `./documentation` directory.
- **Installation Guide**: See the `INSTALL` file for detailed compilation instructions

A detailed manual for the library should be available from the project's web page: [www.openphase.rub.de](https://www.openphase.rub.de)

## Testing

To test the library, compile and run the benchmarks:

```bash
make benchmarks
cd benchmarks
```

Then navigate to a specific benchmark directory and follow the instructions in the README file to run the test.

## Examples

The `examples/` directory contains numerous simulation examples demonstrating various capabilities of OpenPhase:

### Solidification and Melting
- **SolidificationFe**: Iron solidification
- **SolidificationNiAl**: Nickel-Aluminum solidification
- **SolidificationFeC**: Iron-Carbon solidification
- **SolidificationMgAl**: Magnesium-Aluminum solidification
- **MeltingFe**: Iron melting process
- **EutecticII**: Eutectic solidification

### Grain Growth
- **NormalGG**: Normal grain growth
- **FacetedGG**: Faceted grain growth

### Phase Transformations
- **Pearlite**: Pearlite formation
- **PrecipitationNiTi**: Precipitation in Nickel-Titanium
- **Superalloys**: Superalloy phase transformations

### Sintering
- **SolidPhaseSintering**: Solid-state sintering
- **LiquidPhaseSintering**: Liquid-phase sintering

### Fluid Dynamics and Flow
- **Flow**: Basic flow simulation
- **ThermalCompressibleFlow**: Compressible thermal flow (Poiseuille, cylinder, packed bed)
- **Wetting**: Liquid droplet wetting on solid surfaces

### Reactive Flows
- **ReactiveFlow**: Reactive flow simulations including:
  - Flame0D: Zero-dimensional flame
  - FlameSpeed1D: One-dimensional flame speed
  - FlameWall2D: Two-dimensional flame-wall interaction
  - FlamePackedBed2D: Flame propagation in packed beds

### Heat Transfer
- **HeatDiffusion**: Heat diffusion processes
- **LatentHeat**: Latent heat effects

Each example includes source code and input files to help you get started.

## Authors

OpenPhase is developed by researchers at ICAMS (Interdisciplinary Centre for Advanced Materials Simulation) at Ruhr University Bochum, in collaboration with OpenPhase Solutions GmbH.

For a complete list of contributors, see the `AUTHORS` file.

## Publications

OpenPhase has been used in numerous research publications. For a comprehensive list of featured and related publications, visit the [Publications page](https://openphase.rub.de/publications.html).

## License

This project is licensed under the GNU General Public License v3.0. See the `LICENSE` file for details.

## Contact

- **Email**: [info@openphase.rub.de](mailto:info@openphase.rub.de)
- **Project Website**: [www.openphase.rub.de](https://www.openphase.rub.de)
- **Repository**: [Github Repository](https://github.com/ICAMS/OpenPhase)

## Contributing

If you would like to contribute, fork the repository, work on your changes, and contact the developer team at [info@openphase.rub.de](mailto:info@openphase.rub.de) for any contributions.

## Acknowledgments

OpenPhase has been developed with contributions from researchers at Ruhr University Bochum and OpenPhase Solutions GmbH. We thank all contributors and users for their support and feedback.

---

**Note**: This is the Academic version of OpenPhase (v2.0). For the latest version and commercial licensing options, visit [www.openphase-solutions.com](https://openphase-solutions.com/).
