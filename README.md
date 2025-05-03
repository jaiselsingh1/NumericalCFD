<div align="center">
  <img src="./numerical_cfd_banna.png?raw=true" alt="banner image titled numerical cfd">
</div>

# A Collection of Numerical Computational Fluid Dynamics Codes

This repository contains a collection of Julia and Python codes implementing various numerical methods for solving problems in computational fluid dynamics (CFD).  The codebase covers a range of topics, including incompressible flow, supersonic wave drag, and ordinary differential equation (ODE) integration.

## Key Features

* **Incompressible Flow Solver (Python & MATLAB):**  Solves for the potential flow around a NACA 0012 airfoil using the Successive Over-Relaxation (SOR) method.  Includes comparison with experimental data.  The MATLAB version also generates contour plots of the potential field.
* **Supersonic Wave Drag Calculation (Julia & C++):** Computes the wave drag for a parabolic airfoil profile at various Mach numbers using a numerical method and compares the results to an analytical solution.  The C++ code offers two approaches: a sliding window method and a full-domain method, demonstrating differences in computational efficiency and accuracy.
* **ODE Integration Methods (C++):** Implements Euler-Cauchy, Improved Euler-Cauchy, and Runge-Kutta 4th order methods for solving ordinary differential equations. Includes error analysis and visualization.
* **Laval Nozzle Flow Solver (Julia):**  A numerical solver for flow through a Laval nozzle, comparing numerical results against an exact solution. Includes visualizations of the solution.


## Technologies Used

* **Programming Languages:** Julia, Python, C++
* **Libraries:**
    * Julia: `LinearAlgebra`, `Plots`
    * Python: `numpy`, `matplotlib`
    * C++: `matplotplusplus`, `Eigen`
* **Build System (C++):** CMake


## Prerequisites

**For all codes:**

* A suitable text editor or IDE.

**For Julia codes:**

* Julia 1.11.1 or later (check your version with `julia --version`).
* Install necessary Julia packages:  You can install the required packages using the `Project.toml` file included in the repository. Navigate to the repository directory in your terminal and run `] instantiate` within the Julia REPL.

**For Python codes:**

* Python 3.7 or later.
* Install necessary Python packages: `pip install numpy matplotlib`

**For C++ codes:**

* A C++ compiler (e.g., g++).
* CMake (for building the projects).
* Install external libraries:  The `build.sh` scripts in the `ode_integration` and `supersonic_wave_drag` directories handle the installation of `matplotplusplus` and `Eigen`.  Follow the instructions within those scripts.  Note that these scripts are designed for Linux/WSL environments; modification might be necessary for other operating systems.


## Installation Guide

**1. Cloning the Repository:**

Clone the repository to your local machine using Git:

```bash
git clone <repository_url>
```

**2. Julia Codes:**

* Navigate to the Julia code directory using the command line.
* Open the Julia REPL within the directory.
* Run `] instantiate` to install the required packages listed in `Project.toml`.
* Run the Julia scripts (e.g., `julia laplaceeq.jl`, `julia supersonicwavedrag.jl`, `julia criticalflow_laval.jl`, `julia NACA0012.jl`) from your terminal.

**3. Python Codes:**

* Install the required packages using pip: `pip install numpy matplotlib`.
* Run the Python scripts (e.g., `python incompressible_flow.py`, `python laval_nozzle.py`) from your terminal.

**4. C++ Codes:**

* Navigate to the `ode_integration` or `supersonic_wave_drag` directory.
* Make sure you have a C++ compiler and CMake installed.
* Run the `build.sh` script in the respective directory.  This script will install external libraries and compile the code. This script assumes a Linux or WSL environment and uses `pacman` for package management. Adaptations may be needed for other systems (e.g., using `apt-get` on Debian-based systems or manually installing the dependencies and compiling with CMake).
* The compiled executable will be located in the `build` directory.  Run it from the terminal (e.g., `./build/hw -c config.inp`). Results will be in the `/results` directory.

## Usage Examples

**Julia Laval Nozzle Solver:**

To run the Laval nozzle solver with an outlet velocity of -sqrt(1.5), run:

```julia
julia criticalflow_laval.jl
```

The script will generate plots comparing the numerical and exact solutions.

**Python Incompressible Flow Solver:**

To run the incompressible flow solver with a grid spacing of 0.02, run:

```bash
python incompressible_flow.py
```

The script will output convergence information and generate a contour plot of the potential flow (MATLAB version).


**C++ ODE Solver:**

After building the `ode_integration` project (following the installation steps for C++), run the executable:

```bash
./build/hw -c config.inp
```

The `config.inp` file (which needs to be created) specifies the ODE parameters and integration method.  The script will produce plots showing the solution and error analysis.  Refer to `ode_integration/source/main.cpp` for details on creating the input file.

**C++ Supersonic Wave Drag Solver:**

After building the `supersonic_wave_drag` project, run the executable similarly to the ODE solver:

```bash
./build/hw -c config.inp 
```

This will output the computed wave drag coefficients for various Mach numbers (refer to `supersonic_wave_drag/source/main.cpp` for details and results).


## Project Structure

```
NumericalCFD/
├── ode_integration/       # C++ ODE solver
│   ├── build.sh           # build script
│   ├── CMakeLists.txt    # CMake build file
│   ├── include/          # Header files
│   │   ├── eulerCauchy.h
│   │   ├── improvedEulerCauchy.h
│   │   └── rungeKutta4.h
│   └── source/           # Source code
│       ├── eulerCauchy.cpp
│       ├── improvedEulerCauchy.cpp
│       ├── main.cpp
│       └── rungeKutta4.cpp
├── supersonic_wave_drag/ # C++ supersonic wave drag solver
│   ├── build.sh
│   ├── CMakeLists.txt
│   ├── include/
│   │   └── solver.h
│   └── source/
│       ├── main.cpp
│       └── solver.cpp
├── incompressible_flow.m  # MATLAB incompressible flow solver
├── incompressible_flow.py # Python incompressible flow solver
├── criticalflow_laval.jl # Julia Laval nozzle flow solver
├── laplaceeq.jl          # Julia Laplace equation solver
├── supersonicwavedrag.jl # Julia supersonic wave drag solver
├── NACA0012.jl           # Julia NACA 0012 airfoil functions
├── Project.toml          # Julia project file
├── Manifest.toml         # Julia manifest file
└── README.md             # This file
```

## Configuration

The C++ projects (`ode_integration` and `supersonic_wave_drag`) use a `config.inp` file for configuration.  You'll need to create this file.  The file format will depend on the specific executable (see the respective `main.cpp` files for details on the required parameters and their format).


## Contributing Guidelines

No explicit contributing guidelines were found within the provided codebase.


## License Information

No license information was explicitly provided within the codebase.  Please clarify the license under which this code is distributed.
