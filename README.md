# Parallel Computation of Independent Spanning Trees in Bubble Sort Networks

This repository contains the implementation and analysis of parallel computation of Independent Spanning Trees (ISTs) in Bubble Sort Networks, developed as part of the CS-3006: Parallel & Distributed Computing course. The project leverages OpenMPI and OpenMP to distribute workloads and optimize performance for large \( n \), addressing factorial complexity (\( n! \)) in permutation enumeration.

## Project Structure

- `Speedup_Graphs/`: Directory for storing speedup and performance graphs.
- `.gitignore`: File to exclude unnecessary files from version control.
- `callgraph.png`: Generated call graph image from profiling data.
- `commands.txt`: Text file containing useful commands.
- `hosts.txt`: Host configuration file for MPI execution.
- `Makefile`: Makefile for compiling, running, profiling, and cleaning the project.
- `PDC_Report.pdf`: Compiled PDF report of the project.
- `README.md`: This file.
- `V1.cpp`: Source code for Version 1 implementation.
- `V2.cpp`: Source code for Version 2 implementation.

## Prerequisites

- **Compiler**: MPI-enabled C++ compiler (`mpic++`).
- **Libraries**: OpenMPI and OpenMP.
- **Tools**: 
  - `gprof` for profiling.
  - `gprof2dot.py` for converting profiling data to a call graph (download using `wget https://raw.githubusercontent.com/jrfonseca/gprof2dot/master/gprof2dot.py`).
  - `dot` (from Graphviz) for rendering the call graph.
- **Operating System**: Ubuntu 24.04 (or compatible).

## Installation

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. Install dependencies:
   ```bash
   sudo apt update
   sudo apt install openmpi-bin openmpi-common libopenmpi-dev gprof graphviz
   ```

3. Download `gprof2dot.py`:
   ```bash
   wget https://raw.githubusercontent.com/jrfonseca/gprof2dot/master/gprof2dot.py
   ```

4. Ensure `hosts.txt` is configured with the hostnames or IP addresses of machines for MPI execution.

## Usage

### Compilation and Execution
Run the default target to compile and execute the program:
```bash
make all
```
This uses `mpirun` with 1 process and outputs performance metrics (e.g., elapsed time, CPU usage, memory usage).

### Profiling
To profile the program and generate a call graph:
```bash
make profile
```
- Compiles with profiling enabled (`-g -pg`).
- Runs the program and generates `gmon.out` profiling data.
- Converts the data to `analysis.txt` using `gprof`.
- Generates `callgraph.dot` and renders it as `callgraph.png` using `gprof2dot.py` and `dot`.
- Outputs performance metrics similar to the default target.

### Cleaning
Remove generated files:
```bash
make clean
```

## Report
The project report is available as `PDC_Report.pdf`, generated using LaTeX. It details the implementation, optimization strategies, and performance analysis for both Version 1 and Version 2.


## Acknowledgments
- Aneeq Ahmed Malik (22i-1167)
- Abdullah Mehmood (22i-0973)
- Kalbe Raza (22i-0794)
