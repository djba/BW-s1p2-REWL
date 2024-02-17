CUDA Fortran programs for the simulation of the spin-$1/2$ Baxter-Wu model using replica exchange Wang-Landau sampling

All the necessary files are provided in the gpu directory. The source code is structured as follows:
- Makefile
- parameters.cuf - in this file we configure the simulation parameters.
- host-functions.f90 - here we define the functions which run on the CPU.
- device-functions.cuf - in this file we define the subroutines with device and global attributes.
- bw-s1p2-gpu.cuf - the main program file

The following variables can be set by modifying the parameters.cuf file accordingly:
- L - system size
- stage0 - the stage from which we compute the heat capacity maximum and test its fluctuations
- replica_step - After this number of steps we try replica exchanges between walkers from neighbouring windows.
- flatness_step - After this number of steps we check the histogram for flatness.
- deltaE - the width of energy windows (the overlapping region was set to 50% of the width)
- NR - number of replicas for each energy window
- BLOCK_SIZE - number of threads per block
- eps_min - parameter for histogram flatness test

The cpu directory contains the source code for the standard CPU-based Wang-Landau sampling (running on a single core). 

Makefiles are provided for both the gpu and cpu code.
Prior running the Makefile please update the PATH variable according to your NVIDIA HPC SDK installatiion which contains the pgf90 compiler (set_path.sh script)

source set_path.sh

In order to run the simulation change directory to 'data' and modify script.sh (set the desired number of consecutive runs and set the system size which will appear in the output file name in order to match the value from parameters.cuf). 
