#!/bin/bash
#SBATCH --account=sc3260
#SBATCH --partition=mic
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --output=bench-slurm.out

setpkgs -a intel_cluster_studio_compiler
setpkgs -e openmpi_intel
setpkgs -a openmpi_intel14

# rebuild executable if required
make >>errfile-slurm 2>>err2-slurm

# run using 8 MPI processes on system size of 1000
./simulate-slurm.sh "$1" 4000 2>&1 >>err2-slurm | \
  tee -a errfile-slurm | grep -Po "^[0-9]+\\.[0-9]+$" | head -n1
