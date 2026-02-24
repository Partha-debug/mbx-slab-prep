#!/bin/bash
#SBATCH --job-name=mbx_64w
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --output=run_%j.out
#SBATCH --error=run_%j.err

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 FFTW/3.3.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Set `LAMMPS_EXE` to your LAMMPS binary before submitting, for example:
#   export LAMMPS_EXE=/path/to/lmp
# If not set, fall back to `lmp` on your PATH.
LAMMPS_EXE=${LAMMPS_EXE:-lmp}

cd $SLURM_SUBMIT_DIR

mpirun -np $SLURM_NTASKS $LAMMPS_EXE -in in.equil.lammps
