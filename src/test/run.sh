#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=2
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=14
#SBATCH --time=2:00:00
module load gcc/10.2.0 intel-mkl openmpi-ucx-gpu
ulimit -s unlimited 
export OMP_NUM_THREADS=14
export OMP_STACKSIZE=1024m
srun --export=all -N 1 -n 2 -c 14 ../gfortran/bin/cccd ccc.in

