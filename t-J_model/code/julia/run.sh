#!/bin/bash -l
#SBATCH -J qpd-tk
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 1000
#SBATCH --time=1:00:00
#SBATCH -A G73-29
#SBATCH -p topola
#SBATCH --output="./reports/tasker.out"

module load apps/julia/1.1.1

srun -n 1 julia tasker.jl
