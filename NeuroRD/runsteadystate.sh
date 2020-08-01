#!/bin/bash
# Wall clock limit:
#SBATCH --time=8:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3600M
#SBATCH --ntasks=8

TLARGE=12000000

date

echo "mpirun -np 8 python runmodel_parallel.py ${TLARGE} 0.01 0 0 1 0 0 0 0 0 0 0 None 1000.0 8"
mpirun -np 8 python runmodel_parallel.py ${TLARGE} 0.01 0 0 1 0 0 0 0 0 0 0 None 1000.0 8

