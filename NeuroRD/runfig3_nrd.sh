#!/bin/bash
# Wall clock limit:
#SBATCH --time=8:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3400M
#SBATCH --ntasks=8
TASK_ID=0

TLARGE=12000000
TSHORT=500000
TSHORTB=1100000
CAFLUX=2500.0
GLUFLUX=0.0
LFLUX=0.0
ACHFLUX=0.0

module load Java/1.8.0_131
date

echo "mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 1 1.0 300000.0 $CAFLUX 0.05 0.05 0.05 1 100000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8"
mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 1 1.0 300000.0 $CAFLUX 0.05 0.05 0.05 1 100000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8
