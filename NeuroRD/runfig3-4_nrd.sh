#!/bin/bash
# Wall clock limit:
#SBATCH --time=36:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=3600M
#SBATCH --ntasks=8

#Fig 4 (and Fig 5 control): Columns 1-2 (LFS and 4xHFS with full fluxes)
#Fig 5:                     Columns 3-9 (4xHFS without Ca, without beta-adr ligand, without Glu, without Glu&ACh; LFS without beta (not shown in Fig), without Glu, without Glu&ACh) 
TASKS=(    0      1      1      1      1      1      0      0      0)
CAFLUXES=( 1900.0 1900.0 0.0    1900.0 1900.0 1900.0 1900.0 1900.0 1900.0)
LFLUXES=(  10.0   10.0   10.0   0.0    10.0   10.0   10.0   10.0   10.0 )
GLUFLUXES=(20.0   20.0   20.0   20.0   0.0    0.0    0.0    10.0   0.0   )
ACHFLUXES=(20.0   20.0   20.0   20.0   20.0   0.0    20.0   20.0   0.0   )

TLARGE=12000000
TSHORT=500000
TSHORTB=1100000

for iflux in `seq 0 8`
do
    CAFLUX=${CAFLUXES[iflux]}
    LFLUX=${LFLUXES[iflux]}
    GLUFLUX=${GLUFLUXES[iflux]}
    ACHFLUX=${ACHFLUXES[iflux]}
    TASK_ID=${TASKS[iflux]}

    if [ "$TASK_ID" == "0" ]
    then 
      echo "mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 900 5.0 3 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8 #LFS"
      mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 900 5.0 3 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8 #LFS
    else
      echo "mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 100 100.0 3 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8 #4*HFS-3s (4s)"
      mpirun -np 8 python runmodel_parallel.py ${TSHORT} 0.01 40000 100 100.0 3 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 tstop${TLARGE}_tol0.01_onset0.0_n0_freq1.0_dur0.0_flux0.0_Lflux0.0_Gluflux0.0_AChflux0.0_Ntrains0_trainT0.0_8seeds.mat 10.0 8 #4*HFS-3s (4s)
    fi
done

