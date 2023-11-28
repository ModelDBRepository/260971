
#Fig 3 (and Fig 4 control): Columns 1-2 (LFS and 4xHFS with full fluxes)
#Fig 4:                     Columns 3-9 (4xHFS without Ca, without beta-adr ligand, without Glu, without Glu&ACh; LFS without beta (not shown in Fig), without Glu, without Glu&ACh) 
TASKS=(    0      1      1      1      1      1      0      0      0)
CAFLUXES=( 1900.0 1900.0 0.0    1900.0 1900.0 1900.0 1900.0 1900.0 1900.0)
LFLUXES=(  10.0   10.0   10.0   0.0    10.0   10.0   10.0   10.0   10.0 )
GLUFLUXES=(20.0   20.0   20.0   20.0   0.0    0.0    0.0    10.0   0.0   )
ACHFLUXES=(20.0   20.0   20.0   20.0   20.0   0.0    20.0   20.0   0.0   )

TSHORT=5000000
ONSET=4040000

for iflux in `seq 0 8`
do
    CAFLUX=${CAFLUXES[iflux]}
    LFLUX=${LFLUXES[iflux]}
    GLUFLUX=${GLUFLUXES[iflux]}
    ACHFLUX=${ACHFLUXES[iflux]}
    TASK_ID=${TASKS[iflux]}

    if [ "$TASK_ID" == "0" ]
    then 
      echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None #LFS"
      python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None
    else
      echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None #4xHFS-3s"
      python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None
    fi
done

CAFLUX=${CAFLUXES[1]}
LFLUX=${LFLUXES[1]}
GLUFLUX=${GLUFLUXES[1]}
ACHFLUX=${ACHFLUXES[1]}

#For Fig. 3, we need high-resolution Ca data. Better to use these instead of the scripts for iflux=0 and iflux=1 above
echo "python model_nrn_altered_noU_noninterp.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None #LFS"
python model_nrn_altered_noU_noninterp.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None

echo "python model_nrn_altered_noU_noninterp.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None #4xHFS-3s"
python model_nrn_altered_noU_noninterp.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None
