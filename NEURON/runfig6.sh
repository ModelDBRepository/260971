
TSHORT=5000000
ONSET=4040000

TRAINISIS=(-200.0 -180.0 -160.0 -140.0 -120.0 -100.0 -80.0 -60.0 -50.0 -40.0 -30.0 -25.0 -20.0 -15.0 -10.0 -5.0 0.0 5.0 10.0 15.0 20.0 25.0 30.0 40.0 50.0 60.0 80.0 100.0 120.0 140.0 160.0 180.0 200.0)
LFLUXES=(  0.0 20.0 20.0 0.0 )
GLUFLUXES=(0.0 0.0  0.0  0.0 )
ACHFLUXES=(0.0 20.0 0.0  20.0)

ICELL=5
CACOEFF=1.0
ECON=0.0004
WNMDA=7.0
LOCATION=apic250-300

for ilflux in 0 1 2 3
do
  for iISI in 0 5 8 10 11 12 13 14 15 16 17 18 19 20 21 22 24 27 32
  do
      LFLUX=${LFLUXES[ilflux]}
      GLUFLUX=${GLUFLUXES[ilflux]}
      ACHFLUX=${ACHFLUXES[ilflux]}
      TRAINISI=${TRAINISIS[iISI]}

      if [ -f nrn_tstop5000000_tol1e-06_onset4040000.0_n120_freq1.0_dur3.0_CaCoeff1.0_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${TRAINISI}_icell${ICELL}_Econ${ECON}_wNMDA${WNMDA}_${LOCATION}.mat ]
      then
        echo "nrn_tstop5000000_tol1e-06_onset4040000.0_n120_freq1.0_dur3.0_CaCoeff1.0_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_Ntrains1_trainT100000.0_pair${TRAINISI}_icell${ICELL}_Econ${ECON}_wNMDA${WNMDA}_${LOCATION}.mat exists"
      else
        echo "python model_nrn_paired.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION None #1Hz-paired-4AP-somatic-burst"
        python model_nrn_paired.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION None #1Hz-paired-4AP-somatic-burst
      fi
  done
done
