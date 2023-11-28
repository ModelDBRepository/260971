
TSHORT=5000000
ONSET=4040000

TRAINISIS=(-200.0 -180.0 -160.0 -140.0 -120.0 -100.0 -80.0 -60.0 -50.0 -40.0 -30.0 -25.0 -20.0 -15.0 -10.0 -5.0 0.0 5.0 10.0 15.0 20.0 25.0 30.0 40.0 50.0 60.0 80.0 100.0 120.0 140.0 160.0 180.0 200.0 -210.0 -220.0 -230.0 -2.5 2.5)
LFLUXES=(  0.0 0.05 0.05 0.0 )
GLUFLUXES=(0.0 0.0  0.0  0.0 )
ACHFLUXES=(0.0 0.05 0.0  0.05)


ICELL=1
CACOEFF=1.0
ECON=0.001
WNMDA=3.2
LOCATION=apic250-300
NSYN=10

for ilflux in 1 0 2 3 #The main simulations for 4 different conditions (ACh on/off; NA on/off)
do
  for iISI in 0 5 8 10 11 12 13 14 15 16 17 18 19 20 21 22 24 27 32 #make it `seq 0 37 ` to have high-resolution ISIs
  do
      LFLUX=${LFLUXES[ilflux]}
      GLUFLUX=${GLUFLUXES[ilflux]}
      ACHFLUX=${ACHFLUXES[ilflux]}
      TRAINISI=${TRAINISIS[iISI]}

      if [ -f nrn_tstop5000000_tol1e-06_onset4040000.0_n120_freq1.0_dur3.0_CaCoeff1.0_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_contnm3560000.0_600000.0_Ntrains1_trainT100000.0_pair${TRAINISI}_icell${ICELL}_pulseamp5.0_Nsyn10_Ninputs1_Econ${ECON}_wNMDA${WNMDA}_${LOCATION}.mat ]
      then
        echo "nrn_tstop5000000_tol1e-06_onset4040000.0_n120_freq1.0_dur3.0_CaCoeff1.0_Lflux${LFLUX}_Gluflux${GLUFLUX}_AChflux${ACHFLUX}_contnm3560000.0_600000.0_Ntrains1_trainT100000.0_pair${TRAINISI}_icell${ICELL}_pulseamp5.0_Nsyn10_Ninputs1_Econ${ECON}_wNMDA${WNMDA}_${LOCATION}.mat exists, skipping"
      else
        echo "python model_nrn_paired_contnm_var.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION $NSYN None 3560000 600000 #1Hz-paired-4AP-somatic-burst"
        python model_nrn_paired_contnm_var.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION $NSYN None 3560000 600000 #1Hz-paired-4AP-somatic-burst
      fi
  done
done

for ilflux in 1 #Redo some simulations for ACh on and NA on where the Ca2+ transients are saved in high resolution. This is needed for panels E-G.
do
  for iISI in 6 10 20
  do
    LFLUX=${LFLUXES[ilflux]}
    GLUFLUX=${GLUFLUXES[ilflux]}
    ACHFLUX=${ACHFLUXES[ilflux]}
    TRAINISI=${TRAINISIS[iISI]}
    echo "python model_nrn_paired_contnm_var_noninterp.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION $NSYN None 3560000 600000 #1Hz-paired-4AP-somatic-burst"
    python model_nrn_paired_contnm_var_noninterp.py ${TSHORT} 1e-6 $ONSET 120 1.0 3 $CACOEFF $LFLUX $GLUFLUX $ACHFLUX 1 100000 $TRAINISI $ICELL $ECON $WNMDA $LOCATION $NSYN None 3560000 600000 #1Hz-paired-4AP-somatic-burst
  done
done
