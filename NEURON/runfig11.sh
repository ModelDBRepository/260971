
TSHORT=5000000
ONSET=4040000

#              (PKC activation rate altered) (S831 phosph. by PKC blocked)   (GluR2 insertion rate)
#              (Fig. 1D-G)                   (Fig. 1J)                       (Fig. 1A)             
ALTEREDS=(      411 411 411 411 411           166-169-181-184-218-221-233-236 385-387-389          )
ALTEREDCOEFFS=( 0.1 0.3 1.0 3.0 10.0          0.0-0.0-0.0-0.0-0.0-0.0-0.0-0.0 22.4-22.4-22.4       )

CAFLUX=1900
GLUFLUX=20.0
ACHFLUX=20.0
LFLUX=10.0
#control 4xHFS, used in Fig. 11A,B,I,J
echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None"
python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None

#Fig. 11A: 4xHFS with altered (old) GluR2 membrane insertion rated
echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None Ca 1.0 ${ALTEREDS[6]} ${ALTEREDCOEFFS[6]}"
python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None  Ca 1.0 ${ALTEREDS[6]} ${ALTEREDCOEFFS[6]}

#Fig. 11B: Old vs. new CaM activation model
echo "python model_nrn_oldCaM_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None CaM 0.55"
python model_nrn_oldCaM_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None CaM 0.55

#Fig. 11C: Old vs. new CaM activation model
echo "python simsteadystate_li2020.py 0.0 50.0 0.0"
python simsteadystate_li2020.py 0.0 50.0 0.0
echo "python simsteadystate_oldCaM_li2020.py 0.0 50.0 0.0"
python simsteadystate_oldCaM_li2020.py 0.0 50.0 0.0

CAFLUX=1900
#Fig. 11D--G: PKCp activation, rate fitted to data from LFS experiments
for ialtered in 0 1 2 3 4
do
  echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5.0 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None Ca 1.0 ${ALTEREDS[ialtered]} ${ALTEREDCOEFFS[ialtered]}"
  python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5.0 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None Ca 1.0 ${ALTEREDS[ialtered]} ${ALTEREDCOEFFS[ialtered]}
done

#Fig. 11H: Fit the PKA-cAMP binding rate
#1) Calculate the PKA activation time series using the original model
echo "python model_nrn_testPKA_withdiss.py 22000 1e-08 800.0 1 1.0 16000.0 0.64"
python model_nrn_testPKA_withdiss.py 22000 1e-08 800.0 1 1.0 16000.0 0.64

#2) Calculate the PKA activation time series using the single-step PKA activation model with different reaction rates and see which fits best
echo "python fit_cAMP_withdiss_1d.py 0.64"
python fit_cAMP_withdiss_1d.py 0.64

#Fig. 11I: Old vs. new PKA activation model
echo "python model_nrn_oldPKA_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None"
python model_nrn_oldPKA_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None

#Fig. 11J: PKC does not phosphorylate S831
echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None Ca 1.0 ${ALTEREDS[5]} ${ALTEREDCOEFFS[5]}"
python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None  Ca 1.0 ${ALTEREDS[5]} ${ALTEREDCOEFFS[5]}

