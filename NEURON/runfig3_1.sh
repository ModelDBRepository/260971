
TSHORT=5000000
ONSET=4040000

GLUR1S=(0.0 2.0 0.7)
GLUR2S=(2.0 0.0 1.3)

for iratio in 0 1 2
do
  GLUR1=${GLUR1S[iratio]}
  GLUR2=${GLUR2S[iratio]}

  CAFLUX=1900.0
  LFLUX=10.0
  GLUFLUX=20.0
  ACHFLUX=20.0

  echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None GluR1-GluR1_memb-GluR2-GluR2_memb ${GLUR1}-${GLUR1}-${GLUR2}-${GLUR2} #LFS"
  python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 900 5 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 1 100000 None GluR1-GluR1_memb-GluR2-GluR2_memb ${GLUR1}-${GLUR1}-${GLUR2}-${GLUR2}

  echo "python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None GluR1-GluR1_memb-GluR2-GluR2_memb ${GLUR1}-${GLUR1}-${GLUR2}-${GLUR2} #4xHFS-3s"
  python model_nrn_altered_noU.py ${TSHORT} 1e-6 $ONSET 100 100 3.0 $CAFLUX $LFLUX $GLUFLUX $ACHFLUX 4 4000 None GluR1-GluR1_memb-GluR2-GluR2_memb ${GLUR1}-${GLUR1}-${GLUR2}-${GLUR2}

done
