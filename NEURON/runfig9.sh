
TSHORT=5000000
ONSET=4040000

IMEASS=(0 1 2 3 4 5 6 7 8 9 10 7 8)
LABELS=(fewer fewer fewer fewer fewer fewer fewer manyb manyb fewer fewer fewer fewer)
SEEDS=(1 1 1 1 1 1 1 30 30 1 1 1 1)
NSAMPS=(1000 1000 1000 1000 1000 1000 1000 500 500 1000 1000 1000 1000)

for ITASK in 0 1 2 3 4 5 6 7 8 9 10 11 12
do
  DATA_ID=${IMEASS[ITASK]}
  LABEL=${LABELS[ITASK]}
  SEED=${SEEDS[ITASK]}
  NSAMP=${NSAMPS[ITASK]}
  echo "python fitter_${LABEL}_check_given.py $DATA_ID $SEED $NSAMP 1 1 0"
  python fitter_${LABEL}_check_given.py $DATA_ID $SEED $NSAMP 1 1 0
done
