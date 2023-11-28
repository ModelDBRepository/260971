
#This will take approximately 3 x 5 hours
for flux in 0.0 0.005 0.05
do
  echo "python simsteadystate_flexible.py $flux"
  python simsteadystate_flexible.py $flux
done

