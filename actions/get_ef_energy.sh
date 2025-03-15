echo 'EF,energy' > energy.csv
for i in $(cat list_ef); do echo ${i}, $(grep '  without' OUTCAR_$i |tail -n 1| awk '{print $NF}'); done  >> energy.csv

