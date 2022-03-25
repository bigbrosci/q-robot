#!/usr/bin/env bash 
E0=$(grep "f  =" OUTCAR |awk '{print  $10}'| paste -sd+ | bc)
E1=$(echo "scale = 3 ; $E0/2000" |bc)
echo $E1

