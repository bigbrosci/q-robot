#!/usr/bin/env bash 
# get the TS image from multiple NEB jobs  for the next freq calcs.
for i in *; do 
cd $i  
ts=$(ta| sort -k 8| head -n 1 | awk '{print $1}') 
cp $ts ../../2_freq/$i; 
cd -
done
