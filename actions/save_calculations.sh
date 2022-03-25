#!/usr/bin/env bash 
# To save the current  calculations
if [ $# -eq 1 ]; then 
  for i in CONTCAR POSCAR  OUTCAR DOSCAR XDATCAR OSZICAR vasprun.xml ; do 
    mv "$i" "$i"-$1
  done
  cp CONTCAR-$1 POSCAR 
else  
  num=$(ls OUTCAR* | wc -l)
  for i in CONTCAR POSCAR  OUTCAR DOSCAR XDATCAR OSZICAR vasprun.xml ; do
    mv "$i" "$i"-$num
  done  
  cp CONTCAR-$num POSCAR 
fi 
