#!/usr/bin/env bash 
#modify the NCORE parameter in INCAR file
num="$1"
sed -i '/NCOR/d' INCAR 
echo "NCORE = $num " >> INCAR
