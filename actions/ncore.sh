#!/usr/bin/env bash 
num="$1"
sed -i '/NCOR/d' INCAR 
echo "NCORE = $num " >> INCAR
