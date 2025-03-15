#!/usr/bin/env bash 
# for file in IBZKPT KPOINTS faster; do find  . -name ${file}* | xargs -I {} rm {} ; done
file="$1"
find  . -name ${file}* | xargs -I {} rm {} 
