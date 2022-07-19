#!/usr/bin/env bash 
folder=$1
if [ ! -d "$1" ]; then 
mkdir -p  "$1"
find . -maxdepth 1  -type f | xargs -I {} mv {} "$1" -f 
else
echo "$1" is exist! Watch out!
exit
fi 

