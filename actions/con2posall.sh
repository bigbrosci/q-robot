#!/usr/bin/env bash 
for i in $(find . -name CONTCAR); do 
    path=$(dirname $i)
    mv $i $path/POSCAR
done 
