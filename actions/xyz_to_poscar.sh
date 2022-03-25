#!/usr/bin/env bash 

echo "Converting $1 to POSCAR"

sed -n 1p $1 > temp.xyz 
echo 'Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 pbc="T T T" ' >>  temp.xyz
sed -n '3,$p' $1  >> temp.xyz 
ase gui temp.xyz -o POSCAR
rm temp.xyz 

echo 'Done, For Visualization use command: ase gui POSCAR'

