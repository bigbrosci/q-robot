#!/usr/bin/env bash 
# Written by Qiang 
# To get the Magnetization of specific atoms
# To use it: 
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method A : If you want to get the information of the 14th atom 
###########   bash get_mag.sh 14
# Method B : If you want to get the information of atoms: 2 4 6 8 10
###########   for i in 2 4 6 8 10 ; do bash get_mag.sh $i ; done 
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get the line number of the last magnetization (x)
Nline=$(grep -nr 'magnetization (x)'  OUTCAR  | tail -n 1 | awk '{print $1}' |sed  's/://g')

# Get  the total atoms number 
Natom=$(sed -n 7p CONTCAR  |tr ' ' '\n'  |  sed '/^$/d' | paste -sd+ | bc)

# Start line of Magnetization part 
Nstart=$(($Nline+4))

# End line of Magnetization part 
Nend=$(($Nstart+$Natom-1))

# Extract the Magnetization part 
sed -n "$Nstart, $Nend p" OUTCAR > mag-temp 

# Read and print the specific atom's magnetization 
sed -n "$1 p" mag-temp
