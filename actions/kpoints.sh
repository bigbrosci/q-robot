#!/bin/bash 

##########################################
# Rodrigo García-Muelas @ iciq-Tgn       #
# Jan 10, 2013                           #
# For automatic generation of Γ-centered #
# Monkhorst-Pack grid of KPOINTS         #
##########################################

# Input ##################################   
# 1 2 3 # number of divisions            #
# 4     # title of first line (optional) #
##########################################

if [ -z $4 ] ; then  # If title if empty

cat >KPOINTS<<!
K-POINTS
 0
Gamma
  $1 $2 $3
  0 0 0
!

else                 # If title if specified

cat >KPOINTS<<!
$4
 0
Gamma
  $1 $2 $3
  0 0 0 
!

fi
