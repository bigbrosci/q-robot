#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This script is from Tang Gang.I only modified it a little bit.

1) Find Fermi Level from OUTCAR
2) Find the top of the valence band and the bottom of the conduction band from EIGENVAL  

 Usage: chmod +x gbandedges
        gbandedges [path] [ef]
        ef: trial Fermi energy
        ef read from OUTCAR if not given
 Last modified on Jun.21,2014
'''
from sys import argv, exit
from os import path

ief = 0
if len(argv) == 1:
    dir = '.'
elif len(argv) == 2:
    try:
        ef = float(argv[1])
        dir = '.'
        ief = 1
    except: 
        dir = argv[1].rstrip('/')
elif len(argv) == 3:
    dir = argv[1].rstrip('/')
    ef = float(argv[2])
    ief = 1
else:
    print( "Usage: gbandedges [path] [Ef]")
    print( "       ef: trial Fermi energy")
    print( "       ef read from OUTCAR if not given")
    exit(0)

### Get E_fermi 
fn_outcar = dir + '/OUTCAR'
if not ief:
    if path.exists(fn_outcar):
        for line in open(fn_outcar, 'r'):
            if line.find('E-fermi') > -1:
                ef=float(line.split()[2])
    else:
        print( 'OUTCAR does not exist!')    
        exit(0)

### Read EIGENVAL file
fn_eig = dir + '/EIGENVAL'
if not path.exists(fn_eig):
    print( 'EIGENVAL does not exist!'    )
    exit(0)   
f_eig = open(fn_eig)
lines = f_eig.readlines()
f_eig.close

nspin   = int(lines[0].split()[3])   # Get ISPIN  
nkpoint = int(lines[5].split()[1]) # Get the KPOINTS numbers
nband   = int(lines[5].split()[2])   # Get the band numbers

eup=ef+99.0
edown=ef-99.0
k=[]
for i in range(nkpoint):
    k.append(lines[7 + i * (nband + 2)].split()[0:3])   # For each KPOINT, there are two extralines
    for j in range(nband):
        for s in range(nspin):
            band_energy = float(lines[8 + i * (nband + 2) + j].split()[s + 1])
            d = band_energy - ef 
            if d > 0 :
                if band_energy < eup:   # Find the lowest unoccupied orbital: LUMO
                    eup = band_energy
                    kup = i
                    bup = j
                    sup = s
            if d < 0:          
                if band_energy > edown: # Find the highest occupied orbital: HOMO
                    edown = band_energy
                    kdown = i
                    bdown = j
                    sdown = s

print('%-8s  %6s    %7s  %5s  %4s' % ('', 'Energy', 'KPOINTS', 'BAND', 'SPIN') )
print('%-8s  %5.3f  %7s  %7s  %4s' % ('CBM', eup, kup+1, bup+1, sup+1) )
print('%-8s  %5.3f  %7s  %7s  %4s' % ('VBM', edown, kdown+1, bdown+1, sdown+1) )
print('%-8s  %5.3f' % ('Efermi', ef)) 
print('%-8s  %5.3f' % ('Band-gap', eup-edown))
