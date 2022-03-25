#!/usr/bin/env python3

from lattice import *

try:
    lines = read_car('CONTCAR')[0]
except FileNotFoundError:
    lines = read_car('POSCAR')[0]
 
a, b, c, A, V = get_abc(lines)

print('Length_a\tLength_b\tLength_c\tArea/A^2\tVolume/A^3\t')
print('%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t' %(a, b, c, A, V))
