#!/usr/bin/env python3
'''Read the bader outputs (VTST) and calculate the charges on the atoms '''
import sys
from lattice import * 
from potcar import *
atom_s = sys.argv[1:]

lines = read_car('POSCAR')[0]
atom_list = get_atom_list(lines, atom_s)
ele_list  = [get_ele_name(lines, i) for i in atom_list]
dict_potcar = get_multiple_potcar_infor('POTCAR')[0]


try: 
    with open ('ACF.dat', 'r') as infile:
        lines_bader = infile.readlines()
        charge_lines = [line for line in lines_bader if line.rstrip().split()[0].isdigit()]
except:
    print('Can not find ACF.dat file, run Command below and  then rerun this.\n\n')
    print('chgsum.pl AECCAR0 AECCAR2 && bader CHGCAR -ref CHGCAR_sum')
#
charge_sum = []
title = 'Elment\tNo.\tCHARGE\tZVAL\tZVAL-CHARGE'
print(title)
charge_sum.append(title+'\n')
for num, atom in enumerate(atom_list):
    ele_atom = get_ele_name(lines, atom)
    charge  = float(charge_lines[atom -1].strip().split()[4])
    try:
        zval     = dict_potcar.get(ele_atom).get('ZVAL')
    except:
        ele_atom = ele_atom + '_sv'
        zval     = dict_potcar.get(ele_atom).get('ZVAL')
    charge_atom = '%s\t%s\t%.4f\t%s\t%.4f' %(ele_atom, atom, charge, zval, float(zval) - float(charge) )
    print(charge_atom)
    charge_sum.append(charge_atom+'\n')

file_out = open('bader_charges.dat', 'w')
file_out.writelines(charge_sum)
file_out.close()
