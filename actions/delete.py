#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Delete atoms from POSCAR '''
import sys
from lattice import *

if len(sys.argv[:]) < 3: 
    print('\nCommand Usgae:  delete.py file_to_be_deleted  atoms_to_be_deleted\n')
    print('For Example  :  delete.py POSCAR C H O 10 12 14' )
    print('\nYou are going to delete all C H O atoms and the 12th and 14th atoms from POSCAR file' )
    exit()
else:
    file_in = sys.argv[1]
    atom_s = sys.argv[2:]

lines = read_car(file_in)[0]
atom_list = get_atom_list(lines, atom_s)
atom_list =sorted(atom_list, reverse=True)
print(atom_list)

'''Important:
If we want to delete atoms: 1 2 3 4 5 6 7 
The delete sequence would be 7 6 5 4 3 2 1.
That is the reason we sort the list in Descending Order by using sorted() function
'''

#print atom_list
out_name = file_in + '_deleted'
for i in atom_list:
    lines = delete_one_atom(lines, i)
file_out = open(out_name, 'w')
file_out.writelines(line for line in lines)
file_out.close()

print('\nThe output file is:\t%s' %(out_name))
