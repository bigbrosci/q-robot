#!/usr/bin/env python3
''' Sort the POSCAR coordinates in z directions
 New version on 08-05-2018
 Author: Qiang 
 Python version : >= 2.6  
 To use it:  python sortcar.py file_to_be_sorted 
 for example:  python sortcar.py POSCAR 
'''
from collections import defaultdict
import numpy as np
import sys 
from lattice import * 

in_file = sys.argv[1]

def get_elements(ele, lines):
    coord_total = []
    my_list = []
    tf_dict = {}
    list_range = dict_car2.get(ele)
    
    for j in list_range:
        coord_list = [float(i) for i in lines[j+8].strip().split()[0:3]]
        coord_list.append(j)
        tf = '  '.join(lines[j+8].strip().split()[3:])
        tf_dict.update({j:tf})
        my_list.append(tuple(coord_list))              
    dtype = [('x', float), ('y', float), ('z', float), ('label', int)]
    data =  np.array(my_list, dtype=dtype)
    data2 = np.sort(data,order='z') 
    return data2, tf_dict

## Generate the Sorted POSCAR file

lines, dict_car1, dict_car2 = read_car(in_file)
out_name = in_file + '_sorted'
file_out = open(out_name, 'w')
file_out.writelines(lines[0:9]) # write the head part

# Write the coordinate part
for i in dict_car1.keys(): 
    data2, tf_dict = get_elements(i, lines)
    for j in data2:
        j_list = list(j)
        j_list[-1] = tf_dict.get(j_list[-1])
        for k in j_list:
            file_out.write('%s  ' %(k))
        file_out.write('\n')    

print('\n The output file is : %s \n ' %(out_name))
