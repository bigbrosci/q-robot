#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lattice import *
import numpy as np

def get_atom_orbital_list(lines, argv_list):
    '''
    To split the input arguments into atom_list and orbital_list, especially for dos_extract.py
    '''
    dict_car1, dict_car2 = get_dicts(lines)
    ele_list = [i.split('-')[0] for i in dict_car2.keys()]
    atom_list_raw = []
    orbital_list_raw = []
    
    for i in argv_list:
        if '-' not in i:   # for example: atoms 1-5 
            if str(i).isdigit():
                atom_list_raw.append(str(i))
            elif i in ele_list: # For Example: Ru C H O                                
                atom_list_raw.append(i)
            else:
                orbital_list_raw.append(i)
# for one input item, If 1) there is no '-' in it, 2) it is not in element list, so It will be the orbital name 
        else:
            num_start = int(i.split('-')[0])
            num_end   = int(i.split('-')[1])+1
            for j in range(num_start, num_end):
                if j not in atom_list_raw:
                    atom_list_raw.append(j)
    return atom_list_raw, orbital_list_raw

def get_orbital_dict():
    ###  Get sprcific orbital dos information for one specific atom ###
    orbital_name = ['s','py','pz','px','dxy','dyz','dz2','dxz','dx2','f1','f2','f3','f4','f5','f6','f7']
    orbital_num = range(1,17)
    #dict_orbital =  {orbital_name[i]:orbital_num[i] for i in range(0, len(orbital_name))}
    dict_orbital_1 = dict(zip(orbital_name, orbital_num))
    
    
    orbital_1 = ['p', 'd', 'f', 't2g', 'eg']
    orbital_2  =[ ['px','py','pz'], ['dxy','dyz','dxz','dz2','dx2'], 
    ['f1','f2','f3','f4','f5','f6','f7'], ['dxy','dyz','dxz'],  ['dz2','dx2']]
       
    dict_orbital_2 = dict(zip(orbital_1, orbital_2))    
    return dict_orbital_1, dict_orbital_2

def get_orbital_list(orbital_list_raw):
    dict_orbital_1, dict_orbital_2 = get_orbital_dict()
    ### Convert input information to orbital names ###
    orbital_list = []
    for i in orbital_list_raw:
        if i in dict_orbital_2.keys():
            for j in dict_orbital_2.get(i):
                orbital_list.append(j)
        elif i in dict_orbital_1.keys():
            orbital_list.append(i)
    print('You select orbitals:\t %s' %(orbital_list))
    return orbital_list

### READ DOSCAR ###
def read_doscar():
    f = open("DOSCAR", 'r')
    lines = f.readlines()
    f.close()
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    line_infor = lines[index].strip().split()
    if len(line_infor) == 5:
        E_max, E_min, nedos, e_fermi = [float(i) for i in line_infor[:-1]]
        nedos = int(nedos)
    elif len(line_infor) == 4:
        E_max = float(line_infor[0])
        E_min = float(line_infor[1][0:-5])
        nedos = int(line_infor[1][-5:])
        e_fermi = float(line_infor[2])
    ISPIN = None     
    index = 8
    if len(lines[index].strip().split()) == 5:
        ISPIN = 2
    else: 
        ISPIN = 1
    information = [natoms, E_max, E_min, nedos, e_fermi, ISPIN]    
    return lines, information

def write_dos0(lines, information):
    natoms, E_max, E_min, nedos, e_fermi, ISPIN = information[:]
    fdos = open("DOS0", 'w')
    if ISPIN == 1:
        fdos.write('%15s,%15s,%15s \n'  %('Energy', 'DOS', 'Integrated DOS'))
    elif ISPIN == 2:
        fdos.write('%10s,%10s,%10s,%10s,%10s \n'  %('Energy', 'DOS_up', 'DOS_down', 'Integrated DOS_up', 'Integrated DOS_down' ))
    
    for n in range(nedos):
        index = n + 6
        line_index = lines[index].strip().split()
        e = float(line_index[0]) - e_fermi
        fdos.write(('%10.5f,' + ','.join(line_index[1:]) + '\n') %(e))
    fdos.close()

def get_energy_list(lines, information):
    natoms, E_max, E_min, nedos, e_fermi, ISPIN = information[:]
    energy_list = []
    for n in range(nedos):
        index = n + 6
        line_index = lines[index].strip().split()
        e = float(line_index[0]) - e_fermi    
        energy_list.append(e)
    energy = np.array(energy_list)
    return energy    

###  Get sprcific atoms dos information ###
def get_single_atom(num, lines, information): 
    natoms, E_max, E_min, nedos, e_fermi, ISPIN = information[:]
    atom_dos = [] 
    index = int(num) * (nedos + 1 )  + 5
     
    for n in range(nedos):
        index +=1 
        dos = [float(i) for i in lines[index].strip().split()]
        dos[0] = dos[0] - e_fermi 
        atom_dos.append(dos)
    atom_dos = np.array(atom_dos)    
    return atom_dos


def get_single_orbital(num, orbital, lines, information):
    dos_atom = get_single_atom(num, lines, information)
#    print(dos_atom)
    dict_orbital_1 = get_orbital_dict()[0]
    n_orbi = dos_atom.shape[-1] # n_orbi is the number of orbitals +1, this can be used to check ISPIN 
        
    spin_1_list = [2,5,10,17]
    spin_2_list = [3,9,19,33]
    dos_position = dict_orbital_1.get(orbital)
    
   # energy = dos_atom[:,0]
    if n_orbi in spin_2_list:   
        dos_up = dos_atom[:, (dos_position * 2 - 1)]
        dos_down = dos_atom[:,(dos_position * 2)]
        atom_orbital = np.array([dos_up, dos_down])
    elif n_orbi in spin_1_list: 
        dos_position = dict_orbital_1.get(orbital)
        energy = dos_atom[:,0]  # Check in the future
        dos = dos_atom[:,dos_position]
        atom_orbital = np.array([dos])
    return atom_orbital 

