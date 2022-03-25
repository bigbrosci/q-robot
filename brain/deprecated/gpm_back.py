#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Python modules for operations using Gaussian and VASP
GPM stands for Gaussian Python Modules
'''
from __future__ import print_function
from lattice import *
import numpy as np 
import math
from itertools import combinations
import itertools
import pmutt
from pmutt.statmech.rot import RigidRotor
from pmutt.statmech.trans import FreeTrans
from pmutt.statmech.vib import HarmonicVib
from pmutt.statmech import elec
import ase 
from ase.visualize import view
import pmutt.io.gaussian as gr
from pmutt.statmech import StatMech, presets
import numpy as np


from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os, sys
from rdkit.Chem import Lipinski,rdMolDescriptors
from smiles2xyz import *


#### Part-1 Analyze the input files and generate the coords
def read_file(file_read):
    '''Return the lines and name of the file '''
    file_in = open(file_read, 'r')
    lines = file_in.readlines()
    file_in.close()
    name = file_in.name.split('.')[0]
    return name, lines


def get_coord_gjf(lines):   
    '''Get Coords from gjf and com file ''' 
    ### Get the start line number of coordinates part
    
    def calc_items(line):
        ''' See how many items in one line'''
        return len(line.strip().split())

    for num_line, line in enumerate(lines):
        if calc_items(line) == 0:
            if calc_items(lines[num_line+2]) == 0:
                if calc_items(lines[num_line+3]) == 2:
                    coord_num = num_line + 4 
                    break
    ### Get the end line number of coordinates part  
    for num_line, line in enumerate(lines[coord_num:len(lines)]):
        if len(line.strip().split()) == 0:
            coords = lines[coord_num:num_line+coord_num]
            break
    return coords

def get_coord_xyz(lines):
    '''Get Coords from xyz file ''' 
    line_start = 2
    for num, line in enumerate(lines[2:]):
        if len(line.strip()) == 0:
            line_end = num + line_start
            break
        else:
            line_end = len(lines)

    coord = lines[line_start:line_end]
    return coord

def get_coord_log(lines):
    '''Get Coords from log file '''
    look_up_start = 'Coordinates (Angstroms)'
    look_up_end   = '---------------------------------------------------------------------'
    
    start_l = []
    end_l = []
    for num, line in enumerate(lines):
        if look_up_start in line:
            start_l.append(num)
        if look_up_end in line:
            end_l.append(num)
          
    start_num = start_l[-1] + 3 
    end_num = end_l[end_l.index(start_l[-1]+2)+1]
    
    coords = []
    for line in lines[start_num:end_num]:
        infor = line.strip().split()
        name = dict_element.get(infor[1])
        xyz  = '\t'.join(infor[-3:])
        coord = name + '\t' + xyz + '\n'
        coords.append(coord)
    return coords

def file_analyzer(file_in):
    '''Get Coords from any input files '''
    
    def read_file(file_read):  
        '''Return the lines and name of the file '''
        file_in = open(file_read, 'r')
        lines = file_in.readlines()
        file_in.close()
        name = file_in.name.split('.')[0]
        return name, lines
    
    suffix = file_in.split('.')[-1]
    lines = read_file(file_in)[1]
    if suffix in ['com', 'gjf']:
        coords = get_coord_gjf(lines)
    elif suffix == 'xyz':
        coords = get_coord_xyz(lines)
    elif suffix == 'log':
        coords = get_coord_log(lines)
    return coords

##### Part-2 Analyze the coords and Get inforamtions
def xyz_analyzer(coords):
    '''Note: there will be two dictionaries:
    1) one is based on the atoms: {A1:Coord_A1, A2: Coord_A2, ...}
    2) and the other one is on element: {ele_1: [coord_1, coord_2, ...], ele_2: [coord_1, coord_2 ..] ..}
    '''    
    atom_sequence_list  = [] ## ['C1', 'C2', 'H3']
    atom_coord_list     = []  ## contains all the coordinates of each atoms in the gjf file
    ele_list            = []    ## contains all the elements in the gjf file, no duplicated ones
    dict_ele_coord      = {}  ## convert the coordinates to a dictionary based on element 
    dict_atom_coord     = {} ## dictionary based on atoms
    for num, line in enumerate(coords):
        infor       = line.strip().split()
        ele         = infor[0]
        atom_coord   = np.array([float(i) for i in infor[1:]])
        atom_sequence_list.append(ele+'-'+str(num))
        atom_coord_list.append(atom_coord)
        if ele not in ele_list:
            ele_list.append(ele)
            dict_ele_coord.update({ele:[atom_coord]})
        else:
            dict_ele_coord.get(ele).append(atom_coord)
    dict_atom_coord = dict(zip(atom_sequence_list, atom_coord_list))
    return dict_ele_coord, dict_atom_coord

def get_ele_list(coords, ele):
    '''Get all the atom sequence numbers of one element and the coordinates '''
    '''The sequence number should +1 when use get_coord_atom functions '''
    ele_list = []
    coord_list = []
    total_ele  = []
    for num, line in enumerate(coords):
        infor = line.strip().strip().split()
        ele_line = infor[0]
        total_ele.append(ele_line)
        ele_coord = np.array([float(i) for i in infor[1:]])
        if ele_line == ele:
            ele_list.append(num)
            coord_list.append(ele_coord)
    if ele not in total_ele:
        print('Not found %s' %(ele))
    else:
        dict_ele = dict(zip(ele_list, coord_list))    
        return ele_list, coord_list, dict_ele    
            
def get_dis(v1, v2):
    '''Get thre distance between two points in cartesian coordinates '''
    distance = np.linalg.norm(v1 - v2)
    return distance

def get_angle(u,v):
    ''' Get the angle between two vectors'''
    dot_cross = np.dot(u, v)
    c = dot_cross / np.linalg.norm(u) / np.linalg.norm(v)
    angle = np.arccos(np.clip(c, -1, 1))
    return angle

def get_dihedral(p):
    '''From https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python '''
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def get_normal_vector(p1, p2, p3):
    ''' Get the normal vector that is perpendicular to the plane'''
    ## These two vectors are in the plane
    v1 = p3 - p1
    v2 = p2 - p1
    ## the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    #a, b, c = cp
    #
    ## This evaluates a * x3 + b * y3 + c * z3 which equals d
    #d = np.dot(cp, p3)
    #
    #print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return cp

def get_atom_list_g09(coords, atom_s):
    '''     Get the atom number list such as [10, 11, 12,13] '''
    total_atoms = len(coords)
    atom_list = []
    def atom_list_append(num):
        if num not in atom_list:
            atom_list.append(num)
            
    for i in atom_s:
        if i.isdigit():
            atom_list_append(int(i))
        else:
            if '-' in i:
                n_start = int(i.split('-')[0])
                if i.split('-')[1].isdigit():
                    n_end = int(i.split('-')[1])
                else: # for the case: 18-, from 18 to the end.
                    n_end = total_atoms
                for m in range(n_start, n_end + 1):
                    atom_list_append(m)
            else:
                for num, line in enumerate(coords):
                    if line.strip().split()[0] == i:
                        atom_list_append(num + 1)                   
    return atom_list

def get_coord_atom(coords, atom_num):
    '''Get the xyz coordinates of the atom via the sequence number '''
    infor = coords[atom_num - 1].strip().split()
    atom_coord  = np.array([float(i) for i in infor[1:]])
    return atom_coord    

def del_one_element(coords,ele):
    '''Get the sequence number of ele in the coords and return the new coords without this element '''
    ele_list = []
    coord_new = coords[:]
    for num, line in enumerate(coords):
        infor =line.strip().split()
        if infor[0] == ele:
            ele_list.append(num)
            coord_new_remove(line)
    return ele_list, coord_new


def shift_coords(coords, vector):
    '''Shift all atoms with the same vector '''
    new_coords  = []
    for line in coords:
        infor = line.strip().split()
        ele_coord = np.array([float(i) for i in infor[1:]]) 
        new_coord = ele_coord + vector
        new_coords.append('%s %s\n' %(infor[0], ' '.join([str(round(i, 5)) for i in new_coord])))
    return new_coords

def shift_coords_all(coords, vectors):
    '''Shift all atoms with differnt vectors '''
    new_coords = []
    for num, line in enumerate(coords):
        infor = line.strip().split() 
        ele_coord = np.array([float(i) for i in infor[1:]])
        new_coord = ele_coord + vectors[num] 
        new_coords.append('%s %s\n' %(infor[0], '  '.join([str(round(i, 5)) for i in new_coord])))
    return new_coords
 
### Rotation
def rotate_vector(axis, theta, v):
    """ 
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d 
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d 
    eulerangles = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    new_v = np.dot(eulerangles, v)
    return new_v
    
def rotate_one_atom(coords, coord_1, coord_2, coord_3, theta):
    if np.array_equal(coord_1,coord_2):
        ''' Default axis along the z'''
        axis = np.array([0, 0, 1])    
    else:
        axis  = coord_2  - coord_1 
        
    vector_3 = coord_3 - coord_1    
    angle_213= get_angle(vector_3, axis)
    d_13 = get_dis(coord_1, coord_3)
    d_12 = get_dis(coord_1, coord_2)
    
    if np.array_equal(coord_3,coord_1) or np.array_equal(coord_3,coord_2):
        v_rotated = coord_3
    else:  
    # get the coordination of the perpendicular point 
        if np.array_equal(coord_1, coord_2): 
            point_perpen = coord_1
        else: 
            point_perpen = coord_1 + d_13 * math.cos(angle_213) / d_12 * axis
    ## shifted coordination for the atom we want to rotate 
        v = coord_3 - point_perpen
    # rotation_part
        rotated_v = rotate_vector(axis, theta, v)
    # shift the rotated coordinate back by the vector from  perpendicular point 
        v_rotated = rotated_v + point_perpen
    return v_rotated

def rotate_atoms(coords, coord_1, coord_2, atom_list, theta):   
    '''Along the axis(atom_A, atom_B), rotate atoms by theta angle (in degree).
Generall Command: rotate.py atom_A atom_B atom_1 atom_2 ... Angle    
For example: To rotate the C H O atoms 15 degree along the axis formed by No.2 and No.10 atoms.
rotate.py 2 10 C H O 15  
This function is used to analyze the arguments: 2 10 C H O 15  in the command above.
And return the basic information for rotations in the next step.
'''
    for i in atom_list:
        ele_atom = coords[i-1].strip().split()[0]
        P3 = get_coord_atom(coords, i)
        new_coord_i = rotate_one_atom(coords, coord_1, coord_2, P3, theta)
        coords[i-1] = ele_atom + ' ' + '  '.join([str(round(j,5)) for j in new_coord_i]) + '\n'
    return coords

###  Part3 Save files
def save_xyz(coords, name):
    out_name = name + '.xyz'
    file_out = open(out_name, 'w')
    file_out.write(str(len(coords)) + '\n')
    file_out.write(name + '\n')
    file_out.writelines(coords)
    file_out.close()

def save_xyz_by_ele(dict_ele_coord):
    #anchor = dict_coord.get('H')[-1]  ### The atom you want to put at the origin point (0,0,0)
    file_name = 'new_geo.xyz'
    atom_num = 0
    for key, values in dict_ele_coord.items():
        atom_num += len(values) 
    file_out = open(file_name, 'w')
    file_out.write(str(atom_num) + '\n')
    file_out.write(file_name + '\n')
    for key, value in dict_ele_coord.items():
        if len(value) > 1: 
            for atom_coord in value:
                file_out.write('%s %s\n' %(key, ' '.join([str(round(j,5)) for j in atom_coord])))
#                new_coord = one_coord - anchor  # For shift atoms
#                file_out.write('%s %s\n' %(key, ' '.join([str(j) for j in new_coord])))
        else: 
            file_out.write('%s %s\n' %(key, ' '.join([str(round(j,5)) for j in value[0]])))
#            new_coord = value[0] - anchor        
#            file_out.write('%s %s\n' %(key, ' '.join([str(round(j,4)) for j in new_coord]))) 
    file_out.close()

def save_xyz_by_atom(dict_atom_coord):
    out_name = 'new_geo.xyz'
    atom_num = len(dict_atom_coord.keys())
    file_out = open(out_name, 'w')
    file_out.write(str(atom_num) + '\n')
    file_out.write(out_name + '\n')
    for key, value in dict_atom_coord.items():
        file_out.write('%s %s\n' %(key.split('-')[0], ' '.join([str(round(j,5)) for j in value])))
#        new_coord = value - anchor        
 #       file_out.write('%s %s\n' %(key.split('-')[0], ' '.join([str(j) for j in new_coord]))) 
    file_out.close()
    
### Update the gif file 
def save_gjf(coords, name):
    mem        = '10GB'
    functional = 'm062x'
    basis_low      = '6-31g'
    basis_high     = '6-311++g(d,p)'
    task1      = 'opt' 
    task2      = 'freq'
    charge     =  0
    multiplicity = 1
    ### Overwrite the gjf file 
    chk_name = name.split('/')[-1]
    out_name = name + '.gjf'
    file_out = open(out_name, 'w')
    file_out.write('%schk=%s\n' %('%', chk_name))
    file_out.write('%smem=%s\n' %('%', mem))
    file_out.write('# %s/%s %s\n' %(functional, basis_low, task1))
    file_out.write(' \n%s \n \n'     %(chk_name))
    file_out.write('%s %s\n' %(charge, multiplicity))
    file_out.writelines(coords)
    file_out.write('\n--Link1--\n')
    file_out.write('%schk=%s\n' %('%', chk_name))
    file_out.write('%smem=%s\n' %('%', mem))
    file_out.write('#p %s/%s Guess=Read Geom=Checkpoint %s %s\n' %(functional, basis_high, task1, task2))
    file_out.write(' \n%s \n \n'     %(chk_name))
    file_out.write('%s %s\n \n' %(charge, multiplicity)) 
    file_out.close()

#### To POSCAR
def get_mass_center(coords):
    '''Get the mass center of the molecule '''
    sum_mass_coord = np.array([0.0,0.0,0.0]) 
    sum_mass       = 0.0
    for line in coords:
        infor = line.strip().split()
        atom_ele   = infor[0]
        atom_coord = np.array([float(i) for i in infor[1:]])
        atom_mass  = float(atomic_mass.get(atom_ele))
        sum_mass += atom_mass
        sum_mass_coord += atom_mass * atom_coord 
    return sum_mass_coord / sum_mass

def box_center_coords(coords):
    '''Shift the molecule by using the mass center, this will do a lot of help when 
    convert the xyz to POSCAR to avoid the PBC problem. 
    1) the default box size for VASP calculations is: 16 * 16 * 16 \AA^3
    '''
    mass_center = get_mass_center(coords)
    vector      = np.array([8.0, 8.0, 8.0]) - mass_center
    new_coords  = shift_coords(coords, vector)
    return new_coords

def save_poscar(coords, name):
    '''Save the coords to POSCAR file 
    1) the previous atom sequence information in xyz or gjf file  will be lost since VASP sort the atoms via element.
    '''
    new_coords = box_center_coords(coords)
    dict_ele_coord = xyz_analyzer(new_coords)[0]
    print(dict_ele_coord.keys())
    out_name = 'POSCAR_' + name
    file_out = open(out_name, 'w')
    file_out.write(name + '\n' + '1.0 \n')
    file_out.write('%s 0.0 0.0\n0.0 %s 0.0 \n0.0 0.0 %s\n' %('16', '16', '16'))
    file_out.write(' '.join(dict_ele_coord.keys()) + '\n')
    file_out.write(' '.join([str(len(i)) for i in dict_ele_coord.values()]) + '\n')
    file_out.write('Selective dynamics\nCartesian\n')
    for i in dict_ele_coord.values():
        for j in i:
            file_out.write(' '.join([str(k) for k in j]) + ' T  T  T\n')
    file_out.close()
 
def save_poscar_from_g09(f_in):
    coords = file_analyzer(f_in)
    dict_ele_coord = xyz_analyzer(coords)[0]
    out_name = 'POSCAR_' + f_in.split('.')[0]
    file_out = open(out_name, 'w')
    file_out.write(f_in + '\n' + '1.0 \n')
    
    dire_a = dict_ele_coord.get('Tv')
    
    #dire_a = np.array([[16,0,0],[0, 16,0], [0,0,16]])
    
    lattice_a = ' '.join([str(i) for i in dire_a[0]]) + '\n'
    lattice_b = ' '.join([str(i) for i in dire_a[1]]) + '\n'
    lattice_c = ' '.join([str(i) for i in dire_a[2]]) + '\n'
    
    file_out.write('%s %s %s' %(lattice_a, lattice_b, lattice_c))
    
    del dict_ele_coord['Tv']
    
    file_out.write(' '.join(dict_ele_coord.keys()) + '\n')
    file_out.write(' '.join([str(len(i)) for i in dict_ele_coord.values()]) + '\n')
    file_out.write('Selective dynamics\nCartesian\n')
    for i in dict_ele_coord.values():
        for j in i:
            file_out.write(' '.join([str(k) for k in j]) + ' T  T  T\n')
    file_out.close()
    
def find_elements(coords, anchor, ele, threshold):
    list_g = []
    coord_anchor = get_coord_atom(coords, anchor)
    for num, line in enumerate(coords):
        infor = line.strip().split()
        coord_i = np.array([float(i) for i in infor[1:]])
        if infor[0] == ele:
            dis_i = get_dis(coord_anchor, coord_i)
            if  0.1 < dis_i <= threshold:
                list_g.append(num)
    return list_g
            
#hatree_to_ev   = 27.2116     # eV
#hatree_to_kcal = 627.509   # kcal/mol
#hatree_to_kj   = 2625.50  # kJ/mol
#cal_to_j       = 4.184 

###Get information from Log file
def hatree_to_ev(energy):
    '''Convert energies in hartree to other units '''
    hatree_to_ev   = 27.2116     # eV
    return hatree_to_ev * energy

def hatree_to_kcal(energy):
    '''Convert energies in hartree to other units '''
    hatree_to_kcal = 627.509   # kcal/mol
    return hatree_to_kcal * energy    

def hatree_to_kj(energy):
    '''Convert energies in hartree to other units '''
    hatree_to_kj   = 2625.50  # kJ/mol 
    return hatree_to_kj * energy    
    
def get_infor_from_log(log_file):    
    '''Get some useful information from G09 log file ''' 
    file_in = open(log_file, 'r')
    lines = file_in.readlines()
    file_in.close()     
        
    dict_search = {
    'E_state':        'The electronic state',
    'freq'   :        'Frequencies --',            
    'Stoichiometry':  'Stoichiometry',
    'spin'   :        'Multiplicity',
    'zpe_c'  :        'Zero-point correction=',                      
    'etot_c' :        'Thermal correction to Energy=',              
    'enthalpy_c'   :  'Thermal correction to Enthalpy=',             
    'gibbis_c'     :  'Thermal correction to Gibbs Free Energy=',   
    'Ezpe'   :        'Sum of electronic and zero-point Energies=',
    'Etot'   :        'Sum of electronic and thermal Energies=',
    'Enthalpy'     :  'Sum of electronic and thermal Enthalpies=',
    'Gibbis'       :  'Sum of electronic and thermal Free Energies=' ,
    }
    
    dict_log  = {}
    
    freq_line = 0
    e_dft = 0
    Stoichiometry = None 
    E_state       = None
    frequencies = []
    
    for line_num, line in enumerate(lines): 
        infor = line.rstrip().split()
        if dict_search.get('zpe_c') in line: 
            freq_line = line_num
        elif 'SCF Done' in line:
            try:
                e_dft = float(infor[4])
            except ValueError:
                e_dft = 0
        elif dict_search.get('Stoichiometry') in line:
            Stoichiometry = infor[-1].split('(')[0]
        elif dict_search.get('E_state') in line:
            E_state = infor[-1]
        elif dict_search.get('spin') in line:
            spin = (int(infor[-1]) - 1 ) / 2  ## 2s+1 = Multiplicity, in pmutt, the S is used
        elif  dict_search.get('freq') in line:
            freqs = [float(i) for i in line.strip().split('--')[1].split()]
            frequencies.extend(freqs) 
    
    dict_log.update({'E_state' :  E_state})
    dict_log.update({'e_dft' :  e_dft})
    dict_log.update({'Stoichiometry'  :  Stoichiometry})
    dict_log.update({'spin': spin})
    
    if freq_line > 0:   # If ! > 0 means the job is not a frequency calculation     
        dict_log.update({'freq' :   frequencies})  
        if min(frequencies) < 0:
            dict_log.update({'Stability': 'NoneStable'})
        else:
            dict_log.update({'Stability': 'Stable'})
        
        item_list =  ['zpe_c', 'etot_c', 'enthalpy_c', 'gibbis_c', 'Ezpe', 'Etot', 'Enthalpy', 'Gibbis'] # Thermo infor    
        for num, line in enumerate(lines[freq_line:freq_line+8]):
            E_line = float(line.strip().split('=')[1].split()[0])
            dict_log.update({item_list[num]:E_line})      
            
    return dict_log 

def freq_correct(file_in, scale):
    '''Correct the xyz coordinate based on the vibration mode of imaginary frequencies '''
    look_up = 'Frequencies --'
    frequencies = []
    dict_freq = {}
    with open(file_in, 'r') as f_in:
        lines = f_in.readlines()
        coords = get_coord_log(lines)
        for num, line in enumerate(lines):
            if look_up in line:
                 freqs = [float(i) for i in line.strip().split('--')[1].split()]
                 dict_freq.update({num:freqs})
                 frequencies.extend(freqs)
    
    l_f  = min(frequencies)  # Get the smallest imaginary frequncy
    n_f = list(dict_freq.keys())
    
    for key, val in dict_freq.items():
        if l_f in val:
            f_p = val.index(l_f)
            l_s = key + 5
            l_e = n_f[n_f.index(key)+1] - 2
    
    coords_correct = []
    for line in lines[l_s:l_e]:
        infor = line.strip().split()
        coords_correct.append([float(i) for i in infor[f_p*3+2:f_p*3+5]])
    
    vectors = scale * np.array(coords_correct)
    new_coords = shift_coords_all(coords, vectors)
    
    return new_coords


def get_HG_dict(log_file, spin=0):
    
#    ''' 0.970 is the scale factor for M062x/6-311++G(d,p) vibration frequncies. ''' 
#    converter = 27.211386246  ## From hartree to eV
#    mol = ase.io.read(log_file.replace('log', 'xyz'))
#    
#    zpe = gr.read_zpe(log_file)   # zero point energy in hartree 
#    E_zpe = gr.read_electronic_and_zpe(log_file)  # ZPE corrected energy in hartree
#    E_ele = E_zpe - zpe   # Energy without ZPE correction in hartree
#    
#    E_ele_pmutt = E_ele * converter  # Energy without any ZPE corrections in eV 
#    E_zpe_pmutt = (E_ele + zpe * 0.970 )  * converter  # Energy with scaled ZPE correction in eV
#    
#    Freq   = [ float(i) * 0.970 for i in gr.read_frequencies(log_file)]  
#    
#    Mass   = gr.read_molecular_mass(log_file)
#    rot_symm = gr.read_rot_symmetry_num(log_file)
#    rot_tem  = gr.read_rotational_temperatures(log_file)
#
#    mol_trans = FreeTrans(n_degrees=3,atoms=mol)
#    mol_rot = RigidRotor(symmetrynumber=rot_symm, atoms=mol)
#    mol_vibs = HarmonicVib(Freq)
#    mol_elec = elec.GroundStateElec(potentialenergy=E_ele_pmutt, spin=spin)
#    mol_statmech = StatMech(trans_model=mol_trans, 
#                            vib_model=mol_vibs,
#                            rot_model=mol_rot,
#                            elec_model=mol_elec)
#    thermo_dict = {}

    for Tem in range(50,1501,10):
        Cp_pmutt = mol_statmech.get_Cp('eV/K', T=Tem)
        Cv_pmutt = mol_statmech.get_Cv('eV/K', T=Tem)
        S_pmutt = mol_statmech.get_S('eV/K', T= Tem)
        H_pmutt = mol_statmech.get_H('eV', T= Tem)
        G_pmutt = mol_statmech.get_G('eV', T= Tem)
        thermo_dict.update({Tem:[E_ele_pmutt, E_zpe_pmutt, H_pmutt, G_pmutt, Cp_pmutt, Cv_pmutt, S_pmutt]}) 

    return thermo_dict

def get_HG_dict_T(log_file, T, spin = 0):
    ''' 0.970 is the scale factor for M062x/6-311++G(d,p) vibration frequncies. ''' 
    converter = 27.211386246  ## From hartree to eV
    mol = ase.io.read(log_file.replace('log', 'xyz'))
    
    zpe = gr.read_zpe(log_file)   # zero point energy in hartree 
    E_zpe = gr.read_electronic_and_zpe(log_file)  # ZPE corrected energy in hartree
    E_ele = E_zpe - zpe   # Energy without ZPE correction in hartree
    
    E_ele_pmutt = E_ele * converter  # Energy without any ZPE corrections in eV 
    E_zpe_pmutt = (E_ele + zpe * 0.970 )  * converter  # Energy with scaled ZPE correction in eV
    
    Freq   = [ float(i) * 0.970 for i in gr.read_frequencies(log_file)]  
    
    Mass   = gr.read_molecular_mass(log_file)
    rot_symm = gr.read_rot_symmetry_num(log_file)
    rot_tem  = gr.read_rotational_temperatures(log_file)

    mol_trans = FreeTrans(n_degrees=3,atoms=mol)
    mol_rot = RigidRotor(symmetrynumber=rot_symm, atoms=mol)
    mol_vibs = HarmonicVib(Freq)
    mol_elec = elec.GroundStateElec(potentialenergy=E_ele_pmutt, spin=spin)
    mol_statmech = StatMech(trans_model=mol_trans, 
                            vib_model=mol_vibs,
                            rot_model=mol_rot,
                            elec_model=mol_elec)
    
    Tem = float(T)
    Cp_pmutt = mol_statmech.get_Cp('eV/K', T=Tem)
    Cv_pmutt = mol_statmech.get_Cv('eV/K', T=Tem)
    S_pmutt = mol_statmech.get_S('eV/K', T= Tem)
    H_pmutt = mol_statmech.get_H('eV', T= Tem)
    G_pmutt = mol_statmech.get_G('eV', T= Tem)
    
    
    key_list = ['E_elec','E_zpe', 'H', 'G', 'Cp', 'Cv', 'S']
    value_list = [E_ele_pmutt, E_zpe_pmutt, H_pmutt, G_pmutt, Cp_pmutt, Cv_pmutt, S_pmutt]
    
    thermo_dict = dict(zip(key_list, value_list))
    
#    thermo_dict = {}
#    thermo_dict.update({'E_elec':E_ele_pmutt})
#    thermo_dict.update({'E_zpe':E_zpe_pmutt})
#    thermo_dict.update({'H':H_pmutt})
#    thermo_dict.update({'G':G_pmutt})
#    thermo_dict.update({'Cp':Cp_pmutt})
#    thermo_dict.update({'S':S_pmutt})
#    thermo_dict.update({Tem:[E_ele_pmutt, E_zpe_pmutt, H_pmutt, G_pmutt, Cp_pmutt, Cv_pmutt, S_pmutt]}) 
    
    return thermo_dict
    

def get_HG_dict_T_G4(log_file, T, spin = 0):
    
    converter = 27.211386246   
    ''' 0.970 is the scale factor for M062x/6-311++G(d,p) vibration frequncies. ''' 
    mol = ase.io.read(log_file.replace('log', 'xyz'))
    

    f = open(log_file, 'r')
    lines = f.readlines()
    f.close()
    
    list_energy = []
    for line in lines:
        infor = line.strip().split()
        if 'E(ZPE)' in line:
            list_energy.append(infor[1])
        if 'G4(0 K)=' in line:
            list_energy.append(infor[2])
        if 'G4 Enthalpy'  in line:
            list_energy.append(infor [2])
            list_energy.append(infor[-1])
    #        print(ZPE, E_zpe, Enthalpy, Gibbs) 
    
    zpe = float(list_energy[0])   # zero point energy in hartree 
    E_zpe = float(list_energy[1]) # ZPE corrected energy in hartree
    E_ele = E_zpe - zpe   # Energy without ZPE correction in hartree
    
    E_ele_pmutt = E_ele * converter  # Energy without any ZPE corrections in eV 
    E_zpe_pmutt = E_zpe  * converter  # Energy with scaled ZPE correction in eV
    
    Freq   = [ float(i) * 0.965 for i in gr.read_frequencies(log_file)]  
    
    Mass   = gr.read_molecular_mass(log_file)
    rot_symm = gr.read_rot_symmetry_num(log_file)
    rot_tem  = gr.read_rotational_temperatures(log_file)

    mol_trans = FreeTrans(n_degrees=3,atoms=mol)
    mol_rot = RigidRotor(symmetrynumber=rot_symm, atoms=mol)
    mol_vibs = HarmonicVib(Freq)
    mol_elec = elec.GroundStateElec(potentialenergy=E_ele_pmutt, spin=spin)
    mol_statmech = StatMech(trans_model=mol_trans, 
                            vib_model=mol_vibs,
                            rot_model=mol_rot,
                            elec_model=mol_elec)
    
    Tem = float(T)
    Cp_pmutt = mol_statmech.get_Cp('eV/K', T=Tem)
    Cv_pmutt = mol_statmech.get_Cv('eV/K', T=Tem)
    S_pmutt = mol_statmech.get_S('eV/K', T= Tem)
    H_pmutt = mol_statmech.get_H('eV', T= Tem)
    G_pmutt = mol_statmech.get_G('eV', T= Tem)
    
    
    key_list = ['E_elec','E_zpe', 'H', 'G', 'Cp', 'Cv', 'S']
    value_list = [E_ele_pmutt, E_zpe_pmutt, H_pmutt, G_pmutt, Cp_pmutt, Cv_pmutt, S_pmutt]
    
    thermo_dict = dict(zip(key_list, value_list))
    
    return thermo_dict

def get_dict_infor(log_file):
    ''' Return lots of information in a dictionary format''' 
    f = open(log_file, 'r')
    lines = f.readlines()
    f.close()
    
    dict_infor = get_infor_from_log(log_file)
    name = log_file.split('.')[0]
    coords = get_coord_log(lines) 
    save_xyz(coords, name)
    xyz_file = name + '.xyz'
    smiles, formula  = get_smiles_m1(xyz_file)
    if '=' in smiles and ']' in smiles:
        smiles, formula = get_smiles_m2(xyz_file)
    
    dict_infor.update({'log_file':log_file})
    dict_infor.update({'name':name})
#    dict_infor.update({'coords':coords}) 
    dict_infor.update({'smiles':smiles})
    dict_infor.update({'formula':formula})
    
    dict_ele_coord = xyz_analyzer(coords)[0]
    for key, value in dict_ele_coord.items():
        '''Number of atoms for elements'''
        new_key = 'N_' + key   # Add N_ to avoid the mixing of the H (element) and H (Enthalpy)
        dict_infor.update({new_key:len(value)})           
    
    for line in lines:    
        infor_line = line.rstrip().split()
        if 'G4(0 K)=' in line:
            dict_infor['Ezpe'] = infor_line[2]
            dict_infor.update({'method':'G4'})
        if 'G4 Enthalpy'  in line:
            dict_infor['Enthalpy'] = infor_line[2]
            
    return  dict_infor

def add_thermo_to_dict_infor(dict_infor):  # T
    '''Calculate the formation of enthalpies at 0 K and 298.15 K ''' 
    T = 298.15
    hartree_to_eV = 27.211386246
    hartree_to_kcal = 627.5095
    
    
    dict_exp = {#'''Enthalpy formation energies at 0K and 298.15K, from Hf https://atct.anl.gov/Thermochemical%20Data/version%201.118/index.php''' 
        'C': (170.027485659656, 171.338432122371),
        'H': (51.6333652007648, 52.1027724665392),
        'O': (58.9971319311663, 59.5671606118547)
    }

    if dict_infor.get('method')  == 'G4' :
        dict_g09 = {   # g4 in hartree
                'C':(-37.834168, -37.831808),
                'H':(-0.50142,-0.49906),
                'O':(-75.045501,-75.043141)
                }
    else:
        dict_g09 = {   # M06-2X in hartree
            'C': (-37.840428,-37.838068),
            'H': (-0.498195, -0.495835),
            'O': (-75.06104, -75.05868)
            }
        
    def calc_delta_f_H_0(dict_hg_T):
        
        try: 
            N_C = int(dict_infor.get('N_C'))
        except:
            N_C = 0
            
        try: 
            N_O = int(dict_infor.get('N_O'))
        except:
            N_O = 0
            
        try: 
            N_H = int(dict_infor.get('N_H'))
        except:
            N_H = 0
            
       
        E_ZPE = dict_hg_T.get('E_zpe') / hartree_to_eV
        E_H = dict_hg_T.get('H')  / hartree_to_eV
        
        Delta_H_formation_0 =   ( N_C * dict_exp.get('C')[0] + N_O * dict_exp.get('O')[0] + N_H * dict_exp.get('H')[0] ) - hartree_to_kcal * (N_C * dict_g09.get('C')[0] + N_O * dict_g09.get('O')[0] +  N_H * dict_g09.get('H')[0] - E_ZPE) 
        Delta_H_formation_298 = ( N_C * dict_exp.get('C')[1] + N_O * dict_exp.get('O')[1] + N_H * dict_exp.get('H')[1] ) - hartree_to_kcal * (N_C * dict_g09.get('C')[1] + N_O * dict_g09.get('O')[1] +  N_H * dict_g09.get('H')[1] - E_H) 
        dict_formation = {}
        dict_formation.update({'Hf_0': Delta_H_formation_0}) 
        dict_formation.update({'Hf_T': Delta_H_formation_298})        
        
        return  dict_formation
    

    spin = dict_infor.get('spin')
    log_file = dict_infor.get('log_file')
    if dict_infor.get('method') == 'G4' :
        dict_hg_T = get_HG_dict_T_G4(log_file, 298.15, spin=spin)
    else:
        dict_hg_T = get_HG_dict_T(log_file, 298.15, spin=spin)
    dict_infor.update(dict_hg_T)
    
    dict_formation = calc_delta_f_H_0(dict_hg_T)
    dict_infor.update(dict_formation)
        
    return  dict_infor
    
#def get_HG_dict(log_file, spin=0):
#    converter = 27.2216 
#    mol = ase.io.read(log_file.replace('log', 'xyz'))
#    zpe = gr.read_zpe(log_file)
#    E_ele_pmutt = gr.read_electronic_and_zpe(log_file) - zpe
#    E_elec = E_ele_pmutt * converter  # From hartree to eV
#    Freq   = gr.read_frequencies(log_file)
#    Mass   = gr.read_molecular_mass(log_file)
#    rot_symm = gr.read_rot_symmetry_num(log_file)
#    rot_tem  = gr.read_rotational_temperatures(log_file)
#
#    mol_trans = FreeTrans(n_degrees=3,atoms=mol)
#    mol_rot = RigidRotor(symmetrynumber=rot_symm, atoms=mol)
#    mol_vibs = HarmonicVib(Freq)
#    mol_elec = elec.GroundStateElec(potentialenergy=E_elec, spin=spin)
#    mol_statmech = StatMech(trans_model=mol_trans, 
#                            vib_model=mol_vibs,
#                            rot_model=mol_rot,
#                            elec_model=mol_elec)
#    thermo_dict = {}
#    E_zpe_pmutt = gr.read_electronic_and_zpe(log_file)
#    E_ele_pmutt = gr.read_electronic_and_zpe(log_file) - zpe
#    for Tem in range(50,1501,10):
#        Cp_pmutt = mol_statmech.get_Cp('eV/K', T=Tem)
#        Cv_pmutt = mol_statmech.get_Cv('eV/K', T=Tem)
#        S_pmutt = mol_statmech.get_S('eV/K', T= Tem)
#        H_pmutt = mol_statmech.get_H('eV', T= Tem)
#        G_pmutt = mol_statmech.get_G('eV', T= Tem)
#        thermo_dict.update({Tem:[E_ele_pmutt, E_zpe_pmutt, H_pmutt, G_pmutt, Cp_pmutt, Cv_pmutt, S_pmutt]})
#    return thermo_dict

#def get_HG(log_file, Tem, spin=0):
#    converter = 27.2216 
#    zpe = gr.read_zpe(log_file)
#    E_elec = (gr.read_electronic_and_zpe(log_file) - zpe)
#    
#    E_elec = E_elec * converter  # From hartree to eV
#    Freq   = gr.read_frequencies(log_file)
#    Mass   = gr.read_molecular_mass(log_file)
#    rot_symm = gr.read_rot_symmetry_num(log_file)
#    rot_tem  = gr.read_rotational_temperatures(log_file)
#
#    mol_trans = FreeTrans(n_degrees=3,atoms=mol)
#    mol_rot = RigidRotor(symmetrynumber=rot_symm, atoms=mol)
#    mol_vibs = HarmonicVib(Freq)
#    mol_elec = elec.GroundStateElec(potentialenergy=E_elec, spin=spin)
#    mol_statmech = StatMech(trans_model=mol_trans, 
#                            vib_model=mol_vibs,
#                            rot_model=mol_rot,
#                            elec_model=mol_elec)
#    
#    S_pmutt = mol_statmech.get_S('eV/K', T= Tem)
#    H_pmutt = mol_statmech.get_H('eV', T= Tem)
#    G_pmutt = mol_statmech.get_G('eV', T= Tem)
##    return E_elec, zpe, H_pmutt/converter, G_pmutt/converter
#    return H_pmutt/converter, G_pmutt/converter
    
    #from pmutt.statmech import StatMech, presets
    #test = StatMech(name='test', symmetrynumber=1, atoms = mol, vib_wavenumbers=Freq, potentialenergy=E_elec, **presets['idealgas'])
    #print(test.get_H('kcal/mol', T=Tem))
    #print(help(pmutt.constants.R))
    
#    from ase.build import molecule
#    from ase.calculators.emt import EMT
#    from ase.optimize import QuasiNewton
#    from ase.vibrations import Vibrations
#    from ase.thermochemistry import IdealGasThermo
#    
#    atoms = mol
#    vib_energies = np.array(Freq)/8065.54429
#    thermo = IdealGasThermo(vib_energies=vib_energies,
#                            potentialenergy=E_elec,
#                            atoms=atoms,
#                            geometry='nonlinear',
#                            symmetrynumber=1, spin=0)
##    ZPE = thermo.get_ZPE_correction()
##    S_ase = thermo.get_entropy(temperature=Tem, pressure=101325.)
#    G_ase = thermo.get_gibbs_energy(temperature=Tem, pressure=101325.)
#    H_ase = thermo.get_enthalpy(temperature=Tem)
#    
#    print(H_ase/converter, G_ase/converter)
#    print(H_pmutt/converter, G_pmutt/converter)
  
##compare(1.0)  
#compare(27.2216)  


#### Analyze the xyz file to get the position information
def get_dis_all(coords):
    '''Calculate all the possible bonds in the structure''' 

    max_dis_cc_single = 1.60
    max_dis_cc_double = 1.45
    max_dis_co_single =  1.50    
    max_dis_co_double =  1.30    
    max_dis_ch_single = 1.20
    max_dis_oh_single = 1.10
    max_dis_oo_hbond  = 3.5
    max_dis_hbond     = 2.4
    
    def get_bond_pairs(coords, A, B, dis_AB_max, dis_AB_min=0):
        '''Get all possible A-B bond pairs'''
        A_list, A_coord_list, dict_ele_A = get_ele_list(coords, A)
        B_list, B_coord_list, dict_ele_B = get_ele_list(coords, B)
        AB_total = []

        if A == B:
            pairs_AB = combinations(A_list, 2) 
        else:    
            pairs_AB = list(itertools.product(A_list, B_list))
            
        for pair in pairs_AB:
            v_a = dict_ele_A.get(pair[0])
            v_b = dict_ele_B.get(pair[1])
            dis_AB = get_dis(v_a, v_b)
            if dis_AB_min  < dis_AB <= dis_AB_max :
                AB_total.append(list(pair))
        return AB_total
    
    cc_single = get_bond_pairs(coords, 'C', 'C', max_dis_cc_single, max_dis_cc_double) 
    cc_double = get_bond_pairs(coords, 'C', 'C', max_dis_cc_double) 
    co_single = get_bond_pairs(coords, 'C', 'O', max_dis_co_single, max_dis_co_double) 
    ch_single = get_bond_pairs(coords, 'C', 'H', max_dis_ch_single)
    oh_single = get_bond_pairs(coords, 'O', 'H', max_dis_oh_single)
    oo_hbond  = get_bond_pairs(coords, 'O', 'O', max_dis_oo_hbond)
    oh_long_bond    = get_bond_pairs(coords, 'O', 'H', max_dis_hbond, max_dis_oh_single)

    return cc_single, cc_double, co_single, ch_single, oh_single, oo_hbond, oh_long_bond
        
def analyze_rings(cc_double):
    '''This function works on the C=C bond list, and Count the number of rings and return the atoms in the ring '''
    ring_list = []
    if len(cc_double)%6 == 0 :  
        n_ring = int(len(cc_double)/6)
#        print('%s rings found' %(n_ring))
    pairs = cc_double[:] 
    #pairs = cc_double.copy() # Python2 does not support this
    for i in range(0, n_ring):        
        l_ring_i = []
        l_ring_i.extend(pairs[0])
        for num, pair in enumerate(pairs): 
            if any(ele in pair for ele in l_ring_i):
                l_ring_i.extend(pair)
                l_ring_i = list(set(l_ring_i))
        ring_list.append(l_ring_i)
        
        for i in l_ring_i:
            for pair in pairs:
                if i in pair:
                    pairs.remove(pair)
    return ring_list

def analyze_c1(cc_single, cc_double):
    '''Obtain all the C in Ar which connects to the aliphatic Chain '''
    ring_list = analyze_rings(cc_double)
    
    single_c = [i for pair in cc_single for i in pair]
    ring_c   = [i for pair in ring_list for i in pair]
    
    list_c1 = []
    for i in single_c: 
        if single_c.count(i) == 1:  # C1 is not form two single bonds
            if i in ring_c:       # C1 also in ring
               list_c1.append(i)
    return list_c1           
    
def analyze_chain(cc_single, cc_double):
    ring_list = analyze_rings(cc_double)
    ring_c   = [i for pair in ring_list for i in pair]
    list_c1 = analyze_c1(cc_single, cc_double)
        
    list_c_alpha = []
    
    for c1 in list_c1: 
        for pair in cc_single:
            if c1 in pair:
                c_alpha = [i for i in pair if i is not c1][0]
                if c_alpha not in ring_c:
                    list_c_alpha.append([c1, c_alpha])
                    cc_single.remove(pair)
    
    list_c_beta = []
    for one_c_alpha in list_c_alpha:
        c_alpha = one_c_alpha[1]
        for pair in cc_single:
            if c_alpha in pair:
                c_beta = [i for i in pair if i is not c_alpha][0]
                if c_beta not in ring_c:
                    one_c_beta = one_c_alpha[:]
                    one_c_beta.append(c_beta)
                    list_c_beta.append(one_c_beta)
                    cc_single.remove(pair)
                
    list_c_gamma = []
    for one_c_beta in list_c_beta:
        c_beta = one_c_beta[2]
        for pair in cc_single:
            if c_beta in pair:
                c_gamma = [i for i in pair if i is not c_beta][0]
                if c_gamma not in ring_c:
                    one_c_gamma = one_c_beta[:]
                    one_c_gamma.append(c_gamma)
                    list_c_gamma.append(one_c_gamma)
                    cc_single.remove(pair)
                
#    print(list_c1)            
#    print(list_c_alpha)
#    print(list_c_beta)        
#    print(list_c_gamma)        
    return  list_c1, list_c_alpha, list_c_beta, list_c_gamma           

def analyze_O_chain(cc_single, cc_double, co_single):
    ring_list = analyze_rings(cc_double)
    ring_c   =  [i for pair in ring_list for i in pair]
    
    list_ar = []  ### the c in the chain connects to the O
    list_aro = []   ### the o connects to the Ar
    list_aroc = []
    list_aroar = []
    for c_ring in ring_c:
        for pair in co_single:
            if c_ring in pair: 
                list_ar.append(c_ring)
                list_aro.append(pair)
                co_single.remove(pair)

    for one_aro in list_aro:
        o =  one_aro[1]               
        for pair in co_single:
            if o in pair:
                one_aroc = one_aro[:]
                one_aroc.append(pair[0])
                if pair[0] in ring_c: 
                    list_aroar.append(one_aroc)
                else:                    
                    list_aroc.append(one_aroc)
    print(list_ar, list_aro, list_aroc, list_aroar)                    
    return list_ar, list_aro, list_aroc, list_aroar
        
    
def analyze_terminal(ch_single, co_single):
    '''Obtain all terminal groups such as CH3 and CH2OH'''
    list_ch3 = []
    list_ch2 = []
    ch_list = [i for pair in ch_single for i in pair]
#    co_list = [i for pair in co_single for i in pair]
    for ele in ch_list:
        if ch_list.count(ele) == 3:
            ele_ch3 = []
            ele_ch3.append(ele)
            ele_ch3.extend(pair[1] for pair in ch_single if ele in pair)
            if ele_ch3 not in list_ch3:
                list_ch3.append(ele_ch3)
        elif ch_list.count(ele) == 2:
            ele_ch2 = []
            ele_ch2.append(ele)
            ele_ch3.extend(pair[1] for pair in ch_single if ele in pair)
            if ele_ch2 not in list_ch2:
                list_ch2.append(ele_ch2)

    list_ch2oh = []
    for one_ch2 in list_ch2:
        c = one_ch2[0]
        for pair in co_single:
            if c in pair: 
                one_ch2oh = one_ch2[:]
                one_ch2oh.append(pair[1])
                for pair_oh in oh_single:
                    if pair[1] in pair_oh:
                        one_ch2oh.append(pair_oh[1])
                        
    return list_ch3, list_ch2oh
    
def analyze_och3(ch_single, co_single):
    '''Get the OCH3 groups '''
    list_ch3 = analyze_terminal(ch_single, co_single)[0]
    list_och3 = []
    for ele_ch3 in list_ch3:
        ele_c = ele_ch3[0]
        for pair_co in co_single:
            if ele_c in pair_co:
                one_och3 = []
                one_och3.append(pair_co[1])
                one_och3.extend(ele_ch3)
                list_och3.append(one_och3)
    return list_och3

def analyze_cch3(cc_single, cc_double, ch_single, co_single):    
    '''Calculate all the C connects to CH3 groups'''
    ring_list = analyze_rings(cc_double)
    ring_c   = [i for pair in ring_list for i in pair]
    
    list_ch3 = analyze_terminal(ch_single, co_single)[0]
    
    list_cch3 = []   # C--CH3, the C is from aliphatic chain
    list_arch3 = []  # Ar--Ch3, the C will be C-alpha in the simplified models
        
    for ele_ch3 in list_ch3:
        one_cch3 = [] # Stores the atoms in one C-CH3 group.
        ele_c = ele_ch3[0]
        for pair in cc_single:
            if ele_c in pair:
                c_to_ch3 = [i for i in pair if i is not ele_c][0]
                one_cch3.append(c_to_ch3)
                one_cch3.extend(ele_ch3)
                if c_to_ch3 in ring_c:
                    list_arch3.append(one_cch3)
                else:
                    list_cch3.append(one_cch3)
    return list_cch3, list_arch3

def analyze_coch3(ch_single, co_single, cc_double):
    '''Get the Ar-OCH3 groups or C-OCH3, note: C-OCH3 might can not be observed bacause all the OCH3 connects to Ar '''
    list_och3 = analyze_och3(ch_single, co_single)
    ring_list = analyze_rings(cc_double)
    ring_c   = [i for pair in ring_list for i in pair]
    
    list_aroch3 = []
    list_coch3  = [] 
    for ele_och3 in list_och3:
        ele_o = ele_och3[0]
        for pair_co in co_single:
            if ele_o in pair_co:
                one_aroch3 = []
                one_coch3  = []
                c_possibile = [i for i in pair_co if i is not ele_o][0]
                if c_possibile not in ele_och3: 
                    if c_possibile in ring_c:
                        one_aroch3.append(pair_co[0])
                        one_aroch3.extend(ele_och3)
                        list_aroch3.append(one_aroch3)
                    else:
                        one_coch3.append(pair_co[0])    
                        one_coch3.extend(ele_och3)
                        list_coch3.append(one_coch3)
    return list_aroch3, list_coch3

def analyze_hbond(oh_single, oo_hbond, co_single):
    list_o_donor = [pair[0] for pair in oh_single]
    list_o_acceptor = []
    
    case_22 = []
    case_21 = []
    for pair_oo in oo_hbond:
        o1, o2 = pair_oo
        if o1 in list_o_donor and o2 in list_o_donor:
            case_22.append(pair)
        elif o1 in list_o_donor or o2 in list_o_donor:
            case_21.append(pair)
#    for donor_o in list_o_donor:
#        acceptors = []
#        for pair_oo in oo_hbond:
#            if donor_o in pair_oo:
#                acceptor_o = [i for i in pair_oo if i is not donor_o][0]
#                acceptors.append(acceptor_o)

def calc_hbond_num1(coords):
    
    ele_list = [line.rstrip().split()[0] for line in coords]
    if 'O' in ele_list and 'C' in ele_list:
        cc_single, cc_double, co_single, ch_single, oh_single, oo_hbond, oh_long_bond = get_dis_all(coords)    
        pair_hbonds = []
        for pair_oo in oo_hbond:
            o1, o2 = pair_oo
            for pair_oh in oh_single:
                if o1 == pair_oh[0]:
                    h1 = pair_oh[1]
                    if [o2, h1] in oh_long_bond:
                        pair_hbonds.append([o2,h1])
                elif o2 == pair_oh[0]:
                    h2 = pair_oh[1]
                    if [o1, h2] in oh_long_bond:
                        pair_hbonds.append([o1,h2])
        return len(pair_hbonds)
    else:
        return 0
#    print(pair_hbonds)
    
def calc_hbond_num2(coords):
    ele_list = [line.rstrip().split()[0] for line in coords]
    if 'O' in ele_list:
        cc_single, cc_double, co_single, ch_single, oh_single, oo_hbond, oh_long_bond = get_dis_all(coords)    
        num_hbonds = 0
        h_list = [pair[1] for pair in oh_long_bond]
        h_list_oh = [pair[1] for pair in oh_single]
        for h in h_list:
            if h in h_list_oh:
                num_hbonds += 1
        return num_hbonds  
    else:
        return 0     

#################
def analyze_ring_positions(cc_single, cc_double, co_single, one_ring):
    list_c1 = analyze_c1(cc_single, cc_double)
    list_aroch3 = analyze_coch3(ch_single, co_single, cc_double)[0]
    list_c345 = [i[0] for i in list_aroch3]    
    
    for c in one_ring:
        if c in list_c1:
            c1 = c
            
    if len(list_c345) == 3: 
        pairs_345 = list(itertools.product(list_345, list_345))
        cc_double_345 = []
        for pair in pairs_345:
            if pair in cc_double:
                cc_double_345.append(pair)
        c_345 = [i for pair in cc_double_345 for i in pair]
        for i in c_345:
            if c_345.count(i) == 2:
                c4 = i 
                
#cc_single, cc_double, co_single, ch_single, oh_single, oo_hbond = get_dis_all(coords)

#list_c1, list_c_alpha, list_c_beta, list_c_gamma  = analyze_chain(cc_single, cc_double)
#analyze_O_chain(cc_single, cc_double, co_single)
#print(list_c_beta)
#list_ch3, list_ch2 = analyze_ch3(ch_single) 
#list_cch3, list_arch3 = analyze_cch3(cc_single, cc_double, ch_single)
#list_och3 = analyze_och3(ch_single, co_single)
#list_aroch3, list_coch3= analyze_coch3(ch_single, co_single, cc_double)

#print(list_c_beta)
#def check_55(ring_1, ring_2, cc_total):
#    list_c = list(itertools.product(ring_1, ring_2))
#    c5_atoms = []
#    for i in list_c:       
#        if i in cc_total or i[::-1] in cc_total:
#            print('found 5-5 bond')
#            c5_atoms.extend(i)
#    return c5_atoms        
#            
#c5_atoms = check_55(ring_list[0], ring_list[1], cc_total)
#print(c5_atoms)
#            


        
