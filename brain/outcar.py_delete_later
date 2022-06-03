#!/usr/bin/env python3

''' Get the basic information from OUTCAR
From OUTCAR, we can get pelnty of informations
# Basic information 
1) the vasp version
2) the INCAR paramters
3) KPOINTS numbers
4) the atom positions, cell size ('VOLUME and BASIS-vectors are now')
5) POTCAR information  (look_pot    = 'POTCAR')
6) time for each eletronic and ion steps.  (LOOP, LOOP+)
7) Is the optimization converge or not (reached required accuracy )
# Calculated properties 
0) Energy ('energy(sigma->0)')
1) fermi and vacuum energies (workfunction, dos)
2) frequency calculations (IBRION = 5, POTIM = 0.015, NFREE = 2)
3) electron occupations (band structures)
4) magnetization (LORBIT = 11, look_mag = 'magnetization (x)')
5) optimization steps and the forces on each atom.
6) VDW version, parameters ( DFTD3 )
7) NEB vtst 
8) DFT + U 

Remember: In more cases, bash commands such as grep are more handful than using python scripts. 
So: many of the functions in this script are just for fun!!!

'''
import os
import numpy as np


#if os.path.isfile('OUTCAR'):
#    f = open('OUTCAR')
#    lines_o = f.readlines()
#    f.close()
#else:
#    print('No OUTCAR in Current Path. Bye!')
#    exit()    

## The following look_up parameters are listed by their sequence in the OUTCAR file

look_pot            = 'POTCAR'
look_incar_start    = 'Startparameter'
look_separate        = '------------------------------'*2  # the first one after  look_incar_start
look_vectors        = 'VOLUME and BASIS-vectors are now'
look_kpoints        = 'irreducible k-points'
look_position       = 'POSITION      '  ## or  'TOTAL-FORCE (eV/Angst)'
look_iteration      = 'Iteration'
look_energy         = 'energy(sigma->0) ='
look_time_ele       = 'LOOP'
look_time_ion       = 'LOOP+'
look_fermi          =  'E-fermi'
look_vacuum         = 'vacuum level'
look_freq           = 'f  ='
look_freq_i         = 'f/i='
look_converge       = 'reached required accuracy'
look_magnetization  = 'magnetization (x)'
look_vdw            = 'IVDW'

line_pot            =  []
line_incar_start    =  []
line_separate       =  []
line_vectors        =  []
line_kpoints        =  []
line_position       =  []
line_iteration      =  []
line_energy         =  []
line_time_ele       =  []
line_time_ion       =  []
line_fermi          =  []
line_vacuum         =  []
line_freq           =  []
line_freq_i         =  []
line_converge       =  []
line_magnetization  =  []
line_vdw            =  []


def get_vasp_version(lines_o):
    vasp_version        = lines_o[0].strip().split()[0]
    return vasp_version

def get_dict_line(lines_o):
    
    look_list = [look_pot, look_incar_start, look_separate, look_vectors, look_kpoints, look_position, look_iteration, look_energy,
                 look_time_ele, look_time_ion, look_fermi, look_vacuum, look_freq, look_freq_i, look_converge, look_magnetization, look_vdw]
    line_list = [line_pot, line_incar_start, line_separate, line_vectors, line_kpoints, line_position, line_iteration, line_energy,
                 line_time_ele, line_time_ion, line_fermi, line_vacuum, line_freq, line_freq_i, line_converge, line_magnetization, line_vdw]
    
    dict_line = dict(zip(look_list, line_list))
    
    for num, line in enumerate(lines_o):
        for k, v in dict_line.items():
            if k in line: 
                v.append(num)            
    return dict_line 
        
def get_potcars(lines_o, look = look_pot):
    dict_line = get_dict_line(lines_o)
    infor_l = []
    for i in dict_line.get(look_pot):
        infor_l.append(lines_o[i].strip().split(':')[1].strip())    
    potcars  = infor_l[0:len(infor_l)/2]
    return potcars
#potcars = get_potcars()
#print(potcars)      
    
def get_incar(lines_o):
    dict_line = get_dict_line(lines_o)
    def analyze_incar(lines):
        dict_incar = {}
        for line in lines:
            if '=' in line: 
                if ';' in line:
                    line_ele = line.split(';')
                    for i in line_ele:
                        item_ele = [y.strip() for y in i.split('=')]
                        try:
                            dict_incar[item_ele[0]] = item_ele[1].split()[0]
                        except IndexError:
                            pass
                else:
                    item_ele = [y.strip() for y in line.rstrip().split('=')]
                    if 'LDAU' in line:
                        dict_incar[item_ele[0].split()[-1]] = item_ele[1]
                    elif 'POMASS' in line or 'ZVAL' in line or 'RWIGS' in line:
                        dict_incar[item_ele[0]] = item_ele[1]
                    elif 'DFIELD' in line:
                        pass
                    else:    
                        dict_incar[item_ele[0]] =  item_ele[1].split()[0].strip()
        return dict_incar
    
    incar_start_num = dict_line.get(look_incar_start)[0]
    separate_num    = dict_line.get(look_separate)
    for num, line_num in enumerate(separate_num):
        if   15 < incar_start_num - line_num <30:
            incar_end_num = separate_num[num + 1]
    incar_lines = lines_o[incar_start_num : incar_end_num]        
    dict_incar  = analyze_incar(incar_lines)
    return dict_incar
#dict_incar = get_incar(lines_o)
#print(dict_incar)

def get_volume_vectors(lines_o, look_vectors=look_vectors):
    dict_line = get_dict_line(lines_o)
    volume_lines = dict_line.get(look_vectors)[-1]
    volume       = float(lines_o[volume_lines+3].strip().split(':')[1].strip())
    va           = [float(i) for i in lines_o[volume_lines+5].split()[0:3] ]
    vb           = [float(i) for i in lines_o[volume_lines+6].split()[0:3] ]
    vc           = [float(i) for i in lines_o[volume_lines+7].split()[0:3] ]
    length_abc   = [float(i) for i in lines_o[volume_lines+10].split()[0:3] ]
    vector       = np.transpose(np.array([va, vb, vc]))
    return vector, length_abc, volume
#print(get_volume_vectors(look_vectors))

def get_kpoints(lines_o, look_kpoints=look_kpoints):
    dict_line = get_dict_line(lines_o)
    line_kpoints = dict_line.get(look_kpoints)[0]
    kpoints      = lines_o[line_kpoints].split()[1]
    return  kpoints
#print(get_kpoints())

def get_position(lines_o, look_position = look_position):
    dict_line = get_dict_line(lines_o)
    line_position_start  =  dict_line.get(look_position)[-1]
    separate_num    = dict_line.get(look_separate)
    for num, line_num in enumerate(separate_num):
        if line_num == line_position_start + 1:
            line_position_end = separate_num[num + 1]
    position_lines = lines_o[line_position_start+2: line_position_end]
    cartesian_corrd = []
    for i in position_lines:
        cartesian_corrd.append(i.split()[0:3])
    return cartesian_corrd

def get_iteration_infor(lines_o):
    dict_line = get_dict_line(lines_o)
    lines_iteration = dict_line.get(look_iteration)
    lines_time_ele_raw  = dict_line.get(look_time_ele)
    lines_time_ion  = dict_line.get(look_time_ion)
    lines_time_ele  = [ i for i in lines_time_ele_raw if i not in lines_time_ion]
    #print(lines_iteration)
    
    for num, line_num in enumerate(lines_iteration):
        line_ele = lines_o[line_num].split()
        ionic_step = int(line_ele[2].replace('(',''))
        ele_step   = int(line_ele[3].replace(')',''))
        ele_time    = float(lines_o[lines_time_ele[num]].split()[-1])
        print(ionic_step, ele_step, ele_time) 

def get_fermi(lines_o):
    dict_line = get_dict_line(lines_o)
    line_fermi = dict_line.get(look_fermi)[-1]
    e_fermi    = lines_o[line_fermi].split()[2]
    return e_fermi

def get_vacuum(lines_o):
    dict_line = get_dict_line(lines_o)
    line_vacuum = dict_line.get(look_vacuum)[-1]     
    e_vacuums   = lines_o[line_vacuum].split()[-2:]
    vacuum_up, vacuum_dn = e_vacuums
    return vacuum_up, vacuum_dn

def get_freq(lines_o):
    '''Imaginary Frequencies are not concluded '''
    dict_line = get_dict_line(lines_o)
    lines_freq = dict_line.get(look_freq)
    nu  = []
    zpe = []
    for i in lines_freq:
        line_ele = lines_o[i].split()
        nu.append(float(line_ele[7]))
        zpe.append(float(line_ele[9]))
    return nu, zpe
        
def get_freq_i(lines_o):
    '''Only consider the Imaginary Frequencies'''
    dict_line = get_dict_line(lines_o)
    lines_freq_i = dict_line.get(look_freq_i)
    nu     = []
    zpe = []
    for i in lines_freq_i:
        line_ele = lines_o[i].split()
        nu.append(float(line_ele[6]))
        zpe.append(float(line_ele[8]))
    return nu, zpe

def converge_or_not(lines_o):
    dict_line = get_dict_line(lines_o)
    line_converge  =  dict_line.get(look_converge)
    is_converge_or_not = False
    if len(line_converge) == 1:
        is_converge_or_not = True
    return is_converge_or_not   
#print converge_or_not()    
    
def get_mag(lines_o):
    dict_line = get_dict_line(lines_o)
    line_mag_start = dict_line.get(look_magnetization)[-1] + 4 
    line_mag_end = 1
    for num, line in enumerate(lines_o[line_mag_start:]):
        if 'tot  ' in line:
            line_mag_end = num - 1       
    lines_mag = lines_o[line_mag_start: line_mag_start + line_mag_end] 
    dict_mag = {}
    for line in lines_mag:
        line_ele = line.split()
        dict_mag[int(line_ele[0])] = [float(i) for i in line_ele[1:]]
    '''Note: the key is integers and the value is a list for magnetization of s p d   '''    
    return dict_mag    
#print(get_mag())
    
def get_vdw(lines_o):
    dict_line = get_dict_line(lines_o)
    line_vdw  = dict_line.get(look_vdw)[-1]
    vdw = lines_o[line_vdw-1].rstrip()
    return vdw
#get_vdw()    
    
def get_energy(lines_o):
    dict_line = get_dict_line(lines_o)
    line_energy = dict_line.get(look_energy)[-1]
    energy = lines_o[line_energy].rstrip().split()[-1]
    return float(energy)
