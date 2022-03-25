#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os 

if not os.path.isfile('vasprun.xml'):
    print('Error: No vasprun_xml file found!')
    exit()

f_in = open('vasprun.xml', 'r')
lines_vasprun = f_in.readlines()
f_in.close()

look_version    = 'name="version"'
look_incar      = 'incar>'
look_parameter  = 'parameters>'
look_atominfor  = 'atominfo>'
look_structure  = 'structure'
look_nedos      = 'NEDOS'
look_fermi      = 'efermi'
look_kpoints    = 'kpoints>'

line_version    = []
line_incar      = []
line_parameter  = []
line_atominfor  = []
line_structure  = []
line_nedos      = []
line_fermi      = []
line_kpoints    = []
#########################################3
look_list = [look_nedos, look_fermi]
line_list = [line_nedos, line_fermi]

dict_line = dict(zip(look_list, line_list))
### Get all the line numbers for by using the key words in the look_list
for num, line in enumerate(lines_vasprun):
    for k, v in dict_line.items():
        if k in line: 
            v.append(num)

#########################################            
def get_version():
    line_num_version = dict_line.get(look_version)[-1]
    version          = lines_vasprun[line_num_version].split()[2].split('>')[1]
    return version

def get_nedos():
    line_num_dos = dict_line.get(look_nedos)[-1]
    nedos        = lines_vasprun[line_num_dos].split()[3].split('<')[0]   # line.strip().split()[3][:-4]
    return nedos

def get_fermi():
    line_num_fermi  = dict_line.get(look_fermi)[-1]
    e_fermi     =  lines_vasprun[line_num_fermi].split()[2]
    return float(e_fermi)

def get_kpoints():
    line_kpoints_start,line_kpoints_end = dict_line.get(look_kpoints)
    