#!/usr/bin/env python3 
import ase
import ase.io
#import ase.io.espresso as ep
from ase.thermochemistry import IdealGasThermo
import os, subprocess
import numpy as np
from ase import Atoms
from scipy import constants  as sc
import pandas as pdpd


import matplotlib.pyplot as plt

# RytoeV = sc.physical_constants['Rydberg constant times hc in eV'][0]

'''Constants'''
Ry_to_eV = 13.605662285137 
eV_to_cm = 8065.54429
aBohr = 0.52917720859


def get_lines(file_in):
    f_in = open(file_in,'r')
    f_lines = f_in.readlines()
    f_in.close()
    return f_lines

def get_abc_from_xyz(xyz_file):
    '''Get the crystal lattice parameters from the extended xyz format '''
    lines = get_lines(xyz_file)
    cell = lines[1].split('"')[1].split()
    cell_a = '  '.join(cell[:3]) 
    cell_b = '  '.join(cell[3:6]) 
    cell_c = '  '.join(cell[6:]) 
    return cell_a, cell_b, cell_c

def get_coord_xyz(file_in):
    '''Read xyz file and get the coordinates'''
    lines = get_lines(file_in)
    '''Get Coords from xyz file ''' 
    line_start = 2
    for num, line in enumerate(lines[2:]):
        if len(line.strip()) == 0:
            line_end = num
            break
        else:
            line_end = len(lines)

    coords = lines[line_start:line_end+line_start]
    return coords

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
    
def collect_freq(freq_file='dynmat.mold'):
    ''''read the QE freq output and return the frequncy lists using cm-1 and eV units'''
    lines = get_lines(freq_file)
    n_start= 2
    n_end = 0
    for num, line in enumerate(lines):
        if line.rstrip() == '[FREQ]':
            n_start = num + 1
        if line.rstrip() == '[FR-COORD]':
            n_end = num
    freq_list = [float(i.rstrip()) for i in lines[n_start:n_end]] 
    freq_list = [i for i in freq_list if i!= 0] 
    
    vibs = np.array(freq_list)
    print(vibs)
    vib_energies = vibs / eV_to_cm  # convert cm-1 to eV
    return vibs, vib_energies  

def read_qeout(qe_out):    
    '''Get some useful information from QE out file ''' 
    lines = get_lines(qe_out)
    
    '''Items to search the QE output file '''
    dict_search = {
    'Natom':            'number of atoms/cell',
    'cell':             'celldm(1)',
    'lattice':          'crystal axes: (cart. coord. in units of alat)',
    'job':              'End of BFGS Geometry Optimization',         
    'potentialenergy':  '!    total energy',
    'spin'   :          'total magnetization ',
    'coord_start_n':    'ATOMIC_POSITIONS',                   
#    'coord_start_n':    'Begin final coordinates',                      
#    'coord_end_n' :     'End final coordinates',              
    }
    '''Initialize the Items '''
    Natom = 0 
    lca = 0 # lattice constant in a direction http://elcorto.github.io/pwtools/written/background/pwscf.html
    n_lattice = 0
    job = 'UnDone'
    potentialenergies = []
    spin = [0]
    coord_start_n   =  0
    coord_end_n     = 0
    
    '''Search the output file '''
    for line_num, line in enumerate(lines): 
        infor = line.rstrip().split()
        if dict_search.get('Natom') in line:
            Natom = int(infor[-1])
        elif dict_search.get('cell') in line:
            # n_cell = line_num
            lca = float(infor[1])
        elif dict_search.get('lattice') in line:
            n_lattice = line_num + 1 
        elif dict_search.get('potentialenergy') in line:   
            potentialenergies.append(float(infor[4]) * Ry_to_eV)
        elif dict_search.get('spin') in line:
            spin.append(infor[3])
        elif dict_search.get('job') in line:
            job = 'Done'
#            print('\nCongrates, Job Is Done!!!\n')
        elif dict_search.get('coord_start_n') in line:
            coord_start_n = line_num + 1
            coord_end_n = coord_start_n + Natom 
#        elif dict_search.get('coord_end_n') in line:
#            coord_end_n = line_num
    if job == 'UnDone':
        print('Your job is not normmaly terminated. Check it before using the values.')
    
    a1,a2,a3 = (i.rstrip().split()[3:6] for i in lines[n_lattice:n_lattice+3])
    a1 = np.array([float(i) for i in a1]) * lca * aBohr
    a2 = np.array([float(i) for i in a2]) * lca * aBohr
    a3 = np.array([float(i) for i in a3]) * lca * aBohr
    abc = np.array([a1,a2,a3])
    
    '''Summarize the search results and save it to a dictionary '''      
    dict_log  = {}
    dict_log.update({'latttice'     :abc})
    dict_log.update({'coords' : lines[coord_start_n:coord_end_n]})
    dict_log.update({'job'          :  job})
    dict_log.update({'potentialenergy' : potentialenergies[-1]})
    dict_log.update({'spin'         :  spin[-1]})
            
    return dict_log 

def save_xyz(qe_out):
    '''Save QE optmized geometry to xyz file'''
    out_name = qe_out.split('.')[0] + '.xyz'
    dict_log = read_qeout(qe_out)
    abc = dict_log.get('latttice')
    coords = dict_log.get('coords')
    geo = Atoms()
    geo.set_cell(abc)
    geo.set_pbc(((True,True,True)))
    for i in coords:
        infor = i.rstrip().split()
        ele = infor[0]
        position_ele = np.array([float(i) for i in infor[1:]])
        geo.append(ele)
        geo.positions[-1] = position_ele
    geo.center()
    ase.io.write(out_name, geo, format='xyz')
    ase.io.write('POSCAR', geo, format='vasp', vasp5=True)
    return geo
    
    
def make_nscf_qein(qe_in):
    '''Read the QE input and out put to make input for NSCF and projwfc jobs'''
    lines = get_lines(qe_in)
    qe_out = qe_in.replace('in', 'out')
    coord_line = 0
    
    prefix = None
    for line_num, line in enumerate(lines):
        if 'calculation' in line:
            lines[line_num] = line.replace('relax','nscf')
        elif './temp/' in line:
            lines[line_num] = line.replace('./temp/', '../temp/')
        elif '2 2 1' in line:
            lines[line_num] = line.replace('2','3')
        elif 'ATOMIC_POSITIONS' in line:
            coord_line = line_num
        elif 'prefix' in line:
            prefix = line.rstrip().split('=')[1].strip().replace('\'','')
    dict_log = read_qeout(qe_out) 
    coords = dict_log.get('coords')
    new_lines = lines[:coord_line] + coords
    
    f_out = open('nscf.in','w')
    f_out.writelines(new_lines)
    f_out.close()
    
    projwfc = open('/THFS/home/iciq-lq/bin/qerobot/input_temps/temp_projwfc.in')
    new_projwfc = projwfc.read().replace('PREFIX',prefix)
    with open('projwfc.in', 'w') as f:
        f.write(new_projwfc)
        

def get_gibbs_qe(qe_out, geometry= 'nonlinear', SN=1, SP = 0):
    T = 298.15 
    P = 101325.0 
    '''Thermochemistry of QE calculations '''
    vibs, vib_energies =  collect_freq('dynmat.mold')
    atoms = save_xyz(qe_out)
    if len(atoms.positions) ==0:
        atoms = save_xyz('../'+qe_out)
#    atoms = ase.io.read(qe_out.replace('out,xyz'))  # ase gui  XXX.out@-1 -o out.xyz
        qe_out_infor = read_qeout(qe_out)
    else:
        qe_out_infor = read_qeout(qe_out)
    potentialenergy = qe_out_infor.get('potentialenergy')
    
    # SN = 1 # symmetrynumber
    # SP = 0 # Spin, # singlet:0, doublet: 0.5, triplet: 1 
    
    thermo = IdealGasThermo(vib_energies=vib_energies,
                            potentialenergy=potentialenergy,
                            atoms=atoms,
                            geometry=geometry, #'nonlinear', ‘monatomic’, ‘linear’, or ‘nonlinear’
                            symmetrynumber=SN, 
                            spin=SP) 
    
    G = thermo.get_gibbs_energy(temperature=T, pressure=P)    
    
    return G            


def find_dos(atom,orbital):
    '''RETURN the file name of wanted dos '''
    path = '.'
    dos_files = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if 'pdos.pdos_atm' in name:
                dos_files.append(name)
    
    atom_key = 'atm#'+str(atom)+'('
    atom_list = []
    for dos_file in dos_files:
        if atom_key in dos_file:
            atom_list.append(dos_file)
    wanted_dos_file = None 
    for dos_file in atom_list:
        dos_orbital = dos_file.split('#')[-1].replace('(','').replace(')','')
        if dos_orbital == orbital:
            wanted_dos_file  =  dos_file

    if wanted_dos_file == None:
        print('Can not find the PDOS file for Orbital \t %s \t of Atom %s.' %(orbital, int(atom)))
    return wanted_dos_file

def save_3d(natom):
    '''Save the dz2 orbital of atom '''
    dos_file = find_dos(natom,'3d')
    dos_data = pd.read_csv(dos_file,delim_whitespace=True)
    e = np.array(dos_data['#']) # Energy
    
    
    lines = get_lines(dos_file)
    N_orbitals = lines[0].count('(E)')
    if N_orbitals > 6:
        
        dz2up= dos_data.iloc[:,3] # dz2 up 
        dz2dw = dos_data.iloc[:,4] # dz2 down
        dxzup= dos_data.iloc[:,5] # dxz up 
        dxzdw = dos_data.iloc[:,6] # dxz down    
        dyzup= dos_data.iloc[:,7] # dyz up 
        dyzdw = dos_data.iloc[:,8] # dyz down  
        dx2up= dos_data.iloc[:,9] # dx2-y2 up 
        dx2dw = dos_data.iloc[:,10] # dx2-y2down        
        dxyup= dos_data.iloc[:,11] # dxy up 
        dxydw = dos_data.iloc[:,12] # dxy down     
        
        dos_dict = {'Energy' : e, 
                    'dz2_up':dz2up, 'dz2_dw':dz2dw, 
                    'dxz_up': dxzup, 'dxz_dw': dxzdw,
                    'dyz_up': dyzup, 'dyz_dw': dyzdw,
                    'dx2-y2_up': dx2up, 'dx2-y2_dw': dx2dw,
                    'dxy_up': dxyup, 'dxy_dw': dxydw,
                    }
        
        data_out = pd.DataFrame(dos_dict)
        file_out = 'dos_3d_' + str(natom) + '.csv'
        data_out.to_csv(file_out)

    else: 
        dz2 = dos_data.iloc[:,2] # dz2   
        dxz = dos_data.iloc[:,3] # dxz  
        dyz = dos_data.iloc[:,4] # dyz    
        dx2 = dos_data.iloc[:,5] #  dx2-y2
        dxy = dos_data.iloc[:,6] # dxy   
        
        dos_dict = {'Energy' : e, 
                    'dz2':dz2, 
                    'dxz': dxz, 
                    'dyz': dyz, 
                    'dx2-y2': dx2, 
                    'dxy': dxy,
                    }
        
        data_out = pd.DataFrame(dos_dict)
        file_out = 'dos_3d_' + str(natom) + '.csv'
        data_out.to_csv(file_out)
    

def plot_dos(e,dos):
    '''Not finish yet '''
    # dup = np.array(dos_data['ldosup(E)']) # Total Up
    # ddw  = np.array(dos_data['ldosdw(E)']) #Total Down
    # fig, ax = plt.subplots()
    # ax.plot(e, dz2up, '--',color = 'r')
    # ax.plot(e, dz2dw, '--',color = 'b')
    # ax.set_xlabel('Energy/ eV')
    # ax.set_ylabel('Density of States')
    # ax.set_title(r'DOS')
    
    # # Tweak spacing to prevent clipping of ylabel
    # fig.tight_layout()
    # plt.show()
