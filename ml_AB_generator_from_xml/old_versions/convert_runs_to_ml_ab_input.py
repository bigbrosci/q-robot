import os
from ase.io import read

# Last edit: gxy @ 2/27/2023 at 2:43 PM Eastern

# functions:
#
# read_in_vasp_calc(location)
#     parse in a vasp calculation from <location>
#     returns data in an array
# 
# convert_data_to_str(num, data)
#     convert data (parsed by read_in_vasp_calc) into a long string
#     returns data in a string
# 
# create_ML_AB_input(calc_data)
#     uses an array of calculation data (parsed by read_in_vasp_calc) to build "ML_AB_merged"
#     No return type
#

block_sep_0 = '**************************************************\n'
block_sep_1 = '--------------------------------------------------\n'
block_sep_2 = '==================================================\n'

def read_in_vasp_calc(location):
    original_cwd = os.getcwd()
    os.chdir(location)

    # Read in 1) system name, 2) number of atom types, 3) number of atoms 
    #         4) atom types, 5) atom numbers, 6) primitive lattice vectors
    # All from the CONTCAR file
    # I will assume the CONTCAR file is written in vasp5 format
    with open('CONTCAR', 'r') as f:
        contcar_header = list(f)[:7]
    
    # 1st line of CONTCAR is system name
    system_name = contcar_header[0].strip()
    
    # 2nd line is the scaling factor. Almost always we use 1
    univ_scaling = float(contcar_header[1].strip())
    
    # 3rd~5th lines are the cell vectors, scaled down by the scaling factor
    cell_vectors = [[float(x)*univ_scaling for x in line.strip().split()] for line in contcar_header[2:5]]
    
    # Atom types and numbers of each atom
    atom_symbols = contcar_header[5].strip().split()
    n_atoms_in_type = [int(x) for x in contcar_header[6].strip().split()]
    
    # Read in atomic positions and forces from OUTCAR
    geom = read('OUTCAR@-1')
    atom_positions = geom.get_positions()
    atom_forces = geom.get_forces(apply_constraint=False)
    
    # Read in total energy and cell stresses directly from OUTCAR
    # Note that the converged "TOTEN" is used by VASP's MD energy conservation and FF training
    # Cell stresses are in the order of "XX, YY, ZZ, XY, YZ, and ZX"
    # We'll take the last one for each
    toten = []
    stress = []
    with open('OUTCAR', 'r') as f:
        for line in f:
            if 'free  energy   TOTEN' in line:
                toten.append(float(line.strip().split()[4]))
            if 'in kB' in line:
                stress.append([float(x) for x in line.strip().split()[2:]])
    
    """
    # Debug
    print('system_name', system_name)
    print('univ_scaling', univ_scaling)
    print('cell_vectors', cell_vectors)
    print('atom_symbols', atom_symbols)
    print('n_atoms_in_type', n_atoms_in_type)
    print('atom_positions', atom_positions)
    print('total energy', toten[-1])
    print('atom_forces', atom_forces)
    print('stress', stress[-1])
    """
    
    data = [system_name, len(atom_symbols), sum(n_atoms_in_type), atom_symbols, n_atoms_in_type,
            cell_vectors, atom_positions, toten[-1], atom_forces, stress[-1]]
    os.chdir(original_cwd)
    return data

def convert_data_to_str(num, data):
    conv = ''

    # Configuration number
    conv += block_sep_0
    conv += 'Configuration num. %s'%(int(num))
    conv += '\n'

    # System name
    conv += block_sep_2
    conv += 'System name\n'
    conv += block_sep_1
    conv += data[0]
    conv += '\n'

    # Number of atom types
    conv += block_sep_2
    conv += 'The number of atom types\n'
    conv += block_sep_1
    conv += str(int(data[1]))
    conv += '\n'

    # Number of atoms
    conv += block_sep_2
    conv += 'The number of atoms\n'
    conv += block_sep_1
    conv += str(int(data[2]))
    conv += '\n'

    # Atom types and numbers
    conv += block_sep_0
    conv += 'Atom types and atom numbers\n'
    conv += block_sep_1
    for i in range(len(data[3])):
        conv += '%s %s\n'%(data[3][i], data[4][i])

    # Skip CTIFOR block

    # Lattice vectors
    conv += block_sep_2
    conv += 'Primitive lattice vectors (ang.)\n'
    conv += block_sep_1
    for vect in data[5]:
        conv += '%.16E %.16E %.16E\n'%tuple(vect)

    # Atomic positions
    conv += block_sep_2
    conv += 'Atomic positions (ang.)\n'
    conv += block_sep_1
    for vect in data[6]:
        conv += '%.16E %.16E %.16E\n'%tuple(vect)

    # Total energy
    conv += block_sep_2
    conv += 'Total energy (eV)\n'
    conv += block_sep_1
    conv += '%.10E\n'%(data[7])

    # Forces
    conv += block_sep_2
    conv += 'Forces (eV ang.^-1)\n'
    conv += block_sep_1
    for vect in data[8]:
        conv += '%.16E %.16E %.16E\n'%tuple(vect)

    # Stress
    conv += block_sep_2
    conv += 'Stress (kbar)\n'
    conv += block_sep_1
    conv += 'XX YY ZZ\n'
    conv += block_sep_1
    conv += '%.16E %.16E %.16E\n'%tuple(data[9][:3])
    conv += block_sep_1
    conv += 'XY YZ ZX\n'
    conv += block_sep_1
    conv += '%.16E %.16E %.16E\n'%tuple(data[9][3:])

    return conv

def create_ML_AB_input(calc_data):
    # Some basic information about the dataset
    n_configs = len(calc_data)
    max_num_types = 0
    atom_types = []
    max_atoms_in_system = 0
    max_atoms_per_type = 0
    
    # Loop over calc_data to find basic information about the dataset
    for data in calc_data:
        # Check if number of types is more than current max
        print(data[1])
        if data[1] > max_num_types:
            max_num_types = data[1]

        # Check if total number of atoms is more than current max
        if data[2] > max_atoms_in_system:
            print(max_atoms_in_system)
            max_atoms_in_system = data[2]

        # Check atom types and max atoms per type
        for i in range(len(data[3])):
            symbol = data[3][i]
            number = data[4][i]
            if symbol not in atom_types:
                atom_types.append(symbol)
            if number > max_atoms_per_type:
                max_atoms_per_type = number

    # Put together the ML_AB file
    with open('ML_AB_merged', 'w') as f:
        # Title
        f.write('1.0 Version\n')

        # Number of training structures
        f.write(block_sep_0)
        f.write('The number of configurations\n')
        f.write(block_sep_1)
        f.write('%s\n'%(n_configs))

        # Max number of atom type
        f.write(block_sep_0)
        f.write('The maximum number of atom type\n')
        f.write(block_sep_1)
        f.write('%s\n'%(max_num_types))

        # Atom types in datafile
        f.write(block_sep_0)
        f.write('The atom types in the data file\n')
        f.write(block_sep_1)
        for i in range(int(len(atom_types)/3)):
            f.write('%s %s %s\n'%tuple(atom_types[i*3:(i+1)*3]))
        for atom_type in atom_types[(i+1)*3:]:
            f.write('%s '%(atom_type))
        f.write('\n')

        # Max number of atoms per system
        f.write(block_sep_0)
        f.write('The maximum number of atoms per system\n')
        f.write(block_sep_1)
        f.write('%s\n'%(max_atoms_in_system))

        # Max number atoms per type
        f.write(block_sep_0)
        f.write('The maximum number of atoms per atom type\n')
        f.write(block_sep_1)
        f.write('%s\n'%(max_atoms_per_type))

        # Reference atomic energy. I'm going to leave this block as all 0
        f.write(block_sep_0)
        f.write('Reference atomic energy (eV)\n')
        f.write(block_sep_1)
        for i in range(int(len(atom_types)/3)):
            f.write('%.16E %.16E %.16E\n'%(0, 0, 0))
        for atom_type in atom_types[(i+1)*3:]:
            f.write('%.16E '%(0))
        f.write('\n')

        # Atomic masses
        # Please look this up from the VASP POTCAR files
        # Note that the masses are slightly different from those listed elsewhere
        # E.g. mass of H = 1 while it is usually listed as 1.008 elsewhere
        f.write(block_sep_0)
        f.write('Atomic mass\n')
        f.write(block_sep_1)
        
        # Number of basis sets per atom type
        # This is calculated by ML_ISTART = 3, so I will leave this as dummy vars
        f.write(block_sep_0)
        f.write('The numbers of basis sets per atom type\n')
        f.write(block_sep_1)
        for i in range(int(len(atom_types)/3)):
            f.write('%s %s %s\n'%(1, 1, 1))
        for atom_type in atom_types[(i+1)*3:]:
            f.write('%s '%(1))
        f.write('\n')

        # Basis sets for each atom type
        for atom_type in atom_types:
            f.write(block_sep_0)
            f.write('Basis set for %s\n'%(atom_type))
            f.write(block_sep_1)
            f.write('1 1\n')

        # Loop over calc_data and write dataset to file
        for i in range(len(calc_data)):
            f.write(convert_data_to_str(i+1, calc_data[i]))

    return

location = '.'
data = read_in_vasp_calc(location)
create_ML_AB_input(data)
