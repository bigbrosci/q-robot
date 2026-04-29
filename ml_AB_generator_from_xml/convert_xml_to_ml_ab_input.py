import os
from ase.io import read
import xml.etree.ElementTree as ET
from ase.io.vasp import read_vasp_xml
from lxml import etree

# Last edit: gxy @ 2/27/2023 at 2:43 PM Eastern
# Last edit: qli @ 01/11/2025 at 2:43 PM Eastern

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
    

    geometries = read_vasp_xml('vasprun.xml', index=slice(None)) 
    
    
    # Initialize a list to store coordinates for all structures
    atom_positions_all = []
    forces_all = []
    
    # Iterate over all structures
    for geom in geometries:
        # Extract atomic positions for the current structure
        positions = geom.get_positions()
        atom_positions_all.append(positions)

        forces = geom.get_forces(apply_constraint=False)
        forces_all.append(forces)
    

    # Parse the vasprun.xml file using lxml
    parser = etree.XMLParser(recover=True)  # Recover from malformed XML
    tree = etree.parse('./vasprun.xml', parser)
    root = tree.getroot()
    
    # Initialize  lists to store the total energy and stress for all configurations
    toten_all = []
    stress_all = []   # in kB
    
    # Iterate through all <energy> tags
    for energy in root.findall(".//calculation/energy"):
        e_fr_energy = energy.find("i[@name='e_fr_energy']")
        if e_fr_energy is not None:
            toten_all.append(float(e_fr_energy.text))
    print(toten_all)
    
    # Find all <varray name="stress"> blocks
    for stress in root.findall(".//calculation/varray[@name='stress']"):
        rows = stress.findall("v")
        if rows and len(rows) == 3:  # Ensure there are 3 rows (XX, XY, XZ; YX, YY, YZ; ZX, ZY, ZZ)
            # Extract stress values from the rows
            stress_matrix = [[float(value) for value in row.text.split()] for row in rows]
    
            # Arrange stress in the order: XX, YY, ZZ, XY, YZ, ZX
            stress_all.append([
                stress_matrix[0][0],  # XX
                stress_matrix[1][1],  # YY
                stress_matrix[2][2],  # ZZ
                stress_matrix[0][1],  # XY
                stress_matrix[1][2],  # YZ
                stress_matrix[2][0],  # ZX
            ])



    # Debug
    print('system_name', system_name)
    print('univ_scaling', univ_scaling)
    print('cell_vectors', cell_vectors)
    print('atom_symbols', atom_symbols)
    print('n_atoms_in_type', n_atoms_in_type)
    print('atom_positions', len(atom_positions_all))
    print('total energy', len(toten_all))
    print('atom_forces', len(forces_all))
    print('stress', len(stress_all))

    
    data = [] 
    for i in range(len(toten_all)): 
        data_i = [system_name, len(atom_symbols), sum(n_atoms_in_type), atom_symbols, n_atoms_in_type,
                cell_vectors, atom_positions_all[i], toten_all[i], forces_all[i], stress_all[i]]
        data.append(data_i)
    # data = [system_name, len(atom_symbols), sum(n_atoms_in_type), atom_symbols, n_atoms_in_type,
    #         cell_vectors, atom_positions, toten[-1], atom_forces, stress[-1]]
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
    conv += block_sep_2
    conv += 'CTIFOR\n'
    conv += block_sep_1
    conv += '2.000000000000000E-002\n'    

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
        if data[1] > max_num_types:
            max_num_types = data[1]

        # Check if total number of atoms is more than current max
        if data[2] > max_atoms_in_system:
            # print(max_atoms_in_system)
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
    print('The output file is named as: ML_AB_merged' )
    return


# location = sys.argv[1]
# calc_data = read_in_vasp_calc(location) 
# create_ML_AB_input(calc_data)


# Main script
if __name__ == "__main__":
    # Get all folders in the current directory
    current_path = os.getcwd()
    folders = [d for d in os.listdir(current_path) if os.path.isdir(os.path.join(current_path, d))]

    # Initialize an empty list to store the merged calc_data
    merged_calc_data = []

    # Process each folder
    for folder in folders:
        folder_path = os.path.join(current_path, folder)
        try:
            # Read in calculation data for the current folder
            calc_data = read_in_vasp_calc(folder_path)
            print(f"Processed {folder}: {len(calc_data)} items")
            # Append the data to the merged list
            merged_calc_data.extend(calc_data)
        except Exception as e:
            print(f"Error processing folder {folder}: {e}")

    # Pass the merged data to create_ML_AB_input
    create_ML_AB_input(merged_calc_data)
    print("All folders processed!")