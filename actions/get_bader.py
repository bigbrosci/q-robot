#!/usr/bin/env python3 

def read_acf(acf_file):
    """
    Read ACF.dat file and extract Bader charges for each atom.
    """
    with open(acf_file, 'r') as file:
        lines = file.readlines()
        charge_data = lines[2:-4]
        bader_charges = [float(line.split()[4]) for line in charge_data]
    return bader_charges

def read_potcar_with_zval(potcar_file):
    """
    Read POTCAR file and extract elements along with their ZVAL values.
    """
    elements_zval = {}
    with open(potcar_file, 'r') as file:
        current_element = ""
        for line in file:
            if 'VRHFIN' in line:
                current_element = line.split('=')[1].split(':')[0].strip()
            if 'ZVAL' in line:
                zval = float(line.split(';')[1].split('=')[1].strip().split()[0])
                elements_zval[current_element] = zval
    return elements_zval

def read_poscar(poscar_file):
    """
    Read POSCAR file and determine the number of each type of atom.
    """
    with open(poscar_file, 'r') as file:
        lines = file.readlines()
        atom_counts = [int(x) for x in lines[6].split()]
    return atom_counts

def calculate_bader_charge(acf_file, potcar_file, poscar_file):
    """
    Calculate the Bader charge for each atom.
    """
    bader_charges_raw = read_acf(acf_file)
    elements_zval = read_potcar_with_zval(potcar_file)
    atom_counts = read_poscar(poscar_file)

    # Calculate the adjusted Bader charges
    output = []
    atom_index = 1
    for i, count in enumerate(atom_counts):
        element = list(elements_zval.keys())[i]
        zval = elements_zval[element]
        for _ in range(count):
            adjusted_charge = zval - bader_charges_raw[atom_index - 1]
            output.append((atom_index, element, adjusted_charge))
            atom_index += 1

    return output

# To use the script, replace the file paths with your actual file paths
print('Index,Element,Charge')
output = calculate_bader_charge("./ACF.dat", "./POTCAR", "./POSCAR")
for  i in output:
    list_i = [str(j) for j in list(i)]
    print(','.join(list_i))
