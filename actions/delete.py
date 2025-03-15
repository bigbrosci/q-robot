#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Delete atoms from POSCAR using ASE'''

import sys
from ase.io import read, write
from ase import Atoms

def get_atoms_to_delete(args, atoms):
    """
    Analyze the command-line arguments and return the indices of atoms to be deleted.
    Arguments can be either element names (e.g., 'C', 'H', 'O') or atom indices (e.g., '10', '12', '14').
    Duplicates in the arguments will be removed.
    """
    elements_to_delete = []
    atom_indices_to_delete = []

    # Parse the arguments
    for arg in args:
        try:
            # Try to parse as an index (1-based)
            index = int(arg)
            atom_indices_to_delete.append(index - 1)  # Convert to 0-based index for ASE
        except ValueError:
            # If not an index, treat it as an element name
            elements_to_delete.append(arg)

    # Remove duplicates by converting to sets
    elements_to_delete = list(set(elements_to_delete))
    atom_indices_to_delete = list(set(atom_indices_to_delete))

    # Identify atoms to delete based on elements or indices
    atoms_to_delete = set()

    # If elements are provided, find indices of those atoms
    if elements_to_delete:
        for i, atom in enumerate(atoms):
            if atom.symbol in elements_to_delete:
                atoms_to_delete.add(i)

    # If atom indices are provided, add those indices to the set
    if atom_indices_to_delete:
        for idx in atom_indices_to_delete:
            if atoms[idx].symbol in elements_to_delete or not elements_to_delete:
                atoms_to_delete.add(idx)

    return atoms_to_delete


def delete_atoms_and_save(file_in, atoms_to_delete, atoms):
    """
    Delete the specified atoms and save the modified structure to a new file.
    """
    # Remove atoms from the Atoms object
    atoms = atoms[[i for i in range(len(atoms)) if i not in atoms_to_delete]]

    # Write the modified atoms to a new file
    out_name = file_in.replace("POSCAR", "POSCAR_deleted")
    write(out_name, atoms)

    print(f'\nThe output file is: {out_name}')


if len(sys.argv) < 3:
    print('\nCommand Usage: delete.py POSCAR element_or_index1 element_or_index2 ...')
    print('Example: delete.py POSCAR C H O 10 12 14')
    print('This will delete all C, H, O atoms and the 12th and 14th atoms from the POSCAR file')
    exit()

file_in = sys.argv[1]
args = sys.argv[2:]

# Read POSCAR file using ASE
atoms = read(file_in)

# Get atoms to delete based on the command-line arguments
atoms_to_delete = get_atoms_to_delete(args, atoms)

# Delete atoms and save the new POSCAR file
delete_atoms_and_save(file_in, atoms_to_delete, atoms)

