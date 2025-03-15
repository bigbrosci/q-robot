#!/usr/bin/env python3 
from ase.io import read, write
from ase.constraints import FixAtoms

# Load the POSCAR file
atoms = read('POSCAR')

# Identify Ru atoms
ru_indices = [atom.index for atom in atoms if atom.symbol == 'Ru']

# Apply the constraint to fix the Ru atoms
constraint = FixAtoms(indices=ru_indices)
atoms.set_constraint(constraint)

# Save the updated POSCAR with the constraints applied
write('POSCAR_relax', atoms, format='vasp')

print("Modified POSCAR with fixed Ru atoms saved as 'POSCAR_relax'.")

