from ase.io import read, write
from ase.build import surface, sort
from ase.constraints import FixAtoms
import numpy as np

# Read bulk structure from CONTCAR
bulk = read("CONTCAR")

# Get unique element order from CONTCAR
unique_elements = list(dict.fromkeys(bulk.get_chemical_symbols()))

# Define the (111) surface with 4 layers and 30 Å vacuum (will be corrected later)
miller_indices = (1, 1, 1)
slab = surface(bulk, miller_indices, layers=4, vacuum=30.0)

# Use ASE's sort function to maintain element order from CONTCAR
slab = sort(slab, tags=[unique_elements.index(sym) for sym in slab.get_chemical_symbols()])

# Adjust the slab so the lowest atom is at z = 0.1 to avoid visualization issues
z_min = np.min(slab.positions[:, 2])
slab.positions[:, 2] -= z_min - 0.1

# Set vacuum to **15 Å** instead of ASE's default 30 Å
slab.cell[2, 2] = np.max(slab.positions[:, 2]) + 15.0

# **Fix the bottom two layers**
z_coords = slab.positions[:, 2]
sorted_indices = np.argsort(z_coords)
num_fixed = len(slab) // 2  # Bottom half will be fixed

# Apply constraints
constraint = FixAtoms(indices=sorted_indices[:num_fixed])
slab.set_constraint(constraint)

# Write final POSCAR file
write("POSCAR_111", slab, format="vasp")

