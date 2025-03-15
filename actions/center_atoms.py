#!/usr/bin/env python3 
from ase.io import read, write

# Read the POSCAR file containing a cluster in a box
cluster = read('POSCAR', format='vasp')

# Center the cluster in the box
#cluster.center(about=(0.5, 0.5, 0.5))
cluster.center()

# Wrap atoms into the unit cell
cluster.wrap()

# Save the centered structure as 'POSCAR_centered'
write('POSCAR_centered', cluster, format='vasp', vasp5=True)
