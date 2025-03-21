from ase import Atoms, io
import numpy as np
from cluster import * 

atoms = io.read('./POSCAR')  # Reads POSCAR, XYZ, or CIF
adsorption_results = classify_N2_adsorption(atoms)
for category, ru_site_str in adsorption_results:
    print(f"{category:<12},{ru_site_str}")


