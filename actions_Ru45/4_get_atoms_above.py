#!/usr/bin/env python3
"""
Script to identify atoms above the center of mass in the z-direction.
Reads a POSCAR file using ASE and saves atom indices to atoms_above.txt
"""

from ase.io import read
import numpy as np

def get_atoms_above_center_of_mass(poscar_path, output_file='atoms_above.txt'):
    """
    Read POSCAR file, find center of mass, and identify atoms above it in z-direction.
    
    Parameters:
    -----------
    poscar_path : str
        Path to the POSCAR file
    output_file : str
        Output file to save atom indices (default: atoms_above.txt)
    """
    
    # Read the structure from POSCAR
    atoms = read(poscar_path)
    
    # Get atomic masses and positions
    masses = atoms.get_masses()
    positions = atoms.get_positions()
    
    # Calculate center of mass
    total_mass = np.sum(masses)
    center_of_mass = np.sum(positions * masses[:, np.newaxis], axis=0) / total_mass
    com_z = center_of_mass[2]
    
    # Find atoms above the center of mass in z-direction
    atoms_above = []
    for i, (pos, mass) in enumerate(zip(positions, masses)):
        if pos[2] > com_z:
            atoms_above.append(i + 1)
    
    # Save to file
    with open(output_file, 'w') as f:
        for idx in atoms_above:
            f.write(f"{idx}\n")
    
    print(f"Center of mass z-coordinate: {com_z:.6f}")
    print(f"Number of atoms above center of mass: {len(atoms_above)}")
    print(f"Atom indices above center of mass: {atoms_above}")
    print(f"Results saved to {output_file}")

if __name__ == '__main__':
    # Read POSCAR from slab subdirectory and save output to slab folder
    poscar_path = 'slab/POSCAR'
    output_file = 'slab/atoms_above.txt'
    get_atoms_above_center_of_mass(poscar_path, output_file)
