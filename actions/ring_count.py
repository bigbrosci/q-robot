# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 14:12:47 2024

@author: ring test
"""

import numpy as np
from ase import Atoms
from ase.io import read
from ase.neighborlist import NeighborList
from math import degrees, acos

def angle_between(v1, v2):
    # Calculate angle between two vectors in degrees
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = degrees(acos(np.clip(cos_angle, -1.0, 1.0)))
    return angle

def find_hexagonal_rings(atoms, cutoff=1.6, min_angle=110, max_angle=130):
    # Define cutoff radius for neighbor list
    cutoffs = [cutoff] * len(atoms)
    neighbor_list = NeighborList(cutoffs, self_interaction=False, bothways=True)
    neighbor_list.update(atoms)

    hexagonal_rings = set()
    num_atoms = len(atoms)

    for start_atom in range(num_atoms):
        neighbors_A = neighbor_list.get_neighbors(start_atom)[0]
        for B in neighbors_A:
            if B == start_atom:
                continue
            neighbors_B = neighbor_list.get_neighbors(B)[0]
            for C in neighbors_B:
                if C in (start_atom, B):
                    continue
                v1 = atoms[B].position - atoms[start_atom].position
                v2 = atoms[C].position - atoms[B].position
                angle_ABC = angle_between(v1, v2)
                if min_angle <= angle_ABC <= max_angle:
                    neighbors_C = neighbor_list.get_neighbors(C)[0]
                    for D in neighbors_C:
                        if D in (start_atom, B, C):
                            continue
                        v1 = atoms[C].position - atoms[B].position
                        v2 = atoms[D].position - atoms[C].position
                        angle_BCD = angle_between(v1, v2)
                        if min_angle <= angle_BCD <= max_angle:
                            neighbors_D = neighbor_list.get_neighbors(D)[0]
                            for E in neighbors_D:
                                if E in (start_atom, B, C, D):
                                    continue
                                v1 = atoms[D].position - atoms[C].position
                                v2 = atoms[E].position - atoms[D].position
                                angle_CDE = angle_between(v1, v2)
                                if min_angle <= angle_CDE <= max_angle:
                                    neighbors_E = neighbor_list.get_neighbors(E)[0]
                                    for F in neighbors_E:
                                        if F in (start_atom, B, C, D, E):
                                            continue
                                        if F in neighbors_A:
                                            v1 = atoms[E].position - atoms[D].position
                                            v2 = atoms[F].position - atoms[E].position
                                            angle_DEF = angle_between(v1, v2)
                                            if min_angle <= angle_DEF <= max_angle:
                                                v1 = atoms[F].position - atoms[E].position
                                                v2 = atoms[start_atom].position - atoms[F].position
                                                angle_EFA = angle_between(v1, v2)
                                                if min_angle <= angle_EFA <= max_angle:
                                                    ring = sorted([start_atom, B, C, D, E, F])
                                                    hexagonal_rings.add(tuple(ring))

    return [list(ring) for ring in hexagonal_rings]

# Read the XYZ file
atoms = read('./graphene.xyz')

# Find hexagonal rings
hexagonal_rings = find_hexagonal_rings(atoms)

# Sort the rings by the first element of each ring
hexagonal_rings_sorted = sorted(hexagonal_rings, key=lambda x: x[0])

# Print the results
print(f"Number of unique hexagonal rings: {len(hexagonal_rings_sorted)}")
for i, ring in enumerate(hexagonal_rings_sorted):
    print(f"Ring {i+1}: Atoms {ring}")

# Optionally, save the results to a file
with open('./hexagonal_rings.txt', 'w') as f:
    f.write(f"Number of unique hexagonal rings: {len(hexagonal_rings_sorted)}\n")
    for i, ring in enumerate(hexagonal_rings_sorted):
        f.write(f"Ring {i+1}: Atoms {ring}\n")
