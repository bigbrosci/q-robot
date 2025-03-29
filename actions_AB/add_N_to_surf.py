#!/usr/bin/env python3
"""
Merged Adsorption Candidate Generation and Simplification Script

This script performs the following:
  1) Identifies all hollow sites on the top layer of a 4-layer bimetallic (111) slab.
  2) For each candidate hollow site, it:
       - Places an N atom on a copy of the original slab.
       - Analyzes the candidate site by computing the minimum in-plane (xy) distances
         from the N position to the atoms in the second and third layers. The site is
         classified as "hcp" if the minimum distance to any second-layer atom is below a
         tolerance (0.8 Å) or "fcc" if the minimum distance to any third-layer atom is below
         that tolerance. If both apply, the smaller distance determines the type.
       - Computes a local composition label from the three nearest top-layer metal atoms.
       - Determines the underlying metal atom (from the third layer for fcc or second for hcp)
         that is closest (in xy) to the N atom.
  3) For each accepted candidate (non-duplicate):
       - Creates (if needed) a folder named "<comp_label>_<site_type>_<underlying>" (with a suffix if more than one candidate falls in that group).
       - Saves the unsimplified candidate structure in that folder as "POSCAR".
       - Also creates a simplified version—by extracting only metal atoms above a threshold (threshold = (z_max/3)+1.0, where z_max is computed from metal atoms only, ignoring N), applying a FixAtoms constraint to metal atoms, and reducing the cell’s vacuum in z by 4 Å (without scaling atomic positions)—and saves it as "POSCAR_sim".
  4) Duplicate candidates (determined by an N-atom position difference less than 0.2 Å) are saved only in a folder named "hold" (using the same two-file naming scheme).

Usage:
  python merged_ads_simplify.py --input POSCAR --nperlayer 16

Adjust parameters as needed.
"""

import sys, os, argparse, numpy as np
from ase import io, Atom
from ase.constraints import FixAtoms
from scipy.spatial import Voronoi, ConvexHull

# ----------------------------
# Candidate Generation Functions
# ----------------------------

def get_layers(atoms, tol=0.5):
    """
    Identify the top, second, and third layers by sorting atoms by z (highest first)
    and grouping atoms that are within 'tol' of each other.
    Returns: top_layer, second_layer, third_layer (each as an Atoms object).
    """
    positions = atoms.get_positions()
    z_coords = positions[:, 2]
    sorted_indices = np.argsort(z_coords)[::-1]
    layers = []
    current_layer = []
    current_ref = z_coords[sorted_indices[0]]
    for idx in sorted_indices:
        if abs(z_coords[idx] - current_ref) < tol:
            current_layer.append(idx)
        else:
            layers.append(current_layer)
            current_layer = [idx]
            current_ref = z_coords[idx]
    layers.append(current_layer)
    if len(layers) < 3:
        sys.exit("Not enough layers found. Need at least 3 layers.")
    top_layer = atoms[layers[0]]
    second_layer = atoms[layers[1]]
    third_layer = atoms[layers[2]]
    return top_layer, second_layer, third_layer

def get_hollow_sites(top_layer):
    """
    Generate candidate adsorption sites via Voronoi tessellation on the x-y coordinates of the top layer.
    Returns a list of candidate 3D positions (numpy arrays) with the z coordinate set to (z_top + 1.3 Å).
    """
    positions = top_layer.get_positions()[:, :2]
    vor = Voronoi(positions)
    candidate_sites_2d = vor.vertices
    hull = ConvexHull(positions)
    inside_sites = []
    for v in candidate_sites_2d:
        if all(np.dot(eq[:-1], v) + eq[-1] <= 1e-6 for eq in hull.equations):
            inside_sites.append(v)
    z_top = np.max(top_layer.get_positions()[:, 2])
    candidate_sites = [np.array([v[0], v[1], z_top + 1.3]) for v in inside_sites]
    return candidate_sites

def get_site_neighbors_index(top_layer, candidate, cutoff=3.0):
    """
    For a given candidate hollow site (x-y), return indices of top-layer atoms within a cutoff.
    """
    positions = top_layer.get_positions()[:, :2]
    candidate_xy = candidate[:2]
    dists = np.linalg.norm(positions - candidate_xy, axis=1)
    indices = np.where(dists < cutoff)[0]
    return indices

def shift_slab_if_boundary(atoms, neighbor_indices, boundary_tol=0.2):
    """
    If any neighbor (from neighbor_indices) has fractional x-y coordinates near 0 or 1,
    shift the slab so that their average fractional coordinate becomes 0.5.
    """
    cell = atoms.get_cell()
    frac_coords = atoms.get_scaled_positions()
    neighbor_frac = frac_coords[neighbor_indices][:, :2]
    if (np.any(neighbor_frac < boundary_tol) or np.any(neighbor_frac > (1-boundary_tol))):
        avg_frac = np.mean(neighbor_frac, axis=0)
        delta_frac = np.array([0.5, 0.5]) - avg_frac
        shift_cart = np.dot(delta_frac, cell)
        atoms.translate(shift_cart)
    return atoms

def count_site_composition(top_layer, neighbor_indices):
    """
    Count the chemical species among the three top-layer metal atoms forming the adsorption site.
    Returns a label string (e.g., "Ag2Au1") sorted alphabetically.
    """
    symbols = top_layer.get_chemical_symbols()
    comp = {}
    for i in neighbor_indices:
        sym = symbols[i]
        comp[sym] = comp.get(sym, 0) + 1
    label = "".join(f"{k}{comp[k]}" for k in sorted(comp.keys()))
    return label

def analyze_site(n_position, second_layer, third_layer, tol=0.8):
    """
    Analyze the N adsorption site by computing the minimum in-plane (x-y) distance
    from the N position to atoms in the second and third layers.
    Returns:
      - "hcp" if min distance to any second-layer atom is less than tol,
      - "fcc" if min distance to any third-layer atom is less than tol,
      - If both conditions are met, returns the type corresponding to the smaller distance.
      - Returns None if neither condition is met.
    """
    n_xy = np.array(n_position)[:2]
    second_positions = second_layer.get_positions()[:, :2]
    third_positions = third_layer.get_positions()[:, :2]
    d2 = np.min(np.linalg.norm(second_positions - n_xy, axis=1))
    d3 = np.min(np.linalg.norm(third_positions - n_xy, axis=1))
    if d2 < tol and d3 < tol:
        return "hcp" if d2 < d3 else "fcc"
    elif d2 < tol:
        return "hcp"
    elif d3 < tol:
        return "fcc"
    else:
        return None

def get_underlying_atom(n_position, second_layer, third_layer, site_type):
    """
    For a given N position and site type, return the chemical symbol (e.g., "Ag" or "Au")
    of the metal atom in the target layer (third for fcc, second for hcp) that is closest (x-y) to N.
    """
    n_xy = np.array(n_position)[:2]
    layer = third_layer if site_type.lower() == "fcc" else second_layer
    positions = layer.get_positions()[:, :2]
    idx = np.argmin(np.linalg.norm(positions - n_xy, axis=1))
    return layer.get_chemical_symbols()[idx]

# ----------------------------
# Simplification Function
# ----------------------------

def simplify_structure(candidate_atoms):
    """
    Simplify a candidate structure by:
      - Extracting only metal atoms above a threshold and including any N atoms.
      - The threshold is defined as: threshold = (z_max / 3) + 1.0, where z_max is the maximum metal z (ignoring N).
      - Applying a FixAtoms constraint to the metal atoms.
      - Reducing the cell's z-dimension by 4 Å without altering atomic positions.
    Returns a new Atoms object (the simplified structure).
    """
    metal_indices = [i for i, atom in enumerate(candidate_atoms) if atom.symbol.upper() != "N"]
    if not metal_indices:
        sys.exit("No metal atoms found in candidate structure during simplification!")
    metal_z = np.array([candidate_atoms[i].position[2] for i in metal_indices])
    z_max = np.max(metal_z)
    threshold = (z_max / 3.0) + 1.0
    top_metal_indices = [i for i in metal_indices if candidate_atoms[i].position[2] >= threshold]
    n_indices = [i for i, atom in enumerate(candidate_atoms) if atom.symbol.upper() == "N"]
    new_indices = sorted(set(top_metal_indices) | set(n_indices))
    simplified = candidate_atoms[new_indices]
    fixed_indices = [i for i, atom in enumerate(simplified) if atom.symbol.upper() != "N"]
    if fixed_indices:
        simplified.set_constraint(FixAtoms(indices=fixed_indices))
    cell = simplified.get_cell().copy()
    if cell[2,2] < 4.0:
        sys.exit("Cell z-dimension too small to reduce vacuum by 4 Å.")
    cell[2,2] -= 4.0
    simplified.set_cell(cell, scale_atoms=False)
    return simplified

# ----------------------------
# Main Program
# ----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Merged adsorption candidate generation and simplification script."
    )
    parser.add_argument("--input", type=str, required=True,
                        help="Input slab file (e.g., POSCAR)")
    parser.add_argument("--nperlayer", type=int, default=16,
                        help="Number of metal atoms per layer (default: 16 for a (4x4) surface)")
    args = parser.parse_args()

    # Read the original slab.
    slab = io.read(args.input)

    # Identify layers (top, second, third).
    top_layer, second_layer, third_layer = get_layers(slab, tol=0.5)

    # Generate candidate hollow sites from the top layer.
    candidate_sites = get_hollow_sites(top_layer)
    if not candidate_sites:
        sys.exit("No candidate hollow sites found on the top layer.")
    print(f"Found {len(candidate_sites)} candidate hollow sites.")

    accepted_n_positions = []
    dup_tol = 0.2  # duplicate tolerance in Å
    accepted_folder_counts = {}  # folder name -> count
    duplicate_global_count = 0
    os.makedirs("hold", exist_ok=True)

    for i, candidate in enumerate(candidate_sites):
        # Create candidate structure by adding an N atom at the candidate site.
        candidate_structure = slab.copy()
        n_atom = Atom("N", position=candidate.tolist())
        candidate_structure.append(n_atom)

        # Analyze candidate site.
        site_type = analyze_site(candidate, second_layer, third_layer, tol=1.0)
        if site_type is None:
            print(f"Candidate site {i} did not match fcc or hcp criteria; skipping.")
            continue

        # Identify the three nearest top-layer atoms for composition labeling.
        neighbor_indices = get_site_neighbors_index(top_layer, candidate, cutoff=3.0)
        if len(neighbor_indices) < 3:
            print(f"Candidate site {i} skipped: fewer than 3 neighbors found.")
            continue
        if len(neighbor_indices) > 3:
            positions = top_layer.get_positions()[:, :2]
            dists = np.linalg.norm(positions[neighbor_indices] - candidate[:2], axis=1)
            sorted_idx = np.argsort(dists)[:3]
            neighbor_indices = np.array(neighbor_indices)[sorted_idx]

        # Shift candidate structure if needed.
        site_symbols = [top_layer.get_chemical_symbols()[j] for j in neighbor_indices]
        if "B" in site_symbols:
            candidate_structure = shift_slab_if_boundary(candidate_structure, neighbor_indices, boundary_tol=0.2)

        # Determine composition label.
        comp_label = count_site_composition(top_layer, neighbor_indices)
        # Determine underlying metal atom.
        underlying = get_underlying_atom(candidate, second_layer, third_layer, site_type)
        candidate_name = f"POSCAR_{comp_label}_{site_type}_{underlying}"

        # Check for duplicate N positions.
        duplicate_found = False
        for pos in accepted_n_positions:
            if np.linalg.norm(candidate - pos) < dup_tol:
                duplicate_found = True
                break

        if duplicate_found:
            duplicate_global_count += 1
            dup_folder = os.path.join("hold", f"{candidate_name}_{duplicate_global_count}")
            os.makedirs(dup_folder, exist_ok=True)
            # Save unsimplified candidate.
            dup_unsimp = os.path.join(dup_folder, "POSCAR")
            io.write(dup_unsimp, candidate_structure)
            # Save simplified candidate.
            dup_simp = os.path.join(dup_folder, "POSCAR_sim")
            io.write(dup_simp, simplify_structure(candidate_structure))
#            print(f"Candidate site {i} is a duplicate; saved in 'hold' as {dup_folder}.")
            continue
        else:
            accepted_n_positions.append(candidate.copy())
            folder_name = candidate_name.replace('POSCAR_', '')
            os.makedirs(folder_name, exist_ok=True)
            if folder_name not in accepted_folder_counts:
                accepted_folder_counts[folder_name] = 1
                unsimp_path = os.path.join(folder_name, "POSCAR")
                simp_path = os.path.join(folder_name, "POSCAR_sim")
            else:
                accepted_folder_counts[folder_name] += 1
                unsimp_path = os.path.join(folder_name, f"POSCAR_{accepted_folder_counts[folder_name]}")
                simp_path = os.path.join(folder_name, f"POSCAR_sim_{accepted_folder_counts[folder_name]}")
            io.write(unsimp_path, candidate_structure)
            io.write(simp_path, simplify_structure(candidate_structure))
            print(f"Candidate site {i} accepted; saved unsimplified as {unsimp_path} and simplified as {simp_path}.")

if __name__ == "__main__":
    main()

