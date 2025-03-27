import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

# Step 1: Read POSCAR_111
atoms = read("POSCAR_111", format="vasp")

# Step 2: Identify atomic layers with a threshold
z_coords = atoms.positions[:, 2]
z_threshold = 0.5  # Adjust if needed to distinguish layers properly

# Group atoms into layers based on z-coordinates
unique_z = []
for z in sorted(z_coords):
    if not unique_z or abs(z - unique_z[-1]) > z_threshold:
        unique_z.append(z)
unique_z = np.array(unique_z)  # Convert to array for indexing

if len(unique_z) < 3:
    raise ValueError("Not enough distinct layers in the slab to remove top and bottom layers.")

# Step 3: Remove only the topmost and bottommost layers
middle_layers = unique_z[1:-1]  # Exclude the first (bottom) and last (top) layers

# Select atoms that belong to the middle layers
middle_indices = [i for i, atom in enumerate(atoms) if any(abs(atom.position[2] - z) < z_threshold for z in middle_layers)]
atoms_middle = atoms[middle_indices]

# Step 4: Save the modified structure
write("POSCAR_mid", atoms_middle, format="vasp")
print(f"Saved POSCAR_mid with {len(atoms_middle)} atoms, removing the top and bottom layers.")

atoms_new = atoms_middle

# Step 2: Identify atomic layers with a threshold
z_coords = atoms_new.positions[:, 2]
z_threshold = 0.5  # Adjust if needed to distinguish layers properly

# Group atoms into layers based on z-coordinates
unique_z = []
for z in sorted(z_coords):
    if not unique_z or abs(z - unique_z[-1]) > z_threshold:
        unique_z.append(z)
unique_z = np.array(unique_z)  # Convert to array for indexing

if len(unique_z) < 4:
    raise ValueError("Not enough layers in the slab to keep the bottom 4 layers.")

# Step 3: Keep only the bottom 4 layers
bottom_layers = unique_z[:4]  # Select the bottom 4 layers

# Select atoms that belong to the bottom 4 layers
bottom_indices = [i for i, atom in enumerate(atoms_new) if any(abs(atom.position[2] - z) < z_threshold for z in bottom_layers)]
atoms_bottom = atoms_new[bottom_indices]

# Step 5: Fix the bottom 2 layers, relax the top 2 layers
bottom_2_layers = bottom_layers[:2]
fixed_indices = [i for i, atom in enumerate(atoms_bottom) if np.any(np.abs(atoms_bottom.positions[i, 2] - bottom_2_layers) < z_threshold)]

# Debugging: Print the z-coordinates of atoms being fixed
# print("Fixed indices and their z-coordinates:")
# for idx in fixed_indices:
    # print(f"Index: {idx}, z: {atoms_bottom.positions[idx, 2]}")

# Apply the constraint
atoms_bottom.set_constraint(FixAtoms(indices=fixed_indices))

# Step 6: Move the slab so the lowest layer starts at z = 0.1
z_min = np.min(atoms_bottom.positions[:, 2])
atoms_bottom.positions[:, 2] -= z_min - 0.1

# Step 7: Set the vacuum to 15 Å
atoms_bottom.cell[2, 2] = np.max(atoms_bottom.positions[:, 2]) + 15.0

# Step 4: Save the modified structure
write("POSCAR_bottom4", atoms_bottom, format="vasp")
print(f"Saved POSCAR_bottom4 with {len(atoms_bottom)} atoms, keeping only the bottom 4 layers.")
