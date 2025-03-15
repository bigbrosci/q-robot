#!/usr/bin/env python3 
from ase.io import read, write
from ase.geometry import get_distances

def modify_structure(poscar_path):
    # Read the POSCAR file
    structure = read(poscar_path)
    
    # Find indices of all H atoms and the single N atom
    h_indices = [atom.index for atom in structure if atom.symbol == 'H']
    n_index = [atom.index for atom in structure if atom.symbol == 'N'][0]  # Assuming there is only one N atom
    
    # Get the position of the N atom
    n_position = structure.positions[n_index:n_index+1]  # Ensure 2D array for distance calculation
    
    # Initialize a list to collect H atoms to be deleted
    h_to_delete = []

    # Calculate distances and identify H atoms to delete
    for h_index in h_indices:
        # Get the position of the H atom in a 2D array
        h_position = structure.positions[h_index:h_index+1]
        # Calculate distance
        distance = structure.get_distance(n_index, h_index, mic=True) 
        # Access the single distance value correctly
        if distance > 1.5:  # Access the first element of the 2D distance array
            h_to_delete.append(h_index)
    
    # Delete H atoms whose distance to N is greater than 1.5 Å
    if h_to_delete:
        structure = structure[[atom.index for atom in structure if atom.index not in h_to_delete]]
    
    # Save the modified structure
    write('POSCAR_deleted', structure, format='vasp', vasp5=True)

    return f"Deleted {len(h_to_delete)} H atoms. Modified structure saved as 'POSCAR_deleted'."

# Example usage
poscar_path = 'POSCAR'  # Specify the path to your POSCAR file
result = modify_structure(poscar_path)
print(result)

