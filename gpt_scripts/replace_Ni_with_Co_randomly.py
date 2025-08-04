from ase.io import read, write
import numpy as np

def replace_random_atoms(filename='POSCAR', n_replace=120, symbol_from='Ni', symbol_to='Co', output='POSCAR_NiCo'):
    # Step 1: Load structure
    atoms = read(filename)
    
    # Step 2: Get indices of atoms to replace
    ni_indices = [i for i, atom in enumerate(atoms) if atom.symbol == symbol_from]
    if len(ni_indices) < n_replace:
        raise ValueError(f"Only {len(ni_indices)} '{symbol_from}' atoms found, but {n_replace} requested.")
    
    replace_indices = np.random.choice(ni_indices, size=n_replace, replace=False)

    # Step 3: Replace symbols
    for i in replace_indices:
        atoms[i].symbol = symbol_to

    # Step 4: Save new structure
    write(output, atoms, format='vasp')
    print(f"✅ {n_replace} '{symbol_from}' atoms replaced with '{symbol_to}'.")
    print(f"✅ New structure saved as {output}")

# Example usage
if __name__ == '__main__':
    replace_random_atoms(
        filename='POSCAR',
        n_replace=120,
        symbol_from='Ni',
        symbol_to='Co',
        output='POSCAR_NiCo'
    )

