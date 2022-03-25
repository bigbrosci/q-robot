#!/usr/bin/env python3

from smiles2xyz import *
#from chemml.chem import Molecule
import sys 
import os 

ring_out = sys.argv[1]
get_smiles_from_ring(ring_out)
analyze_ring_radical_smiles('rad_ring.out')
analyze_ring_molecule_smiles('mol_ring.out')

#for i in $(cat rad_final_rdkit.out); do python get_final_rdkit.py $i  >> final_rad_rdkit.out ; done
#for i in $(cat mol_final_rdkit.out); do python get_final_rdkit.py $i  >> final_mol_rdkit.out ; done

