import sys, os
import openbabel, pybel
from pybel import *

def get_xyz_from_openbabel(smiles):
    '''use openbabel to convert the SMILES to xyz, mop, and mol files '''
    mol = readstring("smi", smiles)
    mol.make3D(forcefield='mmff94',steps=100)  #print(pybel.forcefields) 
    name = smiles_file_name(smiles) 
    mol.write(format='xyz', filename=name+'.xyz', overwrite=True)
    mol.write(format='mpc', filename=name+'.mop', overwrite=True)
    mol.write(format='mol', filename=name+'.mol', overwrite=True)

smiles = sys.argv[1]
get_xyz_from_openbabel(smiles)
