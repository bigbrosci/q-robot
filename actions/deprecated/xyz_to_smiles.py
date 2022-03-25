#!/usr/bin/env python3
import os, sys
import openbabel, pybel
from pybel import *
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from xyz2mol import *
from smiles2xyz import *  

xyz_file = sys.argv[1]

smiles, formula  = get_smiles_m1(xyz_file)
if '=' in smiles and ']' in smiles:
    smiles, formula = get_smiles_m2(xyz_file)

print(xyz_file, '\t', formula, '\t', smiles, smiles_file_name(smiles))
    
