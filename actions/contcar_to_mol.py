#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 16:46:38 2021
@author: qli
"""
import os, sys
from openbabel import openbabel
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import AllChem

def canonical_hack(mol):
    try:
        new_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)    
        m = Chem.MolFromSmiles(new_smiles)
        new_smiles = Chem.MolToSmiles(m, isomericSmiles=False)    
        return new_smiles
    except:
        return mol 

def get_smiles_from_contar(contcar):
    '''Method using the openbabel and RDkit '''
    mol = next(pybel.readfile('vasp', contcar))
    formula = mol.formula
    smiles = mol.write('smi').split()[0]
    mol = Chem.MolFromSmiles(smiles)
    new_smiles = canonical_hack(mol)
    return new_smiles

#
f_in = sys.argv[1]
full_smiles = get_smiles_from_contar(f_in)
surface_smiles= full_smiles.split('.')[-1]
print('Full_SMILES:\t', full_smiles, '\n'*2)
print('Surface_SMILES:\t', surface_smiles)

