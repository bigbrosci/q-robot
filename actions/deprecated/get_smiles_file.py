#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:07:49 2019

@author: qli
"""
import sys 
from smiles2xyz import get_smiles_file
xyz_file = sys.argv[1]
smiles = get_smiles_file(xyz_file)
print(smiles)
