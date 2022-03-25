#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To Get xyz file from the log output 
"""
import sys, os
from gpm import * 
from smiles2xyz import * 
from optparse import OptionParser

parser = OptionParser()  
parser.add_option("-f",   
                  type = 'string', dest="file", default="XXX",  
                  help="File to be read")  

parser.add_option("-s", 
                  type = "string", dest="smiles", default="C",
                  help="SMILES to be converted") 

(options, args) = parser.parse_args()
file_in = options.file
smiles  = options.smiles

if os.path.isfile(file_in):
    coords = file_analyzer(file_in)  
    name = file_in.split('.')[0]
else:
#    coords = get_coords_from_chemml(smiles)
    get_coords_from_rdkit(smiles)
    name = smiles_file_name(smiles)
#    get_xyz_from_openbabel(smiles)



save_xyz(coords, name)
#save_poscar(coords, name)

