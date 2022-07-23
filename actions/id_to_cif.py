#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 23:13:01 2022

@author: qli
"""
from mp_api import MPRester
import pymatgen as pmg
#from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.inputs import Poscar

mpr = MPRester(api_key="hHyhesqTdnCNr67O79wV1xrKFL7gXe4i")
data = mpr.get_structure_by_material_id("mp-19770")
w = CifWriter(data)
w.write_file('test1.cif')

sd_flags = [[1,1,1]  for i in data.sites]
p = Poscar(data, selective_dynamics= sd_flags)
p.write_file('POSCAR')
