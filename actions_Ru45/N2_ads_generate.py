# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 23:57:20 2025

@author: lqlhz
"""

import sys,os,math
import numpy as np
q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
if q_robot_path not in sys.path:
    sys.path.append(q_robot_path)
from cluster import *

atoms = read('./POSCAR')
add_N2_top_sites(atoms, n_ru_distance=2.0, n_n_distance=1.19)
add_N2_bridge_sites(atoms, n_ru_distance=2.15, n_n_distance=1.19)
add_hollow_sites(atoms, n_ru_distance=1.95, n_n_distance=1.19, n1_height=1.5)
