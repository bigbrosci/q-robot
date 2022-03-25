#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 08:58:43 2019

@author: qli
"""
import sys 
from gpm import *

coords_file = sys.argv[1]

coords = file_analyzer(coords_file)
n_hbonds1 = calc_hbond_num1(coords)
n_hbonds2 = calc_hbond_num2(coords)
print(n_hbonds1, n_hbonds2)
