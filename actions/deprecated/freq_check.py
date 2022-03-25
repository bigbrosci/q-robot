#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To check if there is imaginary frequencies in the log file
"""
import sys 
from gpm import *

log_file = sys.argv[1]
dict_log = get_infor_from_log(log_file)
print(dict_log.get('Stability'))
