#!/usr/bin/env python3
import sys 
from gpm import *

log_file = sys.argv[1]
dict_infor = get_infor_from_log(log_file)
for key, val in dict_infor.items():
    if isinstance(val, float):
#    if isinstance(val, (int, float)):
        print(key, val)