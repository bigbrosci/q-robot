#!/usr/bin/env python3 
from qerobot import *
import sys
log = sys.argv[1]
dict_log = read_qeout(log)
print(dict_log.get('potentialenergy'), '\t eV')
