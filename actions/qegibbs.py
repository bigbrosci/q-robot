#!/usr/bin/env python3 
import sys
from qerobot import * 

infor = sys.argv[:]

qe_out, geometry, SN, SP  = infor[1:]
SN = int(SN)
SP = float(SP)
qe_gibbs = get_gibbs_qe(qe_out, geometry, SN, SP)
print(qe_gibbs)

