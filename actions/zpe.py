#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from outcar import * 

f = open('OUTCAR')
lines_o = f.readlines()
f.close()

dict_incar = get_incar(lines_o)
if dict_incar.get('IBRION') not in ['5', '6', '7', '8']:
    print('This is not frequency calculation!')
    exit()    
        
def check_freq():
    nwrite = dict_incar.get('NWRITE')
    nu_list, zpe_list = get_freq(lines_o)
    #if nwrite == 3:
    #    '''For NWRITE = 3, each frequency is written twice in OUTCAR.   '''
    #    nu_list = list(set(nu_list)) 
    #    zpe_list = list(set(zpe_list))
    nu_list = list(set(nu_list)) 
    zpe_list = list(set(zpe_list))
    return nu_list, zpe_list

nu_list, zpe_list = check_freq()
E_zpe = sum(zpe_list)/2000  # calculate ZPE 1) Sum(hv)/2  2) convert meV to eV
print(E_zpe)
#if os.path.isfile('OUTCAR'):
#    f = open('OUTCAR')
#    lines_o = f.readlines()
#    f.close()
#else:
