#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from outcar import * 


def check_freq_job():
    freq_yes_or_no = 1
    dict_incar = get_incar()
    if dict_incar.get('IBRION') not in ['5', '6', '7', '8']:
        freq_yes_or_no = 0 
    return freq_yes_or_no

if check_freq_job() == 1:
    nu_list, zpe_list = get_freq()
    nu_list = list(set(nu_list)) # remove the duplicated wavenumbers if NWRITE = 3 
    zpe_list = list(set(zpe_list)) # remove the duplicated energies if NWRITE = 3 
    E_zpe = sum(zpe_list)/2000
    print('ZPE:\t', E_zpe, 'eV')


#def check_freq():
#    nwrite = dict_incar.get('NWRITE')
#    nu_list, zpe_list = get_freq(lines_o)
#    #if nwrite == 3:
#    #    '''For NWRITE = 3, each frequency is written twice in OUTCAR.   '''
#    #    nu_list = list(set(nu_list)) 
#    #    zpe_list = list(set(zpe_list))
#    #print(nu_list, zpe_list)
#    nu_list = list(set(nu_list)) 
#    zpe_list = list(set(zpe_list))
#    #print(nu_list, zpe_list)
#    return nu_list, zpe_list

#nu_list, zpe_list = check_freq()
#E_zpe = sum(zpe_list)/2000  # calculate ZPE 1) Sum(hv)/2  2) convert meV to eV
#print(E_zpe)
#if os.path.isfile('OUTCAR'):
#    f = open('OUTCAR')
#    lines_o = f.readlines()
#    f.close()
#else:
