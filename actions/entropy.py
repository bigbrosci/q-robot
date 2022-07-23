#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Writen By Qiang Li from ICIQ in 10-03-2018
Updated 07-03-2019
This script is used to calculate
1) the vibration contributions to the entropy under specific temperature
2) zero point energies

Only vibrations are considered,No Electron contribitions.
Translation and Rotation contributions on the surface are simplified.

Ref: Atkin's Physical chemistry, 10th Edition, Statistical thermodynamics: the vibrational Contribution, Page: 642

To use it: python entropy.py 273.15  
273.15 is the Temperature, of course you can use other T values in K)

'''
import sys, os
import math 
from outcar import *
from scipy import constants as con

f = open('OUTCAR')
lines_o = f.readlines()
f.close()

#if len(sys.argv) == 2:
#    Tem = float(sys.argv[1])  # Temperature

if len(sys.argv) == 1:
    Tem = 298.15
#    sum_pf, TS = cal_entropy(nu_list_simplified, Tem)

elif len(sys.argv) == 2:
    Tem = float(sys.argv[1])  # Temperature
#    sum_pf, TS = cal_entropy(nu_list_simplified, Tem)
else:
    print('Command to use: entropy 293  (293 is the temperature in K)')    

h_p = con.h   # 6.62606957E-34 # J*s Plank Constant
k_b = con.k   # 1.38064852E-23 # m²*kg*s⁻²*K⁻¹ Boltzman Constant 
R_gas = con.R # 8.3144598      # J*mol⁻¹*K⁻¹ Gas Constant 
l_s = con.c   # 299792458      # light speed m * s ⁻¹ 
beta = 1/(k_b * Tem)


dict_incar = get_incar()
if dict_incar.get('IBRION') not in ['5', '6', '7', '8']:
    print('This is not frequency calculation!')
    exit()    


def check_freq():
    nwrite = dict_incar.get('NWRITE')
    nu_list, zpe_list = get_freq()
#    if nwrite == 3:
#        '''For NWRITE = 3, each frequency is written twice in OUTCAR.   '''
#        nu_list = list(set(nu_list)) 
#        zpe_list = list(set(zpe_list))
    set(nu_list) 
    set(zpe_list)
    return nu_list, zpe_list

def get_pf(nu): 
    '''calculate partition function for one vibration mode'''
    x_i = h_p * float(nu) * l_s * beta
    pf_l  = x_i / (math.exp(x_i) - 1)   # Left part in the entropy equation 
    pf_r  = math.log(1 - math.exp(-x_i)) 
    pf    = pf_l - pf_r 
    entropy =  R_gas * pf 
    return entropy

def simplify(nu_list):
    '''Generally the modes below 200 are translations and rotations.'''
    reduce_l = nu_list
    for i, nu in enumerate(nu_list):
        if nu <= 200:
            reduce_l[i] = 200
    return reduce_l        

def cal_entropy(nu_list, Tem):  
    '''  Calculate the Entropy S:
        1) i * 100: convert cm-1 to m-1
        2) /1000/96.486 : convert J * K⁻¹ * mol⁻¹ to eV * K⁻¹ from http://web.utk.edu/~rcompton/constants '''
    sum_pf =  sum(get_pf(i*100) for i in nu_list)  / 1000 / 96.485
    TS = Tem * sum_pf  # entropy contribution: T * S     
    return sum_pf, TS

nu_list, zpe_list = check_freq()
E_zpe = sum(zpe_list)/2000  # calculate ZPE 1) Sum(hv)/2  2) convert meV to eV

#nu_list_simplified = simplify(nu_list)
#sum_pf, TS = cal_entropy(nu_list_simplified, Tem)
#### Unsimplified TS
#For  Gas
sum_pf, TS = cal_entropy(nu_list[:-3], Tem)
#For Surface
sum_pf, TS = cal_entropy(nu_list, Tem)

print('Temperature (K):\t%s\tS (eV/K):\t%6.7f\tTS (eV):\t%6.4f \tE_ZPE (eV):\t%6.4f' %(Tem, sum_pf, TS, E_zpe))


