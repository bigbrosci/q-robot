#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, math 
import numpy as np 
from read_vasprun import * 
from kpoints import *
from read_vasprun import * 


def read_eigenval_pbe(data_file =  'EIGENVAL'): 
    '''The EIGENVAL file will be used to plot the band structures '''
    f_ei = open(data_file, 'r')
    lines_ei = f_ei.readlines()
    f_ei.close()
    num_elec, num_kp, num_band    = [int(i) for i in lines_ei[5].rstrip().split()]
    return lines_ei, num_kp, num_band    

def read_eigenval_hse():
    '''Read the EIGENVAL file from HSE calculation, find the None-weight 0 points.'''
    lines_ei_hse, num_kp_hse, num_band_hse = read_eigenval_pbe()
    num_to_be_delete  = 0   
    for i in range(0, num_kp_hse):
        if  abs(float(lines_ei_hse[7 + i * (num_band_hse + 2) ].split()[3])) >  0 : 
            num_to_be_delete += 1
    return lines_ei_hse, num_kp_hse, num_band_hse, num_to_be_delete 

def save_eigenval_hse():
    '''    Generate the new HSE-EIGENCAL file
    1) write the header
    2) write the new number line of the KPOINTS 
    3) write the data with weight 0
    4) the output can be read by using read_eigenval_pbe() and write_band_energy() functions for plotting.
     '''
    lines_ei_hse, num_kp_hse, num_band_hse, num_to_be_delete   = read_eigenval_hse()
    num_elec_hse = int(lines_ei_hse[5].rstrip().split()[0])

    file_out = open('EIGENVAL-HSE', 'w')
    # Write Header Part
    for i in range(0,5):
        file_out.write(lines_ei_hse[i].rstrip() + '\n')
    file_out.write('%s  %s %s \n' %(num_elec_hse, num_kp_hse - num_to_be_delete,  num_band_hse ) ) 
    # write the weight 0 part
    start_num = 6 + (num_band_hse + 2) * num_to_be_delete
    file_out.writelines(lines_ei_hse[start_num : ])
    file_out.close()
    print 'New EIGENVAL for HSE is named as EIGENVAL-HSE'

def write_band_energy(data_file =  'EIGENVAL'):
    '''Save the band data for plotting. x-axis are calculated from the kpoints.py part.'''
    f_out = open('band-ei.dat', 'w')
    lines_ei, num_kp, num_band = read_eigenval_pbe(data_file)
    k_d = k_distance()
    e_f = get_fermi()
    for i in range(0, num_kp):
        f_out.write('%s  '  %(k_d[i]))
        for j in range(0, num_band):  
           k_i_band_j = float(lines_ei[8 + j + i * (num_band + 2) ].rstrip().split()[1]) - e_f   
           f_out.write('%s ' %(k_i_band_j))
        f_out.write(' \n')
    f_out.close()        




