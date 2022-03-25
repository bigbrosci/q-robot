#!/usr/bin/env python3

import openbabel, pybel
from gpm import *
import sys , os
import subprocess
from subprocess import Popen, PIPE

dire = sys.argv[1]

dict_dir_files = {}
gjf_files = []
log_files = []
for path, subdirs, files in os.walk(dire):
    files_xyz = [i for i in files if '.xyz' in i and 'mop' not in i]
    if len(files_xyz) > 0 :
        dict_dir_files.update({path:files_xyz})


file_check = open('file_check.txt', 'w')        
for key, value in dict_dir_files.items():        
    for xyz_file in value:
        name_old = xyz_file.split('.')[0]
        mop_file = name_old + '_mop.xyz'
        coords_old = file_analyzer(key+'/'+xyz_file)
        coords_new = file_analyzer(key+'/'+mop_file)
        hbond_old = calc_hbond_num1(coords_old)
        hbond_new = calc_hbond_num1(coords_new)
        
        file_check.write('%s\tHbond_old:\t %s \t Hbond_new:\t %s\t' %(key+'/'+xyz_file, hbond_old, hbond_new))
        if hbond_old < hbond_new:
            file_check.write('Check \n')
        else:
            file_check.write('Good \n')

file_check.close()
#name = xyz_file.split('.')[0]

#file_check = open('file_check.txt', 'w')
#
#if xyz_file.count('O') >= 2:     
#    coords = file_analyzer(xyz_file)
#    hbond_num_g09 = calc_hbond_num1(coords)
#
#    ## Generate the mopac input file: XXX.mop  
#    ms = name.count('-r')/2
#    if '-r' in name:
#        line0 = 'PM6-D3H4 GNORM=1.0 SCFCRT=0.0001 EF CHARGE=0  MS=%s \n'  %(str(ms))
#    else:
#        line0 = 'PM6-D3H4 GNORM=1.0 SCFCRT=0.0001 EF CHARGE=0 MS=%s \n'  %(str(ms))
#        
#    mopac_file = open(name+'.mop', 'w')
#    mopac_file.write(line0)
#    mopac_file.write('\n\n')
#    for coord in coords:
#        infor = coord.split()
#        mopac_file.write('%s\t%s\t1\t%s\t1\t%s\t1\t\n' %(infor[0], infor[1], infor[2], infor[3]))
#    mopac_file.close()   
#    
#    ### Run MOPAC    
##    mopac_run = '/opt/mopac/MOPAC2016.exe '+ name + '.mop'
##    subprocess.call(mopac_run, shell=True)
#    
#    mopac_in = name + '.mop'
#    process = subprocess.Popen(['/opt/mopac/MOPAC2016.exe', mopac_in], stdout=PIPE, stderr=PIPE)
#    stdout, stderr = process.communicate()
#    
#    ### Use Openbabel to convert the mopac out to xyz file
#    mols = pybel.readfile('mopout', name+'.out')
#    for mol in mols:
#        mol.write(format='xyz', filename=name+'_mop.xyz', overwrite=True)    
#    
#    ### Calculate the H bond numbers    
#    path = os.getcwd()
#    coords_mop = file_analyzer(path+'/'+name+'_mop.xyz')
#    hbond_num_mop = calc_hbond_num1(coords_mop)
#    
#    ### Compare the old and new H bond numbers
#    file_check.write('%s\tHbond_old:\t %s \t Hbond_new:\t %s\t' %(name, hbond_num_g09, hbond_num_mop))
#    if hbond_num_g09 < hbond_num_mop:
#        file_check.write('Check \n')
#    else:
#        file_check.write('Good \n') 
#        
#    print('%s\tHbond_old:\t %s \t Hbond_new:\t %s' %(name, hbond_num_g09, hbond_num_mop))    
#file_check.close()
