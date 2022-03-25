#!/usr/bin/env python3
from job import * 
from job import check_one_job, get_path
import sys 

if len(sys.argv[:]) == 1:
    path = '.'
    try: 
        check_one_job(path)
    except IOError:
        print('\nPlease make sure that the current directory has vasp input and output files\n')
elif len(sys.argv[:]) == 2:
    '''Can be used to get the path or check job convergence '''
    argu_in = sys.argv[1]
    if argu_in.isdigit() :
        if int(argu_in) > 100000:
            ''' Assume that you want to get job path '''
            print(get_path(argu_in))
        if int(argu_in) in [24, 48, 72, 96, 120, 144, 168, 192]:
            ''' Assume that you want to creat the vasp_run_script in bsc'''
            write_script_bsc(argu_in)
    else:
        '''Assume that you want to check the calculation in the subfolders'''
        try:
            '''Here, argu_in is the path you want to check your jobs.'''
            check_one_job(argu_in)
        except IOError:
            print('\nPlease make sure that %s has vasp input and output files\n' %(sys.argv[1]))    
            
elif len(sys.argv[:]) == 3:            
    script, argu_in_1, argu_in_2 = sys.argv[:]
    if argu_in_1.isdigit():
        write_script_bsc(argu_in_1, argu_in_2)
        if  argu_in_2.isdigit(): 
            argu_in_2 = int(argu_in_2)
            write_script_tekla(argu_in_1, argu_in_2)  
    else:
        write_script_tekla(argu_in_1, argu_in_2)  

elif len(sys.argv[:]) == 4:            
    queue, ncores, job = sys.argv[1:]
    ncores = int(ncores)
    write_script_tekla(queue, ncores, job)              
#    write_script_palma(ncores, job)
#    write_script_bsc(ncores, job)
