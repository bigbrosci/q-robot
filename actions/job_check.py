#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:07:52 2020
@author: qli

Change the user_id from "qli" to your own in the get_job_list() function

"""
import subprocess
from subprocess import Popen, PIPE
import os
import sys 

path = os.path.expanduser('~') + '/job_infor'
id_file = path +'/states.csv'
job_done_file = path +  '/job_done.csv'

def check_files():
    if not os.path.isdir(path):
       os.mkdir(path)
       open(id_file,'w').close()
       open(job_done_file,'w').close()
    else:
       if not os.path.exists(id_file):
           open(id_file,'w').close()
       if not os.path.exists(states_file):
           open(job_done_file,'w').close()

def get_job_list(user_id='qli'):
    '''Get the job list in the queue ''' 
    list_j = [] # list of the job_ID
    process = Popen(['squeue', '-lu',  user_id], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    list_out =  stdout.decode('utf8').split('\n')[2:]
    for i in range(0, len(list_out)-1):
        list_j.append(list_out[i].split()[0])
    return list_j

def get_dir_slurm(job_id):
    '''get the job_path of one job'''
    job_dir = None
#    job_state = None
    command = 'scontrol show job ' + job_id
    process1 = Popen(command, shell = True,  stdout=PIPE, stderr=PIPE)
    stdout1, stderr1 = process1.communicate()
    for i in stdout1.decode('utf8').split('\n'):
        if 'WorkDir' in i:
            #job_dir.append(i.split('=')[1])
            job_dir = i.split('=')[1]
#        if 'JobState' in i:
#            job_state = i.split('=')[1]
    return job_dir #, job_state

def get_recorded_id():
    '''Get a dictionary contains the ID and path infromation '''
    id_dict = {}
    with open(id_file, 'r') as file_in:
        lines = file_in.readlines()
        for line in lines:
            infor = line.rstrip().split(',')
            key = infor[0]
            value = infor[1]
            id_dict[key] = value
    return id_dict


def update_record():
    '''1) update the job information to a file and save it to ~/job_infor/states.csv'''
#    id_file = os.path.expanduser('~') + '/job_infor/states.csv'
    id_dict = get_recorded_id() 
    list_j = get_job_list()
    with open(id_file, 'a+') as file_in:
        for i in list_j:
            if i not in id_dict.keys():
                new_dir = get_dir_slurm(i)
                file_in.write('%s,%s\n' %(i, new_dir))

def check_jobs():
    '''check the job_ID in the states.csv, if it not in the squeque, analyze the outputs and do something.
    1) if the job is done, save the job information to job_done.csv file
    2) if the job is cancelled by any reasons:
    2.1 if CONTCAR is empty, submit the job directly. (In most cases, the job is preempreted.)
    2.2 if the CONTCAR is not empty, backup the calculation and resubmut the job.
     '''
    id_dict = get_recorded_id()
    list_j = get_job_list()
    jd = open(job_done_file, 'a+')
    for job_id, job_dir in id_dict.items():
        if job_id not in list_j: # jobs is in the queue or not in the queue(finished, terminated/cancelled due to wrong inputs, preemption, etc)
            outcar = job_dir + '/OUTCAR'
            contcar = job_dir + '/CONTCAR'
            if os.path.exists(outcar):
                with open(outcar) as f_out:
                    if 'Voluntary'  in f_out.readlines()[-1]:
                        jd.write('%s,%s\n' %(job_id, job_dir))
                        subprocess.call(['sed -i "/' + job_id + '/d" ' + id_file], shell = True)
                        print('%s is done\n \t %s' %(job_id,job_dir))
                    else:
                        if os.path.getsize(contcar) == 0: # CONTCAR is empty, job is killed in less than 1 ionic step
                            subprocess.call(['qsub run_vasp_dar'], shell = True)
                        else:
                            subprocess.call(['save_calculations.sh  && qsub run_vasp_dar'], shell = True)
            else: 
                outcar = job_dir + '4_freq/OUTCAR'
                
                if os.path.exists(outcar):
                    if 'Voluntary'  in f_out.readlines()[-1]:
                        jd.write('%s,%s\n' %(job_id, job_dir))
                        subprocess.call(['sed -i "/' + job_id + '/d" ' + id_file], shell = True)
                        print('%s is done\n \t %s' %(job_id,job_dir))
    jd.close()

update_record()    
check_jobs()
update_record()    
