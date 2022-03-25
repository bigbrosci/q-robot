#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:07:52 2020

@author: qli
"""

from subprocess import Popen, PIPE
import os
import sys 

def get_id_slurm(user_id):
    '''Get the job list in the queue ''' 
    list_j = [] # list for the job_ID
    process = Popen(['squeue', '-lu',  user_id], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    list_out =  stdout.split('\n')[2:]
    for i in range(0, len(list_out)-1):
        list_j.append(list_out[i].split()[0])
    return list_j

def get_dir_slurm(job_id):
    '''get the job_path of one job'''
    job_dir = None
    command = 'scontrol show job ' + job_id
    process1 = Popen(command, shell = True,  stdout=PIPE, stderr=PIPE)
    stdout1, stderr1 = process1.communicate()
    for i in stdout1.split('\n'):
        if 'WorkDir' in i:
            #job_dir.append(i.split('=')[1])
            job_dir = i.split('=')[1]
    return job_dir


def update_record(job_ID):
    '''1) update the job information to a file and 2) print the job and its path'''
    report = os.path.expanduser('~') + '/bin/job_list.txt'
    id_dict = {}
    with open(report, 'a+') as file_in:
        lines = file_in.readlines()
        for line in lines:
            infor = line.rstrip().split(':')
            key = infor[0]
            value = infor[1]
            id_dict[key] = value
        for i in list_j:
            if i not in id_dict.keys():
                new_dir = get_dir_slurm(i)
                id_dict.update({i:new_dir})
                file_in.write('%s:%s\n' %(i, new_dir))
    print(job_ID, '\t', id_dict.get(job_ID))

job_ID = sys.argv[1]
update_record(job_ID)    