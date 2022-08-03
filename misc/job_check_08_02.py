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
id_file = path +'/job_ID.csv'
job_done_file = path +  '/job_done.csv'
job_problem_file = path + '/job_problem.csv'

def check_files():
    '''create folder named: job_infor, and create files in it: job_ID.csv, job_done.csv, job_problem.csv'''
    if not os.path.isdir(path):
       os.mkdir(path)
       open(id_file,'w').close()
       open(job_done_file,'w').close()
       open(job_problem_file,'w').close()
    else:
       if not os.path.exists(id_file):
           open(id_file,'w').close()
       if not os.path.exists(job_done_file):
           open(job_done_file,'w').close()
       if not os.path.exists(job_problem_file):
           open(job_problem_file,'w').close()

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
    command = 'scontrol show job ' + job_id
    process1 = Popen(command, shell = True,  stdout=PIPE, stderr=PIPE)
    stdout1, stderr1 = process1.communicate()
    for i in stdout1.decode('utf8').split('\n'):
        if 'WorkDir' in i:
            job_dir = i.split('=')[1]
    return job_dir 

def get_recorded_id():
    '''Get a dictionary contains the ID and path infromation '''
    id_dict = {}
    with open(id_file, 'r') as file_in:
        lines = file_in.readlines()
        for line in lines:
            infor = line.rstrip().split(',')
#            key = infor[0] ,  value = infor[1] ,  id_dict[key] = value
            id_dict[infor[0]] = infor[1]
    return id_dict


def update_record():
    '''update the job information to a file and save it to ~/job_infor/job_ID.csv'''
    id_dict = get_recorded_id() # the recorded job ID list
    list_j = get_job_list()     # All jobs in the queue
    with open(id_file, 'a+') as file_in:
        for i in list_j:
            if i not in id_dict.keys(): # if the job ID is not recoorded in job_ID.csv file
                new_dir = get_dir_slurm(i)
                file_in.write('%s,%s\n' %(i, new_dir))

def check_ibrion_dimer_from_INCAR(path):
    '''Check the dimer job by reading INCAR file in the path. Return string 'Y' or 'N' which stand  for Yes and  No '''
    ibrion = 2
    dimer_job = 'N'  
    with open(path + '/INCAR', 'r') as incar_in: # Check if it is a dimer job or not. It can also be done via reading OUTCAR
        ibrion = 2
        dimer_job = 'N'
        lines = incar_in.readlines()
        for line in lines:
            if 'IBRION' in line:
                ibrion = int(line.rstrip().split()[-1])
        if ibrion == 44:
            dimer_job = 'Y'
    return dimer_job


def check_ibrion_dimer_from_OUTCAR(path):
    '''Check the improved dimer job by reading OUTCAR in the path. Return string 'Y' or 'N' which stand  for Yes and  No  '''
    dimer_job = 'N'  #     
    outcar = path + '/OUTCAR'
    read_out = open(outcar, 'r')
    lines_out = read_out.readlines()
    read_out.close()    

    # Check if the job is Improved Dimer or not
    ibrion = 2
    dimer_job = 'N'    
    
    for line in lines_out:
        if 'IBRION' in line:
            ibrion = int(line.rstrip().split()[2])
    if ibrion == 44:
        dimer_job = 'Y'
        
    return dimer_job
    
def check_dimer(path):
    '''Check Dimer calculation, Return string 'G' or 'B' which stand  for Good and  Bad '''
    
    curv = [0] #curvature list
    good_or_bad = 'G'
    dimer_job = check_ibrion_dimer_from_OUTCAR(path)
    
    if dimer_job == 'Y':  # Check the curvature in the dimer job. if the curvature
        with open(path + '/OUTCAR','r') as file_in:
            lines = file_in.readlines()
            for line in lines:
                if 'curvature along the dimer direction' in line:
                    if '***' in line:  # if *** in the curvature, the dimer job goes wrong absolutely.
                        good_or_bad = 'B'
                        break
                    else:
                        infor = line.rstrip().split()[-1]
                        curv.append(float(infor))

            if min(curv) <= -30:  # Too negative values might mean that the job is going wrong.
                good_or_bad = 'B'
            if curv[-1] > 1 :     # Good TS always have nagative curvature
                good_or_bad = 'B'
            for cur in curv:
                if cur > 10:      # once a positive curvature is found, the job goes wrong.
                    good_or_bad = 'B'
    return dimer_job, good_or_bad



### save functions.
def save_jp(job_id, job_dir, id_file):
    jp = open(job_problem_file, 'a+')
    '''jp: job with problems, remove the id from id file and save it to job_problem file '''
    jp.write('%s,%s\n' %(job_id, job_dir))
    subprocess.call(['sed -i "/' + job_id + '/d" ' + id_file], shell = True)
    print('%s has problem\n \t %s' %(job_id,job_dir))
    jp.close()

def save_jd(job_id, job_dir, id_file):
    jd = open(job_done_file, 'a+')
    jd.write('%s,%s\n' %(job_id, job_dir))
    subprocess.call(['sed -i "/' + job_id + '/d" ' + id_file], shell = True)
    print('%s is done\n \t %s' %(job_id,job_dir))
    jd.close()
  
def get_nsw_from_OUTCAR(outcar):
    '''Get the NSW in the INCAR file and ionic steps finished'''
    NSW_in = 200  #NSW in the INCAR file
    NSW_out = 0   #Ionic steps finished in the OUTCAR
    with open(outcar, 'r') as f_out:
        lines = f_out.readlines() 
        for line in lines:
            if  'NSW' in line:
                infor = line.rstrip().split()
                NSW_in = int(infor[2])
            if 'Iteration' in line:
                infor = line.rstrip().replace('-', '').replace('(', ' ').split()
                NSW_out = int(infor[1])
    return NSW_in, NSW_out

def check_jobs():
    '''check the job_ID in the states.csv, if it not in the squeque, analyze the outputs and do something.
    1) if the job is done, save the job information to job_done.csv file
    2) if the job is cancelled by any reasons:
    2.1 if CONTCAR is empty, submit the job directly. (In most cases, the job is preempreted.)
    2.2 if the CONTCAR is not empty, backup the calculation and resubmut the job.
     '''
    id_dict = get_recorded_id()
    list_j = get_job_list()
          
    def resubmit(states): 
        ''' Resubit the jobs'''
        if states == 0: # job is killed/terminated in less than 1 ionic step
            print('resubmit', job_dir)
        else:
            print('backup and resubmit', job_dir) 
            subprocess.call(['save_calculations.sh '], shell = True)  # back up the calculations using a bash script.
        subprocess.call(['sbatch run_vasp_cav'], shell = True)
        subprocess.call(['sed -i "/' + job_id + '/d" ' + id_file], shell = True)  # remove the old job information from the id_file
         
            
    def check_outcar(job_id, job_dir, id_file):
        os.chdir(job_dir)
        NSW_in, NSW_out  = get_nsw_from_OUTCAR('./OUTCAR')
        states = 0  # if job is termintated in less than (or equal) 1 ionic step, states = 0 ; for other cases: states = 1.
        with open('OUTCAR', 'r') as f_out:
            lines = f_out.readlines() 
            
            if NSW_in > NSW_out and  'Voluntary'  in lines[-1]: # job is done with NSW smaller than the values in INCAR file
                save_jd(job_id, job_dir, id_file)
           
            else: # Job is terninated or aborted by some reasons
                if NSW_out > 1 :  # less than (equal) 1 ionic step.
                    states = 1
                resubmit(states)
        os.chdir(path)


    for job_id, job_dir in id_dict.items():
        outcar = job_dir + '/OUTCAR'
        if job_id not in list_j: 
            '''jobs not in the queue (finished, terminated/cancelled due to wrong inputs, preemption, etc).'''
            if os.path.exists(outcar):
                check_outcar(job_id, job_dir, id_file)
            else: # OUTCAR is not in the path, probably it has already been  checked and removed manually.
                save_jp(job_id, job_dir, id_file) 
        # if job_id in list_j: # jobs in the queue
        else:
            if os.path.exists(outcar):
                NSW_in, NSW_out = get_nsw_from_OUTCAR(outcar)
                if NSW_out >= 1:
                # if os.path.exists(outcar) and  os.path.getsize(contcar) > 0: # the CONTCAR size is used to make sure that more than (equal) 1 ionic step is finished.
                    dimer_check = check_dimer(job_dir) # Check Improved dimer results, if it goes wrong, kill it
                    if dimer_check[0] == 'Y': 
                        if dimer_check[1] == 'B' : # Do something if dimer goes bad.
                            save_jp(job_id, job_dir, id_file)
                            subprocess.call(['scancel ' + job_id], shell = True)
                        else: # The dimer job in running and the curvature looks fine, do nothing.
                            pass
            else: # OUTCAR if not found in the path. The job has not been started yet because the job is still in the queue. 
                pass

check_files()
update_record()    
check_jobs()
update_record()    
