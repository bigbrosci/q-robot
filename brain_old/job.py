#!/usr/bin/env python3
# Write By Qiang on 27-02-2019
from subprocess import Popen, PIPE
from incar import *
import os 
from data import *
   
### Check Job status
# 1) The job is finish or not? 
# 2) is there problem in electronic convergence?

def get_line(file_in, look_up):
    line_num = 0 
    with open(file_in) as data_in:
        lines = data_in.readlines()
        for num, line in enumerate(lines):
            if look_up in line:
                line_num =  num
    return line_num, lines

def cycles_in(path='.'):
    ''' Get the NSW value from INCAR'''
    nelm_in = 60 
    nsw_in = 0 
    file_in = path + '/INCAR'
    with open(file_in) as data_in:
        lines = data_in.readlines()
    for line in lines:
        if 'NSW' in line:
            nsw_in = (line.rstrip().split('=')[1].split()[0])
        elif 'NELM' in line:
            nelm_in = int(line.rstrip().split('=')[1].split()[0])
    return nelm_in, nsw_in 
   
def cycles_osz(path='.'):
    ''' Check how many steps are finished until so far from OSZICAR ''' 
    file_in = path + '/OSZICAR'
    look_up = 'F='
    nelm_osz = 0
    nsw_osz  = 0 
    try: 
        try: 
            with open(file_in) as data_in:
                lines = data_in.readlines()
                for num, line in enumerate(lines):
                    if look_up in line:
                        nelm_osz = int(lines[num-1].rstrip().split()[1])
                        nsw_osz  = int(line.rstrip().split()[0])
        except IndexError:
            print('The first ionic step has not finished, please try later.') 
    except IOError:
        print('The job has not been started.') 
    return nelm_osz, nsw_osz 
   
def converge_and_finish(path='.'):
    file_in = path + '/OUTCAR'
    look_up_1 = 'reached required accuracy'
    look_up_2 = 'Voluntary context '
    converge_or_not = 'No' 
    finish_or_not   = 'No' 
    try:
        with open(file_in) as data_in: 
            for i in data_in.readlines():
                if look_up_1 in i:
                    converge_or_not = 'Yes'
                if look_up_2 in i:
                    finish_or_not   = 'Yes'
    except IOError:
        print('The job has not been started.') 
    check_out = [converge_or_not, finish_or_not]
    return check_out 

def check_one_job(path='.'):  
    '''Check ionic ''' 
    nelm_in,  nsw_in  = cycles_in(path)
    nelm_osz, nsw_osz = cycles_osz(path)
    convergence, finish = converge_and_finish(path)
    print('%s\tNSW_INCAR:\t %s \t NSW_OSZICAR:\t %s \t Converge:\t %s \t Finish:\t %s' %(path, nsw_in, nsw_osz, convergence, finish))
    
    '''Check electronic''' # Only  print the Warning when it seems that there are some problems in electronic calculations
    if nsw_osz >= 5 : 
        if nelm_osz == nelm_in:
            print('Warning: Check your electronic iterations!!!')
            
### Job_ID Control                       
            
def get_id_slurm(user_id):
    list_j = [] # list for the job_ID
    process = Popen(['squeue', '-lu',  user_id], stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    list_out =  stdout.split('\n')[2:]
    for i in range(0, len(list_out)-1):
        list_j.append(list_out[i].split()[0])
    return list_j

def get_dir_slurm(job_id):
    job_dir = None 
    command = 'scontrol show job ' + job_id 
    process1 = Popen(command, shell = True,  stdout=PIPE, stderr=PIPE)
    stdout1, stderr1 = process1.communicate()
    for i in stdout1.split('\n'):
        if 'WorkDir' in i: 
            #job_dir.append(i.split('=')[1])
            job_dir = i.split('=')[1]
    return job_dir

def get_id_qstat(user_id):
    list_j = [] # list for the job_ID
    command = 'qstat -u ' + user_id
    process = Popen(command, shell = True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    list_out =  stdout.split('\n')[2:]
    for i in range(0, len(list_out)-1):
        list_j.append(list_out[i].split()[0])
    return list_j

def get_dir_qstat(job_id):
    job_dir = None 
    command = 'qstat -j ' + job_id 
    process1 = Popen(command, shell = True,  stdout=PIPE, stderr=PIPE)
    stdout1, stderr1 = process1.communicate()
    for i in stdout1.split('\n'):
        if 'workdir' in i: 
            #job_dir.append(i.split('=')[1])
            job_dir = i.split(':')[1]
    return job_dir 

def update_record():
    list_user_tekla = list_id_nl_tekla
    list_user_bsc = dict_id_nl_bsc.values()
    list_user_bsc.append('qiangli')
    list_user_bsc.append('iciq41406') # my user name in La Palma 
    user_id = os.path.expanduser('~').split('/')[-1]    
    list_j = None
    if user_id  in list_user_tekla:        
        list_j = get_id_qstat(user_id)
        get_dir = get_dir_qstat
        report = '/home/' + user_id  + '/bin/Q_robot/reports/job_list.txt'
    elif user_id in list_user_bsc:
        list_j = get_id_slurm(user_id)
        get_dir = get_dir_slurm
        report = '/home/iciq72/' + user_id  + '/bin/Q_robot/reports/job_list.txt'
    with open(report, 'a+') as file_in:
        id_dict = {}
        lines = file_in.readlines()
        for line in lines:
            key = line.split(':')[0].strip()
            value = line.split(':')[1].strip()
            id_dict[key] = value     

        for i in list_j:
            if i not in id_dict.keys():
                id_dict[i] = get_dir(i)
                file_in.write(i + ':' + get_dir(i) + '\n')  
    return id_dict

def get_path(job_id):
    id_dict = update_record()
    possible_list = [int(i) for i in id_dict.keys() if job_id in i]
    if len(possible_list) >= 1:  
        job_id_match  = str(max(possible_list))        
        path = id_dict.get(job_id_match)
    else: 
        path = None
#        path = os.path.abspath('.')
#        print('Can not find the JOB-ID in the database, please check!\nUse the current path instead!')
    return path
    
def get_version_tekla(q_queue):
    vasp_version = '5.4.4'
    vasp_binary = 'vasp_std'
    if q_queue in ['c8m24', 'c4m8']:
        vasp_version = '5.3.3'
        vasp_binary = 'vasp'
    return vasp_version, vasp_binary
    
    
def check_queue(queue, ncores):
    for i in [12, '12', 'c12', 'c12m48', 'c12m48ib']: 
        if queue == i:
            q_queue, NCORE = ['c12m48ib', 6]
            if (ncores % 12) != 0:
                ncores = 12       
    for i in [12, '24', 'c24', 'c24m128', 'c24m128ib']:
        if queue == i:
            q_queue, NCORE = ['c24m128ib', 12]
            if (ncores % 24) != 0:
                ncores = 24    
    for i in [28, '28', 'c28', 'c28m128', 'c28m128ib']:
        if queue == i :
            q_queue, NCORE = ['c28m128ib', 14]
            if (ncores % 28) != 0:
                ncores = 28    
    for i in [8, '8', 'c8', 'c8m24']:
        if queue == i :
            q_queue, NCORE = ['c8m24', 4]
            if (ncores % 8) != 0:
                ncores = 8    
    for i in [4, '4', 'c4', 'c4m8']:
        if queue == i :
            q_queue, NCORE = ['c4m8', 2]
            if (ncores % 4) != 0:
                ncores = 4
#    else:
#        print('Please specify the queue you are going to use')         
    set_ncore(NCORE)
    return q_queue, ncores             

def write_script_tekla(queue, ncores, job='Q-robot'):
    ncores = int(ncores)    
    q_queue, n_core = check_queue(queue, ncores)
    vasp_version, vasp_binary = get_version_tekla(q_queue)
    script = open('run_vasp_tekla', 'w')
    script.write('#!/bin/bash\n')
    script.write('%s\n# Generated by Q-robot \n%s\n'        %('#'*30,'#'*30)) 
    script.write('#$ -N %s \n'                              %(job))  
    script.write('#$ -cwd\n#$ -masterq %s.q\n'              %(q_queue))
    script.write('#$ -pe  %s_mpi  %s\n'                     %(q_queue, n_core))
    script.write('#$ -m ae\n#$ -M qli@iciq.es\n')
    script.write('#$ -o o_$JOB_NAME.$JOB_ID \n#$ -e e_$JOB_NAME.$JOB_ID \n')
    script.write('%s\n# Load Evironment Variables\n%s\n'    %('#'*30,'#'*30))                          
    script.write('. /etc/profile.d/modules.sh\n')
    script.write('module load vasp/%s\n'                    %(vasp_version))
    script.write('%s \n# Running Job\n%s\n'                 %('#'*30, '#'*30))                              
    script.write('export OMP_NUM_THREADS=1\n')
    script.write('echo $JOB_ID:$PWD >>  $HOME/bin/Q_robot/reports/job_list.txt\n')
    script.write('mpirun -np $NSLOTS %s\n' %(vasp_binary))
    script.close()

def write_script_palma(ncores, job='Q-robot'):
    ncores = int(ncores)
    if ncores % 16 != 0:
        ncores = 16
    script = open('run_vasp_palma', 'w')
    script.write('#!/bin/bash\n')
    script.write('%s\n# Generated by Q-robot \n%s\n'        %('#'*30,'#'*30)) 
    script.write('#SBATCH --time=11:29:00 \n')
    script.write('#SBATCH --job-name=%s\n'                  %(job) )
    script.write('#SBATCH --cpus-per-task=1\n#SBATCH --tasks-per-node=16\n') 
    script.write('#SBATCH --ntasks=%s\n'                    %(ncores))   
    script.write('#SBATCH --output=o.%j\n#SBATCH --error=e.%j\n')
    script.write('module purge\nmodule load intel\nmodule load impi/2018.1\n')
    script.write('module load mkl/2018.1\nmodule load vasp/%s\n' %(vasp_version))
    script.write('mpirun %s\n' %(vasp_binary))
    script.close()
    set_ncore('8')  

def write_script_bsc(ncores, job='Q-robot'):
    ncores = int(ncores)
    vasp_version = '5.4.4'
    vasp_binary  = 'vasp_std'
    if ncores % 24 != 0:
        ncores = 24    
    folders = [f for f in os.listdir('.') if os.path.isdir(f)]    
    if len(folders) >=5 :
        vasp_version = '5.4.1'
        vasp_binary = 'vasp_std_vtst'
        images =  len([i  for i in folders if i.isdigit()]) - 2     
        if ncores % 24 != 0:
            ncores = 12 * images
        elif ncores < 6 * images :
            ncores = 12 * images
            
    script = open('run_vasp_bsc', 'w')
    script.write('#!/bin/bash\n')
    script.write('%s\n# Generated by Q-robot \n%s\n'        %('#'*30,'#'*30)) 
    script.write('#SBATCH --time=11:29:00 \n')
    script.write('#SBATCH --job-name=%s\n'                  %(job) )
    script.write('#SBATCH --cpus-per-task=1\n#SBATCH --tasks-per-node=24\n') 
    script.write('#SBATCH --ntasks=%s\n'                    %(ncores))   
    script.write('#SBATCH --output=o.%j\n#SBATCH --error=e.%j\n')
    script.write('#SBATCH --qos=class_a\nmodule load vasp/%s\n' %(vasp_version))
    script.write('srun /apps/VASP/%s/INTEL/IMPI/bin/%s\n' %(vasp_version, vasp_binary))
    script.close()
    set_ncore('12')
