#!/usr/bin/env bash 

def check_dimer():
    '''Check the dimer job. If it goes wrong, return G or B, G and B stand  for Good and  Bad '''
    dimer_job = 'N'
    curv = []
    good_or_bad = 'G'
    with open('INCAR', 'r') as incar_in: # Check if it is a dimer job or not.
        ibrion = 2
        dimer_job = 'N'
        lines = incar_in.readlines()
        for line in lines:
            if 'IBRION' in line:
                ibrion = str(line.rstrip().split()[-1])
        if ibrion == 44:
            dimer_job = 'Y' 
    
    if dimer_job == 'Y':  # Check the curvature in the dimer job. if the curvature
        with open('OUTCAR','r') as file_in:
            lines = file_in.readlines()
            curv = []
            for line in lines:
                if 'curvature along the dimer direction' in line:
                    infor = line.rstrip().split()[-1]
                    if '*' in infor:
                        good_or_bad = 'B'
                        break 
                    else:
                        curv.append(line.rstrip().split()[-1])
    
        curv = [float(i) for i in curv] 
    
    if min(curv) <= -30:
        good_or_bad = 'B'
    if curv[-1] > 1 :
        good_or_bad = 'B'
    for cur in curv:
        if cur > 10:
            good_or_bad = 'B'
    return good_or_bad
