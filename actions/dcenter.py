#!/usr/bin/env python3
#Writen By Qiang 
import numpy as np
import sys
file_in = sys.argv[1]
f_in = open(file_in, 'r')
lines_in = f_in.readlines()
f_in.close()


if len(sys.argv[:]) == 2: 
    file_in = sys.argv[1]
    f_in = open(file_in, 'r')
    lines_in = f_in.readlines()
    f_in.close()
    start = float(lines_in[0].split()[0])
    end   = float(lines_in[-1].split()[0])

elif len(sys.argv[:]) == 4: 
    script, file_in, start, end  = sys.argv
    start = float(start)
    end = float(end)

else:
   print('Command Usage:')
   print('Integerate for the whole range: dcenter.py file')
   print('Integerate for specific  range: dcenter.py file start end')
   exit()
   
def read_dat(file_in):
    f = open(file_in, 'r')
    lines = f.readlines()
    n_column = len(lines[10].strip().split()) 
    ISPIN = 0  
    if n_column == 2:
        ISPIN = 1 
    elif n_column == 3:
        ISPIN = 2
    return ISPIN 

def integer_ele(x,y):
    interval = x[1] - x[0]
    sum_y = y[1] + y[0]
    sum_xy = x[1] * y[1] + x[0] * y[0] 
    ele_lower = 0.5 * interval * sum_y 
    ele_upper = 0.5 * interval * sum_xy 
    return ele_lower, ele_upper 

def integer_array(x,y):
    list_lower = []
    list_upper = []
    for i in range(1, len(x)):
        x_ele = []
        y_ele = []
        x_ele.append(x[i-1])
        x_ele.append(x[i])
        y_ele.append(y[i-1])
        y_ele.append(y[i])
        list_lower.append(integer_ele(x_ele, y_ele)[0])
        list_upper.append(integer_ele(x_ele, y_ele)[1])
    sum_lower = sum(list_lower)  
    sum_upper = sum(list_upper) 
    return  sum_upper / sum_lower, sum_lower 

def get_dat_range(x, start, end):
    l_start = [] 
    l_end = []
    for i, value  in enumerate(x):
        if value >= start:
            l_start.append(i)  # abbriviation for index_start 
        if value >= end:
            l_end.append(i)  # abbriviation for index_start 
    index_s = l_start[0]
    index_e = l_end[0]
    return index_s, index_e   

ISPIN = read_dat(file_in)
x_in = np.loadtxt(file_in,  usecols=0, unpack=True)
index_s, index_e = get_dat_range(x_in, start, end)
x = x_in[index_s:index_e]

if ISPIN == 1: 
    y_in = np.loadtxt(file_in,  usecols=1, unpack=True)
    y = y_in[index_s: index_e]
    d_center = integer_array(x,y)[0]
    num_elec = integer_array(x,y)[1]
    print( 'd-band center is %s' %(d_center)  )
    print( 'electron counting %s' %(num_elec) )
elif ISPIN == 2 : 
    y1_in, y2_in  = np.loadtxt(file_in,  usecols=(1, 2), unpack=True)
    y1 = y1_in[index_s: index_e]
    y2 = y2_in[index_s: index_e]
    d_center1 = integer_array(x, y1)[0]
    d_center2 = integer_array(x, y2)[0]
    num_elec1 = integer_array(x, y1)[1]
    num_elec2 = integer_array(x, y2)[1]
    print(  'd-band center for SPIN-1 is %10.6f ' %(d_center1) )
    print(  'd-band center for SPIN-2 is %10.6f ' %(d_center2) )
    print(  'd-band_average is %10.6f'            %((d_center1 + d_center2) / 2) )
 
    print(  'Electron counting for ISPIN-1 is %10.6f ' %(num_elec1) )
    print(  'Electron counting for ISPIN-2 is %10.6f ' %(num_elec2) )
    print(  'Total Electron is %10.6f ' %((num_elec1 + num_elec2))  )
