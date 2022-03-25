#!/usr/bin/env python3
import os,sys
from ase import Atoms
import ase
import ase.io 
import sys, os
from  ase.io import read, write
import numpy as np 
import matplotlib.pyplot as plt
from scipy  import interpolate

## Get the energies in the reaction coordinates
dir = sorted([i for i in os.listdir('./') if os.path.isdir(i)])[1:9]
x = [int(i) for i in dir]
y  = []
for i in dir:
    file_in = open('./'+i+'/OUTCAR')
    lines = file_in.readlines()
    for line in lines:
        if '  without' in line:
            E = line.rstrip().split()[-1]
    y.append(float(E))
    # model = ase.io.read('./'+i+'/OUTCAR')
    # E = model.get_total_energy()
    # print(E)
y = [i-y[0] for i in y]

### Plotting part

# 300 represents number of points to make between T.min and T.max
xnew = np.linspace(min(x), max(x), 600)  
## M1
# y_smooth = interpolate.splrep(x, y, s=0)
# ynew = interpolate.splev(xnew, y_smooth, der=0)

# M2
xus = interpolate.InterpolatedUnivariateSpline(x,y)
ynew = xus(xnew)

# ## M3 
# rbf = interpolate.Rbf(x,y)
# ynew = rbf(xnew)

# plt.figure()
plt.plot(xnew,ynew)
system =  sys.argv[1]
plt.legend([system])
# plt.plot(xnew,ynew, x, y)
# plt.legend(['Clean', 'origin'])
plt.xlabel('Reaction Coordinates')
plt.ylabel('Potential Energy / eV')
plt.savefig(system+'.png',dpi=500)
plt.show()

