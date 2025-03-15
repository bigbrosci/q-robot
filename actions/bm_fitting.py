#!/usr/bin/env python3 
import numpy as np 
import sys 
def bm_fitting(file_in):
    '''a, E: lattice parameters and corresponding energies in data file:
        2.8664, -201.834
        2.8554, -201.812
        .....
        BM state equation: https://en.wikipedia.org/wiki/Birch%E2%80%93Murnaghan_equation_of_state
        fitting: https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
     '''
    a, E  = np.loadtxt(file_in, usecols=(0,1), delimiter=',', unpack = True)
    x  = (a)**(-2)

    c3, c2, c1, c0 = np.polyfit(x, E, 3)
    #print 'The fitted BM equation is:'
    #print ' y = %.4f * (x**3) + %.4f * (x**2) + %.4f * (x) + %.4f' %(c3,c2,c1,c0)

    # Get the results of c_1 + 2c_2x + 3c_3x^2 = 0
    x1 = (math.sqrt(4*c2**2 - 12*c1*c3) - 2*c2)/(6*c3)
    lp = 1/math.sqrt(x1)
    print('The final lattice parameter is: %s  ' %(lp) )

file_in = sys.argv[1]    
bm_fitting(file_in)
