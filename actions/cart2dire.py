#!/usr/bin/env python3
import numpy as np 
from  math import cos, sin, sqrt
from lattice import get_abc, read_car
lines = read_car('POSCAR')[0]
#vector = get_vectors(lines)[0]
#list_c = [8.6564608095, 1.8879645958, 0.7298792134]
#list_array = np.array(list_c)
#a = list_array / vector
#print(a)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def get_angle(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def get_vectors(lines):
    '''
    Get the lattice vectors to convert the direct to cartesian
    '''
#    scale = float(lines[1].strip().split()[0]) 
    la1 = np.array([ float(i) for i in lines[2].strip().split() ])
    la2 = np.array([ float(i) for i in lines[3].strip().split() ])
    la3 = np.array([ float(i) for i in lines[4].strip().split() ])
    vector = np.array([la1, la2, la3])
    return vector


va,vb,vc = get_vectors(lines)
a,b,c = get_abc(lines)[0:3]
alpha = get_angle(vc,vb)
beta  = get_angle(vc,va)
gamma = get_angle(va,vb)

omega = a * b * c * sqrt( 1 - (cos(alpha))**2 - (cos(beta))**2 - (cos(gamma))**2 + 2 * cos(alpha) * cos(beta) * cos(gamma) )

v_1 = [ 1 / a , -cos(gamma) / ( a * sin(gamma) ), b * c * ( cos(alpha) * cos(gamma) - cos(beta) ) / (omega * sin(gamma)) ]
v_2 = [ 0,      1 / (b * sin(gamma)),             a * c * ( cos(beta) * cos(gamma) - cos(alpha) ) / (omega * sin(gamma)) ]
v_3 = [ 0,      0,                                a * b * sin(gamma) / omega ]

vector_direct = np.array([v_1, v_2, v_3])
coord = np.array([9.4940864321, 5.4814133571, 0.0000000000])
d,e,f = vector_direct * coord
print(sum(d), sum(e), sum(f))

