#!/usr/bin/env python3
'''
Latest Update: 18-01-2019
Written By Qiang 
'''

import numpy as np 
import math 

#### Part-0 Dictionaries
atomic_mass = dict(H=1.01, He=4.00, Li=6.94, Be=9.01, B=10.81, C=12.01,
                   N=14.01, O=16.00, F=19.00, Ne=20.18, Na=22.99, Mg=24.31,
                   Al=26.98, Si=28.09, P=30.97, S=32.07, Cl=35.45, Ar=39.95,
                   K=39.10, Ca=40.08, Sc=44.96, Ti=47.87, V=50.94, Cr=52.00,
                   Mn=54.94, Fe=55.85, Co=58.93, Ni=58.69, Cu=63.55, Zn=65.39,
                   Ga=69.72, Ge=72.61, As=74.92, Se=78.96, Br=79.90, Kr=83.80,
                   Rb=85.47, Sr=87.62, Y=88.91, Zr=91.22, Nb=92.91, Mo=95.94,
                   Tc=98.00, Ru=101.07, Rh=102.91, Pd=106.42, Ag=107.87,
                   Cd=112.41, In=114.82, Sn=118.71, Sb=121.76, Te=127.60,
                   I=126.90, Xe=131.29, Cs=132.91, Ba=137.33, La=138.91,
                   Ce=140.12, Pr=140.91, Nd=144.24, Pm=145.00, Sm=150.36,
                   Eu=151.96, Gd=157.25, Tb=158.93, Dy=162.50, Ho=164.93,
                   Er=167.26, Tm=168.93, Yb=173.04, Lu=174.97, Hf=178.49,
                   Ta=180.95, W=183.84, Re=186.21, Os=190.23, Ir=192.22,
                   Pt=195.08, Au=196.97, Hg=200.59, Tl=204.38, Pb=207.2,
                   Bi=208.98, Po=209.00, At=210.00, Rn=222.00, Fr=223.00,
                   Ra=226.00, Ac=227.00, Th=232.04, Pa=231.04, U=238.03,
                   Np=237.00, Pu=244.00, Am=243.00, Cm=247.00, Bk=247.00,
                   Cf=251.00, Es=252.00, Fm=257.00, Md=258.00, No=259.00,
                   Lr=262.00, Rf=261.00, Db=262.00, Sg=266.00, Bh=264.00, 
                   Hs=269.00, Mt=268.00)

dict_element = {'1':'H','2':'He','3':'Li','4':'Be','5':'B','6':'C','7':'N','8':'O','9':'F','10':'Ne','11':'Na','12':'Mg','13':'Al','14':'Si','15':'P','16':'S','17':'Cl','18':'Ar','19':'K','20':'Ca','21':'Sc','22':'Ti','23':'V','24':'Cr','25':'Mn','26':'Fe','27':'Co','28':'Ni','29':'Cu','30':'Zn','31':'Ga','32':'Ge','33':'As','34':'Se','35':'Br','36':'Kr','37':'Rb','38':'Sr','39':'Y','40':'Zr','41':'Nb','42':'Mo','43':'Tc','44':'Ru','45':'Rh','46':'Pd','47':'Ag','48':'Cd','49':'In','50':'Sn','51':'Sb','52':'Te','53':'I','54':'Xe','55':'Cs','56':'Ba','57':'La','58':'Ce','59':'Pr','60':'Nd','61':'Pm','62':'Sm','63':'Eu','64':'Gd','65':'Tb','66':'Dy','67':'Ho','68':'Er','69':'Tm','70':'Yb','71':'Lu','72':'Hf','73':'Ta','74':'W','75':'Re','76':'Os','77':'Ir','78':'Pt','79':'Au','80':'Hg','81':'Tl','82':'Pb','83':'Bi','84':'Po','85':'At','86':'Rn','87':'Fr','88':'Ra','89':'Ac','90':'Th','91':'Pa','92':'U','93':'Np','94':'Pu','95':'Am','96':'Cm','97':'Bk','98':'Cf','99':'Es','100':'Fm','101':'Md','102':'No','103':'Lr','104':'Rf','105':'Db','106':'Sg','107':'Bh','108':'Hs','109':'Mt',}
dict_element_2 = {v: k for k, v in dict_element.items()}


def get_dicts(lines):
    '''Return a dictionary from POSCAR file '''
    ele_name  = str(lines[5]).rstrip().split()
    ele_num = [int(i) for i in str(lines[6]).rstrip().split()]
    
    num_ele = 0
    dict_car1 = {}  # key: element-sequence, value: number of atoms
    dict_car2 = {}  # key: element-sequence, value: atom range
    for i, k in enumerate(ele_name):
        name = k + '-' + str(i+1)
        dict_car1.update({name:ele_num[i]})
        list_i = []
        for num in range(1, ele_num[i]+1):
            num_ele +=1 
            list_i.append(num_ele)
        dict_car2.update({name:list_i})   
    return dict_car1, dict_car2

def read_car(file_read):
    '''
    Read POSCAR to get all the lines, element list, atom number list 
    
    '''
    f = open(file_read, 'r')
    lines = f.readlines()
    f.close()
    dict_car1, dict_car2 = get_dicts(lines)
    return lines, dict_car1, dict_car2

def is_direct_or_not(lines):
    is_direct = True
    is_select = True
    if lines[7].strip()[0].upper()  == 'S': # # With Selected T T T, coordination starts from line 9
#        start_num = 9
        if lines[8].strip()[0].upper()  == 'C':
            is_direct = False
    else: 
        is_select = False
#        start_num = 8
        print('----------------------------------------------')
        print( 'Pay Attetion! There is no TTT in the file!   ')
        print( '---------------------------------------------')
        if lines[7].strip()[0].upper()  == 'C':
            is_direct = False 
    return is_direct, is_select

def get_vectors(lines):
    '''
    Get the lattice vectors to convert the direct to cartesian
    '''
#    scale = float(lines[1].strip().split()[0]) 
    
###    Not_used_method.
#    for i in np.arange(2,5):
#       line = [float(i) for i in lines[i].strip().split()]
#       a.append(line[0])
#       b.append(line[1])
#       c.append(line[2])
#    vector = np.array([a,b,c]) 
    
    la1 = np.array([ float(i) for i in lines[2].strip().split() ])
    la2 = np.array([ float(i) for i in lines[3].strip().split() ])
    la3 = np.array([ float(i) for i in lines[4].strip().split() ])
    vector = np.transpose(np.array([la1, la2, la3]))
    
    return vector, la1, la2 

def get_abc(lines):
    '''
    Get the lattice vectors to convert the direct to cartesian
    '''
    scale = float(lines[1].strip().split()[0]) 
    la1 = np.array([ float(i) for i in lines[2].strip().split() ])
    la2 = np.array([ float(i) for i in lines[3].strip().split() ])
    la3 = np.array([ float(i) for i in lines[4].strip().split() ])
    a_length = np.linalg.norm(la1) * scale 
    b_length = np.linalg.norm(la2) * scale 
    c_length = np.linalg.norm(la3) * scale
    A = np.cross(la1, la2)[-1]
    V = np.dot(np.cross(la1, la2), la3)
    return a_length, b_length, c_length, A, V


def get_coordinate(lines, ele_indice):
    '''
    Get the atom xyz coordinate from the POSCAR which is in cartesian format
    '''
    ele_indice = int(ele_indice) + 8
    coordinate =  np.array([float(i) for i in lines[ele_indice].strip().split()[0:3]])
    return coordinate

def determinelayers(lines,threshold=0.5):
    layerscount = {}
#    print(lines)
    line_total = sum([int(i) for i in lines[6].strip().split()]) + 8
    z = [float(lines[i].split()[2]) for i in range(9, line_total+1)]
    seq = sorted(z)
    min = seq[0]
    sets = [min]
    for j in range(len(seq)):
        if abs(seq[j]-min) >= threshold:
            min = seq[j]
            sets.append(min)        
    for i in range(1,len(sets)+1):
        layerscount[i] = []            
        for k in range(len(z)):   
            if abs(z[k]-sets[i-1]) <= threshold:
                layerscount[i].append(k+1)
    return layerscount


#def determinelayers(z_cartesian):
#    seq = sorted(z_cartesian)
#    min = seq[0]
#    layerscount = {}
#    sets = [min]
#    for j in range(len(seq)):
#        if abs(seq[j]-min) >= threshold:
#            min = seq[j]
#            sets.append(min)
#
#    for i in range(1,len(sets)+1):
#        layerscount[i] = []            
#        for k in range(len(z_cartesian)):   
#            if abs(z_cartesian[k]-sets[i-1]) <= threshold:
#                layerscount[i].append(k)
#                
#    return layerscount
    



def get_angle(u,v):
    '''
    Get the angle between two vectors
    '''
    dot_cross = np.dot(u, v)
    c = dot_cross / np.linalg.norm(u) / np.linalg.norm(v)
    angle = np.arccos(np.clip(c, -1, 1))
    return angle

def get_distance(lines, a, b):
    '''
    a and b are the cartesian coordinates of atom_A and atom_B
    This functions solves the periodic probelmes by default.
    '''
    la1, la2 = get_vectors(lines)[1:]
    l_x = abs(la1[0])
    l_y = abs(la2[1])

    dx, dy, dz = a - b 
    if dy > l_y / 2:
        a[1] = a[1] - l_y
        a[0] = a[0] - la2[0]
    if dy <= -l_y /2: 
        a[1] = a[1] + l_y 
        a[0] = a[0] + la2[0]
    
    dx1, dy1, dz = a - b 
    if dx1 > l_x / 2:
        a[0] = a[0] - l_x 
    if dx1 <= - l_x / 2: 
        a[0] = a[0] + l_x 
    distance = np.linalg.norm(a - b)

    return distance 

def get_intact_one_atom(lines, A, B):
    '''atom B is the anchoring point'''
    la1, la2 = get_vectors(lines)[1:]
    l_x = abs(la1[0])
    l_y = abs(la2[1])
    a = get_coordinate(lines, A)
    b = get_coordinate(lines, B)
    dx, dy, dz = a - b 
    if dy > l_y / 2:
        a[1] = a[1] - l_y
        a[0] = a[0] - la2[0]
    if dy <= -l_y /2: 
        a[1] = a[1] + l_y 
        a[0] = a[0] + la2[0]
    
    dx1, dy1, dz = a - b 
    if dx1 > l_x / 2:
        a[0] = a[0] - l_x 
    if dx1 <= - l_x / 2: 
        a[0] = a[0] + l_x 
    return a
    
def get_intact_molecule(lines, atom_list, atom_B):
    coord_list = []
    coord_list = [get_intact_one_atom(lines, i, atom_B) for i in atom_list]
    return coord_list

def get_distance_direct(a,b):
    '''print the distance  '''
    vector_ab = a - b  
    distance = np.linalg.norm(vector_ab)
    return distance

def get_rotate_infor(lines, raw_infor):   
    '''Along the axis(atom_A, atom_B), rotate atoms by theta angle (in degree).
Generall Command: rotate.py atom_A atom_B atom_1 atom_2 ... Angle    
For example: To rotate the C H O atoms 15 degree along the axis formed by No.2 and No.10 atoms.
rotate.py 2 10 C H O 15  
This function is used to analyze the arguments: 2 10 C H O 15  in the command above.
And return the basic information for rotations in the next step.
'''
    point_1, point_2 = raw_infor[0:2]
    point_1 = int(point_1)
    point_2 = int(point_2)
    atom_s = raw_infor[2:-1]
    theta  = float(raw_infor[-1]) * math.pi / 180
    coord_1 = get_coordinate(lines, point_1)
    coord_2 = get_coordinate(lines, point_2)
    atom_list = get_atom_list(lines, atom_s)
            
    B = point_1   ## B atom is used as the point to fix the PBC conditions
    if point_1 in atom_list:
        atom_list.remove(point_1)     
        if point_2 in atom_list:
            atom_list.remove(point_2)
            coord_2 =  get_intact_one_atom(lines, point_2, B) 
    else: 
        if point_2 in atom_list:
            atom_list.remove(point_2)
            B = point_2
        else:
            B = atom_list[-1]
            
    if point_1 == point_2:
        axis = np.array([0, 0, 1])    
    else:
        axis  = coord_2  - coord_1 
        
    return atom_list, coord_1, coord_2, axis, theta, B
            
def rotate_one_atom(infor, coord_3):
    '''Only use this function after your solve the PBC problems '''
    atom_list, coord_1, coord_2, axis, theta, B = infor[:]
    vector_3 = coord_3 - coord_1    
    angle_213= get_angle(vector_3, axis)
    d_13 = get_distance_direct(coord_1, coord_3)
    d_12 = get_distance_direct(coord_1, coord_2)
    
# get the coordination of the perpendicular point 
    if np.array_equal(coord_1, coord_2): 
        point_perpen = coord_1
    else: 
        point_perpen = coord_1 + d_13 * math.cos(angle_213) / d_12 * axis
## shifted coordination for the atom we want to rotate 
    v = coord_3 - point_perpen
# rotation_part
    rotated_v = rotate_vector(axis, theta, v)
    #rotated_v = Quaternion(axis=axis,angle=theta).rotate(v)  if pyquaternion is installed
# shift the rotated coordinate back by the vector from  perpendicular point 
    v_rotated = rotated_v + point_perpen
    return v_rotated

def get_atoms_pbc_rot(lines, infor):
    atom_list, coord_1, coord_2, axis, theta, B = infor[:]
    lines_pbc_rot = lines[:]
    for A in atom_list:
        coord_A   = get_intact_one_atom(lines, A, B)   # solve the PBC problem
        coord_rot = list(rotate_one_atom(infor, coord_A))  # rotate atom_A 
        ft        = lines[A + 8].rstrip().split()[3:]
        line_ele  = ' '.join([str(i) for i in coord_rot]) + '  ' + '  '.join(ft) + '\n'
        lines_pbc_rot[A + 8] = line_ele
    return lines_pbc_rot 

#def get_axis(point_1, point_2):
#    coor_1 = get_coordinate('POSCAR', point_1)                 
#    coor_2 = get_coordinate('POSCAR', point_2)                 
#                                                               
#    if point_1 == point_2:                                     
#        print('\n z axis croos the atom %s will be used as default.\n' %(point_1))  
#        axis = np.array([0, 0, 1])  
#    else:                                                      
#        axis  = coor_2  - coor_1                               
#    return coor_1, coor_2, axis 

def get_ele_name(lines, atom):
    dict_car2 = get_dicts(lines)[-1]
    atom = int(atom)
    for key, values in dict_car2.items():
        if atom in values:
            ele_name = key.split('-')[0]
    return ele_name 

def get_atom_list(lines, atom_s):
    '''
    Get the atom number list such as [10, 11, 12,13]
    atom_s: the arguments from command
    atom_selected: the list format of atom_s
    '''
    dict_car1, dict_car2  = get_dicts(lines)    
    total_atoms = sum([int(i) for i in str(lines[6]).strip().split()])
    def get_atom_selected(atom_selected):
        atom_list = []
        def atom_list_append(num):
            if num not in atom_list:
                atom_list.append(num)

        for i in atom_selected:
            if i.isdigit():
                atom_list_append(int(i))
            else:
                if '-' in i:
                    infor = i.strip().split('-')
                    n_start = int(infor[0])
                    if infor[1].isdigit():
                        n_end = int(i.split('-')[1])
                    else: # for the case: 18-, from 18 to the end.
                        n_end = total_atoms
                    for m in range(n_start, n_end + 1):
                        atom_list_append(m)
                else:
                    for k, v in dict_car2.items():
                        if i == k.split('-')[0]:
                            for ele in v:
                                atom_list_append(ele)
        #print('You select atoms:\t %s' %(atom_list))
        return atom_list

    atom_selected = []
    if isinstance(atom_s, str):
        atom_selected.append(atom_s)
    else:
        atom_selected.extend(atom_s)
    
    atom_list = get_atom_selected(atom_selected)
    return atom_list

def get_selected_lines(lines, atom_list):
    '''This function returns four parts:
        1) the ele_list: ['C', 'H', 'O'];
        2) the ene_name: [10,  12,   12];
        3) the coordinates: ['x1 y1 z1 T T T\n', 'x2 y2 z2 T T T\n', ...]
        4) the ele, cordinates: [ ['C', 'x1 y1 z1 T T T\n'], ['C', 'x2 y2 z2 T T T\n'], ... ]
    '''
    atom_ele_list = []            
    dict_car1, dict_car2 = get_dicts(lines) 
    for i in atom_list:  
        ele_name = get_ele_name(lines, i)
        atom_ele_list.append(ele_name)           
    
    coord_list = get_intact_molecule(lines, atom_list, atom_list[0])
    coord_dict = dict(zip(atom_list, coord_list))
    dict_num_ele = dict(zip(atom_list, atom_ele_list)) # {34:'H'}
    
    ele_num = []
    line_atoms = []
    coord_s = []
    
    ele_list = list(set(atom_ele_list))    # remove the duplicated elements
    for ele in ele_list:
        list_ele = []
        for k, v in dict_num_ele.items():
            if v == ele:
                list_ele.append(k)
                line = '%s T T T \n ' %(' '.join([str(m) for m in  coord_dict.get(k)]))
                line_atoms.append(line)
                coord_s.append([ele, line])
        ele_num.append(len(list_ele))
    return ele_list, ele_num, line_atoms, coord_s

def delete_one_atom(lines, atom_delete):
    '''Delte only one atom from the lines:
        1) lines are the data from CONTCAR or POSCAR; 
        2) atom_delete is the atom sequence in the original file.
        3) the output is the lines that can be used to write a now POSCAR.
        ''' 
    dict_car1, dict_car2 = get_dicts(lines)
    ele_name  = lines[5].strip().split()
    lines_deleted = lines[:]

    ## Update the element and number parts in POSCAR file
    deleted_ele = 'A-1000'
    num_deleted = 1000
    for k, v in dict_car2.items():
        if atom_delete in v:
            if dict_car1.get(k) == 1 :
                deleted_ele = k
                num_deleted = int(deleted_ele.split('-')[1])
                del dict_car1[k]
                del ele_name[num_deleted-1]   ## Do not use ele_name.remove(ele) method:  C:2, O:3, H:6, O:1, if we delete the last O and use list.remove method, the second one will be deleted.  
            else: 
                dict_car1[k] = dict_car1.get(k) - 1
              
    lines_deleted[5] = ' '.join(ele_name) + '\n' 
    
    line_num = ''
    for num, ele in enumerate(ele_name):
        if num >= num_deleted-1:   
            '''
            C  H  O N P 
            10 29 1 2 8 
            dict_car1 = {'C-1':10, 'H-2':29, 'O-3':1, 'N-4':2, 'P-5':8}
            if we delete the O atom, element list will be : C H N P
            the sequence number of N in the list will be 3. 
            So the key for N atom is 'N-3' now; but in dict_car1, the key for N is still 'N-4'.
            To avoid this case: we first get the seuqence number of deleted elment O, here is 3:
            In the new element list, if the element after O, the sequence number will + 1 to get the right key in dict_car1.    
            num + 2 = num + 1 + 1, the first 1 is to convert the python default start number 0 to 1, the sencond + 1 is described above. 
            '''   
            k = ele + '-' + str(num+2)
        else:
            k = ele + '-' + str(num+1)
        line_num = line_num + '  ' + str(dict_car1.get(k))
    line_num = line_num + '\n'
    lines_deleted[6] = line_num  

    # Delete the coordinate part.
    '''If you use this function to delete namy atoms such as : 1 2 3 4 5, so delete the 5 first, then, 4, 3, 2, 1 '''
    del lines_deleted[atom_delete + 8]
    
    return lines_deleted

def add_one_atom(lines, atom_add):
    '''lines  are lines from the file you want to copy to. 
    atom_add is a list: ['C', 'x y z T T T\n '] ''' 
    dict_car1, dict_car2 = get_dicts(lines)  
    ele_add, line_added = atom_add[:]  
    ele_name  = lines[5].strip().split() # ['C', 'H', ...]
    ele_num = [int(i) for i in lines[6].strip().split()] # [10, 12, ...]  
#    print(atom_add)
    is_in_or_not = False
    for k, v in dict_car1.items():
        if ele_add == k.split('-')[0]:
            is_in_or_not = True
        
    if  not is_in_or_not:
        lines[5] = '%s  %s\n ' %('  '.join(lines[5].strip().split()), ele_add)
        lines[6] = '%s  1\n ' %('  '.join(lines[6].strip().split()))
        lines.append(line_added)
        lines_add = lines
    else: 
        lines_add = []
        indice = ele_name.index(ele_add)  ## If there are two same elements in POSCAR like: Ti O C  H O, and we want to add another O to the POSCAR, the added atom will be after the first one.
        n_break   = sum(ele_num[0: indice+1]) + 8 + 1 
        ele_num[indice] = ele_num[indice] + 1
        lines[6] = '%s\n' %('  '.join([str(i) for i in ele_num]))
        lines_half_above  = lines[:n_break]
        lines_half_above.append(line_added)
        lines_half_below = lines[n_break:]
        lines_add.extend(lines_half_above)
        lines_add.extend(lines_half_below)
    return lines_add

def switch_atoms(lines, atom_1, atom_2):
    lines_s = lines[:]
    line_a1 = lines[atom_1+8]
    line_a2 = lines[atom_2+8]
    lines_s[atom_1+8] = line_a2
    lines_s[atom_2+8] = line_a1
    return lines_s

def get_vector_T(lines, A1, A2, A3, A4):
    coor_A1 = get_coordinate(lines, A1) 
    coor_A2 = get_coordinate(lines, A2) 
    coor_A3 = get_coordinate(lines, A3) 
    coor_A4 = get_coordinate(lines, A4) 
    
    vector_start = (coor_A1 +  coor_A2) / 2 
    vector_end    = (coor_A3  +  coor_A4)/2
    T = vector_end - vector_start
    return T

def move_one_atom(lines, atom, T):
    coord = get_coordinate(lines, atom)
    new_coord = coord + T
    return new_coord

def shift_atoms(lines, atom_list, T):
    for i in atom_list:
        print('Translate Atom %s-%s.' %(get_ele_name(lines, i), i))
        coord = list(str(k) for k in move_one_atom(lines, i, T))
        ft = lines[i+8].rstrip().split()[3:]
        line_ele = coord + ft  
        lines[i + 8] = ' '.join(line_ele) + '\n' 
    return lines

def rotate_vector(axis, theta, v):
    """ 
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d 
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d 
    eulerangles = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    new_v = np.dot(eulerangles, v)
    return new_v

def break_or_not(distance,threshold):
    distance = float(distance)
    threshold = float(threshold)
    status = None
    if distance >= threshold:
        status = 'Yes'
    else:
        status = 'No'
    return status


def get_pairs(a0, b):
    '''
a0 = [1,2,]
b = [1,2,3,4]
pairs = get_pairs(a,b)
print pairs

[[1, 1], [1, 2], [1, 3], [1, 4], [2, 2], [2, 3], [2, 4]]
    '''
    list_pair = []
    a = a0[:]
    count = len(a)
    while count > 0 :
        first = a[0]
        for j in b:
            pair = [first, j]
            pair.sort()
            if pair not in list_pair:
                list_pair.append(pair)
        del a[0]
        count = len(a)
    return  list_pair  

def bottom(file_in):
    '''This function is used to pull the cetered atoms (by using ASE) back to the bottom. '''
    lines = read_car(file_in)[0]
    coord = [float(line.rstrip().split()[2]) for line in lines[9:]]
    bottom = min(coord)
    #out_put = open(file_in + '_bottomed', 'w')
    out_put = open(file_in, 'w')
    out_put.writelines(i for i in lines[0:9])
    for line in lines[9:]:
        infor = line.rstrip().split()
        infor[2] = str(float(infor[2]) - bottom)
        out_put.write('   '.join(infor) + '\n')
    out_put.close()    


def bm_fitting(file_in):
    '''a, E: lattice parameters and corresponding energies in data file:
        2.8664, -201.834
        2.8554, -2.1.812
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
