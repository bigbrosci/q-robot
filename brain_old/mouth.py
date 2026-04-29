#!/usr/bin/env python3
from lattice import read_car 


#bond_type = ['C-H', 'C-O', 'C-C', 'O-O', 'O-H', 'H-H']
#bond_dis_threshold  = [1.2, 1.5, 1.65, 1.55, 1.05, 0.8]
#bond_threshold = dict(zip(bond_type, bond_dis_threshold))
#print(bond_threshold)           

def print_poscar(file_read, ele_name):
    lines, ele_name, ele_num, dict_car1, dict_car2, dict_car3 = read_car(file_read)
    for i in ele_name:                                                 
        n_start = int(dict_car2.get(i).split('---')[0])                
        n_end   = int(dict_car2.get(i).split('---')[1])                
        print('\nAtoms: %3.0f --- %3.0f are\t %s' %(n_start, n_end, i) )   

## Not used ###
#def print_poscar_total(file_read):
#    print '\nThe elements and their sequential numbers  in POSCAR  are: \n'
#    dict_car2 = read_car('POSCAR')[-1]
#    for i in dict_car2.keys():
#        n_start = int(dict_car2.get(i).split('---')[0]) 
#        n_end   = int(dict_car2.get(i).split('---')[1]) 
#        list_for_usr = [] 
#        for j in range(n_start, n_end):
#            list_for_usr.append(int(j))  
#        print i, str(list_for_usr)  
#    print '---' * 10 + '\n'
        
def print_dire_to_cart():
    print( 'Please read the head part of this script and get more information!')
    print( """
    ###################################
    #                                 #
    #for VASP 5.2 or higher versions  #
    #                                 #
    ###################################
    """ )       

def print_output():
    print('''\nThe out_put file is named as POSCAR_out.\n \n Please check the structures before running the calculations.\n''') 

def print_rotate():
    print(
            '''
            Note-1:
            Along the axis(atom_A, atom_B), rotate atoms by theta angle (in degree).
            Command: rotate.py atom_A atom_B atom_1 atom_2 ... Angle    
            For example: to rotate the C H O atoms 15 degree along the axis formed by No.2 and No.10 atoms.
            rotate.py 2 10 C H O 15
            
            Note-2: 
            This script only works on POSCAR file with Cartesian Coordinates.          
            
            ''')
def print_xps():
    print('''
    \nThe POSCAR for XPS calculation is named as POSCAR_xps\n\nWarning! \nWarning! \nWarning! \n 
    Please check the following issues before you submit jobs 
    
    1) POTCAR 
    2) INCAR:
        i)  DFT+D2 
        ii) DFT+U 
        iii) MAGMOM
    All these items above should be consistent with the NEW POSCAR. 
    
   3)Add the following to your INCAR file
    
    ###########################################################
    ICORELEVEL = 2
    CLNT = 1
    CLN = XXX   # main quantum number of excited core electron
    CLL = XXX   # l quantum number of excited core electron
    CLZ = 1 
    ##########################################################
    
    Update the main and l quantum number parameters.
    
    ''')


def print_dos_extract():
    print('''
This script is used to get the specific DOS information.
For example: the 5th atoms' px orbital 
Part of the codes are copied from the split_dos.py script of VTST.
A lots of ideas were also inspired by this script.  
Author: Qiang Li
Email: lqcata@gmail.com or qli@bigbrosci.com
'''

### Usage Commands ###

'''
In terminal, type: 

python dos_extract_v5.py  atoms orbitals  out_put_file_name  

'''


### Input Instruction for atoms / elements:   ###
'''
Please Type the Atom(s) names or their sequential Numbers  you are intersted in : 

Choose 1) or 2) input style only:

1) Use the element symbols: 
Example A): Ru        # For all atoms of sinlge element
Example B): Ru C H O  # For all atoms of many elements 

2) Use the sequential numbers in the POSCAR or CONTCAR:
Example A): 10           # for single atom: here is the 10th atom 
Example B): 1 2 3 4 5 10 # for many atoms: here are the 1st, 2nd, 3rd, 4th, 5th and 10th  atoms 
Example C): 1-5 10       # the same results as B)

3) Note 1): 
If the elements in your POSCAR or CONTCAR are like this:  Ru C H O C  (There are two C elements)
Please use the 2) input method. Otherwise only the first C atoms will be considered. 

'''
### Input Instruction for orbitals: ###
'''
Please Type the Orbitals you are intersted in:

Choose 1) or 2) input style:

1) Use the orbital symbols:
Example A):  p   # for px, py and pz orbitals 
Example B):  d   # for dxy, dyz, dxz, dz2 and dx2 orbitals
Example C): p d  # sum of  A) and B) 
Example D): t2g  # for dxy, dyz, dxz orbitals
Example E): eg   # for dx2, dz2  orbitals 

2) Use the specific orbital names: 
Example A) px        # to get the px orbital dos only 
Example B) f1         # to get f1 orbital dos only 
Example C) dxy dz2    # to get the sum dos of dxy and dz2 orbitals 
Example D) px dxy f1  # to get the sum dos of px, dxy and f1 orbitals 

3) Note: The d(x**2 - y**2) orbital is named dx2 here. 

'''

### Most complicated Examples ###

'''
1) To get the s, total p, total d,  f1, and f2 orbitals for all Ce H O atoms:  
python dos_extract_v3.py Ce H O  s p d f1 f2  dos.dat 

2) To get the px, total d, f3 and f5 orbitals for Ce and O atoms: 
python dos_extract_v3.py Ce O px d f3 f5 dos.dat 

3) To get the t2g orbitals of Co atoms:

python dos_extract_v3.py Co t2g dos.dat 

4) To get the total d orbitals of 1 2 3 4 10 12 atoms:

python dos_extract_v3.py 1-4 10 12 d dos.dat

5) DOS0 file contains total dos for all orbitals of all atoms 

In above two examples, the output file is named dos.dat, of course, you can use other names you like.
'''
)
        
