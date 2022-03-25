#!/usr/bin/env python3
# Read POSCAR and Add the C6 and R0 parameters to INCAR file
# Written By Qiang 
import subprocess  

## Clean the vdW parameters in INCAR file
subprocess.call('sed -i "/VDW/d" INCAR', shell=True)
 
c6dict = {
'Au': '7.308',
'Ni': '2.6263',
'Ru': '4.1678',
'Cu': '2.740',
'O': '0.70',
'Si': '9.23',
'H': '0.14',
'C': '1.75',
'B': '3.13',
'N': '1.23',
'Br': '12.47',
'S': '5.57',
'Pd': '5.51',
'Pt': '7.00',
'Cl': '5.07',
'Ag': '5.481',
'Co': '2.565',
'Rh': '4.364',
'Ir': '6.163',
'Os': '5.878',
'P': '7.84',
'Fe': '2.60',
    } 
     
r0dict = {
'Au': '1.823',
'Ni': '1.562',
'Ru': '1.639',
'Cu': '1.562',
'O' : '1.342',
'Si': '1.716',
'H': '1.001',
'C': '1.452',    
'B': '1.485',
'N': '1.397',
'Br': '1.749',
'S': '1.683',
'Pd': '1.69',
'Pt': '1.75',
'Cl': '1.639',
'Ag': '1.819',
'Co': '1.349',
'Rh': '1.677',
'Ir': '1.698',
'Os': '1.504',
'P': '1.705',
'Fe':'1.40',
}
    

poscar =  open('POSCAR', 'r')
ele_list =  poscar.readlines()[5].rstrip().split()

incar = open ('INCAR', 'a')
incar.write('\n#VDW parameters:\n')
incar.write('LVDW   =   T\n')       
incar.write('VDW_VERSION= 2\n')
incar.write('VDW_RADIUS = 40\n')
incar.write('VDW_SCALING = 0.75\n')        
incar.write('VDW_C6 = ')
for i in ele_list:
    incar.write('%s ' %(c6dict.get(i)))
incar.write('\n')
incar.write('VDW_R0 = ')
for i in ele_list:
    incar.write('%s ' %(r0dict.get(i)))
incar.write('\n')
    
poscar.close()
incar.close()    
