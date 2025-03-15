#!/usr/bin/eng python3 

file_in = open('CoPc.xyz')
lines = file_in.readlines()
file_in.close()

cell = lines[1].split('"')[1].split()
print(cell)
cell_a = '  '.join(cell[:3]) 
cell_b = '  '.join(cell[3:6]) 
cell_c = '  '.join(cell[6:]) 
print(cell_a, cell_b, cell_c)
