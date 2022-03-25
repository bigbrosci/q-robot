#!/usr/bin/env python3
import re
from journal import * 
from data import *

'''
@article{li2018chirality,
  title={Chirality, Rigidity, and Conjugation: A First-Principles Study of the Key Molecular Aspects of Lignin Depolymerization on Ni-Based Catalysts},
  author={Li, Qiang and L{\'o}pez, N{\'u}ria},
  journal={ACS Catalysis},
  volume={8},
  number={5},
  pages={4230--4240},
  year={2018},
  publisher={ACS Publications},
  doi={10.1021/acscatal.8b00067}
}
'''
#line = 'title  =  {Chirality, Rigidity, and Conjugation: A First-Principles Study of the Key Molecular Aspects of Lignin Depolymerization on Ni-Based Catalysts},'
#title = re.search("title.*\}", bib).group(0)
#print(line.strip().split('=')[1])

###===================
bib_start, bib_end, bib_lines = locate_bib('hhdma.bib')

bib_out = open('hhdma_new.bib', 'w')
for num, start_num in enumerate(bib_start):
    end_num = bib_end[num]
    one_bib = bib_lines[start_num:end_num]
    one_bib_lines = format_one_bib(one_bib)
    bib_out.writelines(one_bib_lines)
#    bib_out.write('}\n')
bib_out.close()    
