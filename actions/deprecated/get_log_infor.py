#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 13:29:22 2019

@author: qli
"""
from gpm import *
from data import *

log_file = sys.argv[1]
 
dict_infor = get_dict_infor(log_file)
#try: 
dict_infor = add_thermo_to_dict_infor(dict_infor)
#except:
#    print('Check your log file %s' %(log_file))


spin = dict_infor.get('spin')


if dict_infor.get('method') == 'G4' :
    try:
        dict_hg_T = get_HG_dicts_G4(log_file, spin)
    except:
        dict_hg_T = {}
else:
    try: 
        dict_hg_T = get_HG_dicts(log_file, spin)
    except:
        dict_hg_T = {}
 
dict_infor.update(dict_hg_T)
#
name = dict_infor.get('name')
save_dict_txt(dict_infor, name+'_thermo')     

