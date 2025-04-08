import sys,os,math
import numpy as np
q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
if q_robot_path not in sys.path:
    sys.path.append(q_robot_path)
from cluster import *

# ### Get active sites

cluster_path = './slab'


l_rhombus, new_planar, new_non_planar = get_active_sites(cluster_path)

for site in l_rhombus:
    final_config = determine_full_configuration(cluster_path, '0.0', site)
    edft = compute_EDFT(final_config, gas_dict)
    print(site, edft)
    
