import sys,os,math
import numpy as np
q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
if q_robot_path not in sys.path:
    sys.path.append(q_robot_path)
from cluster import *

# test_path = '/mnt/c/Users/lqlhz/OneDrive - UMass Lowell/Projects/P1_cluster/Ru_clusters/run_3/Ru30/sum/NH2_ads'
# test_path =  convert_unix_to_windows_path(test_path)
# ### Get active sites

cluster_path = '/mnt/c/Users/lqlhz/OneDrive - UMass Lowell/Projects/P1_cluster/GA_training'
cluster_path = convert_unix_to_windows_path(cluster_path)


l_rhombus, coplanar_rhombus, edge_rhombus = get_active_sites(cluster_path)

for rhomnus_site  in l_rhombus[:1]:
    try: 
        print("Site:", rhomnus_site)
        E_top,  E_bri, E_hollow_NH, E_hollow_N, E_hollow_H = get_energy_one_site(rhomnus_site)
        # print(E_top,  E_bri, E_hollow_NH, E_hollow_N, E_hollow_H)
        E_R0 = E_top ## NH3 adsorption 
        E_R0 = E_R0 + 0.67260574 - -0.398484572
        Ea_R0 = 220.615 * 673 /96485  / 3 # 600K 
        Ear_R0 = Ea_R0 - E_R0
        

        
        E_R1 = E_bri  - E_top  ## NH3 -> NH2
        E_R1 = E_R1  + 0.505576774 + 0.189926157 - 0.67260574
        Ea_R1 = E_R1 * 0.42 + 1.36
        Ear_R1 = Ea_R1 - E_R1
        E_R1  = E_R1 + E_hollow_H

        
        E_R2 = E_hollow_NH  - E_bri ## NH2 -> NH
        E_R2 = E_R2 + 0.257163785 + 0.189926157 - 0.505576774
        Ea_R2 =  E_R2 * 0.88 + 0.82 
        Ear_R2 = Ea_R2 - E_R2
        E_R2  = E_R2 + E_hollow_H

        # 
        E_R3 = E_hollow_N  - E_hollow_NH  ## NH -> N
        E_R3 = E_R3 + -0.01813526  + 0.189926157 - 0.257163785
        Ea_R3 =  E_R2 * 0.47 + 0.81
        Ear_R3 = Ea_R3 - E_R3
        E_R3  = E_R3 + E_hollow_H

        E_R4 =  E_hollow_N * 2 ## N -> 0.5 N2  
        E_R4 = E_R4 + -1.150152398 - -0.01813526 * 2 
        Ea_R4 = - E_R4  * 0.77 + 2.33  # 2.83
        # print(Ea_R4)
        Ear_R4 = Ea_R4 - E_R4
        ### Reverse step does not occur.
        
        E_R5 =  - E_hollow_H * 2  ## H -> 0.5 H2 
        # E_R5 = E_R5 + -0.604615191 - 0.189926157 * 2 
        Ea_R5 = 151.077 * 673 / 96485  / 3  + E_R5 # 600 K 
        Ear_R5 = 151.077 * 673 /96485  / 3    
        print('E_R5:', E_R5, 'Ea_R5:', Ea_R5)
        ### Reverse step does not occur.
        
        Eaf = [Ea_R0, Ea_R1, Ea_R2, Ea_R3, Ea_R4, Ea_R5]
        Ear = [Ear_R0, Ear_R1, Ear_R2, Ear_R3, Ear_R4, Ear_R5]
        Ea = [round(i,2) for i in Eaf]
        print(Ea)
        
        kf0 = [get_k(i) for i in Eaf ]
        kr0 = [get_k(i) for i in Ear]
        
        rates =  get_rate(kf0, kr0)
        print(rates)
        # print(rhomnus_site,  round_scientific(min(rates), 6)) 
    
        # print('Rate_RDS: ', rates[:], '\n')                 
        # print('Input data:')
        # print('NH3 ads: E_R0, #Ea_R0, Ear_R0', [round(i, 2)  for i in [E_R0, Ea_R0, Ear_R0] ])
        # print('NH3 -> NH2: E_R1, Ea_R1, Ear_R1', [round(i, 2)  for i in [E_R1, Ea_R1, Ear_R1] ])
        # print('NH2 -> NH: E_R2, Ea_R2, Ear_R2',  [round(i, 2)  for i in [E_R2, Ea_R2, Ear_R2] ] )
        # print('NH -> N: E_R3, Ea_R3, Ear_R3',  [round(i, 2) for i in  [E_R3, Ea_R3, Ear_R3] ])
        # print('N -> N2: E_R4, Ea_R4, Ear_R4',   [round(i, 2) for i in [E_R4, Ea_R4, Ear_R4] ] )
        # print('H -> H2: E_R5, Ea_R5, Ear_R5',  [round(i, 2) for i in  [E_R5, Ea_R5, Ear_R5] ] )
        # print('Kf0:', kf0)
        # print('Kr0:', kr0)
        # print('MKM running')
    except: 
        pass 
    #     print('Failed for %s', rhomnus_site)
