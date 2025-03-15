#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 18:38:35 2022
1) https://janaf.nist.gov/tables/H-083.html #Database
2) https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html
3) For adsorption: J. Am. Chem. Soc. 2012, 134, 18109−18115
# Sads = 0.70 * S_gas - 3.3 R   
# slope = 0.7 
# intercept = con.gas_constant * 3.3
# TS_ads = TS_gas * slope + intercept


@author: qli
"""
import numpy as np
from ase.io import read, write
from ase.thermochemistry import IdealGasThermo
from scipy import constants as con


out = read('./OUTCAR',  format='vasp-out')
potentialenergy = out.get_potential_energy()


def get_gas_ads_TS():
    nsym = 3 # symmetry number of NH3
    spin = 0 # spin of NH3.  
    tem = 298.15
    pre = 101325
    model = read('./freq/POSCAR')

    vib_energies = []
    with open('./freq/OUTCAR') as f_in:
        lines = f_in.readlines()
        for num, line in enumerate(lines):
            if 'cm-1' in line:
                vib_e = float(line.rstrip().split()[-2])
                vib_energies.append(vib_e)
    
    vib_energies = np.array(vib_energies[:-6])/1000 # For Gas, the last six are translation and rotation
    
    thermo = IdealGasThermo(vib_energies=vib_energies,
                            potentialenergy=potentialenergy,
                            atoms=model,
                            geometry='nonlinear',
                            symmetrynumber=nsym, spin=spin)
    
    zpe = thermo.get_ZPE_correction()

    entropy = thermo.get_entropy(temperature=tem, pressure=pre,verbose=False) 
    TS = tem * entropy
    TS_ads = 0.70 * TS + con.gas_constant * 3.3
    return zpe, TS, TS_ads

zpe, TS, TS_ads = get_gas_ads_TS()
# G = potentialenergy + zpe - TS
# print(G)

# H = thermo.get_enthalpy(temperature=298.15, verbose=False)
# G = thermo.get_gibbs_energy(temperature=0, pressure=101325.,verbose=False)
# print(G)
# entropy = thermo.get_entropy(temperature=298.15, pressure=101325,verbose=False)
# print(entropy*98.485*1000)
# print(H, G, entropy)

# G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)


# for tem in [100, 200, 298.15, 300, 400, 500, 600, 700, 800]:
#     entropy_tem = thermo.get_entropy(temperature=tem, pressure=101325,verbose=False)
#     print(tem, con.Avogadro * con.electron_volt * entropy_tem)
