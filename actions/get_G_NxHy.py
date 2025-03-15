#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 08:14:11 2024
@author: lqlhz
"""

from ase.units import _hplanck, _k, _e, _c, _amu
from ase.io import read
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
import numpy as np
from scipy.constants import h, c, e
import sys 

def wavenumber_to_ev(wavenumber_cm1):
    """ convert wavenumber from cm-1 to eV """
    conversion_factor = (h * c * 100) / e
    energy_ev = wavenumber_cm1 * conversion_factor
    return energy_ev

def get_energy_from_outcar(outcar_path):
    energy = None
    with open(outcar_path, 'r') as file:
        for line in file:
            if '  without' in line:
                energy = line.split()[-1]
    
    if energy:
        return float(energy)
    else:
        return None
    
def extract_wavenumbers_from_outcar(outcar_path):
    wavenumbers = []
    with open(outcar_path, 'r') as file:
        for line in file:
            if ' f  =' in line:
                # The frequency value is usually the 9th element when the line is split by spaces
                wavenumber = float(line.split()[7])
                wavenumbers.append(float(wavenumber))
    # wavenumbers = list(set(wavenumbers))
    wavenumbers = list(dict.fromkeys(wavenumbers))

    # Simplify the wavenumbers directly within this function
    for i, nu in enumerate(wavenumbers):
        if nu <= 100:
            wavenumbers[i] = 100
    return wavenumbers


### Calculate free energies of gas species"
def get_free_energies_species(species, T):  
    outcar = species+'/OUTCAR'
    outcar_freq = species+'/freq/OUTCAR'
    geom = read(species+'/CONTCAR')
    e = get_energy_from_outcar(outcar)
    wavenumbers = extract_wavenumbers_from_outcar(outcar_freq)
    modes = [wavenumber_to_ev(i) for i in wavenumbers]
    
    if species in ['NH_gas', 'CO_gas', 'H2_gas', 'N2_gas', 'O2_gas']:
        geometry = 'linear'
    elif species in ['N_gas', 'H_gas', 'O_gas']:
        geometry = 'monatomic'
    else:
        geometry = 'nonlinear'
        
    
    if species in ['H_gas', 'N_gas', 'NH_gas', 'NH-N_gas', 'NH2-NH_gas']:
        symmetrynumber = 1 
    elif species in ['H2_gas', 'N2_gas', 'NH2_gas', 'NH-NH_gas', 'NH2-N_gas', 'NH2-NH2_gas']:
        symmetrynumber = 2 
    elif species in ['NH3_gas']:        
        symmetrynumber = 3
    else: 
        symmetrynumber = 1 
        
    if species in ['N_gas']:
        spin = 1.5   # 3 unpaired electrons
    elif species in ['NH_gas']:
        spin = 1.0   
    elif species in ['NH2_gas', 'NH-N_gas', 'NH2-NH_gas']:
        spin = 0.5 
    else:
        spin = 0
        
    if 'gas' in species: ### Calculate free energies of gas species"
        thermo = IdealGasThermo(vib_energies = modes, potentialenergy = e, atoms = geom,
                                geometry = geometry, symmetrynumber = symmetrynumber, spin = spin)
        P0 = 1e5 # [Pa]
        g0 = thermo.get_gibbs_energy(temperature = T, pressure = P0, verbose = False) # [eV]
       # S = thermo.get_entropy(temperature = T, pressure = P0, verbose = False)  # eV/K
    
    elif 'surf' in species: ### Calculate free energies of surface species"

        thermo = HarmonicThermo(vib_energies = modes, potentialenergy = e)
        g0 = thermo.get_helmholtz_energy(temperature = T, verbose = False) # [eV}
        #S = thermo.get_entropy(temperature = T, verbose = False)  

    return g0

species = sys.argv[1]

g0 = get_free_energies_species(species, 400)

print(g0)

