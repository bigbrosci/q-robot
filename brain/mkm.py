#!/usr/bin/env python
# coding: utf-8

# # OpenMKM Input and Output
# This notebook describes pmutt's functionality to write OpenMKM CTI and YAML files. We will use the NH3 formation mechanism as a case study.
# 
# ## Topics Covered
# - Read species *ab-initio* data, reactions, lateral interactions, phases, reactor operating conditions, and desired units from a spreadsheet
# - Write the CTI file that can be read by OpenMKM
# - Write a YAML file that can be read by OpenMKM

# ## Input Spreadsheet
# All the data will be imported from the [`./inputs/NH3_Input_data.xlsx`](https://github.com/VlachosGroup/pMuTT/blob/master/docs/source/examples_jupyter/omkm_io/inputs/NH3_Input_Data.xlsx) file. There are several sheets:
# 
# 1. `units` contains the units that types of quantities should be written
# 2. `refs` contains *ab-initio* and experimental data for a handful of gas species to calculate references (optional)
# 3. `species` contains *ab-initio* data for each specie
# 4. `beps` contains Bronsted-Evans-Polanyi relationships for reactions (optional)
# 5. `reactions` contains elementary steps 
# 6. `lateral_interactions` contains lateral interactions between species (optional)
# 7. `phases` contains phases for the species
# 8. `reactor` contains reactor operating conditions and solver tolerances
# 
# The ``refs``, ``beps`` and ``lateral_interactions`` sheets can be deleted and the code written below should still work.

# First, we change the working directory to the location of the Jupyter notebook.

import re
import os
from pathlib import Path

import numpy as np
import pandas as pd
from IPython.display import display

from pmutt import pmutt_list_to_dict
from pmutt.empirical.nasa import Nasa
from pmutt.empirical.references import Reference, References
from pmutt.empirical.shomate import Shomate
from pmutt.io.excel import read_excel
from pmutt.io.omkm import (organize_phases, write_cti, write_thermo_yaml,
                           write_yaml)
from pmutt.mixture.cov import PiecewiseCovEffect
from pmutt.omkm.reaction import BEP, SurfaceReaction
from pmutt.omkm.units import Units


def update_yaml(file_path):
    """
    Update the YAML file in place by performing the following modifications:
      1) Replace "thermo: surface-lateral-interaction" with "thermo: ideal-surface" globally.
      2) In the line that starts with two spaces and 'sticking-coefficient:',
         update the dictionary values so that:
            - 'b' is set to 0.
            - 'Ea' is set to "0.0 kcal/mol".

    Parameters:
      file_path (str): Path to the YAML file (e.g., "thermo.yaml").
    """
    # Read the entire file content.
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Replace all occurrences of the thermo keyword.
    content = content.replace("thermo: surface-lateral-interaction", "thermo: ideal-surface")

    # Replacement function for the sticking-coefficient line.
    def replace_sticking(match):
        inner_content = match.group(1)
        # Replace b: <any number> with b: 0
        inner_content = re.sub(r'b:\s*\d+', "b: 0", inner_content)
        # Replace Ea: "<any value> kcal/mol" with Ea: "0.0 kcal/mol"
        inner_content = re.sub(r'Ea:\s*".*?\s*kcal/mol"', 'Ea: "0.0 kcal/mol"', inner_content)
        return "  sticking-coefficient: {" + inner_content + "}"

    # Use regex to find the line starting with two spaces and 'sticking-coefficient:'
    content = re.sub(r'^  sticking-coefficient:\s*\{(.*?)\}', replace_sticking, content, flags=re.M)

    # Write the modified content back to the file.
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)

    print(f"✅ Updated {file_path}")


vib_dict = {'NH3(T)': [3493.828275, 3483.710224, 3380.976917, 1601.337218, 1593.083158, 1152.474879, 565.120729, 545.201983, 353.276513, 115.51182, 100.062955, 60.913836],
'TS1_NH3(T)': [3536.75787, 3441.705733, 1599.158985, 1484.669619, 786.688583, 636.527507, 448.855255, 410.494574, 138.396713, 96.519331, 50.0],
'NH2(T)': [3505.326237, 3416.622069, 1499.074193, 693.441876, 675.162244, 635.295685, 467.523262, 279.469855, 186.051905],
'Hv1(T)': [1309.815216, 722.994226, 645.912728],
'TS2_NH2(T)': [3367.301946, 1680.23037, 934.078205, 663.303152, 539.449475, 388.599799, 312.382306, 211.618148],
'NH(T)': [3298.323841, 780.751038, 594.237204, 510.881815, 425.919495, 240.365103],
'Hv2(T)': [1309.815216, 722.994226, 645.912728],
'TS3_NH(T)': [1541.246549, 610.817237, 442.82166, 327.065925, 100.302078],
'N(T)': [649.119623, 461.607208, 268.406425],
'Hv3(T)': [1448.228446, 955.85136, 338.524455],
'N_N(T)': [636.454342, 600.178663, 444.676465, 408.028287, 330.614256, 165.870618],
'TS4_N2(T)': [824.001797, 770.689908, 258.340484, 190.421562, 144.687786],
'N2(T)': [1875.97933, 432.118227, 304.804693, 158.135021, 59.049389],
'H(T)': [2414.172895, 1332.346614, 454.863852]}


def prepare_yaml(main_folder):    
    # Find the location of Jupyter notebook
    # Note that normally Python scripts have a __file__ variable but Jupyter notebook doesn't.
    # Using pathlib can overcome this limiation
    # try:
    #     notebook_path = os.path.dirname(__file__)
    # except NameError:
    #     notebook_path = Path().resolve()
        
    # os.chdir(notebook_path)
    
    input_path = os.path.join(main_folder, "inputs", "NH3_Input_Data.xlsx")
    
    
    units_data = read_excel(io=input_path, sheet_name='units')[0]
    units = Units(**units_data)
    
    
    # ### Reading References (optional)
    # Second, we will open the input spreadsheet and read the `refs` sheet.
    try:
        refs_data = read_excel(io=input_path, sheet_name='refs')
    except:
        # If references are not used, skip this section
        print('The "refs" sheet could not be found in {}. Skiping references'.format(input_path))
        refs = None
    else:
        refs = [Reference(**ref_data) for ref_data in refs_data]
        refs = References(references=refs)

    
    # ### Reading Species
    # 
    # Third, we will use the ``refs`` defined before and the ``species`` sheet to convert statistical mechanical data to [``NASA``](https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa.html#pmutt.empirical.nasa.Nasa) objects.
    # Read the species' data
    species_data = read_excel(io=input_path, sheet_name='species')
    # print(species_data)
    for species in species_data:
        if species['name'] in vib_dict.keys():
            species['vib_wavenumbers']  =  vib_dict[species['name']]
    
    # print(species)
    
    # Create NASA polynomials from the species
    species = [Nasa.from_model(references=refs, **ind_species_data)            for ind_species_data in species_data]
    
    
    # ### Adding species from other empirical sources (optional)
    # 
    # Note that OpenMKM also supports [``Shomate``](https://vlachosgroup.github.io/pMuTT/api/empirical/shomate/pmutt.empirical.shomate.Shomate.html#pmutt.empirical.shomate.Shomate) and [``NASA9``](https://vlachosgroup.github.io/pMuTT/api/empirical/nasa/pmutt.empirical.nasa.Nasa9.html) objects. Below, we define a single ``Shomate`` species.
    Ar = Shomate(name='Ar', elements={'Ar': 1}, phase='gas', T_low=298., T_high=1000.,
                 a=np.array([20.78600, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8, -6.19735, 179.999, 0.]))
    
    species.append(Ar)
    
    
    # ### Read reactions
    # 
    # Then, we read the reactions to include.
    # Convert species to dictionary for easier reaction assignment
    species_with_beps = species.copy()
    species_with_beps_dict = pmutt_list_to_dict(species_with_beps)
    
    reactions_data = read_excel(io=input_path, sheet_name='reactions')
    reactions = [SurfaceReaction.from_string(species=species_with_beps_dict, **reaction_data)              for reaction_data in reactions_data]
    
    interactions = None
    # ### Reading Phases
    # 
    # Finally, we read the phases data from Excel and organize it for use in OpenMKM.
    # Read data from Excel sheet about phases
    phases_data = read_excel(io=input_path, sheet_name='phases')
    #print(phases_data)
    phases = organize_phases(phases_data, species=species, reactions=reactions, interactions=interactions)
    
    
    # ## Write Reactor YAML File
    # 
    # The YAML file specifying the reactor configuration can be written using the [``write_yaml``](https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_yaml.html) function. Note that if:
    # - ``units`` is not specified, float values are assumed to be in SI units
    # - ``units`` is specified, float values are consistent with ``unit``'s attributes
    # - you would like a quantity to have particular units, pass the value as a string with the units  (e.g. "10 cm3/s").
    # Path('./outputs').mkdir(exist_ok=True)
    
    outputs_folder = os.path.join(main_folder, "outputs")
    os.makedirs(outputs_folder, exist_ok=True)

    # yaml_path = './outputs/reactor.yaml'
    reactor_data = read_excel(io=input_path, sheet_name='reactor')[0]
    # write_yaml(filename=yaml_path, phases=phases, units=units, **reactor_data)
    
    
    # If you would prefer to return the file as a string instead of writing it, omit the ``filename``.
    # print(write_yaml(phases=phases, units=units, **reactor_data))
    
    # ## Write Thermo/Kinetic YAML File
    # 
    # As of OpenMKM version 0.6.0 onwards, the thermodynamic and kinetic parameters can be written as a YAML file. We recommend using this format over the older CTI format. To generate the Thermo/Kinetic YAML file using pMuTT, use the [``write_thermo_yaml``](https://vlachosgroup.github.io/pMuTT/api/kinetic_models/omkm/pmutt.io.omkm.write_thermo_yaml.html) function
    T = reactor_data['T']
    out_name = os.path.join(main_folder, "outputs", "thermo.yaml")
    write_thermo_yaml(T=T,
                      phases=phases,
                      species=species,
                      reactions=reactions,
                      lateral_interactions=interactions,
                      units=units,
                      filename=out_name)  
    
    update_yaml(out_name)
    
    
def run_mkm(main_folder, T):
        
    """
    This example solves a plug flow reactor problem, where the chemistry is
    surface chemistry. The specific problem simulated is the partial oxidation of
    methane over a platinum catalyst in a packed bed reactor. To avoid needing to solve a
    DAE system, the PFR is approximated as a chain of successive WSRs. See surf_pfr.py
    for a more advanced implementation that solves the DAE system directly.
    
    Requires: cantera >= 3.0
    Keywords: catalysis, reactor network, surface chemistry, plug flow reactor,
              packed bed reactor
    """
    
    import csv
    import cantera as ct
    import sys 
    
    # unit conversion factors to SI
    cm = 0.01
    minute = 60.0
    #######################################################################
    # Input Parameters
    #######################################################################
    
    tc = float(T)  # Temperature in Celsius
    t = tc + 273.15  # convert to Kelvin
    length = 1 * cm  # Catalyst bed length
    area = 1.0 * cm**2  # Catalyst bed area
    cat_area_per_vol = 3000.0 / cm  # Catalyst particle surface area per unit volume
    velocity = 3 * 60 * cm / minute  # gas velocity
    porosity = 1  # Catalyst bed porosity
    residence_time = length / velocity # unit: s
    R = 8.314
    
    # input file containing the surface reaction mechanism
    yaml_file = os.path.join(main_folder, "outputs", 'thermo.yaml')
    
    outfile_pfr_surf_cov = os.path.join(main_folder, "outputs", "pfr_surface_coverage.csv")
    outfile_pfr_gas_x    = os.path.join(main_folder, "outputs", "pfr_gas_x.csv")
    mkm_out = os.path.join(main_folder, "outputs", 'mkm.out')
    
    # The PFR will be simulated by a chain of 'NReactors' stirred reactors.
    NReactors = 201
    dt = 1.0
    
    #####################################################################
    
    # import the gas model and set the initial conditions
    gas = ct.Solution(yaml_file, 'gas')
    gas.TPX = t, ct.one_atm, 'NH3:1, N2:0, H2:0'
    gas_concentrations_inlet = gas.concentrations
    
    # import the surface model
    surf = ct.Interface(yaml_file, 'terrace', [gas])
    surf.TP = t, ct.one_atm
    surf_concentrations_inlet = surf.concentrations
    
    rlen = length/(NReactors-1)
    rvol = area * rlen * porosity
    
    # catalyst area in one reactor
    cat_area = cat_area_per_vol * rvol
    
    mass_flow_rate = velocity * gas.density * area * porosity
    
    # The plug flow reactor is represented by a linear chain of zero-dimensional
    # reactors. The gas at the inlet to the first one has the specified inlet
    # composition, and for all others the inlet composition is fixed at the
    # composition of the reactor immediately upstream. Since in a PFR model there
    # is no diffusion, the upstream reactors are not affected by any downstream
    # reactors, and therefore the problem may be solved by simply marching from
    # the first to last reactor, integrating each one to steady state.
    
    TDY = gas.TDY
    cov = surf.coverages
    
    with open(mkm_out, "w") as fout:
            
        fout.write('    distance       X_NH3        X_N2        X_H2 \n')
        
        # create a new reactor
        gas.TDY = TDY
        r = ct.IdealGasReactor(gas, energy='off')
        r.volume = rvol
        
        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(gas, name='upstream')
        
        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(gas, name='downstream')
        
        # Add the reacting surface to the reactor. The area is set to the desired
        # catalyst area in the reactor.
        rsurf = ct.ReactorSurface(surf, r, A=cat_area)
        
        # The mass flow rate into the reactor will be fixed by using a
        # MassFlowController object.
        m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)
        
        # We need an outlet to the downstream reservoir. This will determine the
        # pressure in the reactor. The value of K will only affect the transient
        # pressure difference.
        v = ct.PressureController(r, downstream, primary=m, K=1e-5)
        
        sim = ct.ReactorNet([r])
        
        outlist_pfr_surf_cov = []
        outlist_pfr_gas_x = []
        
        
        for n in range(NReactors):
            # Set the state of the reservoir to match that of the previous reactor
            gas.TDY = r.thermo.TDY
            upstream.syncState()
            sim.reinitialize()
            sim.advance_to_steady_state(max_steps=10000, residual_threshold=0.00001, atol=0.00001, return_residuals=True)
            dist = n * rlen * 1.0e3  # distance in mm
        
            if n % 10 == 0:
                fout.write('  {0:10f}  {1:10f}  {2:10f}  {3:10f}\n'.format(dist, *r.thermo['NH3', 'N2', 'H2'].X))
        
            # print(rsurf.kinetics.coverages)
            outlist_pfr_surf_cov.append([dist] + list(rsurf.kinetics.coverages))
            outlist_pfr_gas_x.append([dist, r.thermo.P / ct.one_atm] + list(r.thermo.X))    
        
        gas_concentrations_outlet = gas.concentrations
        
        with open(outfile_pfr_surf_cov, 'w', newline="") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Distance (mm)'] + surf.species_names)
            writer.writerows(outlist_pfr_surf_cov)
        
        with open(outfile_pfr_gas_x, 'w', newline="") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Distance (mm)', "Pressure (atm)"] + gas.species_names)
            writer.writerows(outlist_pfr_gas_x)
        
        
        gas_reactant_id = 1 # The order in yaml file
        # tof = (gas_concentrations_outlet[1] - gas_concentrations_outlet[0]) * 1e3 / residence_time / (1.3e-10 * 1e4 * cat_area_per_vol)
        # assuming the first one is active site.
        tof = (gas_concentrations_inlet[1] - gas_concentrations_outlet[1]) * 1e3 / residence_time / (surf_concentrations_inlet[0] * 1e3 * cat_area_per_vol)
        conversion = (gas_concentrations_inlet[1] - gas_concentrations_outlet[1]) / gas_concentrations_inlet[1]
        
        
        
        fout.write(str(rsurf.kinetics.forward_rate_constants) + "\n")
        fout.write(str(rsurf.kinetics.delta_enthalpy) + "\n")
        fout.write("Turnover frequency: " + str(tof) + " mol s-1 mol_cat-1\n")
        fout.write("Conversion " + str(conversion) + "\n")
