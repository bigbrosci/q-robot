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
T = float(sys.argv[1])
#######################################################################
# Input Parameters
#######################################################################

tc = T  # Temperature in Celsius
length = 1 * cm  # Catalyst bed length
area = 1.0 * cm**2  # Catalyst bed area
cat_area_per_vol = 3000.0 / cm  # Catalyst particle surface area per unit volume
velocity = 3 * 60 * cm / minute  # gas velocity
porosity = 1  # Catalyst bed porosity
residence_time = length / velocity # unit: s
R = 8.314

# input file containing the surface reaction mechanism
yaml_file = 'thermo.yaml'

outfile_pfr_surf_cov = 'pfr_surface_coverage.csv'
outfile_pfr_gas_x = 'pfr_gas_x.csv'

# The PFR will be simulated by a chain of 'NReactors' stirred reactors.
NReactors = 201
dt = 1.0

#####################################################################

t = tc + 273.15  # convert to Kelvin

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

print('    distance       X_NH3        X_N2        X_H2')

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
        print('  {0:10f}  {1:10f}  {2:10f}  {3:10f}'.format(
            dist, *r.thermo['NH3', 'N2', 'H2'].X))

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


print("Results saved")

gas_reactant_id = 1 # The order in yaml file
# tof = (gas_concentrations_outlet[1] - gas_concentrations_outlet[0]) * 1e3 / residence_time / (1.3e-10 * 1e4 * cat_area_per_vol)
# assuming the first one is active site.
tof = (gas_concentrations_inlet[1] - gas_concentrations_outlet[1]) * 1e3 / residence_time / (surf_concentrations_inlet[0] * 1e3 * cat_area_per_vol)
conversion = (gas_concentrations_inlet[1] - gas_concentrations_outlet[1]) / gas_concentrations_inlet[1]

print(rsurf.kinetics.forward_rate_constants)
print(rsurf.kinetics.delta_enthalpy)

print("Turnover frequency: ", tof, "mol s-1 mol_cat-1")
print("Conversion ", conversion)

