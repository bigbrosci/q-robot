import matplotlib.pyplot as plt
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.outputs import Locpot, Outcar
import numpy as np
# Plot settings
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'
def extract_fermi_energy(outcar_path):
    fermi_energy = None
    with open(outcar_path, 'r') as file:
        for line in file:
            if 'Fermi energy' in line:
                fermi_energy = float(line.split()[2])
    if fermi_energy is None:
        raise ValueError("Fermi energy not found in OUTCAR")
    return fermi_energy

locpot = Locpot.from_file('LOCPOT')

planar_average_potential = locpot.get_average_along_axis(2)

vacuum_level = np.max(planar_average_potential)

outcar_path = 'OUTCAR'
fermi_energy = extract_fermi_energy(outcar_path)

work_function = vacuum_level - fermi_energy

print("Planar Average Potential (sample):", planar_average_potential[:10], "...")
print("Vacuum Level:", vacuum_level, "eV")
print("Fermi Energy:", fermi_energy, "eV")
print("Work Function:", work_function, "eV")

z_length = locpot.structure.lattice.c
num_z_points = len(planar_average_potential)
z_positions = np.linspace(0, z_length, num_z_points)

plt.figure(figsize=(9, 7))
plt.plot(z_positions, planar_average_potential, linewidth=2, label="Planar Average Potential (Ru-Ba)")
plt.axhline(y=vacuum_level, color='#fdbd00', linestyle=':', linewidth=3, label="Vacuum Level")
plt.axhline(y=fermi_energy, color='#2da94f', linestyle=':', linewidth=3, label="Fermi Energy")
plt.axhline(y=work_function, color='#ea4335', linestyle=':', linewidth=3, label="Work Function")
plt.xlabel("z-direction (Å)", fontsize=22)
plt.ylabel("Potential (eV)", fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)
legend = plt.legend(fontsize=20)
legend.get_frame().set_alpha(0)
plt.tight_layout()

plt.savefig("Pot_vs_Z.png", dpi=1000)
plt.show()
