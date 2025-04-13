#!/usr/bin/env python3 
import re
import numpy as np
from ase.thermochemistry import HarmonicThermo

def extract_epot_from_outcar(path="OUTCAR"):
    with open(path, 'r') as f:
        lines = f.readlines()
    for line in reversed(lines):
        if '  without' in line:
            try:
                return float(line.strip().split()[-1])
            except ValueError:
                continue
    raise ValueError("未找到能量行（包含 '  without'）")

def extract_vib_energies_from_outcar(path="OUTCAR"):
    vib_energies = []
    imag_freqs = []

    with open(path, 'r') as f:
        for line in f:
            if "f/i=" in line and "THz" in line:
                # 虚频行
                match = re.search(r'([-]?\d+\.\d+)\s+cm-1', line)
                if match:
                    imag_freqs.append(float(match.group(1)))
            elif "f  =" in line and "cm-1" in line :
                # 实频行，从最后一列提取 meV
                try:
                    energy_meV = float(line.strip().split()[-2])
                    energy_eV = energy_meV / 1000.0
                    vib_energies.append(energy_eV)
                except (IndexError, ValueError):
                    continue

    return vib_energies, imag_freqs

# === 主流程 ===
vib_energies, imag_freqs = extract_vib_energies_from_outcar("OUTCAR")
E_pot = extract_epot_from_outcar("OUTCAR")

thermo = HarmonicThermo(vib_energies=np.array(vib_energies),
                         potentialenergy=E_pot,
                         ignore_imag_modes=True)

T = 298.15
print(f"提取到实频数量: {len(vib_energies)}")
print(f"E_pot (from OUTCAR): {E_pot:.6f} eV")
#print(f"ZPE: {thermo.get_ZPE_correction():.6f} eV")
print(f"Gibbs correction @ {T} K: {thermo.get_helmholtz_energy(T):.6f} eV")

if imag_freqs:
    print(f"\n⚠️ 检测到 {len(imag_freqs)} 个虚频:")
    for i, f in enumerate(imag_freqs):
        print(f"  - Imaginary mode {i+1}: {f:.2f} cm⁻¹")

