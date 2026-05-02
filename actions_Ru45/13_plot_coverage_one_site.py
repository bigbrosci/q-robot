#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grouped bar chart of surface coverages for three EF cases.
- Colors follow the agreed mapping by REAL EF:
    data -0.6 -> real +0.6 -> red  (#d62728)
    data  0.0 -> real  0.0 -> orange (#ff7f0e)
    data +0.6 -> real -0.6 -> blue (#1f77b4)
- Legend label shows REAL EF (sign-inverted from data tag).
Usage:
    python3 13_plot_coverage_one_site.py [<site_dir>]
    If no site_dir provided, processes all sites under mkm_inputs/
"""
import os, re, sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LinearLocator, FormatStrFormatter

plt.rcParams.update({
    "font.family": "Times New Roman",
    "pdf.fonttype": 42,   # editable text in Illustrator
    "ps.fonttype": 42,
    "legend.fontsize": 16,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
})

# --- Config ---
# Folder EF tags (data tags)
ef_list = ["0.6", "0.0", "-0.6"]

# --- Standard EF color codes by REAL EF (mapped from data tags) ---
# data -0.6 -> real +0.6 (red); data 0.0 -> real 0.0 (orange); data +0.6 -> real -0.6 (blue)
color_map = {
    "-0.6": "#d62728",  # data -0.6 -> real +0.6
    "0.0":  "#ff7f0e",  # real 0.0
    "0.6":  "#1f77b4",  # data +0.6 -> real -0.6
}

# Exclude columns with these substrings
exclude_keywords = ["RU", "N_N", "N2", "TS", "Hv1", "Hv2", "Hv3", "Distance (mm)"]

def plot_coverage(base_dir):
    """Generate coverage plot for a single site directory."""
    base_path = Path(base_dir)
    
    # --- Load surface coverage (last row, exclude unwanted species) ---
    coverage_data = {}
    species = []

    for ef in ef_list:
        file_path = base_path / f"EF_{ef}" / "outputs" / "pfr_surface_coverage.csv"
        if not file_path.is_file():
            print(f"  [WARN] Missing: {file_path}")
            continue

        df = pd.read_csv(file_path)
        last_row = df.iloc[-1]

        # Keep only species without excluded keywords
        filtered_species = [col for col in df.columns if all(key not in col for key in exclude_keywords)]
        species = filtered_species  # assume same across EFs
        coverage_data[ef] = [last_row[s] for s in species]

    if not coverage_data:
        print(f"  No coverage data found in {base_dir}")
        return

    # --- Plot grouped bar chart ---
    x = np.arange(len(species))
    bar_width = 0.25

    fig, ax = plt.subplots(figsize=(5, 5))

    for i, ef in enumerate(ef_list):
        if ef in coverage_data:
            real_ef = -float(ef)
            label = f"EF = 0.0 V/Angstrom" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Angstrom"
            offset = (i - (len(ef_list)-1)/2) * bar_width
            ax.bar(x + offset, coverage_data[ef], width=bar_width,
                   label=label, color=color_map.get(ef, "#000000"))

    # --- Formatting ---
    ax.set_xticks(x)

    # Clean labels: (T)->*, drop 'RU', add subscripts like NH3 -> NH$_3$
    species_clean = [s.replace("(T)", "*").replace("RU", "") for s in species]

    def fmt_subscripts(labels):
        return [re.sub(r'([A-Za-z\)\)])(\d+)', r'\1$_{\2}$', s) for s in labels]

    species_fmt = fmt_subscripts(species_clean)
    ax.set_xticklabels(species_fmt)
    ax.set_ylabel("Surface Coverage (theta)", fontsize=20)
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylim(0, 1)
    ax.margins(y=0)
    ax.yaxis.set_major_locator(LinearLocator(5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.axhline(0, color='gray', linestyle='--', lw=0.8)
    ax.legend(frameon=False, fontsize=15)
    plt.tight_layout()

    # --- Save ---
    output_path = base_path / "surface_coverage_grouped_bar.png"
    plt.savefig(str(output_path), dpi=300)
    plt.close()
    print(f"  Saved plot to: {output_path}")


def main():
    # If a specific site is provided as argument, process only that one
    if len(sys.argv) > 1:
        site_dir = sys.argv[1]
        print(f"Processing: {site_dir}")
        plot_coverage(site_dir)
        return
    
    # Otherwise, search for all site directories under mkm_inputs/
    mkm_base = Path("mkm_inputs")
    if not mkm_base.is_dir():
        print(f"[ERROR] '{mkm_base}' directory not found.")
        sys.exit(1)
    
    # Find all site directories (folders that contain EF_* subfolders)
    site_dirs = []
    for item in sorted(mkm_base.iterdir()):
        if item.is_dir() and not item.name.startswith("EF_"):
            # Check if this directory contains EF_* folders
            has_ef = any(sub.name.startswith("EF_") for sub in item.iterdir() if sub.is_dir())
            if has_ef:
                site_dirs.append(item)
    
    if not site_dirs:
        print(f"[ERROR] No site directories found under {mkm_base}/")
        sys.exit(1)
    
    print(f"Found {len(site_dirs)} site(s) to process:\n")
    for site_path in site_dirs:
        print(f"Processing: {site_path}")
        plot_coverage(str(site_path))
        print()
    
    print("All sites processed.")


if __name__ == "__main__":
    main()
