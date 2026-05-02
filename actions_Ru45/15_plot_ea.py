#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15_plot_ea.py

Grouped barchart and trend plots of Ea1..Ea4 over EF = -0.6, 0.0, +0.6.
Reads from mkm_out.csv and processes all sites.

Usage:
  python3 15_plot_ea.py                              # processes all sites from mkm_out.csv
  python3 15_plot_ea.py -s 8_11_17                   # process single site
  python3 15_plot_ea.py -i mkm_outputs/mkm_out.csv

Output:
  - Ea_grouped_<site>.png: grouped barchart for each site
  - Ea1_<site>.png, Ea2_<site>.png, Ea3_<site>.png, Ea4_<site>.png: trend plots
  All saved to mkm_inputs/<site>/ directory
"""

import argparse, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# ---------- style (match your existing script) ----------
def configure_style():
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "legend.fontsize": 16,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
        "axes.labelsize": 20,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
    })

# ---------- utilities ----------
def norm_ef(ef: str) -> str:
    """'-0.6'->'-0.6', '0'->'0.0', '0.0'->'0.0' """
    return f"{float(ef):.1f}"

# ---------- color & label helpers (by REAL EF) ----------
COL_NEG  = "#1f77b4"  # real -0.6
COL_ZERO = "#ff7f0e"  #  0.0
COL_POS  = "#d62728"  # real +0.6

def color_for_data_ef(data_ef_str: str) -> str:
    v = float(data_ef_str)
    real = -v
    if abs(real) < 5e-4:
        return COL_ZERO
    return COL_POS if real > 0 else COL_NEG

def label_for_data_ef(data_ef_str: str) -> str:
    v = float(data_ef_str)
    real = -v
    return "EF = 0.0 V/Å" if abs(real) < 5e-4 else f"EF = {real:+.1f} V/Å"

# ---------- plotting ----------
def plot_grouped_ea(ea_by_ef: dict, site: str, site_dir: Path):
    """
    ea_by_ef: {'-0.6': [Ea1..Ea4], '0.0': [...], '0.6': [...]}
    Save to site_dir/Ea_grouped_<site>.png
    """
    configure_style()

    labels = ["Ea1", "Ea2", "Ea3", "Ea4"]
    x = np.arange(4)

    fig, ax = plt.subplots(figsize=(5, 5))

    ef_order = ["0.6", "0.0", "-0.6"]
    n_ef = len(ef_order)
    group_span = 0.90
    bar_width = group_span / n_ef
    start = -group_span/2 + bar_width/2
    offsets = [start + i*bar_width for i in range(n_ef)]

    for i, ef in enumerate(ef_order):
        vals = ea_by_ef.get(ef, [np.nan]*4)
        ax.bar(x + offsets[i], vals, width=bar_width,
               color=color_for_data_ef(ef), edgecolor="black",
               label=label_for_data_ef(ef))

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=20)
    ax.set_ylabel("Ea (eV)",fontsize=20)
    ax.yaxis.set_major_locator(LinearLocator(5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(frameon=False, loc="best", fontsize=15)

    fig.tight_layout()
    outfile = site_dir / f"Ea_grouped_{site}.png"
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outfile


def plot_ea_trends(ea_by_ef: dict, site: str, site_dir: Path):
    """
    Make four scatter figures: Ea1_<site>.png .. Ea4_<site>.png in site_dir
    X: REAL EF [-0.6, 0.0, +0.6]; Y: Ea (eV)
    """
    configure_style()
    data_order = ["0.6", "0.0", "-0.6"]
    xa = [-float(k) for k in data_order]
    ys = []
    for k_ea in range(4):
        vals = []
        for k in data_order:
            ea_list = ea_by_ef.get(k, [float("nan")]*4)
            val = ea_list[k_ea] if (ea_list[k_ea] is not None) else float("nan")
            vals.append(val)
        ys.append(vals)
    
    labels = ["Ea1","Ea2","Ea3","Ea4"]
    palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, lab in enumerate(labels):
        fig, ax = plt.subplots(figsize=(6.0, 4.5))
        ax.scatter(xa, ys[idx], s=48, color=palette[idx], edgecolors='black', linewidths=0.8)
        ax.set_xlabel("EF (V/Å)")
        ax.set_ylabel("Ea (eV)")
        ax.set_xticks([-0.6, 0.0, 0.6])
        ax.set_xticklabels(["-0.6", "0.0", "+0.6"])
        ax.yaxis.set_major_locator(LinearLocator(5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.tick_params(axis='both', which='major', labelsize=20)
        fig.tight_layout()
        
        out_png = site_dir / f"{lab}_{site}.png"
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"    Saved -> {out_png.resolve()}")


# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Extract activation energies from mkm_out.csv and plot grouped barchart + trends.")
    ap.add_argument("-s", "--site", default=None, help="site id, e.g., 8_11_17 (optional; if omitted, processes all sites in summary)")
    ap.add_argument("-i", "--input", default="mkm_outputs/mkm_out.csv", help="mkm_out.csv file path")
    args = ap.parse_args()

    summary_file = Path(args.input)
    if not summary_file.is_file():
        print(f"[ERROR] Summary file not found: {summary_file}", file=sys.stderr)
        sys.exit(1)
    
    # Read the CSV
    try:
        df_summary = pd.read_csv(summary_file)
    except Exception as e:
        print(f"[ERROR] Failed to read CSV: {e}", file=sys.stderr)
        sys.exit(1)
    
    # If site is provided, process only that one
    if args.site:
        process_site_from_csv(args.site, df_summary)
        return
    
    # Otherwise, process all unique sites in the summary
    sites = df_summary['site'].unique()
    print(f"Found {len(sites)} site(s) to process:\n")
    for site in sorted(sites):
        print(f"Processing: {site}")
        try:
            process_site_from_csv(site, df_summary)
        except Exception as e:
            print(f"  [WARN] Failed: {e}")
        print()
    
    print("All sites processed.")


def process_site_from_csv(site: str, df_summary: pd.DataFrame):
    """Process a single site from GA_prediction_summary.csv."""
    # Filter rows for this site
    df_site = df_summary[df_summary['site'] == site].copy()
    
    if df_site.empty:
        print(f"  [SKIP] No records found for site: {site}")
        return
    
    # Activation energy columns in CSV
    ea_cols = [f"Ea{i}" for i in range(1, 5)]
    if not all(col in df_site.columns for col in ea_cols):
        print(f"  [SKIP] Missing activation energy columns")
        return
    
    # Collect Ea values by EF
    ea_by_ef = {}
    for ef_val in sorted(df_site['EF'].unique()):
        ef_norm = norm_ef(str(ef_val))
        row = df_site[df_site['EF'] == ef_val].iloc[0]
        ea_list = [row[col] for col in ea_cols]
        ea_by_ef[ef_norm] = ea_list
    
    # Ensure site directory exists
    site_dir = Path("mkm_inputs") / site
    site_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate grouped barchart
    outfile = plot_grouped_ea(ea_by_ef, site, site_dir)
    print(f"  Saved grouped chart -> {outfile.resolve()}")
    
    # Generate trend plots
    plot_ea_trends(ea_by_ef, site, site_dir)

if __name__ == "__main__":
    main()
