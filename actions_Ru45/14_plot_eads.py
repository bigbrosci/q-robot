#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# ---------- utilities ----------
def norm_ef(ef: str) -> str:
    return f"{float(ef):.1f}"

# ---------- plotting ----------
def configure_jacs_style():
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "legend.fontsize": 16,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
    })

COLOR_MAP = {
    "0.6": "#1f77b4",
    "0.0":  "#ff7f0e",
    "-0.6":  "#d62728",
}
ORDER = ["0.6", "0.0", "-0.6"]

def fmt_subscripts(labels):
    return [re.sub(r'([A-Za-z\)\]])(\d+)', r'\1$_{\2}$', s) for s in labels]

def plot_grouped_bars(df: pd.DataFrame, site: str, outfile: Path,
                      ylabel="Adsorption energy (eV)", xlabel="Species",
                      width=0.25, figsize=(10, 4)):
    configure_jacs_style()
    ef_cols = [c for c in df.columns if c.startswith("EF_")]
    key_to_col = {f"{float(c.split('_',1)[1]):.1f}": c for c in ef_cols}

    states = df.index.tolist()
    x = np.arange(len(states))
    fig, ax = plt.subplots(figsize=figsize)

    n = len(ORDER)
    for i, key in enumerate(ORDER):
        col = key_to_col.get(key)
        if col is None:
            continue
        offset = (i - (n - 1) / 2) * width
        real_ef = -float(key)
        label = f"EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å"
        ax.bar(x + offset, df[col].values, width=width,
               label=label, color=COLOR_MAP.get(key))

    xticks = [s.replace("_ads", "").replace('E_', '') + '*' for s in states]
    xticks =  fmt_subscripts(xticks)
    ax.set_xticks(x)
    ax.set_xticklabels(xticks, rotation=0)

    ax.set_ylabel(ylabel, fontsize = 20)
    ax.set_xlabel(xlabel, fontsize = 20)
    ax.legend(frameon=False)
    ax.grid(False)
    ax.yaxis.set_major_locator(LinearLocator(5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.tick_params(axis='both', which='major', labelsize=20)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Extract adsorption energies from mkm_out.csv and plot grouped barchart.")
    ap.add_argument("-s", "--site", default=None, help="site id, e.g., 40_43_44 (optional; if omitted, processes all sites in summary)")
    ap.add_argument("-i", "--input", default="mkm_outputs/mkm_out.csv", help="mkm_out.csv file path")
    ap.add_argument("-o", "--out", default=None, help="output image name (PNG)")
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
        process_site_from_csv(args.site, df_summary, args.out)
        return
    
    # Otherwise, process all unique sites in the summary
    sites = df_summary['site'].unique()
    print(f"Found {len(sites)} site(s) to process:\n")
    for site in sorted(sites):
        print(f"Processing: {site}")
        try:
            process_site_from_csv(site, df_summary, None)
        except Exception as e:
            print(f"  [WARN] Failed: {e}")
        print()
    
    print("All sites processed.")


def process_site_from_csv(site: str, df_summary: pd.DataFrame, out_file: str):
    """Process a single site from GA_prediction_summary.csv."""
    # Filter rows for this site
    df_site = df_summary[df_summary['site'] == site].copy()
    
    if df_site.empty:
        print(f"  [SKIP] No records found for site: {site}")
        return
    
    # Adsorption energy columns in CSV (note: no E_N2 in current version)
    preferred = ["E_N", "E_H", "E_NH", "E_NH2", "E_NH3"]
    available = [col for col in preferred if col in df_site.columns]
    
    if not available:
        print(f"  [SKIP] No adsorption energy columns found")
        return
    
    # Create output DataFrame with EF values as columns
    data = {}
    for ef_val in sorted(df_site['EF'].unique()):
        ef_norm = norm_ef(str(ef_val))
        row = df_site[df_site['EF'] == ef_val].iloc[0]
        col = f"EF_{ef_norm}"
        data[col] = [row[sp] for sp in available]
    
    df = pd.DataFrame(data, index=available)
    
    # Ensure site directory exists
    site_dir = Path("mkm_inputs") / site
    site_dir.mkdir(parents=True, exist_ok=True)
    
    # Save table in site directory
    csv_out = site_dir / f"adsorption_by_EF_{site}.csv"
    df.to_csv(csv_out)
    print(f"  Saved table -> {csv_out.resolve()}")

    # Generate plot in site directory
    img_out = Path(out_file) if out_file else site_dir / f"adsorption_barchart_{site}.png"
    plot_grouped_bars(df, site=site, outfile=img_out)
    print(f"  Saved figure -> {img_out.resolve()}")

if __name__ == "__main__":
    main()
