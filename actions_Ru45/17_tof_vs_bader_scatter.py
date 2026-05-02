#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
17_tof_vs_bader_scatter.py

Plot TOF vs. Bader charge for three EFs (real EF = - data EF):
- Read tog_*.csv and bader_*.dat files, align by Index (1-based vs 0-based)
- Create individual EF plots (each EF one figure)
- Create a combined plot containing all three EFs on the same axes
"""

import argparse
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

# ---------- EF color scheme (by REAL EF) ----------
COL_POS  = "#d62728"  # real +0.6  (maps from data EF_-0.6)
COL_ZERO = "#ff7f0e"  # real  0.0
COL_NEG  = "#1f77b4"  # real -0.6  (maps from data EF_+0.6)

def color_for_real_ef(real: float) -> str:
    if abs(real) < 5e-4:
        return COL_ZERO
    return COL_POS if real > 0 else COL_NEG

# ---------- Style ----------
def configure_style():
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "legend.fontsize": 14,
        "axes.labelsize": 20,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
    })

# ---------- IO ----------
def load_tog_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if not {"atom_index", "mean_TOF"}.issubset(df.columns):
        raise ValueError(f"Missing columns in {path}: need 'atom_index' and 'mean_TOF'")
    out = pd.DataFrame({
        "Index": df["atom_index"].astype(int) + 1,   # convert to 1-based to match Bader
        "TOF": pd.to_numeric(df["mean_TOF"], errors="coerce")
    })
    return out.dropna(subset=["TOF"]).reset_index(drop=True)

def load_bader_dat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    def pick(colname):
        tgt = colname.strip().lower()
        for c in df.columns:
            if str(c).strip().lower() == tgt:
                return c
        for c in df.columns:
            if tgt in str(c).strip().lower():
                return c
        raise KeyError(f"Column '{colname}' not found in {path}. Available: {list(df.columns)}")
    col_idx = pick("index")
    col_charge = pick("charge")
    out = pd.DataFrame({
        "Index": pd.to_numeric(df[col_idx], errors="coerce").astype("Int64"),
        "Charge": pd.to_numeric(df[col_charge], errors="coerce")
    }).dropna(subset=["Index", "Charge"])
    out["Index"] = out["Index"].astype(int)
    return out.reset_index(drop=True)

def parse_data_ef_from_name(p: Path) -> float:
    """Extract EF value (e.g., ±0.6) from file name."""
    m = re.search(r"([+-]?\d*\.?\d+)", p.stem)
    if not m:
        raise ValueError(f"Cannot parse EF from file name: {p.name}")
    return float(m.group(1))

# ---------- Plot (single EF) ----------
def scatter_one(df_merge: pd.DataFrame, real_ef: float, outpng: Path):
    configure_style()
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.scatter(df_merge["Charge"].values, df_merge["TOF"].values,
               s=40, color=color_for_real_ef(real_ef), edgecolors="black", linewidths=0.6)
    ax.set_xlabel("Bader charge")
    ax.set_ylabel("TOF")
    title = "EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å"
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {outpng}  (n={len(df_merge)})")

# ---------- Plot (combined all EFs) ----------
def plot_combined(panels: list[tuple[float, pd.DataFrame]], outpng: Path):
    """
    Combine three EF datasets into one figure, each EF as a color-coded scatter.
    """
    configure_style()
    panels = [(r, df) for (r, df) in panels if df is not None and not df.empty]
    if not panels:
        print("[WARN] No data to plot in combined figure.")
        return

    # Sort by EF
    order = {-0.6: 0, 0.0: 1, 0.6: 2}
    panels.sort(key=lambda x: order.get(x[0], 99))

    fig, ax = plt.subplots(figsize=(7, 5))
    for real_ef, df in panels:
        label = "EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å"
        ax.scatter(df["Charge"], df["TOF"], s=40,
                   color=color_for_real_ef(real_ef), edgecolors="black", linewidths=0.6,
                   label=label)

    ax.set_xlabel("Bader charge")
    ax.set_ylabel("TOF")
    ax.legend(frameon=False, loc="upper left")
    ax.set_xlim([-0.3,0.3])
    ax.set_ylim([0,13])
    ax.yaxis.set_major_locator(LinearLocator(5))
    ax.tick_params(axis='both', which='major', labelsize=18)
    fig.tight_layout()
    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined figure -> {outpng}")

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Scatter of TOF vs Bader charge for three EFs (real EF = - data EF).")
    ap.add_argument("--tog", nargs="*", default=["mkm_outputs/tog_-0.6_sim.csv", "mkm_outputs/tog_0.0_sim.csv", "mkm_outputs/tog_+0.6_sim.csv"],
                    help="TOF csv files (default: mkm_outputs/tog_*.csv)")
    ap.add_argument("--bader", nargs="*", default=["bader_-0.6.dat", "bader_0.0.dat", "bader_0.6.dat"],
                    help="Bader dat files.")
    ap.add_argument("--out-prefix", default="tof_vs_bader_EF_", help="Prefix for individual figures.")
    ap.add_argument("--combined-out", default="tof_vs_bader.png", help="Output file for combined plot.")
    args = ap.parse_args()

    # Map data EF → path
    tog_by_data_ef = {parse_data_ef_from_name(Path(f)): Path(f) for f in args.tog if Path(f).is_file()}
    bader_by_data_ef = {parse_data_ef_from_name(Path(f)): Path(f) for f in args.bader if Path(f).is_file()}

    targets = [-0.6, 0.0, 0.6]  # REAL EF
    panels = []

    for real in targets:
        data_value = -real  # real = - data
        key_tog = min(tog_by_data_ef.keys(), key=lambda x: abs(x - data_value))
        key_bdr = min(bader_by_data_ef.keys(), key=lambda x: abs(x - data_value))

        df_tog = load_tog_csv(tog_by_data_ef[key_tog])
        df_bdr = load_bader_dat(bader_by_data_ef[key_bdr])
        df = pd.merge(df_tog, df_bdr, on="Index", how="inner").sort_values("Index").reset_index(drop=True)

        if df.empty:
            print(f"[WARN] No overlap for EF {real:+.1f}")
            continue

        sign = "+" if real > 0 else "-" if real < 0 else ""
        outname = f"{args.out_prefix}{sign}{abs(real):.1f}.png".replace("+-", "-").replace("--", "-")
        scatter_one(df, real_ef=real, outpng=Path(outname))
        panels.append((real, df))

    plot_combined(panels, outpng=Path(args.combined_out))

if __name__ == "__main__":
    main()

