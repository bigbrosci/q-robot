#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grouped bar chart of TOP-3 |DRC| per EF (x = EFs; bars colored by Reaction) with formatted reaction labels.

Change in this version:
- Bars with |DRC| < 0.1 are not plotted (threshold can be adjusted via DRC_DRAW_THRESHOLD).

Usage:
    python3 12_plot_drc_onesite.py [SITE_DIR]
Output:
    <SITE_DIR>/drc_top2_grouped_by_rxn_fmt.png
"""

import os
import sys
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------------- Config ----------------
DRC_DRAW_THRESHOLD = 0.1  # Do not draw bars with |DRC| < this threshold

# ---------------- Style (Times New Roman + MathText) ----------------
plt.rcParams.update({
    "font.family": "Times New Roman",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "legend.fontsize": 12,
    "axes.labelsize": 18,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
})

# ---------------- Safe palette for reactions (avoid EF colors: red/orange/blue) ----------------
SAFE_REACTION_COLORS = [
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#bcbd22",  # olive
    "#7f7f7f",  # gray
    "#2ca02c",  # green
    "#8c6d31",  # mustard/brown
    "#525252",  # dark gray
    "#bdb76b",  # dark khaki
    "#17a773",  # teal-green (not blue)
    "#393b79",  # indigo
    "#637939",  # olive green
    "#843c39",  # dark brick
    "#7b4173",  # plum
]

# Optional overrides: map raw Reaction strings to specific colors
REACTION_COLOR_OVERRIDES: Dict[str, str] = {
    # "N2 + RU(T) <=> N2(T)": "#9467bd",
}

# ---------------- Utilities ----------------
def _detect_ef_folders(site_dir: Path) -> List[str]:
    """Return EF_* folder names sorted so that REAL EF ascends: -0.6, 0.0, +0.6."""
    tags = []
    for name in os.listdir(site_dir):
        if name.startswith("EF_"):
            try:
                _ = float(name.replace("EF_", ""))
                tags.append(name)
            except Exception:
                pass
    # Real EF = - (data EF), so sort by -dataEF
    tags.sort(key=lambda s: -float(s.replace("EF_", "")))
    return tags


def _real_ef_label(ef_folder: str) -> str:
    """EF_-0.6 -> '+0.6'; EF_0.6 -> '-0.6' (returns string with sign)."""
    try:
        data_val = float(ef_folder.replace("EF_", ""))
        return f"{-data_val:+.1f}"
    except Exception:
        return ef_folder



def _col_by_name(df: pd.DataFrame, target: str) -> str:
    """Find column by case-insensitive name; fall back to contains-match."""
    t = target.strip().lower()
    norm = {str(c).strip().lower(): c for c in df.columns}
    if t in norm:
        return norm[t]
    for k, v in norm.items():
        if t in k:
            return v
    raise KeyError(f"Column '{target}' not found. Available: {list(df.columns)}")


def _read_drc_csv(path: Path) -> pd.DataFrame:
    """Read drc_results.csv; drop 'initial' row; keep [ReactionID, Reaction, DRC]."""
    df = pd.read_csv(path)
    first_col = df.columns[0]
    df = df[df[first_col].astype(str).str.strip().str.lower() != "initial"].copy()

    col_drc = _col_by_name(df, "DRC")
    try:
        col_rxn = _col_by_name(df, "Reaction")
    except KeyError:
        col_rxn = df.columns[0]  # fallback to the first column

    keep = [first_col, col_drc, col_rxn]
    df = df[keep].copy()
    df.rename(columns={df.columns[0]: "ReactionID", col_drc: "DRC", col_rxn: "Reaction"}, inplace=True)
    df["DRC"] = pd.to_numeric(df["DRC"], errors="coerce")
    df = df.dropna(subset=["DRC"])
    return df.reset_index(drop=True)


def _top2_by_abs_drc(df: pd.DataFrame) -> pd.DataFrame:
    """Return Top-3 rows by |DRC| (keep signed DRC)."""
    d = df.copy()
    d["absDRC"] = d["DRC"].abs()
    d = d.sort_values("absDRC", ascending=False).head(2).reset_index(drop=True)
    return d[["ReactionID", "Reaction", "DRC"]]


def _shorten(s: str, n: int = 60) -> str:
    if not isinstance(s, str):
        return ""
    return s if len(s) <= n else s[:n] + "…"

def _format_reaction_math(s: str) -> str:
    """
    Format reaction string for mathtext:
      - RU(T), (T) -> *
      - NH3 / NH2 / N2 / H2 -> with subscripts
      - <=> -> \\rightleftharpoons
      - Wrap with $...$
    """
    import re
    if not isinstance(s, str):
        return s
    s = s.replace("RU(T)", "*")
    s = s.replace("(T)", "*")
    s = re.sub(r'(?<![A-Za-z0-9_])NH3(?![A-Za-z0-9_])', r'NH{_3}', s)
    s = re.sub(r'(?<![A-Za-z0-9_])NH2(?![A-Za-z0-9_])', r'NH{_2}', s)
    s = re.sub(r'(?<![A-Za-z0-9_])N2(?![A-Za-z0-9_])',  r'N{_2}',  s)
    s = re.sub(r'(?<![A-Za-z0-9_])H2(?![A-Za-z0-9_])',  r'H{_2}',  s)
    s = s.replace("<=>", r"\rightleftharpoons")
    return f"${s}$"


# ---------------- Main plotting ----------------
def plot_grouped_by_rxn(site_dir: str = ".") -> None:
    site = Path(site_dir)
    if not site.is_dir():
        print(f"[ERROR] '{site_dir}' is not a directory.")
        sys.exit(1)

    ef_folders = _detect_ef_folders(site)
    if not ef_folders:
        ef_folders = ["EF_-0.6", "EF_0.0", "EF_0.6"]  # fallback

    # Collect Top-3 per EF
    per_ef_top2: Dict[str, pd.DataFrame] = {}
    for ef in ef_folders:
        csv_path = site / ef / "outputs" / "drc_results.csv"
        if not csv_path.is_file():
            print(f"[WARN] Missing: {csv_path}")
            continue
        try:
            df = _read_drc_csv(csv_path)
            top2 = _top2_by_abs_drc(df)
            per_ef_top2[ef] = top2
            # console summary
            print(f"\n[{ef}] Top-3 by |DRC|:")
            for i, r in top2.iterrows():
                print(f"  Top{i+1}: DRC={r['DRC']:+.4f} | {r['ReactionID']} | {r['Reaction']}")
        except Exception as e:
            print(f"[WARN] Failed to parse {csv_path}: {e}")

    if not per_ef_top2:
        print("❌ No valid DRC data found.")
        return

    # Panels in ascending REAL EF
    panels = _detect_ef_folders(site) or ef_folders
    m = len(panels)
    ef_labels = [_real_ef_label(ef) for ef in panels]

    # Build reaction→color mapping
    rxn_order = []
    for ef in panels:
        top2 = per_ef_top2.get(ef)
        if top2 is None:
            continue
        for i in range(min(3, len(top2))):
            rxn = str(top2.loc[i, "Reaction"])
            if rxn not in rxn_order:
                rxn_order.append(rxn)

    rxn_color: Dict[str, str] = {}
    rxn_color.update(REACTION_COLOR_OVERRIDES)
    palette_idx = 0
    for rxn in rxn_order:
        if rxn not in rxn_color:
            rxn_color[rxn] = SAFE_REACTION_COLORS[palette_idx % len(SAFE_REACTION_COLORS)]
            palette_idx += 1

    # Data cube: ranks x EFs
    heights = np.full((3, m), np.nan)
    reactions = [["", "", ""] for _ in range(m)]
    for j, ef in enumerate(panels):
        top2 = per_ef_top2.get(ef)
        if top2 is None:
            continue
        for rank in range(min(3, len(top2))):
            val = float(top2.loc[rank, "DRC"])
            # apply threshold: suppress bars with small |DRC|
            if abs(val) < DRC_DRAW_THRESHOLD:
                val = np.nan
            heights[rank, j] = val
            reactions[j][rank] = str(top2.loc[rank, "Reaction"])

    # Plot
    x = np.arange(m)
    n_ranks = 3
    group_span = 0.90
    width = group_span / n_ranks
    start = -group_span/2 + width/2
    offsets = [start + i*width for i in range(n_ranks)]

    fig, ax = plt.subplots(figsize=(max(5, 2.0*m), 5))
    drawn_rxns = set()
    for j in range(m):
        for r in range(3):
            val = heights[r, j]
            if np.isnan(val):
                continue
            rxn = reactions[j][r]
            ax.bar(x[j] + offsets[r], val, width=width,
                   color=rxn_color.get(rxn, "#7f7f7f"), edgecolor="black")
            drawn_rxns.add(rxn)

    # X
    ax.set_xticks(x)
    ax.set_xticklabels([("0.0" if round(float(v), 1) == 0 else f"{round(float(v), 1):+.1f}") for v in ef_labels])
    ax.set_xlabel("EF (V/Å)", fontsize=20)
    ax.set_ylabel("DRC")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.tick_params(axis="both", labelsize=18)

    # Legend (only reactions that actually have at least one drawn bar)
    legend_order = [r for r in rxn_order if r in drawn_rxns]
    if legend_order:
        legend_patches = [Patch(facecolor=rxn_color[r], edgecolor="black",
                                label=_format_reaction_math(_shorten(r))) for r in legend_order]
        ncol = 1 if len(legend_patches) <= 8 else 2
        ax.legend(handles=legend_patches, frameon=False, loc="best", ncol=ncol, fontsize=13)

    fig.tight_layout()
    out = site / "drc_top2_grouped_by_rxn_fmt.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n✅ Saved figure: {out}")
    print(f"   (Bars with |DRC| < {DRC_DRAW_THRESHOLD} were suppressed.)")


def main():
    # If a specific site is provided as argument, process only that one
    if len(sys.argv) > 1:
        site_dir = sys.argv[1]
        plot_grouped_by_rxn(site_dir)
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
        plot_grouped_by_rxn(str(site_path))
        print()
    
    print("✅ All sites processed.")


if __name__ == "__main__":
    main()
