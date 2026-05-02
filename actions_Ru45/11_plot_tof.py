#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read drc_summary_output.csv (columns include: path, TOF, ... where path like '0_11_34/EF_-0.6'),
aggregate per-atom TOFs for each EF, write dict & per-atom mean CSV, and plot heatmaps
(xy/xz/yz) using mean TOF as dq_substrate.

Outputs per EF value (e.g., -0.6, 0.0, 0.6):
  - tof_<ef>_dict.txt     # Python-like dict text: {0: [...], 1: [...], ...}
  - tog_<ef>_sim.csv      # columns: atom_index, mean_TOF

Plots (per EF): tof_xy_EF<ef>, tof_xz_EF<ef>, tof_yz_EF<ef>

NOTE:
- Assumes a function `draw_cell(ref_geom, layer_indices, dq_substrate, basename, plane, ef_str, vmin, vmax)`
  is available (either imported or defined elsewhere in your project).
"""
# === Standard library ===
import os
import re
from pathlib import Path
from pprint import pformat
from typing import Dict, List, Tuple
from collections import defaultdict
from sys import argv

# === Third-party ===
import numpy as np
import pandas as pd
from PIL import Image

# === Matplotlib ===
import matplotlib.pyplot as plt
plt.rc('font', size=18)
from matplotlib import cm, colormaps
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D  # if you use 3D plots
from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator

# === ASE (Atomic Simulation Environment) ===
from ase import Atom, Atoms
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp
from ase.neighborlist import natural_cutoffs
from ase.visualize import view

def configure_jacs_style():
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
    
# ---------------- Files ----------------
CSV = Path("mkm_outputs/mkm_out.csv")
POSCAR_SLAB = Path("slab/POSCAR")

# ---------------- Helpers ----------------
def choose_geometry_file() -> Path:
    """Load POSCAR from slab folder for plotting."""
    if POSCAR_SLAB.exists():
        return POSCAR_SLAB
    raise FileNotFoundError(f"POSCAR not found at {POSCAR_SLAB.resolve()}")

def parse_path(path_str: str) -> Tuple[str, float]:
    """
    'mkm_inputs/8_11_17/EF_0.6' -> ('8_11_17', 0.6)
    Handles paths with mkm_inputs prefix. Site key is the folder before EF_*.
    """
    s = str(path_str).strip()
    parts = s.split("/")
    
    # Find the EF token (should be the last part or near the end)
    ef_token = None
    site_key = None
    
    for i, part in enumerate(parts):
        if part.startswith("EF_"):
            ef_token = part
            # Site key is the part before EF
            if i > 0:
                site_key = parts[i - 1]
            break
    
    if not ef_token or not site_key:
        raise ValueError(f"Bad path format: {s!r} (expected format like 'mkm_inputs/8_11_17/EF_0.6')")
    
    # Extract EF value (handle both EF_-0.6 and EF_0.6)
    ef_str = ef_token.split("EF_")[-1]
    ef_val = float(ef_str)
    
    return site_key, ef_val

def site_to_indices(site_key: str) -> List[int]:
    """'13_29_45' -> [12, 28, 44] (convert from 1-based to 0-based indexing)."""
    sep = "_" if "_" in site_key else " "
    idxs = [int(x) - 1 for x in site_key.split(sep) if x != ""]  # Convert 1-based to 0-based
    if len(idxs) < 1:
        raise ValueError(f"Empty site key: {site_key!r}")
    return idxs

def build_per_atom_dict_for_ef(df: pd.DataFrame, natoms: int, ef_val: float) -> Dict[int, List[float]]:
    """
    For a given EF value, build: {0: [], 1: [], ..., natoms-1: []}
    Append each row's TOF to all atoms in its 'site' triplet.
    """
    tof_dict: Dict[int, List[float]] = {i: [] for i in range(natoms)}
    sub = df[df["EF_value"] == ef_val].copy()
    for _, row in sub.iterrows():
        site_atoms = site_to_indices(row["site"])
        tof_val = float(row["TOF"])
        for a in site_atoms:
            if 0 <= a < natoms:
                tof_dict[a].append(tof_val)
            else:
                raise IndexError(f"Atom {a} out of range 0..{natoms-1} (site={row['site']})")
    return tof_dict

def write_tof_dict_text(tof_dict: Dict[int, List[float]], ef_val: float, output_dir: Path = Path(".")) -> Path:
    """Write dict to 'tof_<ef>_dict.txt' with numeric EF in filename."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fname = f"tof_{ef_val:+.1f}_dict.txt".replace("+0.0", "0.0").replace("-0.0", "0.0")
    filepath = output_dir / fname
    text = pformat(tof_dict, width=120, compact=False)
    filepath.write_text(text + "\n", encoding="utf-8")
    return filepath

def write_mean_csv(tof_dict: Dict[int, List[float]], ef_val: float, output_dir: Path = Path(".")) -> Path:
    """Write 'tog_<ef>_sim.csv' with columns: atom_index, mean_TOF (empty -> 0.0)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for idx in sorted(tof_dict.keys()):
        vals = tof_dict[idx]
        mean_val = float(np.max(vals)) if len(vals) > 0 else 0.0
        rows.append((idx, mean_val))
    df_out = pd.DataFrame(rows, columns=["atom_index", "mean_TOF"])
    fname = f"tog_{ef_val:+.1f}_sim.csv".replace("+0.0", "0.0").replace("-0.0", "0.0")
    filepath = output_dir / fname
    df_out.to_csv(filepath, index=False)
    return filepath

def load_mean_vector(csv_name: str, natoms: int) -> np.ndarray:
    """Read per-atom mean CSV and build a length-natoms vector."""
    df = pd.read_csv(csv_name)
    if not {"atom_index", "mean_TOF"}.issubset(df.columns):
        raise ValueError(f"{csv_name} must contain columns 'atom_index' and 'mean_TOF'")
    arr = np.zeros(natoms, dtype=float)
    idx = pd.to_numeric(df["atom_index"], errors="coerce").astype(int)
    val = pd.to_numeric(df["mean_TOF"], errors="coerce").fillna(0.0).to_numpy()
    bad = (idx < 0) | (idx >= natoms)
    if bad.any():
        raise IndexError(f"{csv_name}: atom_index out of range: {idx[bad].tolist()}")
    arr[idx.to_numpy()] = val
    return arr

def draw_cell(
    s, indices, dq, name, plane, ef, min_bader, max_bader,
    use_discrete=True,          # True=离散分箱；False=连续色带
    palette_base='rainbow',       # 推荐“彩虹”效果：'turbo'/'rainbow'/'gist_rainbow'/'nipy_spectral'/'hsv'
    n_bins=10,                  # 离散分箱数
    reverse=False               # 是否反转色带
):

    import numpy as np
    # --- 采样坐标与高度 ---
    if plane == 'xy':
        xy = [i.position[:2] for i in s]
        z = [i.position[2] for i in s]
    elif plane == 'xz':
        xy = [[i.position[0], i.position[2]] for i in s]
        z = [i.position[1] for i in s]
    elif plane == 'yz':
        xy = [[i.position[1], i.position[2]] for i in s]
        z = [i.position[0] for i in s]
    else:
        raise ValueError("plane 必须是 'xy'/'xz'/'yz' 之一")

    cov_radii = natural_cutoffs(s, mult=1.0)
    atoms_all = list(zip(xy, z, cov_radii, list(dq)))

    # 只取需要显示的原子，并按“高度”排序（确保遮挡正确）
    selected_atoms = [atoms_all[i] for i in indices]
    selected_atoms.sort(key=lambda a: a[1])
    dq_bader = np.asarray([a[3] for a in selected_atoms], dtype=float)

    # --- 盒矢量与画布布局 ---
    cell = s.get_cell()
    cell_width = float(np.linalg.norm(cell[0]))
    cell_height = float(np.linalg.norm(cell[1]))

    fig_height = 10
    fig_width = 10

    top_pad = 0.5
    bot_pad = top_pad * 1.5
    scaled_width = (fig_height - top_pad - bot_pad) * (cell_width / cell_height if cell_height != 0 else 1.0)
    left_pad = (fig_width - scaled_width) / 3
    right_pad = left_pad * 2

    # 视窗范围（按需修改）
    xlims = [3, 17]
    ylims = [3, 17]

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    direct_1 = plane[0] + r' (Å)'
    direct_2 = plane[1] + r' (Å)'
    ax.set(xlim=xlims, ylim=ylims, xlabel=direct_1, ylabel=direct_2)
    ax.set_aspect('equal', adjustable='box')
    fig.subplots_adjust(
        left=left_pad/fig_width,
        right=1-right_pad/fig_width,
        bottom=bot_pad/fig_height,
        top=1-top_pad/fig_height
    )

    # --- 彩虹色图 & 归一化 ---
    base = colormaps[palette_base]
    if reverse:
        base = base.reversed()

    if use_discrete:
        # 离散：从 0~1 均匀采样 n_bins 个颜色
        colors = [base(x) for x in np.linspace(0.0, 1.0, n_bins, endpoint=False)]
        cmap = ListedColormap(colors, name=f"{palette_base}_{n_bins}{'_r' if reverse else ''}")
        # 分箱
        boundaries = np.linspace(min_bader, max_bader, n_bins + 1)
        norm = BoundaryNorm(boundaries, cmap.N, clip=True)
        tick_locs = (boundaries[:-1] + boundaries[1:]) / 2.0
    else:
        # 连续
        cmap = base
        norm = Normalize(vmin=min_bader, vmax=max_bader)
        tick_locs = np.linspace(min_bader, max_bader, 5)

    # --- 画原子（圆形） ---
    circles = [patches.Circle(xy=tuple(a[0]), radius=a[2]) for a in selected_atoms]
    patcoll = PatchCollection(circles, edgecolors='black', linewidths=2, cmap=cmap, norm=norm)
    patcoll.set_array(dq_bader)
    if not use_discrete:
        patcoll.set_clim([min_bader, max_bader])
    ax.add_collection(patcoll)

    # --- 颜色条 ---
    cbar_left   = 1 - right_pad/fig_width * 0.9
    cbar_bottom = bot_pad/fig_height
    cbar_width  = right_pad/fig_width * 0.6
    cbar_height = 1 - top_pad/fig_height - bot_pad/fig_height
    cbar_ax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, cax=cbar_ax)
    cbar.set_label("log10(TOF / mol s$^{-1}$ mol$_{cat}^{-1}$)", fontsize=40, labelpad=32)
    # Force integer log10 ticks on colorbar (only this block changed)
    from matplotlib.ticker import FixedLocator, FormatStrFormatter
    import numpy as np
    tmin = int(np.floor(min_bader))
    tmax = int(np.ceil(max_bader))
    int_ticks = np.array([-18, -15, -12, -9, -6, -3, 0, 3])
    cbar.set_ticks(int_ticks)
    cbar.ax.yaxis.set_major_locator(FixedLocator(int_ticks))
    cbar.ax.yaxis.set_minor_locator(NullLocator())
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    cbar.outline.set_linewidth(3.0)
    cbar.ax.tick_params(axis='y', which='both', width=2.5, length=9, pad=20, labelsize=40, direction='out')
    cbar.ax.tick_params(which='minor', length=0)

    # --- EF 标注 ---
    ef_display = 0.0 if abs(-float(ef)) < 1e-5 else -float(ef)
    ax.text(0.95, 0.95, f'EF = {ef_display} V/Å',
            transform=ax.transAxes, ha='right', va='top', fontsize=40)

    # 坐标轴样式
    ax.xaxis.label.set_size(50)
    ax.yaxis.label.set_size(50)
    ax.set_xticks(np.linspace(xlims[0], xlims[1], 5))
    ax.set_yticks(np.linspace(ylims[0], ylims[1], 5))
    ax.tick_params(axis='both', which='major', labelsize=50)
    ax.tick_params(axis='x', pad=10)
    ax.tick_params(axis='y', pad=5)

    output_file = f'{name}_{ef_display}.png'
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close(fig)

# ---------------- Main ----------------

# ---------- New helpers: draw on existing axes & triptych composer ----------

def _draw_panel_on_ax(ax, s, indices, dq_log10, plane, ef_display, tmin, tmax,
                      use_discrete=True, palette_base='rainbow', n_bins=8, reverse=False):
    """Render one EF panel on given Matplotlib Axes (no colorbar). dq_log10 is log10(TOF)."""
    import numpy as np
    from matplotlib import patches
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm, Normalize
    from matplotlib import colormaps

    # Coordinates
    if plane == 'xy':
        xy = [i.position[:2] for i in s]; z = [i.position[2] for i in s]
    elif plane == 'xz':
        xy = [[i.position[0], i.position[2]] for i in s]; z = [i.position[1] for i in s]
    elif plane == 'yz':
        xy = [[i.position[1], i.position[2]] for i in s]; z = [i.position[0] for i in s]
    else:
        raise ValueError("plane must be 'xy'/'xz'/'yz'")

    cov_radii = natural_cutoffs(s, mult=1.0)
    atoms = list(zip(xy, z, cov_radii, list(dq_log10)))
    selected_atoms = [atoms[i] for i in indices]
    selected_atoms = sorted(selected_atoms, key=lambda atom: atom[1])
    values = np.asarray([a[3] for a in selected_atoms], dtype=float)

    # Axis limits and labels
    xlims = [3, 17]; ylims = [3, 17]
    ax.set(xlim=xlims, ylim=ylims, xlabel=plane[0] + ' (Å)', ylabel=plane[1] + ' (Å)')
    ax.set_aspect('equal', adjustable='box')

    # Colormap/norm
    base = colormaps[palette_base]
    if reverse and hasattr(base, "reversed"):
        base = base.reversed()
    if use_discrete:
        n_bins_eff = max(1, int(tmax - tmin))
        colors = [base(x) for x in np.linspace(0.0, 1.0, n_bins_eff, endpoint=False)]
        cmap = ListedColormap(colors, name=f"{palette_base}_{n_bins_eff}{'_r' if reverse else ''}")
        boundaries = np.arange(tmin, tmax + 1, 1.0)
        norm = BoundaryNorm(boundaries, cmap.N, clip=True)
    else:
        cmap = base
        norm = Normalize(vmin=tmin, vmax=tmax)

    # Draw circles
    pats = [patches.Circle(xy=tuple(a[0]), radius=a[2]) for a in selected_atoms]
    coll = PatchCollection(pats, edgecolors='black', linewidths=2, cmap=cmap, norm=norm)
    coll.set_array(values)
    if not use_discrete:
        coll.set_clim([tmin, tmax])
    ax.add_collection(coll)

    # EF label (real-value display already provided as ef_display)
    ax.text(0.96, 0.96, f'EF = {ef_display:+.1f} V/Å',
            transform=ax.transAxes, ha='right', va='top', fontsize=26)

    # Ticks & labels
    ax.xaxis.label.set_size(28)
    ax.yaxis.label.set_size(28)
    ax.set_xticks(np.linspace(xlims[0], xlims[1], 5))
    ax.set_yticks(np.linspace(ylims[0], ylims[1], 5))
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.tick_params(axis='x', pad=10)
    ax.tick_params(axis='y', pad=5)

    return coll, cmap, norm
def plot_triptych(ref_geom, layers_plot, mean_vectors, ef_values, plane, vmin, vmax,
                  output_dir: Path = Path("."), palette_base='rainbow', use_discrete=True, reverse=False):
    """Compose a 1x3 panel for the chosen plane with a single colorbar on the rightmost (+0.6 real EF)."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import FixedLocator, FormatStrFormatter, NullLocator

    # Convert to log10 range (input vmin/vmax may come from ln)
    LN10 = np.log(10.0)
    tmin = float(vmin / LN10)
    tmax = float(vmax / LN10)
    # Force the displayed color range to at least [-15, 3] so we can show full tick set
    tmin = min(int(np.floor(tmin)), -15)
    tmax = max(int(np.ceil(tmax)), 3)
    if tmax <= tmin:
        tmax = tmin + 1

    # Order panels by REAL EF: [-0.6, 0.0, +0.6]
    order = sorted(ef_values, reverse=True)  # data: +0.6 | 0.0 | -0.6  # data +0.6 first (real -0.6), then 0.0, then -0.6 (real +0.6)

    fig, axes = plt.subplots(1, 3, figsize=(22, 8), constrained_layout=True)
    final_cmap = None; final_norm = None

    for j, ef_data in enumerate(order):
        ax = axes[j]
        vec = mean_vectors[ef_data]
        vec = np.clip(vec.astype(float), 1e-30, None)
        dq_log10 = np.log10(vec)
        ef_display = -float(ef_data)  # real EF
        coll, cmap, norm = _draw_panel_on_ax(ax, ref_geom, layers_plot[0], dq_log10, plane, ef_display,
                                             tmin=tmin, tmax=tmax,
                                             use_discrete=use_discrete, palette_base=palette_base, reverse=reverse)
        if j == 2:
            final_cmap, final_norm = cmap, norm

    
# Single colorbar to the right of last axis
    if final_cmap is not None:
        cax = axes[-1].inset_axes([1.02, 0.0, 0.06, 1.0])
        sm = ScalarMappable(norm=final_norm, cmap=final_cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label("log$_{10}$(TOF / mol s$^{-1}$ mol$_{cat}^{-1}$)", fontsize=24, labelpad=16)

        # Integer ticks: 3, 0, -3, -6, -9, -12, -15 (always show full set)
        desired = [3, 0, -3, -6, -9, -12, -15]
        cbar.set_ticks(desired)
        cbar.ax.yaxis.set_major_locator(FixedLocator(desired))
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        cbar.ax.yaxis.set_minor_locator(NullLocator())
        # thicker colorbar & ticks
        cbar.outline.set_linewidth(3.0)
        cbar.ax.tick_params(axis='y', which='both', width=2.5, length=10, labelsize=18, direction='out')

    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / f"tof_{plane}_triptych.png"

    fig.savefig(str(out), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved combined figure: {out.resolve()}")
def main():
    # Output directory
    output_dir = Path("mkm_outputs")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"[Info] Output directory: {output_dir.resolve()}\n")
    
    # Geometry
    geom_file = choose_geometry_file()
    ref_geom = read(geom_file.as_posix(), format="vasp")
    natoms = len(ref_geom)
    print(f"Geometry: {geom_file.name} (natoms={natoms})")

    # Read CSV
    if not CSV.exists():
        raise FileNotFoundError(f"Missing {CSV.resolve()}")
    raw = pd.read_csv(CSV, sep=None, engine="python")
    cols_lower = {c.lower(): c for c in raw.columns}
    path_col = cols_lower.get("path")
    tof_col = cols_lower.get("tof")
    if not path_col or not tof_col:
        raise ValueError(f"{CSV} must include 'path' and 'TOF' columns. Found: {raw.columns.tolist()}")

    df = raw.rename(columns={path_col: "path", tof_col: "TOF"})[["path", "TOF"]].copy()
    df["TOF"] = pd.to_numeric(df["TOF"], errors="coerce")
    df = df.dropna(subset=["TOF"]).copy()

    # Parse site & EF
    parsed = df["path"].apply(parse_path)
    df["site"] = parsed.apply(lambda x: x[0])
    df["EF_value"] = parsed.apply(lambda x: x[1])

    # Which EF values exist in file (sorted numeric)
    ef_values = sorted(df["EF_value"].dropna().unique().tolist())
    if not ef_values:
        print("No EF entries found in CSV.")
        return
    print(f"Found EF values: {', '.join(f'{e:+.1f}' for e in ef_values)}")

    # Per-EF aggregation → dict + mean CSV
    mean_files = {}
    for ef in ef_values:
        print(f"\nProcessing EF={ef:+.1f} ...")
        tof_dict = build_per_atom_dict_for_ef(df, natoms, ef)

        dict_path = write_tof_dict_text(tof_dict, ef, output_dir)
        csv_path = write_mean_csv(tof_dict, ef, output_dir)
        mean_files[ef] = csv_path

        print(f"  Saved dict: {dict_path.name}")
        print(f"  Saved mean: {csv_path.name}")

    # Build global color scale across EF
    global_max = 0.0
    mean_vectors = {}
    for ef, csv_name in mean_files.items():
        vec = load_mean_vector(str(csv_name), natoms)
        mean_vectors[ef] = vec
        max_here = float(np.nanmax(vec)) if vec.size else 0.0
        if max_here > global_max:
            global_max = max_here
    vmin, vmax = np.log([1e-8, (global_max * 1.02 if global_max > 0 else 1.0)])
    # vmin, vmax = 0.0, (global_max * 1.02 if global_max > 0 else 1.0)
    print(f"\nColor scale: vmin={vmin:.3g}, vmax={vmax:.3g}")

    # Compose per-plane triptych combining three EF panels; single colorbar on rightmost (+0.6 real EF)
    layers_plot = [list(range(natoms))]
    for plane in ("xy", "xz", "yz"):
        plot_triptych(ref_geom, layers_plot, mean_vectors, ef_values, plane, vmin, vmax,
                      output_dir=output_dir, palette_base='rainbow', use_discrete=True, reverse=False)

    print(f"\n[Info] All outputs saved to: {output_dir.resolve()}")
    print("All done.")

# --------------- If this file is run directly ---------------
if __name__ == "__main__":
    # 这里假设 draw_cell 在你的环境中已经可用；若需要占位测试，可临时定义一个空函数：
    configure_jacs_style()
    try:
        draw_cell  # type: ignore  # noqa
    except NameError:
        def draw_cell(ref_geom, layer_indices, dq_substrate, basename, plane, ef_str, vmin, vmax):
            # 占位：实际环境请用你项目里的 draw_cell 实现
            print(f"[draw_cell] {basename}: plane={plane}, EF={ef_str}, "
                  f"vmin={vmin:.3g}, vmax={vmax:.3g}, atoms={len(layer_indices)}")
    main()

