#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""TOF vs Bader vs CN: scatter or PES surface.
   - x: CN (from POSCAR connectivity), integers 3..12
   - y: Bader charge, limited to [-0.4, 0.4]
   - z/color: TOF
   - REAL EF coloring (scatter) or per-EF surfaces (surface mode)
"""

import argparse, sys, os, re, ast
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
import matplotlib.tri as mtri

# ---------------- Colors by REAL EF ----------------
COL_POS  = "#d62728"   # real +0.6
COL_ZERO = "#ff7f0e"   # real  0.0
COL_NEG  = "#1f77b4"   # real -0.6

def color_for_real_ef(real: float) -> str:
    if abs(real) < 5e-4:
        return COL_ZERO
    return COL_POS if real > 0 else COL_NEG


# save as make_cn_list.py
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs
from ase import Atoms

atoms = read("POSCAR")
# 若你只要 Ru 的 CN：filtered = [a for a in atoms if a.symbol=="Ru"]; atoms = Atoms(filtered)
radii = natural_cutoffs(atoms, mult=0.9)
nl = NeighborList(radii, bothways=True, self_interaction=False); nl.update(atoms)
CN = nl.get_connectivity_matrix(sparse=False).sum(axis=0)
with open("cn_list.txt","w") as f:
    f.write(" ".join(str(int(x)) for x in CN))
print("Wrote cn_list.txt with", len(CN), "values")
 


# ---------------- Style ----------------
def configure_style():
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

# ---------------- CN: compute from POSCAR ----------------
def compute_connections_from_qrobot(poscar_path: Path, metal: str = "Ru", mult: float = 0.9):
    q_robot_path = r"C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain"
    if q_robot_path and q_robot_path not in sys.path:
        sys.path.append(q_robot_path)
    try:
        from cluster import get_connection
        from ase.io import read
    except Exception as e:
        return None, f"[INFO] q-robot import failed ({e})."
    try:
        atoms = read(str(poscar_path))
        connections, *_ = get_connection(atoms, metal=metal, mult=mult)
        return connections, None
    except Exception as e:
        return None, f"[INFO] q-robot get_connection failed ({e})."

def compute_connections_from_ase(poscar_path: Path, mult: float = 0.9):
    from ase.io import read
    from ase.data import covalent_radii
    from ase.neighborlist import NeighborList
    atoms = read(str(poscar_path))
    Z = atoms.get_atomic_numbers()
    cutoffs = [mult * covalent_radii[z] if z < len(covalent_radii) else mult*1.2 for z in Z]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    connections = {i: [] for i in range(len(atoms))}
    for i in range(len(atoms)):
        idx, _ = nl.get_neighbors(i)
        connections[i] = [int(j) for j in idx]
    return connections

def compute_cn_map(poscar_dir: Path, poscar_name: str, metal: str, mult: float,
                   q_robot_path: str | None = None) -> dict[int, list[int]]:
    poscar_path = poscar_dir / poscar_name
    if not poscar_path.is_file():
        raise SystemExit(f"[Error] POSCAR not found at: {poscar_path}")
    if q_robot_path and q_robot_path not in sys.path:
        sys.path.append(q_robot_path)
    conn, info = compute_connections_from_qrobot(poscar_path, metal=metal, mult=mult)
    if conn is not None:
        print("[OK] CN via q-robot/cluster.get_connection")
        return conn
    if info: print(info)
    print("[OK] CN via ASE NeighborList fallback")
    return compute_connections_from_ase(poscar_path, mult=mult)

# ---------------- Data IO ----------------
def load_tog_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if not {"atom_index", "mean_TOF"}.issubset(df.columns):
        raise ValueError(f"Missing columns in {path}: need 'atom_index' and 'mean_TOF'")
    out = pd.DataFrame({
        "Index": df["atom_index"].astype(int) + 1,
        "TOF": pd.to_numeric(df["mean_TOF"], errors="coerce")
    })
    return out.dropna(subset=["TOF"]).reset_index(drop=True)

def load_bader_dat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    def pick(colname):
        tgt = colname.strip().lower()
        for c in df.columns:
            if str(c).strip().lower() == tgt: return c
        for c in df.columns:
            if tgt in str(c).strip().lower(): return c
        raise KeyError(f"Column '{colname}' not found in {path}. Available: {list(df.columns)}")
    col_idx = pick("index"); col_charge = pick("charge")
    out = pd.DataFrame({
        "Index": pd.to_numeric(df[col_idx], errors="coerce").astype("Int64"),
        "Charge": pd.to_numeric(df[col_charge], errors="coerce")
    }).dropna(subset=["Index", "Charge"])
    out["Index"] = out["Index"].astype(int)
    return out.reset_index(drop=True)

def parse_data_ef_from_name(p: Path) -> float:
    m = re.search(r"([+-]?\d*\.?\d+)", p.stem)
    if not m: raise ValueError(f"Cannot parse EF from file name: {p.name}")
    return float(m.group(1))

# ---------------- Build per-EF dataframe ----------------
def build_df_for_real_ef(real_ef: float, tog_path: Path, bader_path: Path, cn_map: dict[int, list[int]]) -> pd.DataFrame:
    df_tog = load_tog_csv(tog_path)
    df_bdr = load_bader_dat(bader_path)
    df = pd.merge(df_tog, df_bdr, on="Index", how="inner").sort_values("Index").reset_index(drop=True)
    if df.empty:
        return df
    df["CN"] = df["Index"].astype(int).map(lambda idx1: len(cn_map.get(idx1 - 1, [])))
    df["real_EF"] = real_ef
    # 过滤至指定范围
    df = df[df["CN"].between(3, 12) & df["Charge"].between(-0.4, 0.4)].copy()
    return df[["Index", "CN", "Charge", "TOF", "real_EF"]]

# ---------------- Plot: scatter (existing) ----------------
def plot_scatter_3d(df_all: pd.DataFrame, outpng: Path, per_ef: bool = False):
    configure_style()
    fig = plt.figure(figsize=(7.4, 6.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev=18, azim=35); ax.dist = 9
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        try:
            axis.pane.set_alpha(0.0); axis.pane.set_edgecolor("none")
        except Exception:
            pass

    # 按 TOF 排序 + 轻微抖动，减少遮挡
    for real_ef, dsub in df_all.groupby("real_EF"):
        dsub = dsub.sort_values("TOF")
        cn = dsub["CN"].to_numpy(dtype=float) + np.random.uniform(-0.08, 0.08, size=len(dsub))
        ax.scatter(cn, dsub["Charge"], dsub["TOF"],
                   s=26, alpha=0.9, depthshade=False,
                   color=color_for_real_ef(real_ef),
                   edgecolors=(0,0,0,0.25), linewidths=0.4,
                   label=("EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å"))

    # 轴设置
    ax.set_xlim(3, 12); ax.set_xticks(np.arange(3, 13, 1))
    ax.set_ylim(-0.4, 0.4); ax.set_yticks(np.linspace(-0.4, 0.4, 5))
    ax.set_zlim(df_all["TOF"].min(), df_all["TOF"].max())
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.tick_params(pad=2, length=3, width=0.6)

    ax.set_xlabel("Coordination number (CN)", labelpad=8)
    ax.set_ylabel("Bader charge", labelpad=8)
    ax.set_zlabel("TOF", labelpad=12)
    ax.legend(frameon=False, loc=(0.02, 0.80), fontsize=12)
    fig.subplots_adjust(left=0.10, right=0.88, bottom=0.10, top=0.95)
    fig.savefig(outpng, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"Saved scatter -> {outpng}")

    if per_ef:
        for real_ef, dsub in df_all.groupby("real_EF"):
            fig = plt.figure(figsize=(7.0, 6.0))
            ax = fig.add_subplot(111, projection="3d")
            ax.view_init(elev=18, azim=35); ax.dist = 9; ax.grid(False)
            ax.scatter(dsub["CN"], dsub["Charge"], dsub["TOF"],
                       s=28, alpha=0.95, depthshade=False,
                       color=color_for_real_ef(real_ef),
                       edgecolors=(0,0,0,0.3), linewidths=0.4)
            ax.set_title("EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å")
            ax.set_xlim(3, 12); ax.set_xticks(np.arange(3, 13, 1))
            ax.set_ylim(-0.4, 0.4); ax.set_yticks(np.linspace(-0.4, 0.4, 5))
            ax.set_zlim(dsub["TOF"].min(), dsub["TOF"].max())
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_xlabel("Coordination number (CN)", labelpad=8)
            ax.set_ylabel("Bader charge", labelpad=8)
            ax.set_zlabel("TOF", labelpad=12)
            fig.subplots_adjust(left=0.10, right=0.88, bottom=0.10, top=0.95)
            outname = Path(f"tof_bader_cn_scatter_EF_{'+' if real_ef>0 else '-' if real_ef<0 else ''}{abs(real_ef):.1f}.png").as_posix().replace('+-','-').replace('--','-')
            fig.savefig(outname, dpi=300, bbox_inches="tight"); plt.close(fig)
            print(f"Saved scatter per-EF -> {outname}")

# ---------------- Plot: PES surface (NEW) ----------------
def plot_pes_surfaces(df_all: pd.DataFrame, out3d_prefix: str = "tof_bader_cn_surface_3d_EF_",
                      out2d_prefix: str = "tof_bader_cn_contour_EF_",
                      levels: int = 15):
    """Per-EF surfaces: 3D trisurf + 2D tricontourf (CN,Bader) with TOF as height/color."""
    configure_style()
    for real_ef, dsub in df_all.groupby("real_EF"):
        if len(dsub) < 3:
            print(f"[WARN] Not enough points to make a surface for EF {real_ef:+.1f}. Skipped.")
            continue

        X = dsub["CN"].to_numpy(dtype=float)
        Y = dsub["Charge"].to_numpy(dtype=float)
        Z = dsub["TOF"].to_numpy(dtype=float)
        tri = mtri.Triangulation(X, Y)

        # --- 3D surface ---
        fig = plt.figure(figsize=(7.2, 5.8))
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(elev=25, azim=40); ax.dist = 9; ax.grid(False)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            try: axis.pane.set_alpha(0.0); axis.pane.set_edgecolor("none")
            except Exception: pass

        surf = ax.plot_trisurf(tri, Z, cmap="viridis", linewidth=0.15, antialiased=True, alpha=0.98)
        # 同时叠加原始散点（淡）
        ax.scatter(X, Y, Z, s=16, color="k", alpha=0.35, depthshade=False, linewidths=0)

        ax.set_xlim(3, 12); ax.set_xticks(np.arange(3, 13, 1))
        ax.set_ylim(-0.4, 0.4); ax.set_yticks(np.linspace(-0.4, 0.4, 5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_xlabel("Coordination number (CN)", labelpad=8)
        ax.set_ylabel("Bader charge", labelpad=8)
        ax.set_zlabel("TOF", labelpad=12)

        cbar = fig.colorbar(surf, shrink=0.75, pad=0.05)
        cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.set_title("EF = 0.0 V/Å" if abs(real_ef) < 5e-4 else f"EF = {real_ef:+.1f} V/Å")
        fig.subplots_adjust(left=0.10, right=0.90, bottom=0.10, top=0.90)
        out3d = f"{out3d_prefix}{'+' if real_ef>0 else '-' if real_ef<0 else ''}{abs(real_ef):.1f}.png".replace('+-','-').replace('--','-')
        fig.savefig(out3d, dpi=300, bbox_inches="tight"); plt.close(fig)
        print(f"Saved PES 3D -> {out3d}")

        # --- 2D contour (tricontourf) ---
        fig2, ax2 = plt.subplots(figsize=(6.6, 5.2))
        cf = ax2.tricontourf(tri, Z, levels=levels, cmap="viridis")
        ax2.tricontour(tri, Z, colors="k", linewidths=0.3, levels=levels)
        ax2.scatter(X, Y, s=10, c="k", alpha=0.35)
        ax2.set_xlim(3, 12); ax2.set_xticks(np.arange(3, 13, 1))
        ax2.set_ylim(-0.4, 0.4); ax2.set_yticks(np.linspace(-0.4, 0.4, 5))
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax2.set_xlabel("Coordination number (CN)")
        ax2.set_ylabel("Bader charge")
        cbar2 = fig2.colorbar(cf, shrink=0.85, pad=0.02)
        cbar2.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        fig2.tight_layout()
        out2d = f"{out2d_prefix}{'+' if real_ef>0 else '-' if real_ef<0 else ''}{abs(real_ef):.1f}.png".replace('+-','-').replace('--','-')
        fig2.savefig(out2d, dpi=300, bbox_inches="tight"); plt.close(fig2)
        print(f"Saved PES 2D -> {out2d}")

# ---------------- Main ----------------
def main():
    parser = argparse.ArgumentParser(description="3D scatter or PES surface of TOF(CN, Bader), colored by REAL EF.")
    # Files
    parser.add_argument("--tog", nargs="*", default=["mkm_outputs/tog_-0.6_sim.csv", "mkm_outputs/tog_0.0_sim.csv", "mkm_outputs/tog_+0.6_sim.csv"])
    parser.add_argument("--bader", nargs="*", default=["bader_-0.6.dat", "bader_0.0.dat", "bader_0.6.dat"])
    # CN compute opts
    parser.add_argument("--slab-dir", default="./slab")
    parser.add_argument("--poscar-name", default="POSCAR")
    parser.add_argument("--metal", default="Ru")
    parser.add_argument("--mult", type=float, default=0.9)
    parser.add_argument("--q-robot-path", default=None)
    # Plot mode
    parser.add_argument("--mode", choices=["scatter","surface"], default="surface",
                        help="scatter: 原先的3D散点; surface: 每个EF生成PES三角面+等高图")
    parser.add_argument("--scatter-out", default="tof_bader_cn_3d.png", help="scatter模式合并输出文件名")
    parser.add_argument("--per-ef", action="store_true", help="scatter模式下，同时输出每个EF单独一张")
    args = parser.parse_args()

    # --- CN map ---
    cn_map = compute_cn_map(Path(args.slab_dir), args.poscar_name, metal=args.metal, mult=args.mult,
                            q_robot_path=args.q_robot_path)

    # --- map DATA EF -> path ---
    tog_by_data_ef, bader_by_data_ef = {}, {}
    for f in args.tog:
        p = Path(f); 
        if not p.is_file(): raise SystemExit(f"[Error] Missing TOF file: {p}")
        tog_by_data_ef[parse_data_ef_from_name(p)] = p
    for f in args.bader:
        p = Path(f);
        if not p.is_file(): raise SystemExit(f"[Error] Missing Bader file: {p}")
        bader_by_data_ef[parse_data_ef_from_name(p)] = p

    # --- build merged frames for REAL EFs ---
    frames = []
    for real in [-0.6, 0.0, 0.6]:
        data_value = -real
        key_tog = min(tog_by_data_ef.keys(), key=lambda x: abs(x - data_value))
        key_bdr = min(bader_by_data_ef.keys(), key=lambda x: abs(x - data_value))
        df_one = build_df_for_real_ef(real, tog_by_data_ef[key_tog], bader_by_data_ef[key_bdr], cn_map)
        if df_one is None or df_one.empty:
            print(f"[WARN] No merged rows for EF {real:+.1f}.")
            continue
        frames.append(df_one)
    if not frames:
        raise SystemExit("[Error] No data to plot after merging.")

    df_all = pd.concat(frames, ignore_index=True)

    # --- Plot ---
    if args.mode == "scatter":
        plot_scatter_3d(df_all, outpng=Path(args.scatter_out), per_ef=args.per_ef)
    else:
        # surface: per-EF surfaces (3D + 2D)
        plot_pes_surfaces(df_all)

if __name__ == "__main__":
    main()

