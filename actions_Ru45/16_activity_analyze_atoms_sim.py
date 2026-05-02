#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16_activity_analyze_atoms_sim.py

新增：在“每个电场（数据 EF）下”对相同 TOF 数值进行去重：
- 如果在某个电场列里，多个原子具有完全相同的 mean_TOF，只保留首次出现（按 atom_index 升序的第一行），
  其余重复者在该电场列置为 NaN（从而该电场下不再绘制该原子的柱子）。
- 不会删除该原子在“其他电场”的数据；但如果它在其他电场也与别人重复，也会在对应电场列同样置为 NaN。
- 用法与参数保持不变。

功能（与原脚本一致，除上面新增外）：
- 针对“在 atoms_above.txt 里的原子”和“不在其中的原子”分别绘制一张**分组柱状图**，展示 **每个原子的 TOF**（不求和）。
- 输入（默认在当前目录）：
  - tog_-0.6_sim.csv, tog_0.0_sim.csv, tog_+0.6_sim.csv  （列：atom_index, mean_TOF）
  - atoms_above.txt  （整数，可含空格/逗号/换行，去重后排序）
- 颜色与标注遵循固定的“真实 EF”配色（数据 EF 取反 -> 真实 EF）：
  真实 -0.6 → #1f77b4（蓝）  [数据 EF +0.6]
  真实  0.0 → #ff7f0e（橙）   [数据 EF  0.0]
  真实 +0.6 → #d62728（红）   [数据 EF -0.6]

用法（保持不变，仅文件名变为本脚本）：
  python3 16_activity_analyze_atoms_sim.py
  # 或自定义：
  python3 16_activity_analyze_atoms_sim.py -a atoms_above.txt       -i tog_-0.6_sim.csv tog_0.0_sim.csv tog_+0.6_sim.csv       --out-above tof_atoms_above_grouped.png       --out-below tof_atoms_below_grouped.png
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ---------- EF color scheme (by REAL EF) ----------
COL_POS  = "#d62728"  # real +0.6  (maps from data EF_-0.6)
COL_ZERO = "#ff7f0e"  # real  0.0
COL_NEG  = "#1f77b4"  # real -0.6  (maps from data EF_+0.6)

def color_for_data_ef(ef_str: str) -> str:
    v = float(ef_str)
    real = -v
    if abs(real) < 5e-4:
        return COL_ZERO
    return COL_POS if real > 0 else COL_NEG

def label_for_data_ef(ef_str: str) -> str:
    v = float(ef_str)
    real = -v
    return "EF = 0.0 V/Å" if abs(real) < 5e-4 else f"EF = {real:+.1f} V/Å"

# ---------- IO helpers ----------
def parse_ef_from_filename(p: Path) -> str:
    """tog_+0.6_sim.csv -> '0.6'; tog_-0.6_sim.csv -> '-0.6'"""
    m = re.search(r'tog_([+-]?\d*\.?\d+)_sim\.csv$', p.name)
    if not m:
        raise ValueError(f"Bad filename (need 'tog_±X.X_sim.csv'): {p}")
    return f"{float(m.group(1)):.1f}"

def load_one_csv(path: Path) -> pd.Series:
    """Return Series of mean_TOF indexed by atom_index (int), name=EF string (data EF)."""
    df = pd.read_csv(path)
    if not {'atom_index', 'mean_TOF'}.issubset(df.columns):
        raise ValueError(f"Missing columns in {path}: need 'atom_index' and 'mean_TOF'")
    s = pd.to_numeric(df['mean_TOF'], errors='coerce')
    ser = pd.Series(s.values, index=df['atom_index'].astype(int))
    ser.sort_index(inplace=True)
    ser.name = parse_ef_from_filename(path)  # '-0.6'/'0.0'/'0.6' (data EF)
    return ser

def load_atoms_list(path: Path) -> list[int]:
    """Parse integers from atoms_above.txt; accept spaces, commas, newlines. Return sorted unique list."""
    if not path.is_file():
        raise FileNotFoundError(f"atoms file not found: {path}")
    text = path.read_text(encoding="utf-8", errors="ignore")
    nums = re.findall(r'-?\d+', text)
    uniq = sorted(set(int(n) for n in nums))
    return uniq

# ---------- de-duplicate per EF (data EF) ----------
def deduplicate_per_ef_exact(df: pd.DataFrame) -> pd.DataFrame:
    """
    对每个数据 EF 列，按数值精确匹配去重：保留首次出现的原子，其余重复者置 NaN。
    - 按 index（atom_index）已排序的顺序确定“首次出现”。
    - 仅对该列（该电场）去重，不影响其它电场列。
    """
    df2 = df.copy()
    # 确保按 atom_index 升序
    df2 = df2.sort_index()
    for col in df2.columns:
        ser = df2[col]
        # pandas 的 duplicated 基于精确相等（NaN 不会计为重复）
        dup_mask = ser.duplicated(keep='first')
        # 将重复者置为 NaN（仅该列）
        df2.loc[dup_mask, col] = np.nan
    return df2

# ---------- plotting ----------
def _configure_style():
    plt.rcParams.update({
        "font.family": "Times New Roman",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "legend.fontsize": 12,
        "axes.labelsize": 18,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
    })

def _format_one_decimal(y, _):
    return f"{y:.1f}"

def _plot_grouped(df: pd.DataFrame, ef_cols: list[str], outpng: Path, title: str):
    """
    df: index = atom_index, columns = data EF strings (subset of ef_cols) with mean_TOF values (部分可能为 NaN).
    ef_cols: desired order of data EF columns ["-0.6","0.0","0.6"]
    """
    _configure_style()
    atoms = df.index.tolist()
    n_atoms = len(atoms)
    if n_atoms == 0:
        print(f"[WARN] No atoms to plot for {outpng.name}, skip.")
        return

    # dynamic figure width to avoid crowding
    fw = max(6.5, 0.45 * n_atoms + 2.0)
    fig, ax = plt.subplots(figsize=(fw, 4.6))

    n_ef = len(ef_cols)
    group_span = 0.90
    bar_width = group_span / n_ef
    start = -group_span/2 + bar_width/2
    offsets = [start + i*bar_width for i in range(n_ef)]

    x = np.arange(n_atoms)

    for i, ef in enumerate(ef_cols):
        if ef not in df.columns:
            continue
        # 将 NaN 高度的条形图忽略（matplotlib 对 NaN 会不绘制）
        ax.bar(x + offsets[i], df[ef].values, width=bar_width,
               color=color_for_data_ef(ef), edgecolor="black",
               label=label_for_data_ef(ef))

    ax.set_xticks(x)
    ax.set_xticklabels([str(a) for a in atoms], rotation=45, ha='right')
    ax.set_ylabel("TOF")
    ax.yaxis.set_major_formatter(FuncFormatter(_format_one_decimal))

    # legend once
    ax.legend(frameon=False, loc="best", ncol=1)

    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved -> {outpng}")

def main():
    ap = argparse.ArgumentParser(description="Per-atom grouped TOF bars for atoms_above and atoms_below (no summation), with per-EF duplicate removal.")
    ap.add_argument("-i", "--inputs", nargs='*',
                    default=["mkm_outputs/tog_-0.6_sim.csv", "mkm_outputs/tog_0.0_sim.csv", "mkm_outputs/tog_+0.6_sim.csv"],
                    help="Input CSVs (default: mkm_outputs/tog_*.csv files)")
    ap.add_argument("-a", "--atoms", default="atoms_above.txt", help="Atoms list file (integers).")
    ap.add_argument("--out-above", default="tof_atoms_above_grouped.png", help="Output image for atoms above.")
    ap.add_argument("--out-below", default="tof_atoms_below_grouped.png", help="Output image for atoms not in above.")
    args = ap.parse_args()

    # Load atoms list
    atoms_file = Path(args.atoms)
    above_list = load_atoms_list(atoms_file)
    above_set = set(above_list)

    # Load CSVs
    series = []
    efs = []
    for p in args.inputs:
        path = Path(p)
        if not path.exists():
            raise SystemExit(f"[Error] not found: {path}")
        ser = load_one_csv(path)
        series.append(ser)
        efs.append(ser.name)

    # Combine to DataFrame aligned by atom_index
    df = pd.concat(series, axis=1)
    df.columns = efs  # data EF strings
    # Ensure desired order exists
    des_efs = ["0.6", "0.0", "-0.6"]
    df = df.reindex(columns=[ef for ef in des_efs if ef in df.columns])

    # === 新增：按电场逐列精确去重（仅置 NaN，不删除整个原子行） ===
    before_non_nan = df.notna().sum()
    df = deduplicate_per_ef_exact(df)
    after_non_nan = df.notna().sum()
    # 打印每个电场列的去重数量
    for col in df.columns:
        removed = int(before_non_nan[col] - after_non_nan[col])
        if removed > 0:
            print(f"[Info] Removed {removed} duplicate value(s) in EF {col} (kept first occurrences).")

    # Split above vs below
    idx_all = df.index.astype(int).tolist()
    above_idx = [i for i in idx_all if i in above_set]
    below_idx = [i for i in idx_all if i not in above_set]

    df_above = df.loc[above_idx] if len(above_idx)>0 else df.iloc[0:0]
    df_below = df.loc[below_idx] if len(below_idx)>0 else df.iloc[0:0]

    # Plot both
    _plot_grouped(df_above, [ef for ef in des_efs if ef in df.columns], Path(args.out_above), "Atoms in atoms_above.txt")
    _plot_grouped(df_below, [ef for ef in des_efs if ef in df.columns], Path(args.out_below), "Atoms not in atoms_above.txt")

if __name__ == "__main__":
    main()
