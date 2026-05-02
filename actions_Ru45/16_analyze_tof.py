#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Summarize TOF distribution by electric field (EF) from drc_summary_output.csv
and optionally plot a bar chart of log10(TOF) for EF = -0.6, 0.0, 0.6.

Input format (CSV with header):
path,TOF
0_11_34/EF_-0.6,0.085360003
0_11_34/EF_0.0,0.016856941
0_11_38/EF_-0.6,0.02354535
0_11_38/EF_0.0,0.000871508
0_11_38/EF_0.6,5.89E-08

Usage:
  # 控制台输出各电场的 TOF 最小/最大值
  python drc_tof_summary_with_plot.py -i drc_summary_output.csv

  # 额外生成 log10(TOF) 柱状图（默认名 ef_log10tof_barchart.png）
  python drc_tof_summary_with_plot.py -i drc_summary_output.csv --plot

  # 指定图片文件名；可选使用配色（与前面一致）
  python drc_tof_summary_with_plot.py -i drc_summary_output.csv --plot -p my.png --use-colors
"""
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import re
import matplotlib.pyplot as plt

# 打印顺序（控制台）：0, 0.6, -0.6
ORDER_PRINT = ["0.0", "0.6", "-0.6"]
# 绘图顺序：-0.6, 0.0, 0.6
ORDER_PLOT = ["-0.6", "0.0", "0.6"]

# 可选配色（与前面的颜色保持一致）；仅在 --use-colors 时应用
COLOR_MAP = {
    "-0.6": "#1f77b4",
    "0.0":  "#2ca02c",
    "0.6":  "#d62728",
}

def parse_ef(path_str: str) -> str:
    """Extract EF string from '.../EF_-0.6' -> '-0.6', normalize to one decimal."""
    if not isinstance(path_str, str):
        return np.nan
    m = re.search(r'EF_([+-]?[0-9]*\.?[0-9]+)', path_str)
    if not m:
        return np.nan
    try:
        return f"{float(m.group(1)):.1f}"
    except Exception:
        return np.nan

def load_table(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df['TOF'] = pd.to_numeric(df['TOF'], errors='coerce')
    df['EF'] = df['path'].map(parse_ef)
    df = df.dropna(subset=['EF', 'TOF']).copy()
    return df

def summarize_by_ef(df: pd.DataFrame) -> pd.DataFrame:
    idx_min = df.groupby('EF')['TOF'].idxmin()
    idx_max = df.groupby('EF')['TOF'].idxmax()

    gcount = df.groupby('EF')['TOF'].size().rename('count')
    dmin = df.loc[idx_min, ['EF', 'path', 'TOF']].rename(columns={'path': 'path_min', 'TOF': 'TOF_min'}).set_index('EF')
    dmax = df.loc[idx_max, ['EF', 'path', 'TOF']].rename(columns={'path': 'path_max', 'TOF': 'TOF_max'}).set_index('EF')

    out = gcount.to_frame().join(dmin, how='left').join(dmax, how='left').reset_index()
    out['EF_float'] = out['EF'].astype(float)
    out = out.sort_values('EF_float').drop(columns=['EF_float'])
    return out

def _log10_stats_by_ef(df: pd.DataFrame) -> pd.DataFrame:
    """Return per-EF min/median/max of log10(TOF); drop non-positive TOF rows."""
    d = df[df['TOF'] > 0].copy()
    if d.empty:
        return pd.DataFrame(columns=['EF','min','median','max'])
    d['log10_TOF'] = np.log10(d['TOF'])
    g = d.groupby('EF')['log10_TOF']
    out = pd.DataFrame({
        'min': g.min(),
        'median': g.median(),
        'max': g.max()
    }).reset_index()
    return out

def plot_log10_tof_barchart(df: pd.DataFrame,
                            outfile: str = 'ef_log10tof_barchart.png',
                            order = ORDER_PLOT,
                            color_map = None,
                            title: str = 'log10(TOF) by EF') -> str:
    stats = _log10_stats_by_ef(df)
    if stats.empty:
        print('[Warn] No positive TOF values; skip plotting.')
        return ''
    # keep only requested order that is present
    ef_avail = [ef for ef in order if ef in set(stats['EF'])]
    stats = stats.set_index('EF').loc[ef_avail]

    y = stats['median'].values
    err_low = (stats['median'] - stats['min']).values
    err_high = (stats['max'] - stats['median']).values

    x = np.arange(len(ef_avail))
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(x, y)  # 默认配色；除非用户开启 --use-colors
    ax.errorbar(x, y, yerr=[err_low, err_high], fmt='none', capsize=4, linewidth=1)

    # 应用颜色（可选）
    if color_map is not None:
        for rect, ef in zip(bars, ef_avail):
            if ef in color_map:
                rect.set_color(color_map[ef])

    ax.set_xticks(x)
    ax.set_xticklabels(ef_avail)
    ax.set_xlabel('EF (V/Å)')
    ax.set_ylabel('median log10(TOF)')
    ax.set_title(title)
    ax.grid(False)  # 无网格
    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return outfile

def main():
    ap = argparse.ArgumentParser(description="Summarize TOF min/max by EF and plot log10(TOF) barchart.")
    ap.add_argument('-i', '--input', default='mkm_outputs/mkm_out.csv', help='Input CSV (default: mkm_outputs/mkm_out.csv)')
    ap.add_argument('-o', '--output', default=None, help='Optional output CSV path for the min/max summary')
    ap.add_argument('--plot', action='store_true', help='Generate log10(TOF) barchart for EF = -0.6, 0.0, 0.6')
    ap.add_argument('-p', '--plot-outfile', default='ef_log10tof_barchart.png', help='Output image path for the barchart')
    ap.add_argument('--use-colors', action='store_true', help='Use predefined colors for EF values')
    args = ap.parse_args()

    csv_path = Path(args.input)
    if not csv_path.exists():
        raise SystemExit(f"[Error] Input file not found: {csv_path}")

    df = load_table(csv_path)
    if df.empty:
        raise SystemExit("[Error] No valid rows after parsing EF and TOF.")

    # Print min/max summary
    summary = summarize_by_ef(df)

    wanted = [ef for ef in ORDER_PRINT if ef in set(summary['EF'])]
    others = [ef for ef in summary['EF'].tolist() if ef not in wanted]
    def _num(x):
        try: return float(x)
        except: return np.inf
    others = sorted(set(others), key=_num)

    for ef in wanted + others:
        row = summary.loc[summary['EF'] == ef]
        if row.empty:
            continue
        r = row.iloc[0]
        print(f"EF = {ef}")
        print(f"  min TOF = {r['TOF_min']:.6g}  (path: {r['path_min']})")
        print(f"  max TOF = {r['TOF_max']:.6g}  (path: {r['path_max']})")
        print("")

    # Optional CSV output
    if args.output:
        out_path = Path(args.output)
        summary.to_csv(out_path, index=False)
        print(f"Saved summary CSV -> {out_path}")

    # Optional plot
    if args.plot:
        cmap = COLOR_MAP if args.use_colors else None
        img = plot_log10_tof_barchart(df, outfile=args.plot_outfile, color_map=cmap)
        if img:
            print(f"Saved bar chart -> {img}")

if __name__ == '__main__':
    main()

