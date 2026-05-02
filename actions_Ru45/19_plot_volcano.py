#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Volcano plots with:
- REAL EF assignment by atoms_above hits (>=2 -> invert sign).
- Optional two-line fits (global / per-ef / peak-anchored).
- Deviation marking & CSV export (threshold on |log10(TOF)-fit|).
- Peak-based classification CSV (left / peak / right) + CN columns via --cn-file.
- E_N adsorption energies from --en-file (E_N.csv) into exports as 'E_N_ads'.
- NEW: Bader charges per-atom from --bader-dir (files like bader_0.6.dat, bader_-0.6.dat, bader_0.0.dat).
  Added columns: bader_1, bader_2, bader_3, bader_average (matched by *data* EF from path).

Input rows are filtered to DRC_1_label == 'r_0006' when that column exists.
"""

import argparse, re, os
from pathlib import Path
from typing import List, Set, Tuple, Optional, Dict, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator

# ---------- Visual defaults ----------
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
DESCRIPTORS = ["E_N", "E_H", "E_NH", "E_NH2", "E_NH3"]  # Fallback defaults

# ---------- Descriptor discovery ----------
def discover_descriptors_from_csv(csv_path: Path, energy_pattern: str = r'^E_[A-Za-z0-9_]+$') -> List[str]:
    """
    Discover energy descriptor columns from CSV file.
    Matches columns starting with 'E_' and excludes known non-descriptor columns like 'EF'.
    Returns sorted list of matching columns found in CSV.
    """
    try:
        df = pd.read_csv(csv_path)
        pattern = re.compile(energy_pattern)
        descriptors = sorted([col for col in df.columns if pattern.match(col) and col not in ['EF', 'E_slab']])
        if descriptors:
            print(f"[Info] Discovered {len(descriptors)} descriptors from {csv_path}: {descriptors}")
            return descriptors
        else:
            print(f"[Warn] No descriptors matching pattern found in {csv_path}. Using defaults.")
            return DESCRIPTORS
    except Exception as e:
        print(f"[Warn] Could not load descriptors from {csv_path}: {e}. Using defaults.")
        return DESCRIPTORS

# ---------- Colors by REAL EF ----------
COL_POS  = "#d62728"  # +0.6
COL_ZERO = "#ff7f0e"  #  0.0
COL_NEG  = "#1f77b4"  # -0.6

def legend_label_from_real_ef(real_ef: float) -> str:
    r = round(float(real_ef), 1)
    return "EF = 0.0 V/Å" if r == 0 else f"EF = {r:+.1f} V/Å"

def color_from_real_ef(real_ef: float) -> str:
    real_ef = float(real_ef)
    return COL_POS if real_ef > 0 else (COL_NEG if real_ef < 0 else COL_ZERO)

# ---------- EF tag/site parsing ----------
_EF_TAG_RE = re.compile(r'(EF_[+-]?\d+(?:\.\d+)?)')

def extract_ef_tag(path_str: str) -> str:
    if not isinstance(path_str, str):
        return ""
    m = _EF_TAG_RE.search(path_str)
    return m.group(1) if m else ""

def parse_ef_value_from_tag(tag: str) -> float:
    return float(tag.split('_', 1)[1])

def parse_site_indices(path_str: str) -> List[int]:
    if not isinstance(path_str, str):
        return []
    first = path_str.split('/', 1)[0]
    out = []
    for tok in first.split('_'):
        tok = tok.strip()
        if tok and tok.lstrip('-').isdigit():
            out.append(int(tok))
    return out

def path_sites_str(path_str: str) -> str:
    if not isinstance(path_str, str):
        return ""
    return path_str.split('/', 1)[0]

# ---------- atoms_above & REAL EF ----------
def load_atoms_above(file_path: Path, *, require_exists: bool) -> Set[int]:
    if not file_path.exists():
        if require_exists:
            raise SystemExit(f"[Error] atoms list not found: {file_path}")
        else:
            print(f"[Warn] '{file_path}' not found. Proceeding with empty atoms_above set.")
            return set()
    atoms = set()
    with open(file_path, 'r') as f:
        for line in f:
            t = line.strip()
            if t and t.lstrip('-').isdigit():
                atoms.add(int(t))
    print(f"Loaded {len(atoms)} atom indices from {file_path}")
    return atoms

def count_atoms_hits(path_str: str, atoms_set: Set[int]) -> int:
    sites = parse_site_indices(path_str)
    return sum(1 for s in sites if s in atoms_set)

def compute_real_ef_from_row(path_str: str, atoms_set: Set[int]) -> float:
    tag = extract_ef_tag(path_str)
    if tag == "":
        return np.nan
    data_ef = parse_ef_value_from_tag(tag)
    hits = count_atoms_hits(path_str, atoms_set)
    return -data_ef if hits >= 2 else data_ef

# ---------- CN helpers ----------
def load_cn_list(cn_file: Optional[str]) -> Optional[List[int]]:
    if not cn_file:
        return None
    p = Path(cn_file)
    if not p.exists():
        print(f"[Warn] CN file not found: {cn_file}. CN columns will be empty.")
        return None
    txt = p.read_text().strip().replace(',', ' ')
    vals = [int(x) for x in txt.split()]
    print(f"Loaded {len(vals)} CN values from {cn_file}")
    return vals

def cn_for_site(cn_list: Optional[List[int]], site_indices: List[int]) -> Tuple[Optional[int],Optional[int],Optional[int],Optional[float]]:
    if not cn_list:
        return (None, None, None, None)
    cn_vals = []
    for idx in site_indices[:3]:
        if 0 <= idx < len(cn_list):
            cn_vals.append(int(cn_list[idx]))
        else:
            cn_vals.append(None)
    while len(cn_vals) < 3:
        cn_vals.append(None)
    avg_vals = [v for v in cn_vals if v is not None]
    avg = float(sum(avg_vals))/len(avg_vals) if avg_vals else None
    return (cn_vals[0], cn_vals[1], cn_vals[2], avg)

# ---------- E_N adsorption energy loader ----------
def load_en_lookup(en_file: Optional[str], atoms_set: Set[int]) -> Dict[str, Any]:
    maps = {'map_site_ef': {}, 'map_site_only': {}}
    if not en_file:
        return maps
    p = Path(en_file)
    if not p.exists():
        print(f"[Warn] E_N file not found: {en_file}. Will skip adding E_N_ads.")
        return maps
    df = pd.read_csv(p)
    en_col = None
    for cand in ['E_N', 'E_N_ads', 'EN', 'E_N_eV']:
        if cand in df.columns:
            en_col = cand; break
    if en_col is None:
        raise SystemExit(f"[Error] Could not find adsorption-energy column in {en_file} (looked for E_N/E_N_ads/EN/E_N_eV).")
    has_path = 'path' in df.columns
    has_site = 'site' in df.columns
    has_EF   = 'EF' in df.columns  # treated as EF_real if present
    for _, row in df.iterrows():
        try:
            en_val = float(row[en_col])
        except Exception:
            continue
        site_str = None; ef_real = None
        if has_path and isinstance(row['path'], str) and row['path']:
            site_str = path_sites_str(row['path'])
            ef_real = compute_real_ef_from_row(row['path'], atoms_set)
        if site_str is None and has_site and isinstance(row['site'], str) and row['site']:
            site_str = row['site'].split('/', 1)[0].strip()
        if site_str is None: continue
        if ef_real is None and has_EF:
            try: ef_real = float(row['EF'])
            except Exception: ef_real = None
        if ef_real is not None and np.isfinite(ef_real):
            maps['map_site_ef'][(site_str, float(np.round(ef_real, 1)))] = en_val
        maps['map_site_only'][site_str] = en_val
    print(f"Loaded E_N lookup from {en_file}: "
          f"{len(maps['map_site_ef'])} site+EF entries, {len(maps['map_site_only'])} site-only entries.")
    return maps

def lookup_en_for_row(row: pd.Series, en_maps: Dict[str, Any]) -> Optional[float]:
    if not en_maps:
        return None
    site = path_sites_str(row['path'])
    try:
        ef = float(np.round(float(row['EF_real']), 1))
    except Exception:
        ef = None
    m_se = en_maps.get('map_site_ef', {})
    m_s  = en_maps.get('map_site_only', {})
    if ef is not None and (site, ef) in m_se:
        return float(m_se[(site, ef)])
    if site in m_s:
        return float(m_s[site])
    return None

# ---------- Bader charge loader (by *data* EF from path) ----------
_BADER_FILE_RE = re.compile(r'bader_([+-]?\d+(?:\.\d+)?)\.dat$')

def load_bader_maps(bader_dir: Optional[str]) -> Dict[float, Dict[int, float]]:
    """
    Scan directory for files 'bader_*.dat'. For each EF value in filename,
    read a table with at least columns: Index, Charge (Index is 1-based).
    Returns: { data_EF(float, rounded 1dp) : { index_1based(int) : charge(float) } }
    """
    maps: Dict[float, Dict[int, float]] = {}
    if not bader_dir:
        return maps
    d = Path(bader_dir)
    if not d.exists():
        print(f"[Warn] Bader dir not found: {bader_dir}. Bader columns will be empty.")
        return maps
    for name in os.listdir(d):
        m = _BADER_FILE_RE.match(name)
        if not m: continue
        ef_val = round(float(m.group(1)), 1)  # data EF
        fpath = d / name
        try:
            df = pd.read_csv(fpath)  # expects 'Index,Element,Charge'
        except Exception:
            # try whitespace
            df = pd.read_csv(fpath, delim_whitespace=True)
        if 'Index' not in df.columns or 'Charge' not in df.columns:
            print(f"[Warn] {name}: missing 'Index' or 'Charge' columns; skip.")
            continue
        charges = {}
        for _, row in df.iterrows():
            try:
                idx1 = int(row['Index'])     # 1-based
                ch   = float(row['Charge'])
                charges[idx1] = ch
            except Exception:
                continue
        maps[ef_val] = charges
        print(f"Loaded Bader: {name}  (EF={ef_val:+.1f}, N={len(charges)})")
    return maps

def lookup_bader_for_row(row: pd.Series, bader_maps: Dict[float, Dict[int, float]]) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Use *data EF* from row['EF_tag'] to select Bader file.
    For site indices (0-based), map to 1-based Index and fetch Charge.
    Return (b1, b2, b3, bavg). None if missing.
    """
    tag = row.get('EF_tag', '')
    if not tag:
        return (None, None, None, None)
    try:
        data_ef = round(float(parse_ef_value_from_tag(tag)), 1)
    except Exception:
        return (None, None, None, None)
    charges = bader_maps.get(data_ef, None)
    if charges is None:
        return (None, None, None, None)
    sites = parse_site_indices(row['path'])
    vals: List[Optional[float]] = []
    for idx0 in sites[:3]:
        idx1 = idx0 + 1  # convert to Bader 1-based Index
        vals.append(float(charges[idx1]) if idx1 in charges else None)
    while len(vals) < 3:
        vals.append(None)
    nums = [v for v in vals if v is not None]
    bavg = float(sum(nums))/len(nums) if nums else None
    return (vals[0], vals[1], vals[2], bavg)

# ---------- load table ----------
def load_table(csv_path: Path, atoms_above: Set[int] | None = None, descriptors: List[str] | None = None, tof_csv: Path | None = None) -> pd.DataFrame:
    """
    Load energy descriptors from csv_path (typically GA_prediction_summary.csv).
    Optionally merge with TOF data from tof_csv (typically drc_summary_output.csv).
    """
    df = pd.read_csv(csv_path)
    
    # Check if this CSV has 'path' column (DRC format) or 'site' column (GA format)
    has_path = 'path' in df.columns
    has_site = 'site' in df.columns
    
    # If no TOF column but we have tof_csv, merge with it
    if 'TOF' not in df.columns and tof_csv:
        print(f"[Info] Loading TOF data from: {tof_csv}")
        df_tof = pd.read_csv(tof_csv)
        
        # Filter DRC data
        if 'DRC_1_label' in df_tof.columns:
            before = len(df_tof)
            df_tof = df_tof[df_tof['DRC_1_label'].astype(str) == 'r_0006'].copy()
            print(f"[Info] Filtered by DRC_1_label == 'r_0006': kept {len(df_tof)}/{before} rows.")
        
        df_tof['TOF'] = pd.to_numeric(df_tof['TOF'], errors='coerce')
        df_tof = df_tof.dropna(subset=['TOF']).copy()
        
        # Extract site and EF from path in TOF data
        df_tof['EF_tag'] = df_tof['path'].map(extract_ef_tag)
        df_tof = df_tof[df_tof['EF_tag'] != ""].copy()
        
        # Extract site name from path (handle both "site/EF_X" and "mkm_inputs/site/EF_X" formats)
        def extract_site(path_str):
            if not isinstance(path_str, str):
                return ""
            parts = path_str.split('/')
            # If path starts with "mkm_inputs", site is the second part; otherwise first part
            if parts[0] == 'mkm_inputs' and len(parts) > 1:
                return parts[1]
            return parts[0]
        
        df_tof['site'] = df_tof['path'].map(extract_site)
        
        # Extract EF value as float from tag
        df_tof['EF'] = df_tof['EF_tag'].apply(lambda tag: float(tag.split('_')[1]) if tag else np.nan)
        
        # Merge descriptor data with TOF data on site and EF
        if has_site and 'EF' in df.columns:
            print(f"[Info] Merging {len(df)} descriptor records with {len(df_tof)} TOF records...")
            df = pd.merge(df, df_tof[['site', 'EF', 'TOF', 'path', 'EF_tag']], 
                         on=['site', 'EF'], 
                         how='inner')
            print(f"[Info] After merge: {len(df)} records with both descriptors and TOF.")
        else:
            raise SystemExit("[Error] Cannot merge: descriptor CSV missing 'site' or 'EF' column.")
        
        # Ensure path is in the merged dataframe
        if 'path' not in df.columns:
            df['path'] = df['site'] + '/' + 'EF_' + df['EF'].astype(str)
        if 'EF_tag' not in df.columns:
            df['EF_tag'] = 'EF_' + df['EF'].astype(str)
    
    elif 'TOF' not in df.columns and not tof_csv:
        print(f"[Error] 'TOF' column not found in {csv_path} and no --tof-csv provided.")
        print(f"[Suggestion] Use: python3 19_plot_volcano.py -i GA_prediction_summary.csv --tof-csv drc_summary_output.csv")
        raise SystemExit(1)
    
    # Standard processing for CSV with TOF already present
    if 'TOF' in df.columns:
        df['TOF'] = pd.to_numeric(df['TOF'], errors='coerce')
        df = df.dropna(subset=['TOF']).copy()
        
        if 'DRC_1_label' in df.columns:
            before = len(df)
            df = df[df['DRC_1_label'].astype(str) == 'r_0006'].copy()
            print(f"Filtered by DRC_1_label == 'r_0006': kept {len(df)}/{before} rows.")
        else:
            print("[Warn] Column 'DRC_1_label' not found; no DRC-based filtering applied.")
        
        if 'EF_tag' not in df.columns and 'path' in df.columns:
            df['EF_tag'] = df['path'].map(extract_ef_tag)
        df = df[df.get('EF_tag', '') != ""].copy()
    
    # Use discovered descriptors if provided, otherwise fall back to defaults
    cols_to_process = descriptors if descriptors else DESCRIPTORS
    for col in cols_to_process:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    atoms_above = atoms_above or set()
    if 'EF_real' not in df.columns:
        if 'path' in df.columns:
            df['EF_real'] = df['path'].apply(lambda p: compute_real_ef_from_row(p, atoms_above))
        elif 'EF' in df.columns and 'site' in df.columns:
            # For GA format: construct path-like string for EF_real computation
            df['EF_real'] = df.apply(lambda row: compute_real_ef_from_row(
                f"{row['site']}/EF_{row['EF']}", atoms_above), axis=1)
    
    df = df.dropna(subset=['EF_real']).copy()
    return df

def pretty_xlabel(xcol: str, x_unit: Optional[str] = None) -> str:
    mapping = {
        "E_N":   r"$E_{\mathrm{N}}$",
        "E_H":   r"$E_{\mathrm{H}}$",
        "E_NH":  r"$E_{\mathrm{NH}}$",
        "E_NH2": r"$E_{\mathrm{NH_2}}$",
        "E_NH3": r"$E_{\mathrm{NH_3}}$",
    }
    # Use mapping if available, otherwise generate subscript for E_* columns
    if xcol in mapping:
        base = mapping[xcol]
    elif xcol.startswith("E_"):
        # Generate subscript for unknown E_* columns (e.g., E_TS1 -> E_TS1)
        suffix = xcol[2:]  # remove "E_"
        base = f"$E_{{{suffix}}}$"
    else:
        base = xcol
    return f"{base} {x_unit}" if x_unit else base

# ---------- two-line fits ----------
def segmented_fit_continuous(x: np.ndarray, ylog: np.ndarray) -> Optional[Tuple[float, float, float, float]]:
    x = np.asarray(x, float); ylog = np.asarray(ylog, float)
    if len(x) < 4: return None
    x_unique = np.unique(x)
    if len(x_unique) < 3: return None
    best = None; best_sse = np.inf
    for x0 in x_unique[1:-1]:
        X1 = x - x0
        H = np.maximum(0.0, X1)
        A = np.column_stack([np.ones_like(x), X1, H])  # [a,b,c]
        coef, *_ = np.linalg.lstsq(A, ylog, rcond=None)
        yhat = A @ coef
        sse = float(np.sum((ylog - yhat)**2))
        if sse < best_sse:
            best_sse = sse
            best = (x0, coef[0], coef[1], coef[2])
    return best

def slopes_through_peak(x: np.ndarray, ylog: np.ndarray, x0: float, y0: float) -> Tuple[Optional[float], Optional[float]]:
    x = np.asarray(x, float); ylog = np.asarray(ylog, float)
    left = x <= x0; right = x >= x0
    def fit_side(mask):
        xs = x[mask]; ys = ylog[mask]
        if len(xs) < 2: return None
        dx = xs - x0; dy = ys - y0
        denom = float(np.sum(dx*dx))
        if denom == 0.0: return None
        return float(np.sum(dx*dy) / denom)
    return fit_side(left), fit_side(right)

def yfit_log_from_params(x: np.ndarray, params: Tuple, mode: str) -> np.ndarray:
    x = np.asarray(x, float)
    if mode in ("global", "per-ef"):
        x0, a, b, c = params
        return np.where(x <= x0, a + b*(x - x0), a + (b + c)*(x - x0))
    elif mode == "peak":
        x0, y0, m_left, m_right = params
        return np.where(x <= x0, y0 + m_left*(x - x0), y0 + m_right*(x - x0))
    else:
        raise ValueError("Unsupported mode for yfit")

# ---------- plotting + exports ----------
def scatter_descriptor_vs_tof(
    df: pd.DataFrame,
    xcol: str = 'E_N',
    outfile: str = 'descriptor_vs_tof.png',
    logy: bool = True,
    ymin_log10: float = -4.0,
    alpha: float = 0.9,
    s: float = 40.0,
    figsize=(5, 6),
    dpi: int = 300,
    x_unit: Optional[str] = None,
    show_legend: bool = True,
    fit2_mode: str = "none",          # 'none' | 'global' | 'per-ef' | 'peak'
    annotate_peak: bool = False,
    mark_threshold: float = 1.0,
    export_dev_csv: Optional[str] = None,
    cn_list: Optional[List[int]] = None,
    export_peak_split_csv: Optional[str] = None,
    en_maps: Optional[Dict[str, Any]] = None,
    bader_maps: Optional[Dict[float, Dict[int, float]]] = None
) -> str:
    if xcol not in df.columns:
        raise SystemExit(f"[Error] Column '{xcol}' not found in CSV.")
    d = df.dropna(subset=[xcol, 'TOF', 'EF_real']).copy()
    if logy:
        d = d[d['TOF'] > 0]
    if d.empty:
        raise SystemExit("[Error] No valid data to plot after filtering.")

    ef_vals = sorted(d['EF_real'].unique())
    fig, ax = plt.subplots(figsize=figsize)

    for val in ef_vals:
        sub = d[np.isclose(d['EF_real'], val)]
        ax.scatter(sub[xcol], sub['TOF'], label=legend_label_from_real_ef(val),
                   s=s, alpha=alpha, edgecolor='none', c=color_from_real_ef(val))

    ax.set_xlabel(pretty_xlabel(xcol, x_unit))
    if logy:
        ax.set_yscale('log', base=10)
        ax.set_ylabel(r'$\log_{10}(\mathrm{TOF})$')
        y_min_data = 10 ** (ymin_log10)
        y_max_data = float(d['TOF'].max())
        safe_max = max(y_max_data, y_min_data * 10)
        y_max_npow = 10 ** np.ceil(np.log10(safe_max * 1.05))
        ax.set_ylim(y_min_data, y_max_npow)
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=15))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{np.log10(y):.0f}"))
    else:
        ax.set_ylabel('TOF')

    xvals = d[xcol].values
    ylog = np.log10(d['TOF'].values)
    x_min, x_max = float(xvals.min()), float(xvals.max())
    mode = fit2_mode.lower()

    deviations: List[Dict] = []
    peak_rows: List[Dict] = []

    def row_common_fields(row: pd.Series) -> Dict[str, Any]:
        sites = parse_site_indices(row['path'])
        cn1, cn2, cn3, cn_avg = cn_for_site(cn_list, sites)
        en_val = lookup_en_for_row(row, en_maps) if en_maps is not None else None
        b1, b2, b3, bavg = lookup_bader_for_row(row, bader_maps or {}) if bader_maps is not None else (None, None, None, None)
        return {
            'site': path_sites_str(row['path']),
            'CN_1': cn1, 'CN_2': cn2, 'CN_3': cn3, 'CN_average': cn_avg,
            'bader_1': b1, 'bader_2': b2, 'bader_3': b3, 'bader_average': bavg,
            'TOF': float(row['TOF']),
            'EF_real': float(row['EF_real']),
            'x': float(row[xcol]),
            'log10TOF': float(np.log10(row['TOF'])),
            'E_N_ads': en_val if (en_val is None or np.isfinite(en_val)) else None
        }

    if mode == "global":
        params = segmented_fit_continuous(xvals, ylog)
        if params:
            x0, a, b, c = params
            xs_l = np.linspace(x_min, x0, 120); xs_r = np.linspace(x0, x_max, 120)
            ylog_l = a + b*(xs_l - x0); ylog_r = a + (b + c)*(xs_r - x0)
            ax.plot(xs_l, 10**ylog_l, color="#000000", linewidth=2.0, alpha=0.95, label="Two-line fit")
            ax.plot(xs_r, 10**ylog_r, color="#000000", linewidth=2.0, alpha=0.95)
            if annotate_peak:
                ax.scatter([x0], [10**a], s=70, marker='*', color="#000000", zorder=5)
                #ax.annotate(f"Peak ~ {x0:.2f}", xy=(x0, 10**a), xytext=(5,5), textcoords="offset points", fontsize=12)

            yfit_log = yfit_log_from_params(xvals, params, "global")
            diffs = ylog - yfit_log
            for i, (xi, yi, di) in enumerate(zip(xvals, ylog, diffs)):
                base = row_common_fields(d.iloc[i])
                base['position'] = 'peak' if abs(xi - x0) < 1e-9 else ('left' if xi < x0 else 'right')
                peak_rows.append(base)
                if abs(di) >= mark_threshold:
                    ax.annotate(base['site'], xy=(xi, 10**yi), xytext=(5,6), textcoords="offset points", fontsize=11)
                    deviations.append({
                        'site': base['site'],
                        'side': 'left' if xi <= x0 else 'right',
                        'diff': float(di),
                        'abs_diff': float(abs(di)),
                        'EF_real': float(d.iloc[i]['EF_real']),
                        'x': float(xi),
                        'log10TOF': float(yi),
                        'E_N_ads': base['E_N_ads'],
                        'bader_1': base['bader_1'], 'bader_2': base['bader_2'],
                        'bader_3': base['bader_3'], 'bader_average': base['bader_average'],
                    })

    elif mode == "per-ef":
        params_map: Dict[float, Optional[Tuple[float,float,float,float]]] = {}
        for val in ef_vals:
            sub = d[np.isclose(d['EF_real'], val)]
            if len(sub) < 4:
                params_map[val] = None
                continue
            params_map[val] = segmented_fit_continuous(sub[xcol].values, np.log10(sub['TOF'].values))
            params = params_map[val]
            if params:
                x0, a, b, c = params
                cfit = color_from_real_ef(val)
                xs_l = np.linspace(x_min, x0, 100); xs_r = np.linspace(x0, x_max, 100)
                ylog_l = a + b*(xs_l - x0); ylog_r = a + (b + c)*(xs_r - x0)
                ax.plot(xs_l, 10**ylog_l, color=cfit, linewidth=2.0, alpha=0.95)
                ax.plot(xs_r, 10**ylog_r, color=cfit, linewidth=2.0, alpha=0.95)
                if annotate_peak:
                    ax.scatter([x0], [10**a], s=70, marker='*', color=cfit, zorder=5)
                    ax.annotate(f"{x0:.2f}", xy=(x0, 10**a), xytext=(5,5), textcoords="offset points", fontsize=12, color=cfit)
        for val in ef_vals:
            params = params_map.get(val)
            sub = d[np.isclose(d['EF_real'], val)]
            if params is None or sub.empty:
                continue
            x0, a, b, c = params
            x_sub = sub[xcol].values
            y_sub_log = np.log10(sub['TOF'].values)
            yfit_log = yfit_log_from_params(x_sub, params, "per-ef")
            diffs = y_sub_log - yfit_log
            for (xi, yi, di, idx) in zip(x_sub, y_sub_log, diffs, sub.index.values):
                base = row_common_fields(d.loc[idx])
                base['position'] = 'peak' if abs(xi - x0) < 1e-9 else ('left' if xi < x0 else 'right')
                peak_rows.append(base)
                if abs(di) >= mark_threshold:
                    ax.annotate(base['site'], xy=(xi, 10**yi), xytext=(5,6), textcoords="offset points", fontsize=11)
                    deviations.append({
                        'site': base['site'],
                        'side': 'left' if xi <= x0 else 'right',
                        'diff': float(di), 'abs_diff': float(abs(di)),
                        'EF_real': float(val), 'x': float(xi), 'log10TOF': float(yi),
                        'E_N_ads': base['E_N_ads'],
                        'bader_1': base['bader_1'], 'bader_2': base['bader_2'],
                        'bader_3': base['bader_3'], 'bader_average': base['bader_average'],
                    })

    elif mode == "peak":
        imax = int(np.argmax(d['TOF'].values))
        x0 = float(d.iloc[imax][xcol])
        y0 = float(np.log10(d.iloc[imax]['TOF']))
        m_left, m_right = slopes_through_peak(xvals, ylog, x0, y0)
        if m_left is not None:
            xs = np.linspace(x_min, x0, 120); ylog_l = y0 + m_left*(xs - x0)
            ax.plot(xs, 10**ylog_l, color="#000000", linewidth=2.0, alpha=0.95, label="Two-line (peak-anchored)")
        if m_right is not None:
            xs = np.linspace(x0, x_max, 120); ylog_r = y0 + m_right*(xs - x0)
            ax.plot(xs, 10**ylog_r, color="#000000", linewidth=2.0, alpha=0.95)
        if annotate_peak:
            ax.scatter([x0], [10**y0], s=80, marker='*', color="#000000", zorder=6)
            # ax.annotate(f"Peak (max TOF) ~ {x0:.2f}", xy=(x0, 10**y0), xytext=(6,8), textcoords="offset points", fontsize=12)
        params = (x0, y0, m_left if m_left is not None else 0.0, m_right if m_right is not None else 0.0)
        yfit_log = yfit_log_from_params(xvals, params, "peak")
        diffs = ylog - yfit_log
        for i, (xi, yi, di) in enumerate(zip(xvals, ylog, diffs)):
            base = row_common_fields(d.iloc[i])
            base['position'] = 'peak' if abs(xi - x0) < 1e-9 else ('left' if xi < x0 else 'right')
            peak_rows.append(base)
            if abs(di) >= mark_threshold:
                # ax.annotate(base['site'], xy=(xi, 10**yi), xytext=(5,6), textcoords="offset points", fontsize=11)
                deviations.append({
                    'site': base['site'],
                    'side': 'left' if xi <= x0 else 'right',
                    'diff': float(di), 'abs_diff': float(abs(di)),
                    'EF_real': float(d.iloc[i]['EF_real']), 'x': float(xi), 'log10TOF': float(yi),
                    'E_N_ads': base['E_N_ads'],
                    'bader_1': base['bader_1'], 'bader_2': base['bader_2'],
                    'bader_3': base['bader_3'], 'bader_average': base['bader_average'],
                })

    if show_legend:
        ax.legend(frameon=False, ncols=1, loc='upper right')
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(outfile, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    # deviations CSV
    if deviations:
        dev_path = export_dev_csv if export_dev_csv else str(Path(outfile).with_suffix('').as_posix() + "_deviations.csv")
        pd.DataFrame(deviations).to_csv(dev_path, index=False)
        print(f"[Saved deviations] {dev_path}  (N={len(deviations)})")
    else:
        print("[Info] No points exceeded the deviation threshold; no deviations CSV.")

    # peak-split CSV
    if peak_rows:
        peak_path = export_peak_split_csv if export_peak_split_csv else str(Path(outfile).with_suffix('').as_posix() + "_peak_split.csv")
        cols = ['site','CN_1','CN_2','CN_3','CN_average',
                'bader_1','bader_2','bader_3','bader_average',
                'position','TOF','EF_real','x','log10TOF','E_N_ads']
        pd.DataFrame(peak_rows, columns=cols).to_csv(peak_path, index=False)
        print(f"[Saved peak split] {peak_path}  (N={len(peak_rows)})")

    return outfile

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Volcano scatter with peak classification, CN, E_N adsorption, and Bader charges.")
    ap.add_argument('-i','--input', default='mkm_outputs/mkm_out.csv', help='Input CSV with descriptors and TOF data (default: mkm_outputs/mkm_out.csv)')
    ap.add_argument('--atoms', type=str, default='atoms_above.txt', help="Atom indices file (default: atoms_above.txt)")
    ap.add_argument('--no-atoms', action='store_true', help='Ignore atoms_above (treat as empty set).')
    ap.add_argument('--x', dest='xcol', default='E_N', help='Descriptor for X axis (default: E_N)')
    ap.add_argument('--all', action='store_true', help='Plot all descriptors (combined EF).')
    ap.add_argument('--batch', action='store_true', help='For each descriptor, also make per-EF plots.')
    ap.add_argument('--out-dir', type=str, default='mkm_outputs', help='Output directory (default: mkm_outputs)')
    ap.add_argument('-o','--output', dest='outfile', default=None, help='Output image path for single mode.')
    ap.add_argument('--x-unit', type=str, default=None, help='Append unit to x label, e.g., \"(eV)\".')
    ap.add_argument('--linear', action='store_true', help='Use linear Y-axis (default is log10).')
    ap.add_argument('--ymin-log10', type=float, default=-4.0, help='Lower bound for log10(TOF) (default: -4).')
    ap.add_argument('--dpi', type=int, default=300, help='Figure DPI (default: 300)')

    # fits + deviations + peak split
    ap.add_argument('--fit2', type=str, choices=['none','global','per-ef','peak'], default='peak',
                    help='Two-line volcano fit (default: peak).')
    ap.add_argument('--annotate-peak', action='store_true', help='Annotate apex on the plot.')
    ap.add_argument('--mark-threshold', type=float, default=1.0, help='|log10(TOF)-fit| threshold for annotation/export.')
    ap.add_argument('--export-dev', type=str, default=None, help='Path for deviations CSV (default derived from figure path).')

    # CN & E_N & Bader
    ap.add_argument('--cn-file', type=str, default=None, help='Path to CN list (space/newline-separated).')
    ap.add_argument('--export-peak-split', type=str, default=None, help='Path for peak-split CSV (default derived from figure path).')
    ap.add_argument('--en-file', type=str, default=None, help='Path to E_N.csv (adsorption energies for N).')
    ap.add_argument('--bader-dir', type=str, default=None, help='Directory containing Bader files like bader_0.6.dat, bader_-0.6.dat.')

    args = ap.parse_args()

    csv_path = Path(args.input)
    if not csv_path.exists():
        raise SystemExit(f"[Error] Input file not found: {csv_path}")

    # Discover descriptors from CSV
    descriptors = discover_descriptors_from_csv(csv_path)
    
    # If using defaults (no E_* columns found), check if they actually exist
    if descriptors == DESCRIPTORS:
        missing_cols = [col for col in descriptors if col not in pd.read_csv(csv_path).columns]
        if missing_cols:
            print(f"\n[Error] Descriptor columns not found in CSV: {missing_cols}")
            print(f"        Available columns: {list(pd.read_csv(csv_path).columns[:15])}")
            print(f"[Suggestion] For volcano plots with energy descriptors, use:")
            print(f"        python3 19_plot_volcano.py -i mkm_inputs/GA_prediction_summary.csv")
            raise SystemExit(1)
    
    # Validate user's xcol choice
    if args.xcol not in descriptors:
        raise SystemExit(f"[Error] Selected descriptor '{args.xcol}' not found in CSV. Available: {descriptors}")
    
    # atoms
    atoms_set: Set[int] = set()
    if not args.no_atoms:
        atoms_path = Path(args.atoms)
        user_specified = args.atoms != 'atoms_above.txt'
        atoms_set = load_atoms_above(atoms_path, require_exists=user_specified)
    else:
        print("[Info] Using empty atoms_above set (--no-atoms).")

    cn_list   = load_cn_list(args.cn_file)
    en_maps   = load_en_lookup(args.en_file, atoms_set)
    bader_maps = load_bader_maps(args.bader_dir)

    df = load_table(csv_path, atoms_set, descriptors=descriptors)
    use_logy = not args.linear
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    def df_by_real_ef(in_df: pd.DataFrame, target: float) -> pd.DataFrame:
        return in_df[np.isclose(in_df['EF_real'], target)].copy()

    if args.batch:
        for col in descriptors:
            fig_path = out_dir / f"volcano_{col.lower()}_all.png"
            dev_path = out_dir / f"volcano_{col.lower()}_all_deviations.csv"
            split_path = out_dir / f"volcano_{col.lower()}_all_peak_split.csv"
            scatter_descriptor_vs_tof(df, xcol=col, outfile=str(fig_path),
                                      logy=use_logy, ymin_log10=args.ymin_log10,
                                      dpi=args.dpi, x_unit=args.x_unit,
                                      show_legend=True,
                                      fit2_mode=args.fit2, annotate_peak=args.annotate_peak,
                                      mark_threshold=args.mark_threshold,
                                      export_dev_csv=str(dev_path),
                                      cn_list=cn_list,
                                      export_peak_split_csv=str(split_path),
                                      en_maps=en_maps,
                                      bader_maps=bader_maps)
            print(f"[Saved] {fig_path}")

        for col in descriptors:
            for ef in (-0.6, 0.0, 0.6):
                sub = df_by_real_ef(df, ef)
                if sub.empty:
                    print(f"[Warn] No data for {col} at EF={ef:+.1f} V/Å; skip.")
                    continue
                token = ("n0p6" if ef < 0 else "p0p6") if abs(ef) > 1e-8 else "0p0"
                fig_path = out_dir / f"volcano_{col.lower()}_ef_{token}.png"
                dev_path = out_dir / f"volcano_{col.lower()}_ef_{token}_deviations.csv"
                split_path = out_dir / f"volcano_{col.lower()}_ef_{token}_peak_split.csv"
                scatter_descriptor_vs_tof(sub, xcol=col, outfile=str(fig_path),
                                          logy=use_logy, ymin_log10=args.ymin_log10,
                                          dpi=args.dpi, x_unit=args.x_unit,
                                          show_legend=False,
                                          fit2_mode=("peak" if args.fit2 == "per-ef" else args.fit2),
                                          annotate_peak=args.annotate_peak,
                                          mark_threshold=args.mark_threshold,
                                          export_dev_csv=str(dev_path),
                                          cn_list=cn_list,
                                          export_peak_split_csv=str(split_path),
                                          en_maps=en_maps,
                                          bader_maps=bader_maps)
                print(f"[Saved] {fig_path}")

    elif args.all:
        for col in descriptors:
            fig_path = out_dir / f"volcano_{col.lower()}_all.png"
            dev_path = out_dir / f"volcano_{col.lower()}_all_deviations.csv"
            split_path = out_dir / f"volcano_{col.lower()}_all_peak_split.csv"
            scatter_descriptor_vs_tof(df, xcol=col, outfile=str(fig_path),
                                      logy=use_logy, ymin_log10=args.ymin_log10,
                                      dpi=args.dpi, x_unit=args.x_unit,
                                      show_legend=True,
                                      fit2_mode=args.fit2, annotate_peak=args.annotate_peak,
                                      mark_threshold=args.mark_threshold,
                                      export_dev_csv=str(dev_path),
                                      cn_list=cn_list,
                                      export_peak_split_csv=str(split_path),
                                      en_maps=en_maps,
                                      bader_maps=bader_maps)
            print(f"[Saved] {fig_path}")

    else:
        out = args.outfile or str(out_dir / f"{args.xcol.lower()}_vs_tof.png")
        dev = args.export_dev or str(Path(out).with_suffix('').as_posix() + "_deviations.csv")
        split = args.export_peak_split or str(Path(out).with_suffix('').as_posix() + "_peak_split.csv")
        img = scatter_descriptor_vs_tof(df, xcol=args.xcol, outfile=out,
                                        logy=use_logy, ymin_log10=args.ymin_log10,
                                        dpi=args.dpi, x_unit=args.x_unit,
                                        show_legend=True,
                                        fit2_mode=args.fit2, annotate_peak=args.annotate_peak,
                                        mark_threshold=args.mark_threshold,
                                        export_dev_csv=dev,
                                        cn_list=cn_list,
                                        export_peak_split_csv=split,
                                        en_maps=en_maps,
                                        bader_maps=bader_maps)
        print(f"[Saved] {img}")

if __name__ == '__main__':
    main()
