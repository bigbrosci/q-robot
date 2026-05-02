import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ===== Global font & output settings (Times New Roman + math fonts) =====
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 16,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "axes.labelsize": 22,
    "axes.titlesize": 18,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "legend.fontsize": 18,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ===== Inputs =====
FILE_PATH = "GA_matrix.csv"
ATOMS_FILE = "atoms_above.txt"  # atom indices that indicate "above" surface

# ===== Helpers =====

def classify_site_type(site):
    """Classify by number of '_' in site string: top / bridge / hollow."""
    count = str(site).count('_')
    if count == 0:
        return 'top'
    elif count == 1:
        return 'bridge'
    elif count == 2:
        return 'hollow'
    else:
        return 'unknown'


def parse_site_indices(site):
    """Parse site like '12' or '11_17_28' -> [12] or [11,17,28]."""
    s = str(site).strip()
    if not s:
        return []
    out = []
    for tok in s.split('_'):
        tok = tok.strip()
        if tok and tok.lstrip('-').isdigit():
            out.append(int(tok))
    return out


def load_atoms_above(fp: Path) -> set:
    """Load integer indices from atoms_above.txt."""
    if not fp.exists():
        raise FileNotFoundError(
            f"'{fp}' not found. Please create it in the same folder, e.g.\n"
            f"cat > {fp.name} << 'EOF'\n4\n10\n12\n...\nEOF"
        )
    s = set()
    for line in fp.read_text().splitlines():
        t = line.strip()
        if t and t.lstrip('-').isdigit():
            s.add(int(t))
    return s


def fmt_field(x: float) -> str:
    """Format EF value to match column suffix in GA_matrix.csv."""
    return f"{float(x):.1f}"


def site_has_any_above_atom(site: str, atoms_above: set) -> bool:
    """Keep site if it has AT LEAST ONE atom index in atoms_above.txt."""
    idx = parse_site_indices(site)
    if len(idx) == 0:
        return False
    return any(i in atoms_above for i in idx)


def is_site_above_all(site: str, atoms_above: set) -> bool:
    """
    For REAL EF mapping:
      - If ALL atoms in site are in atoms_above.txt -> treat as 'above' (real EF = -vasp)
      - Otherwise -> treat as 'not-above' (real EF = vasp)
    """
    idx = parse_site_indices(site)
    if len(idx) == 0:
        return False
    return all(i in atoms_above for i in idx)


# ===== Read data =====
df = pd.read_csv(FILE_PATH)

# Required columns check
if 'site' not in df.columns:
    raise ValueError("GA_matrix.csv must contain column 'site'.")

# ===== Load atoms_above and FILTER to above-adsorption only =====
atoms_above = load_atoms_above(Path(ATOMS_FILE))

df['has_above_atom'] = df['site'].apply(lambda s: site_has_any_above_atom(s, atoms_above))
df_above = df[df['has_above_atom']].copy()

if len(df_above) == 0:
    raise RuntimeError(
        "After filtering (site has any atom in atoms_above.txt), no rows remain.\n"
        "Check atoms_above.txt indices and the 'site' column formatting."
    )

# ===== Site types =====
df_above['site_type'] = df_above['site'].apply(classify_site_type)

# ===== Colors by site_type =====
colors = {
    'top':    '#E24A33',
    'bridge': '#348ABD',
    'hollow': '#988ED5',
}

# ===== Real EF logic using atoms_above.txt =====
df_above['is_above_all'] = df_above['site'].apply(lambda s: is_site_above_all(s, atoms_above))

# ===== Selected REAL EF panels (V/Å) =====
REAL_EFS = [-0.6, 0.0, 0.6]

# Validate required source columns exist
needed_vasp_fields = {fmt_field(v) for v in REAL_EFS}
needed_vasp_fields.add(fmt_field(-0.6))
needed_vasp_fields.add(fmt_field(0.6))
missing = [f"E_ads_{f}" for f in sorted(needed_vasp_fields) if f"E_ads_{f}" not in df_above.columns]
if missing:
    raise ValueError(f"Missing columns in GA_matrix.csv: {missing}")

# Build real-EF energy columns (above-only subset)
for ef in REAL_EFS:
    ef_str = fmt_field(ef)
    col_real = f"E_ads_real_{ef_str}"
    if abs(ef) < 1e-12:
        df_above[col_real] = pd.to_numeric(df_above[f"E_ads_{fmt_field(0.0)}"], errors='coerce')
    else:
        # if is_above_all: real EF = -vasp -> vasp = -real
        # else (mixed): real EF = vasp -> vasp = real
        src_above = f"E_ads_{fmt_field(-ef)}"
        src_below = f"E_ads_{fmt_field(ef)}"
        df_above[col_real] = np.where(
            df_above['is_above_all'].values,
            pd.to_numeric(df_above[src_above], errors='coerce'),
            pd.to_numeric(df_above[src_below], errors='coerce')
        )

# ===== Build subplots =====
fig, axes = plt.subplots(1, len(REAL_EFS), figsize=(15, 6))

# ===== Unified bins across panels =====
all_values = pd.concat([
    df_above[f"E_ads_real_{fmt_field(ef)}"].dropna()
    for ef in REAL_EFS
])

bins = np.arange(
    all_values.min() - 0.05,
    all_values.max() + 0.1,
    0.1
)

# ===== Global max y across stacked hist =====
max_y = 0
for ef in REAL_EFS:
    total_hist = np.zeros(len(bins) - 1)
    col_real = f"E_ads_real_{fmt_field(ef)}"
    for site_type in ['top', 'bridge', 'hollow']:
        values = df_above[df_above['site_type'] == site_type][col_real].dropna()
        hist, _ = np.histogram(values, bins=bins)
        total_hist += hist
    max_y = max(max_y, total_hist.max())

# ===== Plot stacked hist per REAL EF =====
for ax, ef in zip(axes, REAL_EFS):
    bottom = np.zeros(len(bins) - 1)
    col_real = f"E_ads_real_{fmt_field(ef)}"

    for site_type in ['top', 'bridge', 'hollow']:
        values = df_above[df_above['site_type'] == site_type][col_real].dropna()
        hist, _ = np.histogram(values, bins=bins)

        ax.bar(
            bins[:-1],
            hist,
            width=0.1,
            align='edge',
            bottom=bottom,
            color=colors[site_type],
            label=site_type,
            edgecolor='black'
        )

        # annotate counts in the middle of each stacked segment
        for j, h in enumerate(hist):
            if h > 0:
                ax.text(
                    bins[j] + 0.05,
                    bottom[j] + h / 2,
                    str(int(h)),
                    ha='center',
                    va='center',
                    fontsize=14
                )
        bottom += hist

    # panel label: REAL EF
    ef_display = 0.0 if abs(ef) < 1e-5 else ef
    ax.text(
        0.02, 0.95,
        f"EF = {ef_display:+.1f} V/Å" if ef_display != 0 else "EF = 0.0 V/Å",
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=18,
        weight='bold'
    )

    # y limits and ticks
    y_max_plot = max_y + 2
    ax.set_ylim(0, y_max_plot)
    ax.set_yticks(np.arange(0, y_max_plot, 4))

    # x limits
    ax.set_xlim(bins[0], bins[-1])

# ===== Layout & labels =====
plt.subplots_adjust(wspace=0.0)
axes[0].legend(loc='upper right')
plt.tight_layout(rect=[0.06, 0.08, 0.98, 0.96])

fig.text(
    0.5,
    0.07,
    r'Distribution of N adsorption energies (eV) at real EF = -0.6, 0.0, and +0.6 V/Å (sites with any above-atom only)',
    ha='center',
    va='center',
    fontsize=22
)

fig.text(
    0.07,
    0.5,
    'Number of sites',
    va='center',
    ha='center',
    rotation='vertical',
    fontsize=24
)

# ===== Save =====
output_path = 'Ads_distribution_realEF_aboveOnly.png'
plt.savefig(output_path, dpi=300)
plt.show()
plt.close()

print(f"Saved figure to {output_path}")
print(f"Loaded {len(atoms_above)} atoms from {ATOMS_FILE}.")
print(f"Filtered to above-only: kept {len(df_above)} / {len(df)} rows (site contains any atom in atoms_above.txt).")
print("Real EF mapping rule for kept rows:")
print("  - if ALL atoms in site are in atoms_above.txt -> real EF = - (VASP parameter)")
print("  - otherwise (mixed site) -> real EF = (VASP parameter)")
