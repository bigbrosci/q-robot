import os
import sys
import platform

if platform.system() == 'Windows':
    q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
else:
    q_robot_path = '/home/qli/bin/q-robot/brain'

sys.path.insert(0, q_robot_path)

from cluster import *
from mkm import *
from pathlib import Path
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"]   = 14


data_path = Path('.')
slab_path = Path('./slab')

atoms_in = read(str(slab_path / 'POSCAR'))

connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites = get_connection(atoms_in, metal='Ru', mult=0.9)
# print(exposed_top_sites)
sites_all = np.array(hollow_sites, dtype=int).tolist() # 0-based
initial_dir = os.getcwd()
# os.chdir(initial_dir)

with open("atoms_above.txt") as f:  #  1-based file
    atoms_above = {int(line.strip()) for line in f if line.strip()}

# Initialize list to collect all data for summary CSV
all_records = []

for site_0 in sites_all:    
    # Convert to 1-based first
    site_1 = [i+1 for i in site_0]   # convert to 1-based for obtain the energies from GA
    
    # only keep sites whose atoms are all in atoms_above
    if not all(a in atoms_above for a in site_1):
        continue
    
    for EF in [-0.6, 0.0, 0.6]:
        print(f"Site-EF: {site_1} {EF}")
        
        try:
            edft, e_slab = compute_EDFT(data_path, site_1, str(EF))
        except (AttributeError, ValueError) as e:
            print(f"Warning: Could not compute EDFT - {str(e)[:100]}")
            continue
        
        base_dir = "mkm_inputs"
        index_folder = "_".join(map(str, site_1))
        ef_folder = f"EF_{EF:.1f}"  # e.g. EF_-0.3
        
        fig_outputs_folder = os.path.join(base_dir, index_folder, ef_folder)

        try:
            # build_refdict_and_plot_rc_double_NH3(EF, adsorption_data=edft, outputs_folder=fig_outputs_folder)
            replacement_dict = generate_replacement_dict(EF, edft, data_path)
            update_excel_with_replacement(site_1, replacement_dict, EF)
            
            # Collect data for summary CSV
            record = {'site': '_'.join(map(str, site_1)), 'EF': EF}
            record.update(replacement_dict)
            all_records.append(record)
            
        except Exception as e:
            print(f"Warning: Error in replacement dict/update - {str(e)[:100]}")
            continue
        
        os.makedirs(fig_outputs_folder, exist_ok=True)

        out_csv = os.path.join(fig_outputs_folder, "Eads_Ea.csv")
        
        with open(out_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["E_type", "Eads_eV", "EDFT_eV"])  # header
            for species, energy_tuple in edft.items():
                # edft returns tuples like (site, Eads, EDFT)
                if isinstance(energy_tuple, tuple) and len(energy_tuple) >= 3:
                    site_id, eads, edft_val = energy_tuple[0], energy_tuple[1], energy_tuple[2]
                    writer.writerow([species, f"{eads:.3f}", f"{edft_val:.3f}"])
                else:
                    writer.writerow([species, str(energy_tuple), ""])
        
        print(f"Saved: {out_csv}")

# Save summary CSV with all records
if all_records:
    df_summary = pd.DataFrame(all_records)
    summary_csv = os.path.join("mkm_inputs", "GA_prediction_summary.csv")
    df_summary.to_csv(summary_csv, index=False)
    print(f"\n✅ Saved summary: {summary_csv}")
    print(f"   Total records: {len(all_records)}")
else:
    print("\n⚠️ No records collected for summary CSV")

