import os
import subprocess
import pandas as pd
from ase import io
from cluster import classify_N2_adsorption  # Your function from previous steps

def get_final_energy(outcar_path):
    """
    Extract the final energy from OUTCAR using shell commands.
    """
    try:
        result = subprocess.run(
            "grep ' without' OUTCAR | tail -n 1 | awk '{print $NF}'",
            shell=True,
            capture_output=True,
            text=True,
            cwd=os.path.dirname(outcar_path)
        )
        energy = float(result.stdout.strip())
        return energy
    except Exception as e:
        print(f"❌ Failed to read energy from {outcar_path}: {e}")
        return None

def analyze_all_folders(base_dir="."):
    data = []

    for root, dirs, files in os.walk(base_dir):
        if "CONTCAR" in files and "OUTCAR" in files:
            contcar_path = os.path.join(root, "CONTCAR")
            outcar_path = os.path.join(root, "OUTCAR")
            folder = os.path.relpath(root, base_dir)

            try:
                atoms = io.read(contcar_path)
            except Exception as e:
                print(f"❌ Failed to read {contcar_path}: {e}")
                continue

            try:
                results = classify_N2_adsorption(atoms)
                energy = get_final_energy(outcar_path)

                for adsorption_type, ru_site_str in results:
                    data.append({
                        "folder": folder,
                        "Ads_config": adsorption_type,
                        "Ads_site": ru_site_str,
                        "E": energy
                    })

            except Exception as e:
                print(f"❌ Failed to analyze {root}: {e}")
                continue
    def save_list_good(df_analyzed, filename="list_good"):
        with open(filename, "w") as f:
            for folder in df_analyzed["folder"]:
                f.write(f"{folder}\n")
        print(f"✅ Saved good folders to {filename}") 

    def save_list_unknown(df, filename="list_unknown"):
        unknown_folders = df[df["Ads_config"] == "unknown"]["folder"].unique()
        with open(filename, "w") as f:
            for folder in unknown_folders:
                f.write(f"{folder}\n")
        print(f"✅ Saved unknown folders to {filename}")
 


    df = pd.DataFrame(data)
    df.to_csv("data.csv", index=False)
    print("✅ Step 1: Saved full data to data.csv")

    # Step 2: Analyze data and find lowest energy per Ads_site
    df_analyzed = df.loc[df.groupby('Ads_site')['E'].idxmin()].reset_index(drop=True)
    df_analyzed.to_csv("data_analyzed.csv", index=False)
    print("✅ Step 2: Saved lowest energy per Ads_site to data_analyzed.csv")

    save_list_good(df_analyzed)
    save_list_unknown(df)
if __name__ == "__main__":
    analyze_all_folders()

