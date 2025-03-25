import os
import subprocess
import pandas as pd
from ase import io
from cluster import classify_N2_adsorption

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
        print(f"Failed to read energy from {outcar_path}: {e}")
        return None

def analyze_all_folders(base_dir="."):
    data = []
    
    for root, dirs, files in os.walk(base_dir):
        if "CONTCAR" in files and "OUTCAR" in files:
            contcar_path = os.path.join(root, "CONTCAR")
            outcar_path = os.path.join(root, "OUTCAR")
            
            try:
                atoms = io.read(contcar_path)
            except Exception as e:
                print(f"Failed to read {contcar_path}: {e}")
                continue

            try:
                results = classify_N2_adsorption(atoms)
                energy = get_final_energy(outcar_path)

                # Handle multiple N2 molecules, if present
                for n1, n2, adsorption_type, ru_site_str in results:
                    data.append({
                        "folder": os.path.relpath(root, base_dir),
                        "Ads_config": adsorption_type,
                        "Ads_site": ru_site_str,
                        "E": energy
                    })
            except Exception as e:
                print(f"Failed to analyze {root}: {e}")
                continue

    df = pd.DataFrame(data)
    df.to_csv("data.csv", index=False)
    print("✅ Analysis complete. Saved to data.csv")

if __name__ == "__main__":
    analyze_all_folders()

