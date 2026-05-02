import os
import pandas as pd
import numpy as np
from pathlib import Path

def summarize_drc(base_dir="mkm_inputs", output_csv="drc_summary_output.csv"):
    """
    汇总 base_dir 下所有 EF_*/outputs/drc_results.csv：
        - 初始行 (label==initial) 的 TOF
        - 非初始行 DRC 求和
        - 按 |DRC| 排序列出 (value, label) 对
    输出列：path, TOF, Sum(DRC), DRC_1_value, DRC_1_label, ...
    """
    records = []

    for root, _, files in os.walk(base_dir):
        # 只处理 'outputs' 目录
        if os.path.basename(root) != "outputs":
            continue
        if "drc_results.csv" not in files:
            continue

        drc_file = os.path.join(root, "drc_results.csv")

        # 提取 ef_folder，确保是 EF_* 格式
        ef_folder  = os.path.basename(os.path.dirname(root))          # EF_-0.6
        if not ef_folder.startswith("EF_"):
            continue
        
        # 构建完整的相对路径：mkm_inputs/intermediate_folder/EF_folder
        relative_to_base = os.path.relpath(os.path.dirname(root), base_dir)
        rel_path = f"{base_dir}/{relative_to_base}"

        try:
            df = pd.read_csv(drc_file, engine="python")
        except Exception as e:
            print(f"[Warning] Could not read {drc_file}: {e}")
            continue

        # ---------- 1. 读取 initial 行 TOF ----------
        initial_row = df[df.iloc[:, 0] == "initial"]
        if initial_row.empty:
            print(f"[Warning] No 'initial' row in {drc_file}")
            continue
        try:
            tof_initial = float(initial_row["TOF"].values[0])
        except Exception:
            print(f"[Warning] Could not parse TOF from {drc_file}")
            continue

        # ---------- 2. 处理 DRC ----------
        drc_rows = df[df.iloc[:, 0] != "initial"]
        if "DRC" not in drc_rows.columns:
            print(f"[Warning] 'DRC' column missing in {drc_file}")
            continue

        drc_sum = drc_rows["DRC"].sum()

        # 按 |DRC| 排序
        drc_sorted = drc_rows.copy()
        drc_sorted["abs_DRC"] = drc_sorted["DRC"].abs()
        drc_sorted = drc_sorted.sort_values(by="abs_DRC", ascending=False)

        drc_values = [(row["DRC"], row.iloc[0]) for _, row in drc_sorted.iterrows()]

        # ---------- 3. 组织记录 ----------
        record = [rel_path, tof_initial, drc_sum]
        for val, label in drc_values:
            record.extend([val, label])

        records.append(record)

    # ---------- 4. 写汇总 CSV ----------
    if records:
        max_len = max(len(r) for r in records)
        cols = ["path", "TOF", "Sum(DRC)"]
        for i in range((max_len - 3) // 2):
            cols.extend([f"DRC_{i+1}_value", f"DRC_{i+1}_label"])

        df_out = pd.DataFrame(records, columns=cols[:max_len])
        df_out.to_csv(output_csv, index=False)
        print(f"✅ DRC summary saved to: {output_csv}")
        return df_out
    else:
        print("⚠️  No valid drc_results.csv found under", base_dir)
        return None


def merge_ga_and_drc(ga_csv="mkm_inputs/GA_prediction_summary.csv", 
                     drc_csv="drc_summary_output.csv",
                     output_csv="mkm_outputs/mkm_out.csv"):
    """
    Merge GA prediction summary with DRC results on (site, EF) keys.
    """
    output_dir = Path(output_csv).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read GA predictions
    print(f"\n[Info] Reading GA predictions from: {ga_csv}")
    df_ga = pd.read_csv(ga_csv)
    print(f"  GA data shape: {df_ga.shape}")
    print(f"  GA columns: {df_ga.columns.tolist()[:10]}...")
    
    # Read DRC summary
    print(f"[Info] Reading DRC summary from: {drc_csv}")
    df_drc = pd.read_csv(drc_csv)
    print(f"  DRC data shape: {df_drc.shape}")
    print(f"  DRC columns: {df_drc.columns.tolist()[:10]}...")
    
    # Extract site and EF from path in drc_csv
    # Path format: mkm_inputs/site/EF_value
    def extract_site_ef(path_str):
        parts = str(path_str).split('/')
        if len(parts) >= 3:
            site = parts[1]
            ef_str = parts[2]  # e.g., "EF_-0.6"
            if ef_str.startswith("EF_"):
                try:
                    ef_val = float(ef_str.split("_")[1])
                    return site, ef_val
                except:
                    return None, None
        return None, None
    
    df_drc[['site', 'EF']] = df_drc['path'].apply(lambda x: pd.Series(extract_site_ef(x)))
    
    print(f"\n[Info] Extracted site and EF from DRC paths")
    print(f"  Sample sites: {df_drc['site'].unique()[:5].tolist()}")
    print(f"  Sample EF values: {df_drc['EF'].unique()}")
    
    # Merge on site and EF
    print(f"\n[Info] Merging on (site, EF)...")
    df_merged = pd.merge(df_ga, df_drc, on=['site', 'EF'], how='left')
    
    print(f"  Merged shape: {df_merged.shape}")
    print(f"  Rows with TOF data: {df_merged['TOF'].notna().sum()}")
    print(f"  Rows without TOF data: {df_merged['TOF'].isna().sum()}")
    
    # Save merged result
    df_merged.to_csv(output_csv, index=False)
    print(f"\n✅ Merged data saved to: {output_csv}")
    print(f"  Total columns: {len(df_merged.columns)}")
    print(f"  Column sample: {df_merged.columns.tolist()[:15]}")
    
    return df_merged


# 示例调用
if __name__ == "__main__":
    # Step 1: Collect DRC summary
    df_drc = summarize_drc(base_dir="mkm_inputs", output_csv="drc_summary_output.csv")
    
    # Step 2: Merge GA predictions with DRC results
    df_merged = merge_ga_and_drc(
        ga_csv="mkm_inputs/GA_prediction_summary.csv",
        drc_csv="drc_summary_output.csv",
        output_csv="mkm_outputs/mkm_out.csv"
    )
    
    print("\n[Done] All tasks completed successfully!")

