import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import os 

# Load the data
input_file_path = 'raw_data.csv'  # Update the path if needed
df = pd.read_csv(input_file_path)

# Convert '0.0' column to numeric, handling any non-numeric values
df['0.0'] = pd.to_numeric(df['0.0'], errors='coerce')

# Get the minimum energy at 0.0 for each site, filtering out groups with all NaN values
df_unique_list = []
for site, group in df.groupby('site'):
    if not group['0.0'].isna().all():  # Only process if group has at least some non-NaN values
        min_idx = group['0.0'].idxmin()
        df_unique_list.append(df.loc[min_idx])

if df_unique_list:
    df_unique = pd.DataFrame(df_unique_list)
else:
    # If no valid groups found, use the original dataframe
    df_unique = df.copy()

# Function to determine site type based on the number of underscores
def determine_site_type(site):
    site = str(site)
    if isinstance(site, str):
        underscore_count = site.count('_')
        if underscore_count == 1:
            return 'bri'
        elif underscore_count == 2:
            return 'hollow'
        elif underscore_count == 3:
            return 'square'
        elif underscore_count == 0:
            return 'top'
        else:
            return 'N/A'

# Add the 'site_type' column
df_unique['site_type'] = df_unique['site'].apply(determine_site_type)
df_unique = df_unique.fillna('nan')

print('Good Datapoints\n')
for data in   df_unique['Species']:
    print(data)
print('Rename the Good Datapoints')    
for species, site in zip(df_unique['Species'], df_unique['site']):
    print(f"mv {species} {site}")



def get_dipole_polarization_from_sorted_csv(data_df, fill_nans_with_fit=True, verbose=True):
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

    def func(x, a, b, c):
        return -0.5 * a * x**2 - b * x + c

    # 电场范围 -0.6 到 0.6（步长0.2）
    x_values = np.round(np.arange(-0.6, 0.61, 0.2), 1)
    eads_cols = [f"{x:.1f}" for x in x_values]
    eads_cols = [str(x) if x != 0 else '0.0' for x in x_values]

    results_list = []
    data_df = data_df.copy()

    for idx, row in data_df.iterrows():
        site = row['site']
        y_values_orig = row[eads_cols].values.astype(float)
        valid_mask = ~np.isnan(y_values_orig)
        x_fit = x_values[valid_mask]
        y_fit = y_values_orig[valid_mask]

        if len(x_fit) < 3:
            if verbose:
                print(f"Skipping site '{site}' (index {idx}) due to insufficient data points.")
            continue

        try:
            # 初次拟合
            params, _ = curve_fit(func, x_fit, y_fit)
            y_pred = func(x_fit, *params)
            deviation = np.abs(y_pred - y_fit)
            outlier_found = False

            if np.any(deviation > 0.2):
                outlier_found = True
                if verbose:
                    print(f"Site '{site}' (index {idx}) has deviations > 0.2 eV:")
                for xi, yi, ypi, dev in zip(x_fit, y_fit, y_pred, deviation):
                    if dev > 0.2:
                        if verbose:
                            print(f"  E={xi:.1f}: actual={yi:.3f}, predicted={ypi:.3f}, Δ={dev:.3f}")
                        data_df.at[idx, f"{xi:.1f}"] = np.nan

                # 重新拟合
                y_clean = data_df.loc[idx, eads_cols].values.astype(float)
                valid_mask = ~np.isnan(y_clean)
                x_fit = x_values[valid_mask]
                y_fit = y_clean[valid_mask]

                if len(x_fit) < 3:
                    if verbose:
                        print(f"After removing outliers, not enough data for site '{site}' (index {idx})")
                    continue

                params, _ = curve_fit(func, x_fit, y_fit)
                y_pred = func(x_fit, *params)

            # 没有异常：填补原始 NaN
            elif fill_nans_with_fit and np.any(np.isnan(y_values_orig)):
                full_pred = func(x_values, *params)
                for j, val in enumerate(y_values_orig):
                    if np.isnan(val):
                        if verbose:
                            print(f"  Filling NaN for site '{site}' at E={x_values[j]:.1f} with predicted value {full_pred[j]:.3f}")
                        data_df.at[idx, f"{x_values[j]:.1f}"] = full_pred[j]

            a, b, c = params
            R2 = r2_score(y_fit, y_pred)
            mae = mean_absolute_error(y_fit, y_pred)
            rmse = np.sqrt(mean_squared_error(y_fit, y_pred))

            results_list.append({
                'site': site,
                'polarizability': a,
                'dipole': b,
                'c': c,
                # 'R2': R2,
                # 'MAE': mae,
                # 'RMSE': rmse
            })

        except Exception as e:
            if verbose:
                print(f"Error fitting site '{site}' (index {idx}): {e}")
            continue

    results_df = pd.DataFrame(results_list)
    return results_df, data_df

results_df, updated_df = get_dipole_polarization_from_sorted_csv(df_unique)

# 保存拟合结果和清洗后表格
results_df.to_csv("fitting_results.csv", index=False)
updated_df.to_csv("sorted_data.csv", index=False)

### Check the cleaned data
import pandas as pd
import matplotlib.pyplot as plt

# 读取 CSV 文件
df = pd.read_csv('sorted_data.csv')  # 替换为你的实际文件名

# 提取电场列名（浮点型）
ef_columns = [col for col in df.columns if col.replace('-', '').replace('.', '').isdigit()]
ef_values = [float(col) for col in ef_columns]

os.makedirs("2_data_analysis", exist_ok=True)

# 为每个 Species 作图
for species, group in df.groupby('Species'):
    plt.figure(figsize=(8, 6))
    for _, row in group.iterrows():
        energies = row[ef_columns].values.astype(float)
        plt.plot(ef_values, energies, marker='o', label=row['site'])  # 可选择是否加 label
        
    plt.xlabel('Electric field (V/Å)', fontsize=14)
    plt.ylabel('Adsorption energy (eV)', fontsize=14)
    plt.title(f'{species}', fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join("2_data_analysis", f"{species}.png"))
    plt.close()

