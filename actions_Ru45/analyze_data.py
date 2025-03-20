import pandas as pd

# Read the original data.csv file
df = pd.read_csv('data.csv')

# For each Ads_site, find the row with the minimum E value
df_analyzed = df.loc[df.groupby('Ads_site')['E'].idxmin()].reset_index(drop=True)

# Save the resulting data to data_analyzed.csv keeping the same format
df_analyzed.to_csv('data_analyzed.csv', index=False)
