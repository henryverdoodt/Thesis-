import pandas as pd
import pyarrow.parquet as pq

# Read in the Parquet file
table = pq.read_table('/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2030.parquet')

# Convert to a pandas dataframe
df = table.to_pandas()

# Write the data to a CSV file
df.to_csv('/Users/henryverdoodt/Documents/CODE/DATA/ENTSO/WINDOFF/PECD-2021.3-country-Offshore-2030.csv', index=False)
