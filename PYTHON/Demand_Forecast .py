import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Read the CSV file into a DataFrame
df1 = pd.read_csv('/Users/henryverdoodt/Documents/CODE/PYTHON/NT2040_Demand_CY1984.csv')
df2 = pd.read_csv('/Users/henryverdoodt/Documents/CODE/PYTHON/GA2040_Demand_CY1984.csv')
#df3 = pd.read_csv('NT2025_Demand_CY1984.csv')
print("OK_1")
'''
# Define a function to convert the date string to a datetime object
def parse_date(date_str):
    date_str = date_str.replace('M', '').replace('D', '').replace('H', '')
    return datetime.strptime(date_str, '%m,%d,%H')

#print(parse_date("M01,D01,H1"))
#print(df3['PATTERN'])
#print(df3['PATTERN'].apply(parse_date))

# Apply the function to the date column and create a new datetime column
df3['Datetime'] = df3['PATTERN'].apply(parse_date)

#Try alternative because issues with apply function
count=0
for i, row in df3.iterrows():
    count +=1
    date_str = row['PATTERN']
    converted_date = parse_date(date_str)
    df3.at[i, 'PATTERN'] = converted_date
print(count)

# Drop the original date column if desired
#df3.drop('PATTERN', axis=1, inplace=True)

# Print the first few rows of the DataFrame to check the new datetime column
print(df3.head())
'''

# Sum the demand from column 5 to 60 to get the total demand for each hour
df1['Total Demand'] = df1.iloc[:, 4:].sum(axis=1)
df2['Total Demand'] = df2.iloc[:, 4:].sum(axis=1)
#print(df['Total Demand'])
print("OK_2")

# Convert the year, month, day, and hour columns to a single datetime column
df1['Date'] = pd.to_datetime(df1[['YEAR', 'MONTH', 'DAY']])
df2['Date'] = pd.to_datetime(df2[['YEAR', 'MONTH', 'DAY']])
print("OK_3")
# Create a line plot of the total demand versus time
plt.plot(df1['Date'], df1['Total Demand'])
plt.plot(df2['Date'], df2['Total Demand'])
plt.xlabel('Time')
plt.ylabel('Electricity Demand (MW)')
plt.title('Electricity Demand in Europe (Year 2041)')
plt.show()

