import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV files into pandas DataFrames
solar_CNRM_North = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_1979_2099_year_North_EU.csv'); 
print(solar_CNRM_North); 
solar_EARTH_North = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_1979_2099_year_North_EU.csv'); 
solar_HadGEM_North = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_1979_2099_year_North_EU.csv'); 

# Calculate the median, upper quartile, and lower quartile curves for each model
median_CNRM_North = solar_CNRM_North.groupby(solar_CNRM_North.index // (2920)).median();  
print(median_CNRM_North); 
upper_quartile_CNRM_North = solar_CNRM_North.groupby(solar_CNRM_North.index // (2920)).quantile(0.75);  
lower_quartile_CNRM_North = solar_CNRM_North.groupby(solar_CNRM_North.index // (2920)).quantile(0.25);  

median_EARTH_North = solar_EARTH_North.groupby(solar_EARTH_North.index // (2920)).median();  
upper_quartile_EARTH_North = solar_EARTH_North.groupby(solar_EARTH_North.index // (2920)).quantile(0.75);  
lower_quartile_EARTH_North = solar_EARTH_North.groupby(solar_EARTH_North.index // (2920)).quantile(0.25);  

median_HadGEM_North = solar_HadGEM_North.groupby(solar_HadGEM_North.index // (2920)).median(); 
upper_quartile_HadGEM_North = solar_HadGEM_North.groupby(solar_HadGEM_North.index // (2920)).quantile(0.75); 
lower_quartile_HadGEM_North = solar_HadGEM_North.groupby(solar_HadGEM_North.index // (2920)).quantile(0.25); 
print(lower_quartile_HadGEM_North); 

# Plotting the curves
plt.figure(figsize=(10, 6)); 

# Model 1
plt.plot(median_CNRM_North.index, median_CNRM_North['data'], label='CNRM', color='blue'); 
plt.fill_between(median_CNRM_North.index, lower_quartile_CNRM_North['data'], upper_quartile_CNRM_North['data'], alpha=0.1, color='blue'); 

# Model 2
plt.plot(median_EARTH_North.index, median_EARTH_North['data'], label='EARTH', color='green'); 
plt.fill_between(median_EARTH_North.index, lower_quartile_EARTH_North['data'], upper_quartile_EARTH_North['data'], alpha=0.1, color='green'); 

# Model 3
plt.plot(median_HadGEM_North.index, median_HadGEM_North['data'], label='HadGEM', color='red'); 
plt.fill_between(median_HadGEM_North.index, lower_quartile_HadGEM_North['data'], upper_quartile_HadGEM_North['data'], alpha=0.1, color='red'); 

plt.xlabel('Year'); 
plt.ylabel('Capacity Factor'); 
plt.title('Evolution of Solar Capacity Factor in North Europe'); 
plt.legend(); 
plt.grid(True); 
plt.xticks(range(0, len(median_CNRM_North.index), 10), range(1979, 2100, 10)); 


# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_FIG/Solar_CF_North_EU_V2.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

#plt.show()
