import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV files into pandas DataFrames
windoff_CNRM_South = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2099_year_South_EU.csv'); 
print(windoff_CNRM_South); 
windoff_EARTH_South = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_1979_2099_year_South_EU.csv'); 
windoff_HadGEM_South = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_1979_2099_year_South_EU.csv'); 

# Calculate the median, upper quartile, and lower quartile curves for each model
median_CNRM_South = windoff_CNRM_South.groupby(windoff_CNRM_South.index // 2920).median();  
print(median_CNRM_South); 
upper_quartile_CNRM_South = windoff_CNRM_South.groupby(windoff_CNRM_South.index // 2920).quantile(0.75);  
lower_quartile_CNRM_South = windoff_CNRM_South.groupby(windoff_CNRM_South.index // 2920).quantile(0.25);  

median_EARTH_South = windoff_EARTH_South.groupby(windoff_EARTH_South.index // 2920).median();  
upper_quartile_EARTH_South = windoff_EARTH_South.groupby(windoff_EARTH_South.index // 2920).quantile(0.75);  
lower_quartile_EARTH_South = windoff_EARTH_South.groupby(windoff_EARTH_South.index // 2920).quantile(0.25);  

median_HadGEM_South = windoff_HadGEM_South.groupby(windoff_HadGEM_South.index // 2920).median(); 
upper_quartile_HadGEM_South = windoff_HadGEM_South.groupby(windoff_HadGEM_South.index // 2920).quantile(0.75); 
lower_quartile_HadGEM_South = windoff_HadGEM_South.groupby(windoff_HadGEM_South.index // 2920).quantile(0.25); 
print(lower_quartile_HadGEM_South); 

# Plotting the curves
plt.figure(figsize=(10, 6)); 

# Model 1
plt.plot(median_CNRM_South.index, median_CNRM_South['data'], label='CNRM', color='blue'); 
plt.fill_between(median_CNRM_South.index, lower_quartile_CNRM_South['data'], upper_quartile_CNRM_South['data'], alpha=0.3, color='blue'); 

# Model 2
plt.plot(median_EARTH_South.index, median_EARTH_South['data'], label='EARTH', color='green'); 
plt.fill_between(median_EARTH_South.index, lower_quartile_EARTH_South['data'], upper_quartile_EARTH_South['data'], alpha=0.3, color='green'); 

# Model 3
plt.plot(median_HadGEM_South.index, median_HadGEM_South['data'], label='HadGEM', color='red'); 
plt.fill_between(median_HadGEM_South.index, lower_quartile_HadGEM_South['data'], upper_quartile_HadGEM_South['data'], alpha=0.3, color='red'); 

plt.xlabel('Year'); 
plt.ylabel('Capacity Factor'); 
plt.title('Evolution of Windoff Capacity Factor in South Europe'); 
plt.legend(); 
plt.grid(True); 
plt.xticks(range(0, len(median_CNRM_South.index), 10), range(1979, 2100, 10)); 


# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_FIG/Windoff_CF_South_EU.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

#plt.show()
