import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


#############################################################################################################################################################
###############################                                              IMPORT DATA                                      ###############################
#############################################################################################################################################################

# WINDOFF SOUTH EU
windoff_ENTSO_25_1986_2015_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_ENTSO_25_1986_2015_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_ENTSO_25_1986_2015_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_ENTSO_25_1986_2015_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_ENTSO_25_1986_2015_spring_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_ENTSO_25_1986_2015_spring_South_EU.csv', delimiter=',', skiprows=1); 
windoff_ENTSO_25_1986_2015_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_ENTSO_25_1986_2015_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_ENTSO_25_1986_2015_autumn_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_ENTSO_25_1986_2015_autumn_South_EU.csv', delimiter=',', skiprows=1); 

windoff_ENTSO_25_1986_2015_summer_South_EU_3H = np.mean(windoff_ENTSO_25_1986_2015_summer_South_EU.reshape(-1, 3), axis=1); 
windoff_ENTSO_25_1986_2015_summer_South_EU_3H = np.mean(windoff_ENTSO_25_1986_2015_summer_South_EU.reshape(-1, 3), axis=1); 
windoff_ENTSO_25_1986_2015_spring_South_EU_3H = np.mean(windoff_ENTSO_25_1986_2015_spring_South_EU.reshape(-1, 3), axis=1); 
windoff_ENTSO_25_1986_2015_summer_South_EU_3H = np.mean(windoff_ENTSO_25_1986_2015_summer_South_EU.reshape(-1, 3), axis=1); 
windoff_ENTSO_25_1986_2015_autumn_South_EU_3H = np.mean(windoff_ENTSO_25_1986_2015_autumn_South_EU.reshape(-1, 3), axis=1); 

windoff_CNRM_1979_2005_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_1979_2005_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_1979_2005_spring_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_spring_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_1979_2005_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_1979_2005_autumn_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_autumn_South_EU.csv', delimiter=',', skiprows=1); 

windoff_CNRM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_2071_2100_spring_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_spring_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_CNRM_2071_2100_autumn_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_autumn_South_EU.csv', delimiter=',', skiprows=1); 

windoff_EARTH_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_EARTH_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_EARTH_2071_2100_spring_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_spring_South_EU.csv', delimiter=',', skiprows=1); 
windoff_EARTH_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_EARTH_2071_2100_autumn_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_autumn_South_EU.csv', delimiter=',', skiprows=1); 

windoff_HadGEM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_HadGEM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_HadGEM_2071_2100_spring_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_spring_South_EU.csv', delimiter=',', skiprows=1); 
windoff_HadGEM_2071_2100_summer_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1); 
windoff_HadGEM_2071_2100_autumn_South_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_autumn_South_EU.csv', delimiter=',', skiprows=1); 

#############################################################################################################################################################
#print(((((np.median(windoff_CNRM_2071_2100_summer_South_EU) + np.median(windoff_EARTH_2071_2100_summer_South_EU) + np.median(windoff_HadGEM_2071_2100_summer_South_EU))/3) / np.median(windoff_CNRM_1979_2005_summer_South_EU)) - 1 ) * 100); 


# Create a figure and subplots
fig, ax = plt.subplots(); 

# Plot box plots 
flier_properties = dict(marker='o', markerfacecolor='white', markeredgecolor='black', markersize=3, linewidth=0.2); 


bp_summer = ax.boxplot([windoff_CNRM_1979_2005_summer_South_EU, windoff_CNRM_2071_2100_summer_South_EU, windoff_EARTH_2071_2100_summer_South_EU, windoff_HadGEM_2071_2100_summer_South_EU], positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_summer = ax.boxplot([windoff_CNRM_1979_2005_summer_South_EU, windoff_CNRM_2071_2100_summer_South_EU, windoff_EARTH_2071_2100_summer_South_EU, windoff_HadGEM_2071_2100_summer_South_EU], positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_spring = ax.boxplot([windoff_CNRM_1979_2005_spring_South_EU, windoff_CNRM_2071_2100_spring_South_EU, windoff_EARTH_2071_2100_spring_South_EU, windoff_HadGEM_2071_2100_spring_South_EU], positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_summer = ax.boxplot([windoff_CNRM_1979_2005_summer_South_EU, windoff_CNRM_2071_2100_summer_South_EU, windoff_EARTH_2071_2100_summer_South_EU, windoff_HadGEM_2071_2100_summer_South_EU], positions=[13, 14, 15, 16], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_autumn = ax.boxplot([windoff_CNRM_1979_2005_autumn_South_EU, windoff_CNRM_2071_2100_autumn_South_EU, windoff_EARTH_2071_2100_autumn_South_EU, windoff_HadGEM_2071_2100_autumn_South_EU], positions=[17, 18, 19, 20], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 


# Calculate the median value
change4 = ((((np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1)) + np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_EARTH_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1)) + np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_HadGEM_2071_2100_summer_South_EU.csv', delimiter=',', skiprows=1)))/3) / np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/windoff_CNRM_1979_2005_summer_South_EU.csv', delimiter=',', skiprows=1))) - 1 ) * 100; 

# Print the median value
print("the changes in CF is equal to for windoff summer value in South EU:", change4); 
 


# Set colors for the boxes
colors = ['lightblue', 'lightsalmon', 'lightcoral', 'rosybrown']; 
for box, color in zip(bp_summer['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_summer['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_spring['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_summer['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_autumn['boxes'], colors):
    box.set(facecolor=color); 

# Set the x-axis ticks and labels
ax.set_xticks([2.5, 6.5, 10.5, 14.5, 18.5]); 
ax.set_xticklabels(['Annual', 'Summer', 'Spring', 'Summer', 'Autumn']); 

# Create custom legend
legend_elements = [Patch(facecolor=colors[0], label='1979-2005 (CNRM)'),
                   Patch(facecolor=colors[1], label='2071-2100 (CNRM)'),
                   Patch(facecolor=colors[2], label='2071-2100 (EARTH)'),
                   Patch(facecolor=colors[3], label='2071-2100 (HadGEM)')]; 
ax.legend(handles=legend_elements, loc='upper right'); 

# Set the title and y-axis label
plt.title('Windoff CF for South Europe'); 

# Save the plot as a PNG image
#plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_FIG/Windoff_CF_for_South_Europe_outliers_COPER_HIST.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
#plt.close()

# Show the plot
plt.show()





