import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


#############################################################################################################################################################
###############################                                              IMPORT DATA                                      ###############################
#############################################################################################################################################################

# SOLAR NORTH EU
solar_ENTSO_25_1986_2015_year_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_ENTSO_25_1986_2015_year_North_EU.csv', delimiter=',', skiprows=1); 
solar_ENTSO_25_1986_2015_winter_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_ENTSO_25_1986_2015_winter_North_EU.csv', delimiter=',', skiprows=1); 
solar_ENTSO_25_1986_2015_spring_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_ENTSO_25_1986_2015_spring_North_EU.csv', delimiter=',', skiprows=1); 
solar_ENTSO_25_1986_2015_summer_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_ENTSO_25_1986_2015_summer_North_EU.csv', delimiter=',', skiprows=1); 
solar_ENTSO_25_1986_2015_autumn_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_ENTSO_25_1986_2015_autumn_North_EU.csv', delimiter=',', skiprows=1); 

solar_ENTSO_25_1986_2015_year_North_EU_3H = np.mean(solar_ENTSO_25_1986_2015_year_North_EU.reshape(-1, 3), axis=1); 
solar_ENTSO_25_1986_2015_winter_North_EU_3H = np.mean(solar_ENTSO_25_1986_2015_winter_North_EU.reshape(-1, 3), axis=1); 
solar_ENTSO_25_1986_2015_spring_North_EU_3H = np.mean(solar_ENTSO_25_1986_2015_spring_North_EU.reshape(-1, 3), axis=1); 
solar_ENTSO_25_1986_2015_summer_North_EU_3H = np.mean(solar_ENTSO_25_1986_2015_summer_North_EU .reshape(-1, 3), axis=1); 
solar_ENTSO_25_1986_2015_autumn_North_EU_3H = np.mean(solar_ENTSO_25_1986_2015_autumn_North_EU.reshape(-1, 3), axis=1); 

solar_CNRM_2071_2100_year_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_2071_2100_year_North_EU.csv', delimiter=',', skiprows=1); 
solar_CNRM_2071_2100_winter_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_2071_2100_winter_North_EU.csv', delimiter=',', skiprows=1); 
solar_CNRM_2071_2100_spring_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_2071_2100_spring_North_EU.csv', delimiter=',', skiprows=1); 
solar_CNRM_2071_2100_summer_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_2071_2100_summer_North_EU.csv', delimiter=',', skiprows=1); 
solar_CNRM_2071_2100_autumn_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_CNRM_2071_2100_autumn_North_EU.csv', delimiter=',', skiprows=1); 

solar_EARTH_2071_2100_year_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_2071_2100_year_North_EU.csv', delimiter=',', skiprows=1); 
solar_EARTH_2071_2100_winter_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_2071_2100_winter_North_EU.csv', delimiter=',', skiprows=1); 
solar_EARTH_2071_2100_spring_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_2071_2100_spring_North_EU.csv', delimiter=',', skiprows=1); 
solar_EARTH_2071_2100_summer_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_2071_2100_summer_North_EU.csv', delimiter=',', skiprows=1); 
solar_EARTH_2071_2100_autumn_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_EARTH_2071_2100_autumn_North_EU.csv', delimiter=',', skiprows=1); 

solar_HadGEM_2071_2100_year_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_2071_2100_year_North_EU.csv', delimiter=',', skiprows=1); 
solar_HadGEM_2071_2100_winter_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_2071_2100_winter_North_EU.csv', delimiter=',', skiprows=1); 
solar_HadGEM_2071_2100_spring_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_2071_2100_spring_North_EU.csv', delimiter=',', skiprows=1); 
solar_HadGEM_2071_2100_summer_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_2071_2100_summer_North_EU.csv', delimiter=',', skiprows=1); 
solar_HadGEM_2071_2100_autumn_North_EU = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT/solar_HadGEM_2071_2100_autumn_North_EU.csv', delimiter=',', skiprows=1); 
#############################################################################################################################################################


# Create a figure and subplots
fig, ax = plt.subplots(); 

# Plot box plots 
flier_properties = dict(marker='o', markerfacecolor='white', markeredgecolor='black', markersize=3, linewidth=0.2); 


bp_year = ax.boxplot([solar_ENTSO_25_1986_2015_year_North_EU_3H, solar_CNRM_2071_2100_year_North_EU, solar_EARTH_2071_2100_year_North_EU, solar_HadGEM_2071_2100_year_North_EU], positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=True, flierprops=flier_properties); 
bp_winter = ax.boxplot([solar_ENTSO_25_1986_2015_winter_North_EU_3H, solar_CNRM_2071_2100_winter_North_EU, solar_EARTH_2071_2100_winter_North_EU, solar_HadGEM_2071_2100_winter_North_EU], positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=True, flierprops=flier_properties); 
bp_spring = ax.boxplot([solar_ENTSO_25_1986_2015_spring_North_EU_3H, solar_CNRM_2071_2100_spring_North_EU, solar_EARTH_2071_2100_spring_North_EU, solar_HadGEM_2071_2100_spring_North_EU], positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=True, flierprops=flier_properties); 
bp_summer = ax.boxplot([solar_ENTSO_25_1986_2015_summer_North_EU_3H, solar_CNRM_2071_2100_summer_North_EU, solar_EARTH_2071_2100_summer_North_EU, solar_HadGEM_2071_2100_summer_North_EU], positions=[13, 14, 15, 16], widths=0.6, patch_artist=True, showfliers=True, flierprops=flier_properties); 
bp_autumn = ax.boxplot([solar_ENTSO_25_1986_2015_autumn_North_EU_3H, solar_CNRM_2071_2100_autumn_North_EU, solar_EARTH_2071_2100_autumn_North_EU, solar_HadGEM_2071_2100_autumn_North_EU], positions=[17, 18, 19, 20], widths=0.6, patch_artist=True, showfliers=True, flierprops=flier_properties); 


# Set colors for the boxes
colors = ['lightblue', 'lightsalmon', 'lightcoral', 'rosybrown']; 
for box, color in zip(bp_year['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_winter['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_spring['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_summer['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_autumn['boxes'], colors):
    box.set(facecolor=color); 

# Set the x-axis ticks and labels
ax.set_xticks([2.5, 6.5, 10.5, 14.5, 18.5]); 
ax.set_xticklabels(['Annual', 'Winter', 'Spring', 'Summer', 'Autumn']); 

# Create custom legend
legend_elements = [Patch(facecolor=colors[0], label='1986-2015 (ENTSO)'),
                   Patch(facecolor=colors[1], label='2071-2100 (CNRM)'),
                   Patch(facecolor=colors[2], label='2071-2100 (EARTH)'),
                   Patch(facecolor=colors[3], label='2071-2100 (HadGEM)')]; 
ax.legend(handles=legend_elements, loc='upper right'); 

# Set the title and y-axis label
plt.title('Solar CF for North Europe'); 

# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_FIG/Solar_CF_for_North_Europe_outliers.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

# Show the plot
#plt.show()





