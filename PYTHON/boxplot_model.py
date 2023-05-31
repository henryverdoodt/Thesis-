import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


#############################################################################################################################################################
###############################                                              IMPORT DATA                                      ###############################
#############################################################################################################################################################

############################### READ ME ###############################
# Ctrl F to change the names of the desired files
# Comment/ uncomment the SOUTH or NORTH PLOT 
# Do not forget to change the the Title AND the name of the saved file
#######################################################################

# REF CASE 1986-2015
Solar_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 

Solar_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1); 

# CNRM RCP8.5 2071-2100
Solar_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 

Solar_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 

# EARTH RCP8.5 2071-2100
Solar_South_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_South_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_South_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 

Solar_North_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_North_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_North_EARTH_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 

# HadGEM RCP8.5 2071-2100
Solar_South_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_South_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_South_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 

Solar_North_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOnshore_North_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
WindOffshore_North_HadGEM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1); 
#############################################################################################################################################################

# Calculate the median value
#change = ((((np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1)) + np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1)) + np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1)))/3) / np.median(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1))) - 1 ) * 100; 
change = ((((np.percentile(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1), 75) + np.percentile(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_EARTH_2071_2099_fixed_network.csv', delimiter=',', skiprows=1), 75) + np.percentile(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_HadGEM_2071_2099_fixed_network.csv', delimiter=',', skiprows=1), 75))/3) / np.percentile(np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_1986_2015_fixed_network.csv', delimiter=',', skiprows=1), 75)) - 1) * 100

# Print the median value
print("Change in installed capacity for solar in South EU:", change); 


# Create a figure and subplots
fig, ax = plt.subplots(); 

# Plot box plots 
flier_properties = dict(marker='o', markerfacecolor='white', markeredgecolor='black', markersize=3, linewidth=0.2); 

### SOUTH EU PLOT ###
'''
# Original data for boxplots
data_solar_south = [Solar_South_CNRM_1986_2015, Solar_South_CNRM_2071_2100, Solar_South_EARTH_2071_2100, Solar_South_HadGEM_2071_2100]
data_windon_south = [WindOnshore_South_CNRM_1986_2015, WindOnshore_South_CNRM_2071_2100, WindOnshore_South_EARTH_2071_2100, WindOnshore_South_HadGEM_2071_2100]
data_windoff_south = [WindOffshore_South_CNRM_1986_2015, WindOffshore_South_CNRM_2071_2100, WindOffshore_South_EARTH_2071_2100, WindOffshore_South_HadGEM_2071_2100]

# Divide all data by 1000
data_solar_south = [[x / 1000 for x in sublist] for sublist in data_solar_south]
data_windon_south = [[x / 1000 for x in sublist] for sublist in data_windon_south]
data_windoff_south = [[x / 1000 for x in sublist] for sublist in data_windoff_south]

# Plot the boxplots with the updated data
bp_solar_south = ax.boxplot(data_solar_south, positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)
bp_windon_south = ax.boxplot(data_windon_south, positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)
bp_windoff_south = ax.boxplot(data_windoff_south, positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)

#bp_solar_south = ax.boxplot([Solar_South_CNRM_1986_2015, Solar_South_CNRM_2071_2100, Solar_South_EARTH_2071_2100, Solar_South_HadGEM_2071_2100], positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windon_south = ax.boxplot([WindOnshore_South_CNRM_1986_2015, WindOnshore_South_CNRM_2071_2100, WindOnshore_South_EARTH_2071_2100, WindOnshore_South_HadGEM_2071_2100], positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windoff_south = ax.boxplot([WindOffshore_South_CNRM_1986_2015, WindOffshore_South_CNRM_2071_2100, WindOffshore_South_EARTH_2071_2100, WindOffshore_South_HadGEM_2071_2100], positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 

# Set colors for the boxes
colors = ['lightblue', 'lightsalmon' , 'lightcoral', 'rosybrown'];    #, 'lightcoral', 'rosybrown']; 
for box, color in zip(bp_solar_south['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windon_south['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windoff_south['boxes'], colors):
    box.set(facecolor=color); 
'''
### NORTH EU PLOT ###

# Original data for boxplots
data_solar_north = [Solar_North_CNRM_1986_2015, Solar_North_CNRM_2071_2100, Solar_North_EARTH_2071_2100, Solar_North_HadGEM_2071_2100]
data_windon_north = [WindOnshore_North_CNRM_1986_2015, WindOnshore_North_CNRM_2071_2100, WindOnshore_North_EARTH_2071_2100, WindOnshore_North_HadGEM_2071_2100]
data_windoff_north = [WindOffshore_North_CNRM_1986_2015, WindOffshore_North_CNRM_2071_2100, WindOffshore_North_EARTH_2071_2100, WindOffshore_North_HadGEM_2071_2100]

# Divide all data by 1000
data_solar_north = [[x / 1000 for x in sublist] for sublist in data_solar_north]
data_windon_north = [[x / 1000 for x in sublist] for sublist in data_windon_north]
data_windoff_north = [[x / 1000 for x in sublist] for sublist in data_windoff_north]

# Plot the boxplots with the updated data
bp_solar_north = ax.boxplot(data_solar_north, positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)
bp_windon_north = ax.boxplot(data_windon_north, positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)
bp_windoff_north = ax.boxplot(data_windoff_north, positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties)


#bp_solar_north = ax.boxplot([Solar_North_CNRM_1986_2015, Solar_North_CNRM_2071_2100, Solar_North_EARTH_2071_2100, Solar_North_HadGEM_2071_2100], positions=[1, 2, 3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windon_north = ax.boxplot([WindOnshore_North_CNRM_1986_2015, WindOnshore_North_CNRM_2071_2100, WindOnshore_North_EARTH_2071_2100, WindOnshore_North_HadGEM_2071_2100], positions=[5, 6, 7, 8], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windoff_north = ax.boxplot([WindOffshore_North_CNRM_1986_2015, WindOffshore_North_CNRM_2071_2100, WindOffshore_North_EARTH_2071_2100, WindOffshore_North_HadGEM_2071_2100], positions=[9, 10, 11, 12], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 


# Set colors for the boxes
colors = ['lightblue', 'lightsalmon' , 'lightcoral', 'rosybrown'];    #, 'lightcoral', 'rosybrown']; 
for box, color in zip(bp_solar_north['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windon_north['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windoff_north['boxes'], colors):
    box.set(facecolor=color); 


# Set the x-axis ticks and labels
ax.set_xticks([2.5, 6.5, 10.5]); 
ax.set_xticklabels(['Solar', 'WindOnshore', 'WindOffshore']); 

# Create custom legend
legend_elements = [Patch(facecolor=colors[0], label='1986-2015 (CNRM)'),
                   Patch(facecolor=colors[1], label='2071-2099 (CNRM)'),
                   Patch(facecolor=colors[2], label='2071-2099 (EARTH)'),
                   Patch(facecolor=colors[3], label='2071-2099 (HadGEM)')]; 
ax.legend(handles=legend_elements, loc='center right'); 

# Set the title and y-axis label
plt.title('Installed Capacity for South Europe (fixed network)'); 
plt.ylabel('RES Capacity (GW)'); 

# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES_FIG/Installed_Capacity_for_North_Europe_fixed_network.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

# Show the plot
#plt.show()





