import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


#############################################################################################################################################################
###############################                                              IMPORT DATA                                      ###############################
#############################################################################################################################################################

# REF CASE 1986-2015
Solar_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 
WindOnshore_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 
WindOffshore_South_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 

Solar_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 
WindOnshore_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 
WindOffshore_North_CNRM_1986_2015 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_CNRM_1986_2015.csv', delimiter=',', skiprows=1); 

# CNRM RCP8.5 2071-2100
Solar_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_South_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 
WindOnshore_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_South_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 
WindOffshore_South_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_South_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 

Solar_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_North_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 
WindOnshore_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_North_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 
WindOffshore_North_CNRM_2071_2100 = np.loadtxt('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_North_CNRM_2071_2100.csv', delimiter=',', skiprows=1); 


#############################################################################################################################################################


# Create a figure and subplots
fig, ax = plt.subplots(); 

# Plot box plots 
flier_properties = dict(marker='o', markerfacecolor='white', markeredgecolor='black', markersize=3, linewidth=0.2); 


bp_solar_south = ax.boxplot([Solar_South_CNRM_1986_2015, Solar_South_CNRM_2071_2100], positions=[1, 2], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_windon_south = ax.boxplot([WindOnshore_South_CNRM_1986_2015, WindOnshore_South_CNRM_2071_2100], positions=[3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
bp_windoff_south = ax.boxplot([WindOffshore_South_CNRM_1986_2015, WindOffshore_South_CNRM_2071_2100], positions=[5, 6], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 

#bp_solar_north = ax.boxplot([Solar_North_CNRM_1986_2015, Solar_North_CNRM_2071_2100], positions=[1, 2], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windon_north = ax.boxplot([WindOnshore_North_CNRM_1986_2015, WindOnshore_North_CNRM_2071_2100], positions=[3, 4], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 
#bp_windoff_north = ax.boxplot([WindOffshore_North_CNRM_1986_2015, WindOffshore_North_CNRM_2071_2100], positions=[5, 6], widths=0.6, patch_artist=True, showfliers=False, flierprops=flier_properties); 


# Set colors for the boxes
colors = ['lightblue', 'lightsalmon' ];   #, 'lightcoral', 'rosybrown']; 
for box, color in zip(bp_solar_south['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windon_south['boxes'], colors):
    box.set(facecolor=color); 
for box, color in zip(bp_windoff_south['boxes'], colors):
    box.set(facecolor=color); 


# Set the x-axis ticks and labels
ax.set_xticks([1.5, 3.5, 5.5]); 
ax.set_xticklabels(['Solar', 'WindOnshore', 'WindOnshore']); 

# Create custom legend
legend_elements = [Patch(facecolor=colors[0], label='1986-2015 (CNRM)'),
                   Patch(facecolor=colors[1], label='2071-2100 (CNRM)')]; 
ax.legend(handles=legend_elements, loc='upper right'); 

# Set the title and y-axis label
plt.title('Installed Capacity for South Europe'); 

# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES_FIG/Installed_Capacity_for_South_Europe.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

# Show the plot
#plt.show()





