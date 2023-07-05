import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#############################################################################################################################################################
###############################                                              IMPORT DATA                                      ###############################
#############################################################################################################################################################

############################### READ ME ###############################
# Ctrl F to change the names of the desired files
# Comment/ uncomment the SOUTH or NORTH PLOT 
# Do not forget to change the the Title AND the name of the saved file
#######################################################################

# Read the CSV files
df_ref = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/Solar_CNRM_1986_2015_fixed_network_V2.csv'); 
df_model1 = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_CNRM_1986_2015_fixed_network_V2.csv'); 
df_model2 = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOnshore_CNRM_1986_2015_fixed_network_V2.csv'); 
df_model3 = pd.read_csv('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES/WindOffshore_CNRM_1986_2015_fixed_network_V2.csv'); 

print(df_ref); 


# Calculate the percentage change for each model
df_change1 = (df_model1 - df_ref) / df_ref * 100; 
df_change2 = (df_model2 - df_ref) / df_ref * 100; 
df_change3 = (df_model3 - df_ref) / df_ref * 100; 

# Calculate the mean and quartiles for each country
mean_change1 = df_change1.mean(); 
mean_change2 = df_change2.mean(); 
mean_change3 = df_change3.mean(); 
q25_change = df_change1.quantile(0.25); 
q75_change = df_change1.quantile(0.75); 

# Plot the data
countries = df_ref.columns; 
x = range(len(countries)); 
width = 0.2; 

fig, ax = plt.subplots(); 
ax.bar(x, mean_change1, width, label='Model 1'); 
ax.bar([i + width for i in x], mean_change2, width, label='Model 2'); 
ax.bar([i + 2*width for i in x], mean_change3, width, label='Model 3'); 

# Plot quartile lines
for i in x:
    ax.plot([i, i], [q25_change[i], q75_change[i]], color='black', linewidth=1); 

ax.axhline(0, color='black', linestyle='--');  # Zero line

ax.set_ylabel('Percentage Change'); 
ax.set_xlabel('Country'); 
ax.set_xticks([i + width for i in x]); 
ax.set_xticklabels(countries); 
ax.legend(); 

# Set the title and y-axis label
plt.title('Europe Solar Power'); 
plt.ylabel('Percentage Change (%)'); 

# Save the plot as a PNG image
plt.savefig('/Users/henryverdoodt/Documents/CODE/DATA/BOXPLOT_RES_FIG/Installed_Solar_Capacity_Changes_test.png', dpi=300);  # Set the desired filename and dpi value

# Close the plot
plt.close()

# Show the plot
#plt.show()


###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################



