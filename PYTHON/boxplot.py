import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Generate random data for box plots
annual_2000 = np.random.normal(loc=2, scale=1, size=6000)
annual_2100 = np.random.normal(loc=3, scale=1, size=7000)
winter_2000 = np.random.normal(loc=2, scale=1, size=6000)
winter_2100 = np.random.normal(loc=3, scale=1, size=7000)
spring_2000 = np.random.normal(loc=2, scale=1, size=6000)
spring_2100 = np.random.normal(loc=3, scale=1, size=7000)
summer_2000 = np.random.normal(loc=0, scale=1, size=5000)
summer_2100 = np.random.normal(loc=1, scale=1, size=4000)
autumn_2000 = np.random.normal(loc=2, scale=1, size=6000)
autumn_2100 = np.random.normal(loc=3, scale=1, size=7000)



# Create a figure and subplots
fig, ax = plt.subplots()

# Plot box plots 
bp_annual = ax.boxplot([annual_2000, annual_2100], positions=[1, 2], widths=0.6, patch_artist=True, showfliers=False)
bp_winter = ax.boxplot([winter_2000, winter_2100], positions=[3, 4], widths=0.6, patch_artist=True, showfliers=False)
bp_spring = ax.boxplot([spring_2000, spring_2100], positions=[5, 6], widths=0.6, patch_artist=True, showfliers=False)
bp_summer = ax.boxplot([summer_2000, summer_2100], positions=[7, 8], widths=0.6, patch_artist=True, showfliers=False)
bp_autumn = ax.boxplot([autumn_2000, autumn_2100], positions=[9, 10], widths=0.6, patch_artist=True, showfliers=False)


# Set colors for the boxes
colors = ['lightblue', 'lightgreen']
for box, color in zip(bp_annual['boxes'], colors):
    box.set(facecolor=color)
for box, color in zip(bp_winter['boxes'], colors):
    box.set(facecolor=color)
for box, color in zip(bp_spring['boxes'], colors):
    box.set(facecolor=color)
for box, color in zip(bp_summer['boxes'], colors):
    box.set(facecolor=color)
for box, color in zip(bp_autumn['boxes'], colors):
    box.set(facecolor=color)

# Set the x-axis ticks and labels
ax.set_xticks([1.5, 3.5, 5.5, 7.5, 9.5])
ax.set_xticklabels(['Annual', 'Winter', 'Spring', 'Summer', 'Autumn'])

# Create custom legend
legend_elements = [Patch(facecolor=colors[0], label='2000'),
                   Patch(facecolor=colors[1], label='2100')]
ax.legend(handles=legend_elements, loc='upper right')

# Set the title and y-axis label
plt.title('Offshore Wind CF for North Europe')
#plt.ylabel('Temperature')

# Show the plot
plt.show()





