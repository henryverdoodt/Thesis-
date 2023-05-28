import matplotlib.pyplot as plt
import numpy as np

# Cost data for each period
periods = ['Period 1', 'Period 2', 'Period 3', 'Period 4', 'Period 5', 'Period 6', 'Period 7']
#periods = ['', '', '', '', '', '', '']
costs = [0, 50, 170, 97, 100, 60, 60]

# Calculate cumulative cost
cumulative_costs = np.cumsum(costs)

# Create the figure and axes
fig, ax = plt.subplots()

# Plot the bar chart for individual periods
ax.bar(periods, costs, align='center', alpha=0.5, label='Montly Cost')

# Plot the line chart for cumulative cost
ax.plot(periods, cumulative_costs, marker='o', color='blue', linestyle='-', linewidth=2, label='Planned Cumulative Cost')

# Set y-axis label
ax.set_ylabel('Cost')

# Set title and legend
ax.set_title('Cost Chart')
ax.legend()

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Display the plot
plt.show()


'''

import matplotlib.pyplot as plt
import numpy as np

# Cost data for each period
periods = ['Period 1', 'Period 2', 'Period 3', 'Period 4', 'Period 5', 'Period 6', 'Period 7']
planned_costs = [100, 150, 120, 200, 180, 160, 220]
actual_costs = [None] * len(periods)

# Calculate total costs
worker_hours_costs = [80, 120, 90, 150, 130, 110, 180]
administrative_costs = [10, 15, 12, 20, 18, 16, 22]
material_costs = [10, 15, 18, 30, 32, 34, 18]
total_costs = [wh + ad + mat for wh, ad, mat in zip(worker_hours_costs, administrative_costs, material_costs)]

# Create the figure and axes
fig, ax = plt.subplots()

# Plot the bar chart for individual periods
ax.bar(periods, planned_costs, align='center', alpha=0.5, label='Planned Cost')

# Plot the line chart for cumulative cost
ax.plot(periods, np.cumsum(planned_costs), marker='o', color='blue', linestyle='-', linewidth=2, label='Cumulative Cost')

# Set y-axis label
ax.set_ylabel('Cost')

# Set title and legend
ax.set_title('Cost Chart')
ax.legend()

# Create the table
table_data = [
    ['Worker-hours', None] + worker_hours_costs,
    ['Administrative', None] + administrative_costs,
    ['Material', None] + material_costs,
    ['Total', None] + total_costs
]
table_columns = ['Cost Type', 'Planned'] + periods

table = plt.table(cellText=table_data, colLabels=table_columns, cellLoc='center', loc='bottom')

# Adjust table position and size
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5)
table.set_fontsize(10)

# Remove table borders
for key, cell in table.get_celld().items():
    cell.set_linewidth(0)

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Display the plot
plt.show()
'''