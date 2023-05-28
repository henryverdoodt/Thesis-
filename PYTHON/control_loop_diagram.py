import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create figure and axes
fig, ax = plt.subplots()

# Define block positions
input_block = patches.Rectangle((0, 0), 1, 1, edgecolor='black', facecolor='white')
process_block = patches.Rectangle((2, 0), 1, 1, edgecolor='black', facecolor='white')
output_block = patches.Rectangle((4, 0), 1, 1, edgecolor='black', facecolor='white')

# Add blocks to the plot
ax.add_patch(input_block)
ax.add_patch(process_block)
ax.add_patch(output_block)

# Add block labels
ax.text(0.5, 0.5, 'Input', horizontalalignment='center', verticalalignment='center')
ax.text(2.5, 0.5, 'Process', horizontalalignment='center', verticalalignment='center')
ax.text(4.5, 0.5, 'Output', horizontalalignment='center', verticalalignment='center')

# Add feedback arrow
arrow = patches.FancyArrowPatch((3, 0.5), (1, 0.5), arrowstyle='->', mutation_scale=10)
ax.add_patch(arrow)

# Add text annotations on the feedback arrow
ax.text(1.7, 0.5, 'Feedback', horizontalalignment='right', verticalalignment='center', rotation=90)
ax.text(2.3, 0.5, 'Block', horizontalalignment='left', verticalalignment='center', rotation=90)

# Set x and y limits
ax.set_xlim(0, 6)
ax.set_ylim(0, 2)

# Remove axis ticks and labels
ax.axis('off')

# Show the plot
plt.show()
