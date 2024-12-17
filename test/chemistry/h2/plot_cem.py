import matplotlib.pyplot as plt
import numpy as np

# Use LaTeX for high-quality typesetting
plt.rcParams.update({
    'text.usetex': True,
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'lines.linewidth': 2.5,
    'savefig.dpi': 600,  # Ensure high resolution
    'legend.fontsize': 12,
    'axes.formatter.useoffset': False,  # Avoid offsets
    'axes.formatter.limits': [-3, 4]  # Use scientific notation outside of these limits
})

# File paths
cem_file = 'postProcessing/probe/0/cem'
temp_file = 'postProcessing/probe/0/T'

# Load data using numpy (efficient and clean)
time_cem, cem = np.loadtxt(cem_file, comments='#', unpack=True)
time_temp, temp = np.loadtxt(temp_file, comments='#', unpack=True)

# Create a figure with adjusted aspect ratio
fig, ax1 = plt.subplots(figsize=(7, 5))  # Better for publication column size

# Plot CEM data on the left axis
color_cem = 'tab:blue'
ax1.set_xlabel(r'\textbf{Time [s]}')
ax1.set_ylabel(r'\textbf{CEM}', color=color_cem)
ax1.plot(time_cem, cem, color=color_cem, linestyle='-', label='CEM')
ax1.tick_params(axis='y', labelcolor=color_cem)

# Add a second Y-axis for Temperature
ax2 = ax1.twinx()
color_temp = 'tab:red'
ax2.set_ylabel(r'\textbf{Temperature} $T \, [\mathrm{K}]$', color=color_temp)
ax2.plot(time_temp, temp, color=color_temp, linestyle='-', label='Temperature')
ax2.tick_params(axis='y', labelcolor=color_temp)

# Add combined legend outside the plot
fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, frameon=False)

# Tweak layout for publication clarity
fig.tight_layout()

# Export the figure as vector graphics (PDF or SVG)
plt.savefig('cem.pdf', format='pdf', bbox_inches='tight')
