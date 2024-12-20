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
cem_file = 'h2/postProcessing/probe/0/cem'
temp_file = 'h2/postProcessing/probe/0/T'

# Load data using numpy (efficient and clean)
time_cem, cem = np.loadtxt(cem_file, comments='#', unpack=True)
time_temp, temp = np.loadtxt(temp_file, comments='#', unpack=True)

# TODO Skip the first row (omit index 0)
time_cem, cem = time_cem[1:], cem[1:]
time_temp, temp = time_temp[1:], temp[1:]

# Create a figure with adjusted aspect ratio
fig, ax1 = plt.subplots(figsize=(7, 5))  # Better for publication column size

# Plot maximum eigen value data on the left axis
color_cem = 'tab:blue'
ax1.set_xlabel(r'$t \, [\mathrm{s}]$')
ax1.set_ylabel(r'$\lambda^+$', color=color_cem)
ax1.plot(time_cem, cem, color=color_cem, linestyle='-', label=r'$\lambda^+$')
ax1.tick_params(axis='y', labelcolor=color_cem)

# TODO Set y-limits for the left y-axis
# ax1.set_ylim(-.1, .1)  # Adjust as needed for $\lambda^+$ values

# Add a second Y-axis for Temperature
ax2 = ax1.twinx()
color_temp = 'tab:red'
ax2.set_ylabel(r' $T \, [\mathrm{K}]$', color=color_temp)
ax2.plot(time_temp, temp, color=color_temp, linestyle='-', label='Temperature')
ax2.tick_params(axis='y', labelcolor=color_temp)

# Add combined legend outside the plot
fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, frameon=False)

# Tweak layout for publication clarity
fig.tight_layout()

# Export the figure as vector graphics (PDF or SVG)
plt.savefig('cem.pdf', format='pdf', bbox_inches='tight')
