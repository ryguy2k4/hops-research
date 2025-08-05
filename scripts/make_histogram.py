import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


# SET OUTPUT
output_folder = "results/histogram"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
outflow_data = pd.read_csv('data/output/outflow_data.csv')

# SCRIPT
fig, ax = plt.subplots(figsize=(11,8))
ax.hist(outflow_data['delta_PA'], label=f"N = {len(outflow_data[~outflow_data['delta_PA'].isna()])}", color='#1D58A7')
ax.legend(loc='upper left', fontsize=32)
ax.set_title("$\Delta$PA Distribution", fontsize=24)
ax.set_xlabel("$\Delta$PA - smallest angle between binary separation and outflow (degrees)", fontsize=16)
ax.set_ylabel("count", fontsize=16)
fig.savefig(os.path.join(output_folder, "histogram.pdf"))
# fig.savefig(os.path.join(output_folder, "histogram.png"), dpi=300, transparent=True, bbox_inches='tight')
