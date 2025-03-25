import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

output_folder = "results"

# create histogram
df = pd.read_csv('data/output/outflow_data.csv')

fig2, ax = plt.subplots()
ax.hist(df['delta_PA'], label=f"N = {len(df)}")
ax.legend(loc='upper left')
ax.set_xlabel("smallest angle between binary separation and outflow")
ax.set_ylabel("count")
histogram_filename = "histogram.pdf"
histogram_path = os.path.join(output_folder, histogram_filename)
fig2.savefig(histogram_path)