import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

output_folder = "results"

# create histogram
df = pd.read_csv('data/output/outflow_data.csv')
angle = np.abs(df['outflow_angle'] - df['separation_angle'])
angle = np.min([angle, 180 - angle], axis=0)
angle = np.abs(angle)
angle

fig2, ax = plt.subplots()
ax.hist(angle, label=f"N = {len(angle)}")
ax.legend(loc='upper left')
ax.set_xlabel("smallest angle between binary separation and outflow")
ax.set_ylabel("count")
histogram_filename = "histogram.pdf"
histogram_path = os.path.join(output_folder, histogram_filename)
fig2.savefig(histogram_path)