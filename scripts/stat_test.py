# Code From John Tobin
# Modified by Ryan Sponzilli

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import os
from labellines import labelLine

# SET OUTPUT
output_folder = "results/stat_test"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
outflow_data = pd.read_csv('data/output/outflow_data.csv')
outflow_data = outflow_data.loc[~outflow_data['delta_PA'].isna()]

# Generate random position angles and project them
num_samples = 100000
pa = np.random.uniform(0, 90, num_samples)                          # binary PA
ofa = np.random.uniform(0, 90, num_samples)                         # random outflow PA offset
inc = np.random.random_sample(num_samples)                          # inclination
R = np.random.uniform(30.0, 300.0, num_samples)                     # binary separation distance
projx = R * np.cos(np.radians(pa))                                  # project binary PA to 2D plane
projy = R * np.sin(np.radians(pa)) * inc                            # project binary PA to 2D plane
projpa_aligned = 90.0 - np.degrees(np.arctan2(projy, projx))        # expected orthogonal outflow PA
projpa_random = np.abs(ofa - np.degrees(np.arctan2(projy, projx)))  # random outflow PA

def makeCumulate(arrayData):
    # sort ascending
    sorted_pas = np.sort(arrayData, axis=None)
    # get length
    number = len(sorted_pas)
    # Compute cumulative fraction (normalized by total number of elements)
    frac_cum = (np.arange(0, number, dtype=np.float32)) / float(number)
    return sorted_pas, frac_cum

def make_mixed_cumulates(pct_rand_list, include_obs=False):
    cumulates = []
    ks_pvalues = []

    # Generate cumulative distributions for observed data
    c0 = makeCumulate(outflow_data['delta_PA'][:])
    if include_obs:
        cumulates.append(c0)
        ks_pvalues.append({"% random": -1, "p": -1})

    for i in pct_rand_list:
        # for i = 1, add 99% from projpa and 1% from projpa_rand
        index = 1000*(100-i)
        # create different percentage-randomized datasets
        d = np.concatenate((projpa_aligned[:index], projpa_random[index:]))
        # create cumulative distribution
        c = makeCumulate(d)
        cumulates.append(c)
        # perform KS test
        ks_pvalues.append({"% random": i, "p": stats.ks_2samp(c0[0], c[0]).pvalue})
    return pd.DataFrame(ks_pvalues), cumulates

# generate cumulative distributions and p-values for simulated data
data1, _ = make_mixed_cumulates(range(0, 101, 1), include_obs=False)
plot_data, cumulates = make_mixed_cumulates([0, 25, 50, 75, 100], include_obs=True)

# Generate cumulative plot
fig, ax = plt.subplots(figsize=(11,8))
# plot each distribution
colors = ['#4477AA', '#AA3377', '#66CCEE', '#CCBB44', '#228833', '#EE6677']
line_labels = ['Observations', '100% Orthogonal', '75% Orthogonal', '50% Orthogonal', '25% Orthogonal', 'Random']
lines = []
for i, item in enumerate(cumulates):
    line, = ax.step(item[0], item[1], color=colors[i], linewidth=2.0)
    lines.append(line)
# Plot Settings
ax.set_title('$\\Delta$PA Cumulative Frequency Distribution', fontsize=28, pad=10)
ax.set_xlabel('$\\Delta$PA - smallest angle between binary separation and outflow (degrees)', fontsize=20, labelpad=15)
ax.set_ylabel('frequency', fontsize=20, labelpad=15)
ax.set_xlim(0, 90.0)
ax.set_ylim(0, 1.0)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
labelLine(lines[0], label=line_labels[0], x=16, align=False, rotation=0, yoffset=-0.018, fontsize=14)
labelLine(lines[1], label=line_labels[1], x=35, align=False, rotation=24, yoffset=0.03, fontsize=14)
labelLine(lines[2], label=line_labels[2], x=40, align=False, rotation=28, yoffset=0.03, fontsize=14)
labelLine(lines[3], label=line_labels[3], x=45, align=False, rotation=30, yoffset=0.035, fontsize=14)
labelLine(lines[4], label=line_labels[4], x=50, align=False, rotation=32, yoffset=0.03, fontsize=14)
labelLine(lines[5], label=line_labels[5], x=55, align=False, rotation=30, yoffset=0.03, fontsize=14)
legend_labels = [f"p={round(p, 5)}" if (p != -1) else "Observations" for p in plot_data['p']]
ax.legend(lines, legend_labels, fontsize=13)
plt.savefig(os.path.join(output_folder, 'DeltaPA_cumulat_deg.pdf'), dpi=200)

# p-value plot
fig2, ax2 = plt.subplots(figsize=(11,8))
ax2.plot(100 - data1['% random'], data1['p'], color=colors[0])
less_than_01 = 100 - data1.loc[data1['p']<0.1].iloc[0]['% random']
less_than_005 = 100 - data1.loc[data1['p']<0.05].iloc[0]['% random']
less_than_001 = 100 - data1.loc[data1['p']<0.01].iloc[0]['% random']
plt.axhline(0.1, ls='dashed', label=f'% Orthogonal = {less_than_01}', color=colors[1])
plt.axhline(0.05, ls='dashdot', label=f'% Orthogonal = {less_than_005}', color=colors[4])
plt.axhline(0.01, ls='dotted', label=f'% Orthogonal = {less_than_001}', color=colors[5])

ax2.set_title("% Orthogonal Outflows vs p-value", fontsize=28, pad=10)
ax2.set_xlabel("% Orthogonal Outflows", fontsize=20, labelpad=15)
ax2.set_ylabel('p-value', fontsize=20, labelpad=15)
ax2.set_xlim(0,100)
ax2.set_xticklabels(ax2.get_xticklabels(), fontsize=16)
ax2.set_yticklabels(ax2.get_yticklabels(), fontsize=16)
ax2.legend(fontsize=13)
plt.savefig(os.path.join(output_folder, 'p-values.pdf'), dpi=200)

# Delta PA Histrogram (nothing special)
fig, ax = plt.subplots(figsize=(11,8))
ax.hist(outflow_data['delta_PA'], label=f"N = {len(outflow_data[~outflow_data['delta_PA'].isna()])}", color='#1D58A7')
ax.legend(loc='upper left', fontsize=32)
ax.set_title("$\Delta$PA Distribution", fontsize=28, pad=10)
ax.set_xlabel("$\Delta$PA - smallest angle between binary separation and outflow (degrees)", fontsize=20, labelpad=15)
ax.set_ylabel("count", fontsize=20, labelpad=15)

ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
fig.savefig(os.path.join(output_folder, "histogram.pdf"))