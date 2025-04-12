# Code From John Tobin
# Modified by Ryan Sponzilli

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import os

def createCumulatPlotCos(datax, datay, datalabel='', title='', xlabel='', filename='test'):
    fig, ax = plt.subplots(figsize=(11,8))
    
    # get labels
    if datalabel == '':
        datalabel = [''] * len(datax)
    
    # Loop through each dataset and plot it as a step function
    colors = ['#1D58A7', '#FCB316', '#006230', '#007E8E', '#5C0E41', '#7D3E13']
    for i in range(len(datax)):
        ax.step(datax[i], datay[i], label=datalabel[i], color=colors[i])
    
    # PLOT
    ax.set_title(title, fontsize=24)
    ax.set_ylabel('Frequency', fontsize=18)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)

    # legend
    ax.legend(fontsize=10)

    plt.savefig(filename + '.pdf', dpi=200)
    plt.savefig(filename + '.png', dpi=300, transparent=True, bbox_inches='tight')


def createCumulatPlotDeg(datax, datay, datalabel='', title='', xlabel='', filename='test'):
    fig, ax = plt.subplots(figsize=(11,8))
    
    # get labels
    if datalabel == '':
        datalabel = [''] * len(datax)
    
    # Loop through each dataset and plot it as a step function
    colors = ['#1D58A7', '#FCB316', '#006230', '#007E8E', '#5C0E41', '#7D3E13']
    for i in range(len(datax)):
        ax.step(datax[i], datay[i], label=datalabel[i], color=colors[i])
        
    # PLOT
    ax.set_title(title, fontsize=24)
    ax.set_ylabel('frequency', fontsize=16)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_xlim(0, 90.0)
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=10)
    plt.savefig(filename + '.pdf', dpi=200)
    plt.savefig(filename + '.png', dpi=300, transparent=True, bbox_inches='tight')

def makeCumulate(arrayData):
    # sort ascending
    sorted_pas = np.sort(arrayData, axis=None)
    # get length
    number = len(sorted_pas)
    # Compute cumulative fraction (normalized by total number of elements)
    frac_cum = (np.arange(0, number, dtype=np.float32)) / float(number)
    
    return sorted_pas, frac_cum

def perform_test(data, output_path):
    # verify output path exists
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Create histogram plot
    fig, ax1 = plt.subplots(figsize=(11, 8))
    cosi = np.cos(data['delta_PA'][:] * np.pi / 180.0)
    bins = np.linspace(0, 1.0, 11)
    ax1.hist(cosi, bins, fill=True, facecolor='gray', alpha=0.4, histtype='bar', ec='black', linewidth=2.0)
    ax1.set_ylabel('N', fontsize=16)
    ax1.set_xlabel('cos($Delta$[Outflow-Binary])', fontsize=16)
    ax1.set_title('Outflow PA vs. Binary PA', fontsize=24)
    plt.savefig(os.path.join(output_path, 'hist.pdf'))

    # Generate cumulative distributions for observed data
    paCumulat_cos, frac_paCumulat_cos = makeCumulate(cosi)
    paCumulat_deg, frac_paCumulat_deg = makeCumulate(data['delta_PA'][:])

    # Generate random position angles and project them
    num_samples = 100000
    pa = np.random.uniform(0, 90, num_samples) # binary PA
    ofa = np.random.uniform(0, 90, num_samples) # random outflow PA offset
    inc = np.random.random_sample(num_samples) # inclination
    R = np.random.uniform(30.0, 300.0, num_samples) # binary separation distance
    projx = R * np.cos(np.radians(pa)) # project binary PA to 2D plane
    projy = R * np.sin(np.radians(pa)) * inc # project binary PA to 2D plane
    projpa = 90.0 - np.degrees(np.arctan2(projy, projx)) # expected orthogonal outflow PA
    projpa_rand = np.abs(ofa - np.degrees(np.arctan2(projy, projx))) # random outflow PA

    # Create different percentage-randomized datasets
    projpa_25p_rand = np.concatenate((projpa[:75000], projpa_rand[75000:]))
    projpa_50p_rand = np.concatenate((projpa[:50000], projpa_rand[50000:]))
    projpa_75p_rand = np.concatenate((projpa[:25000], projpa_rand[25000:]))

    # Compute cumulative distributions for models
    paCumulat_deg_model, frac_paCumulat_deg_model = makeCumulate(projpa)
    paCumulat_cos_model, frac_paCumulat_cos_model = makeCumulate(np.cos(np.radians(projpa)))
    paCumulat_deg_model_25p_rand, frac_paCumulat_deg_model_25p_rand = makeCumulate(projpa_25p_rand)
    paCumulat_cos_model_25p_rand, frac_paCumulat_cos_model_25p_rand = makeCumulate(np.cos(np.radians(projpa_25p_rand)))
    paCumulat_deg_model_50p_rand, frac_paCumulat_deg_model_50p_rand = makeCumulate(projpa_50p_rand)
    paCumulat_cos_model_50p_rand, frac_paCumulat_cos_model_50p_rand = makeCumulate(np.cos(np.radians(projpa_50p_rand)))
    paCumulat_deg_model_75p_rand, frac_paCumulat_deg_model_75p_rand = makeCumulate(projpa_75p_rand)
    paCumulat_cos_model_75p_rand, frac_paCumulat_cos_model_75p_rand = makeCumulate(np.cos(np.radians(projpa_75p_rand)))
    paCumulat_deg_model_rand, frac_paCumulat_deg_model_rand = makeCumulate(projpa_rand)
    paCumulat_cos_model_rand, frac_paCumulat_cos_model_rand = makeCumulate(np.cos(np.radians(projpa_rand)))

    # Perform statistical tests
    f = open(os.path.join(output_path, 'test_results.txt'), 'w')
    ad_tests = [
        ("AD Test (Degrees)", paCumulat_deg, paCumulat_deg_model),
        ("AD Test 25% random (Degrees)", paCumulat_deg, paCumulat_deg_model_25p_rand),
        ("AD Test 50% random (Degrees)", paCumulat_deg, paCumulat_deg_model_50p_rand),
        ("AD Test 75% random (Degrees)", paCumulat_deg, paCumulat_deg_model_75p_rand),
        ("AD Test (Cos)", paCumulat_cos, paCumulat_cos_model)
    ]
    for label, observed, model in ad_tests:
        f.write(f'{label}: {stats.anderson_ksamp([observed, model])}\n')

    ks_tests = [
        ("KS Test (Degrees)", paCumulat_deg, paCumulat_deg_model),
        ("KS Test 25% random (Degrees)", paCumulat_deg, paCumulat_deg_model_25p_rand),
        ("KS Test 50% random (Degrees)", paCumulat_deg, paCumulat_deg_model_50p_rand),
        ("KS Test 75% random (Degrees)", paCumulat_deg, paCumulat_deg_model_75p_rand),
        ("KS Test (Cos)", paCumulat_cos, paCumulat_cos_model)
    ]
    for label, observed, model in ks_tests:
        f.write(f'{label}: {stats.ks_2samp(observed, model)}\n')
    f.close()

    # Generate cumulative plots
    createCumulatPlotDeg([
        paCumulat_deg, paCumulat_deg_model, paCumulat_deg_model_rand,
        paCumulat_deg_model_25p_rand, paCumulat_deg_model_50p_rand, paCumulat_deg_model_75p_rand
    ], [
        frac_paCumulat_deg, frac_paCumulat_deg_model, frac_paCumulat_deg_model_rand,
        frac_paCumulat_deg_model_25p_rand, frac_paCumulat_deg_model_50p_rand, frac_paCumulat_deg_model_75p_rand
    ], datalabel=[
        'Observations', 'Model - Orthogonal Outflow', 'Model - Random Orientation',
        'Model - 75% Orthogonal Outflow, 25% Random',
        'Model - 50% Orthogonal Outflow, 50% Random',
        'Model - 25% Orthogonal Outflow, 75% Random'
    ], title='$\Delta$PA Cumulative Frequency Distribution', xlabel='$\Delta$PA - smallest angle between binary separation and outflow (degrees)', filename=os.path.join(output_path, 'OutflowPA_cumulat_deg'))

    createCumulatPlotCos([
        paCumulat_cos, paCumulat_cos_model, paCumulat_cos_model_rand
    ], [
        frac_paCumulat_cos, frac_paCumulat_cos_model, frac_paCumulat_cos_model_rand
    ], datalabel=[
        'Observations', 'Model - Orthogonal Outflow', 'Model - Random Orientation'
    ], title='Cumulative Frequency Distribution', xlabel='cos($\Delta$PA)', filename=os.path.join(output_path, 'OutflowPA_cumulat_cos'))

# script options
output_folder = "results/stat_test"
filename = 'OutflowPA_hist'
output_path = os.path.join(output_folder, filename)

# verify output path exists
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# load data
data = pd.read_csv('data/output/outflow_data.csv')
data = data.loc[~data['delta_PA'].isna()]

perform_test(data.loc[data['group'] == 'orion'], os.path.join(output_folder, 'orion'))
perform_test(data.loc[data['group'] == 'perseus'], os.path.join(output_folder, 'perseus'))
perform_test(data, os.path.join(output_folder, 'combined'))