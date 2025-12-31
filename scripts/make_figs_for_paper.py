import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import glob
import yaml
from astropy.coordinates import SkyCoord

import sys
from _create_figs import create_m0_map, create_m8_map, mark_sources, getIdx, plot_outflow_and_separation_vectors, create_cont_map

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
IMAGE_DIRECTORY = config["data_dir"]

# SET OUTPUT
output_folder = "results/figs_for_paper"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
outflow_data = pd.read_csv('data/output/outflow_data.csv')
outflow_data2 = outflow_data.groupby('field').first().reset_index()
by_field = pd.read_csv('data/output/data_by_field.csv')
source_info = pd.read_csv("data/output/source_info.csv", index_col='Main')
source_info.index = source_info.index.str.casefold()

cont_files = glob.glob(f'{IMAGE_DIRECTORY}/*/*cont*.fits')
spw39_files = glob.glob(f'{IMAGE_DIRECTORY}/*/*12co*.fits') + glob.glob(f'{IMAGE_DIRECTORY}/*/*spw39*.fits')
if len(cont_files) != len(spw39_files):
    raise ValueError("Number of continuum files does not match number of 12CO/SPW39 files.")

# main function
def make_compound_plots(field_list, output_name, rows, cols, figsize=None, v_pad_frac=0.15):
    if figsize == None:
        figure = plt.figure(figsize=(17/3*cols, 22/5*rows))
    else:
        figure = plt.figure(figsize=figsize)

    # This loops through each source field with an angle measurement and creates
    # a figure with the separation vector and outflow vector overlayed

    data = by_field[by_field['field'].isin(field_list)]

    subplot_num = 1
    
    for i, field in data.sort_values(by='ra', ascending=True).reset_index().iterrows():

        ### READ DATA
        # verify output path exists and
        # skip already existing files if you don't want to overwrite them
        target_name = field['field']

        # verify output path exists and
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        # open images
        spw39_file = spw39_files[next((i for i, s in enumerate(spw39_files) if target_name.casefold() in s), None)]
        cont_file = cont_files[next((i for i, s in enumerate(cont_files) if target_name.casefold() in s), None)]
        hdu_cont = fits.open(cont_file)[0]
        hdu_spw39 = fits.open(spw39_file)[0]
        distance = source_info.loc[target_name.casefold(), 'Dis'].iloc[0]
        target_info = source_info.loc[target_name.casefold()]

        ### MOMENT MAP PANEL
        # M8 Map
        if pd.isna(field['integrated_channels']):
            # Create Figure
            fig = create_m8_map(hdu_spw39, distance=distance, figure=figure, subplot=(rows, cols, subplot_num), multiimage=True, use_offset_labels=True)
            fig.set_title(f"{target_name} M8")   
        
        # M0 Map
        else:
            # Create Figure
            outflow = outflow_data.loc[outflow_data['field'] == target_name].groupby('field').first().reset_index()
            channels = getIdx([outflow.at[0, 'red_channels'], outflow.at[0, 'blue_channels']])
            fig = create_m0_map(hdu_spw39, channels, sigma=3, distance=distance, figure=figure, subplot=(rows, cols, subplot_num), multiimage=True, use_offset_labels=True)
            fig.set_title(f"{target_name} M0")

            # Plot Vectors
            plot_outflow_and_separation_vectors(fig, outflow_data, target_name)

        # Mark Sources
        mark_sources(fig, target_info, fontsize=8, use_short_label=True)

        # Manage Axis Labels
        
        if (subplot_num-1) % cols != cols - 1:
            if not str(output_name).startswith("appendix"):
                fig.colorbar.set_axis_label_text("")
        if (subplot_num-1) % cols != 0:
            fig.axis_labels.hide_y()
        if (rows > 1) & ((subplot_num-1) < (rows*cols - cols)):
            fig.axis_labels.hide_x()

        ### CONTINUUM PANEL
        subplot_num += cols

        distance = target_info.iloc[0]['Dis']
        size = np.array([3, 3]) * u.arcsecond
        center = SkyCoord(target_info[target_info['important']==True]['RA'].mean(), target_info[target_info['important']==True]['Dec'].mean(), unit=u.degree)
        fig1 = create_cont_map(hdu_cont, size=size, center=center, distance=distance, scalebar_au=100, figure=figure, subplot=(rows, cols, subplot_num), multiimage=True, use_offset_labels=True)
        fig1.set_title("Continuum")

        # mark sources
        mark_sources(fig1, target_info[target_info['important']==True], fontsize=8, use_short_label=True)

        # Manage Axis Labels
        if (subplot_num-1) % cols != 0:
            fig1.axis_labels.hide_y()
        if (rows > 1) & ((subplot_num-1) < (rows*cols - cols)):
            fig1.axis_labels.hide_x()

        subplot_num -= (cols - 1)
        if (i+1) % cols == 0:
            subplot_num += cols

    ### SAVE IMAGE
    # plt.subplots_adjust(wspace=wspace, hspace=hspace)
    figure.savefig(os.path.join(output_folder, f"{output_name}.pdf"), bbox_inches='tight')
    figure.show()


# SCRIPT
# define fields in each figure
fig_1 = ['HH270VLA1', 'HOPS-168']
fig_2 = ['HOPS-290', 'HOPS-288']
no_outflows = ['HOPS-28', 'HOPS-163', 'HOPS-242', 'HOPS-255', 'HOPS-248', 'HOPS-357']
app_1 = by_field[~by_field['field'].isin(no_outflows)]['field'].values.tolist()
app_2 = by_field[by_field['field'].isin(no_outflows)]['field'].values.tolist()
app_full = app_1 + app_2

# generate all figures
make_compound_plots(fig_1, "fig-1", rows=2, cols=2, figsize=(17/3*2, 22/5*2.25))
make_compound_plots(fig_2, "fig-2", rows=2, cols=2, figsize=(17/3*2, 22/5*2.25))
make_compound_plots(app_full[0:4], "appendix-1", rows=4, cols=2)
make_compound_plots(app_full[4:8], "appendix-2", rows=4, cols=2)
make_compound_plots(app_full[8:12], "appendix-3", rows=4, cols=2)
make_compound_plots(app_full[12:16], "appendix-4", rows=4, cols=2)
make_compound_plots(app_full[16:20], "appendix-5", rows=4, cols=2)
make_compound_plots(app_full[20:24], "appendix-6", rows=4, cols=2)
make_compound_plots(app_full[24:28], "appendix-7", rows=4, cols=2)
make_compound_plots(app_full[28:32], "appendix-8", rows=4, cols=2)
make_compound_plots(app_full[32:36], "appendix-9", rows=4, cols=2)
make_compound_plots(app_full[36:40], "appendix-10", rows=4, cols=2)
make_compound_plots(app_full[40:44], "appendix-11", rows=4, cols=2)
make_compound_plots(app_full[44:48], "appendix-12", rows=4, cols=2)
make_compound_plots(app_full[48::], "appendix-13", rows=4, cols=2)