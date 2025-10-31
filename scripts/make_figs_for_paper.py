import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import glob
import yaml

import sys
from _create_figs import create_m0_map, create_m8_map, mark_sources, getIdx, plot_outflow_and_separation_vectors, create_cont_map_2

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

# this tedius function was made by chatgpt
# it takes care of positioning the continuum and moment map panels correctly
def get_subplot_params(
    subplot_num, rows=4, cols=2,
    main_width=0.45, mid_width_pad=0.1,
    top_pad=0.03, bottom_pad=0.03, v_pad_frac=0.15
):
    """
    v_pad_frac: vertical padding as a fraction of main subplot height.
    """
    if not (1 <= subplot_num <= rows * cols):
        raise ValueError("subplot_num out of range.")

    # --- Horizontal sizing ---
    available_for_mains = 1.0 - mid_width_pad * (cols - 1)
    main_w = min(main_width, available_for_mains / cols)
    small_w = main_w / 2.0

    total_main_block = cols * main_w + (cols - 1) * mid_width_pad
    left_margin = max(0.0, (1.0 - total_main_block) / 2.0)
    left_main_cols = [left_margin + i * (main_w + mid_width_pad) for i in range(cols)]

    # --- Vertical sizing ---
    num_main_rows = (rows + 1) // 2
    num_small_rows = rows // 2
    total_units = num_main_rows * 1.0 + num_small_rows * 0.5
    main_h = (1.0 - (top_pad + bottom_pad)) / (total_units + (rows - 1) * v_pad_frac)
    small_h = main_h / 2.0
    v_pad = v_pad_frac * main_h  # scale padding with subplot size

    bottoms = []
    current_top = 1.0 - top_pad
    for r in range(rows):
        if r > 0:
            current_top -= v_pad
        h = main_h if (r % 2 == 0) else small_h
        bottom = current_top - h
        bottoms.append(bottom)
        current_top = bottom

    # --- Pick subplot ---
    idx = subplot_num - 1
    row_idx = idx // cols
    col_idx = idx % cols
    is_main = (row_idx % 2 == 0)

    if is_main:
        width, height = main_w, main_h
        left = left_main_cols[col_idx]
    else:
        width, height = small_w, small_h
        left = left_main_cols[col_idx] + (main_w - small_w) / 2.0

    bottom = bottoms[row_idx]
    return [left, bottom, width, height]

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
            fig = create_m8_map(hdu_spw39, distance=distance, figure=figure, subplot=get_subplot_params(subplot_num, rows=rows, cols=cols, v_pad_frac=v_pad_frac), multiimage=True, use_offset_labels=True)
            fig.set_title(f"{target_name} M8")   
        
        # M0 Map
        else:
            # Create Figure
            outflow = outflow_data.loc[outflow_data['field'] == target_name].groupby('field').first().reset_index()
            channels = getIdx([outflow.at[0, 'red_channels'], outflow.at[0, 'blue_channels']])
            fig = create_m0_map(hdu_spw39, channels, sigma=3, distance=distance, figure=figure, subplot=get_subplot_params(subplot_num, rows=rows, cols=cols, v_pad_frac=v_pad_frac), multiimage=True, use_offset_labels=True)
            fig.set_title(f"{target_name} M0")

            # Plot Vectors
            plot_outflow_and_separation_vectors(fig, outflow_data, target_name)

        # Mark Sources
        mark_sources(fig, target_info, fontsize=8, use_short_label=True)

        # Manage Axis Labels
        if not str(output_name).startswith("appendix"):
            if (subplot_num-1) % cols != cols - 1:
                fig.colorbar.set_axis_label_text("")
        if (subplot_num-1) % cols != 0:
            fig.axis_labels.hide_y()
        if (rows > 1) & ((subplot_num-1) < (rows*cols - 2*cols)):
            fig.axis_labels.hide_x()

        ### CONTINUUM PANEL
        subplot_num += cols

        distance = target_info.iloc[0]['Dis']
        size = np.array([1075/distance, 1000/distance]) * u.arcsecond
        fig1 = create_cont_map_2(hdu_cont, size=size, distance=distance, scalebar_au=100, figure=figure, subplot=get_subplot_params(subplot_num, rows=rows, cols=cols, v_pad_frac=v_pad_frac), multiimage=True, use_offset_labels=True)
        fig1.set_title("Continuum")

        # mark sources
        mark_sources(fig1, target_info, fontsize=8, use_short_label=True)

        # Manage Axis Labels
        if (subplot_num-1) % cols != cols - 1:
            fig1.colorbar.set_axis_label_text("")

        subplot_num -= (cols - 1)
        if (i+1) % cols == 0:
            subplot_num += cols

    ### SAVE IMAGE
    # plt.subplots_adjust(wspace=wspace, hspace=hspace)
    figure.savefig(os.path.join(output_folder, f"{output_name}.pdf"), bbox_inches='tight')
    figure.show()


# SCRIPT
# define fields in each figure
fig_1 = ['HOPS-32', 'HOPS-168', 'HOPS-281', 'Per-emb-17']
fig_2 = ['HOPS-290', 'HOPS-288']
no_outflows = ['HOPS-28', 'HOPS-163', 'HOPS-242', 'HOPS-255', 'HOPS-248', 'HOPS-357']
app_1 = by_field[~by_field['field'].isin(no_outflows)]['field'].values
app_2 = by_field[by_field['field'].isin(no_outflows)]['field'].values

# generate all figures
make_compound_plots(fig_1, "fig_1", rows=2, cols=4, figsize=(17/3*3.8, 22/5*1.7))
make_compound_plots(fig_2, "fig_2", rows=2, cols=2, figsize=(17/3*1.6, 22/5*1.7))
make_compound_plots(app_1[0:9], "appendix-1", rows=6, cols=3)
make_compound_plots(app_1[9:18], "appendix-2", rows=6, cols=3)
make_compound_plots(app_1[18:27], "appendix-3", rows=6, cols=3)
make_compound_plots(app_1[27:36], "appendix-4", rows=6, cols=3)
make_compound_plots(app_1[36:45], "appendix-5", rows=6, cols=3)
make_compound_plots(app_2, "appendix-6", rows=4, cols=3)