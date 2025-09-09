import warnings
warnings.filterwarnings("ignore")

import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
import matplotlib.pyplot as plt
import yaml

from _create_figs import create_m0_map, create_m8_map, mark_sources, mark_sources_2, plot_vector, plot_dotted_vector, getIdx, create_cont_map, plot_outflow_and_separation_vectors
from matplotlib.backends.backend_pdf import PdfPages

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
IMAGE_DIRECTORY = config["data_dir"]

# SET OUTPUT
output_folder = "results"
output_pdf = os.path.join(output_folder, "master_reference.pdf")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
cont_files = glob.glob(f'{IMAGE_DIRECTORY}/*/*cont*.fits')
spw39_files = glob.glob(f'{IMAGE_DIRECTORY}/*/*12co*.fits') + glob.glob(f'{IMAGE_DIRECTORY}/*/*spw39*.fits')
outflow_data = pd.read_csv('data/output/outflow_data.csv')
source_info = pd.read_csv("data/output/source_info.csv", index_col='Main')
source_info.index = source_info.index.str.casefold()
if len(cont_files) != len(spw39_files):
    raise ValueError("Number of continuum files does not match number of 12CO/SPW39 files.")

# SCRIPT
with PdfPages(output_pdf) as pdf:
    for cont, spw39 in zip(cont_files, spw39_files):
        figure = plt.figure(figsize=(30,30))
        figure.set_size_inches(17, 22)
        target_name = os.path.dirname(spw39).split('/')[-1]
        target_info = source_info.loc[target_name]

        hdu = fits.open(spw39)[0]
        distance = target_info.iloc[0]['Dis']

        # Continuum Map
        # create figure
        hdu_cont = fits.open(cont)[0]
        fig1 = create_cont_map(hdu_cont, distance=distance, figure=figure, subplot=(2, 2, 1), multiimage=True)
        fig1.set_title(f"{target_name} Continuum")
        fig1.axis_labels.hide_y()
        fig1.colorbar.set_axis_label_text("")
        # add a marker at each source
        mark_sources(fig1, target_info)



        # M8 Map
        # create figure
        fig2 = create_m8_map(hdu, distance=distance, figure=figure, subplot=(2, 2, 2), multiimage=True)
        fig2.set_title(f"{target_name} 12CO M8")
        # add a marker at each source
        mark_sources(fig2, target_info)



        # M0 Maps with vectors
        outflow = outflow_data[~outflow_data['binary_PA'].isna()].loc[outflow_data['field'].apply(str.casefold) == target_name].groupby('field').first().reset_index()
        if len(outflow) == 0:
            print(f"No outflow data for {target_name}. Skipping M0 Maps with vectors.")
            # SAVE
            pdf.savefig(figure, transparent=True)
            figure.clf()
            continue
        channels = getIdx([outflow.at[0, 'red_channels'], outflow.at[0, 'blue_channels']])
        center_origin = np.array([np.mean([outflow['source_a_ra'], outflow['source_b_ra']]), np.mean([outflow['source_a_dec'], outflow['source_b_dec']])])
        separation_angle_north = outflow['binary_PA']


        # M0 Map with average vector
        fig3 = create_m0_map(hdu, channels, sigma=3, distance=distance, figure=figure, subplot=(2, 2, 3), multiimage=True)
        fig3.axis_labels.hide_y()
        fig3.colorbar.set_axis_label_text("")
        fig3.set_title(f"{target_name} 12CO M0")
        # draw separation vector in both directions
        plot_vector(fig3, center_origin, separation_angle_north, color='white', length=0.005)
        plot_vector(fig3, center_origin, separation_angle_north + 180, color='white', length=0.005)

        # M0 Map with both vectors
        fig4 = create_m0_map(hdu, channels, sigma=3, distance=distance, figure=figure, subplot=(2, 2, 4), multiimage=True)
        fig4.set_title(f"{target_name} 12CO M0")
        # draw separation vector in both directions
        plot_vector(fig4, center_origin, separation_angle_north, color='white', length=0.005)
        plot_vector(fig4, center_origin, separation_angle_north + 180, color='white', length=0.005)
        
        # plot each outflow vector
        for j, source in outflow_data[outflow_data['field'].apply(str.casefold) == target_name].reset_index(drop=True).iterrows():
            # define vector origin at the outflow source
            if source['outflow_source'] == 'both':
                outflow_origin = center_origin
            elif source['outflow_source'] == source['source_a']:
                outflow_origin = np.array([source['source_a_ra'], source['source_a_dec']])
            else:
                outflow_origin = np.array([source['source_b_ra'], source['source_b_dec']])

            # plot average vector on fig3
            outflow_angle_north = source['outflow_PA']
            plot_vector(fig3, outflow_origin, outflow_angle_north, color='red', length=0.005)

            # plot double vectors on fig4
            if ~pd.isna(source['red_outflow_PA']):
                plot_vector(fig4, outflow_origin, source['red_outflow_PA'], color='red', length=0.005)
            if ~pd.isna(source['blue_outflow_PA']):
                plot_vector(fig4, outflow_origin, 180 + source['blue_outflow_PA'], color="#c300ff", length=0.005)
        

        # add a marker at each source with legend
        mark_sources(fig3, target_info)
        mark_sources(fig4, target_info)

        plt.subplots_adjust(wspace=0.4, hspace=0.1)

        # SAVE
        pdf.savefig(figure, transparent=True)
        figure.clf()
