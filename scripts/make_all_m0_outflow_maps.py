import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
import yaml
from matplotlib.backends.backend_pdf import PdfPages

from _create_figs import create_m0_map, mark_sources, plot_vector, getIdx, plot_outflow_and_separation_vectors

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
IMAGE_DIRECTORY = config["data_dir"]


# SET OUTPUT
output_folder = "results"
output_pdf = os.path.join(output_folder, "m0_outflow_maps.pdf")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
outflow_data = pd.read_csv('data/output/outflow_data.csv')
source_info = pd.read_csv("data/output/source_info.csv", index_col='Main')
source_info.index = source_info.index.str.casefold()

# SCRIPT
with PdfPages(output_pdf) as pdf:
    for i, field in outflow_data[~outflow_data['binary_PA'].isna()].groupby('field').agg('first').sort_values('source_a_ra').reset_index().iterrows():
        target_name = field['field']

        image_filename = (glob.glob(f'{IMAGE_DIRECTORY}{target_name.casefold()}/*12co*.fits') +
                          glob.glob(f'{IMAGE_DIRECTORY}{target_name.casefold()}/*spw39*.fits'))[0]
        hdu = fits.open(image_filename)[0]
        distance = field['distance']

        channels = getIdx([field['red_channels'], field['blue_channels']])
        fig = create_m0_map(hdu, channels, sigma=3, distance=distance)
        # Set figure size to letter (8.5 x 11 inches)
        fig._figure.set_size_inches(8.5, 11)
        fig.set_title(f"{target_name} 12CO M0")

        # Plot Vectors
        plot_outflow_and_separation_vectors(fig, outflow_data, target_name, delta_pa_label=True)

        target_info = source_info.loc[target_name.casefold()]
        mark_sources(fig, target_info)

        # SAVE
        pdf.savefig(fig._figure, transparent=True)
        fig._figure.clf()