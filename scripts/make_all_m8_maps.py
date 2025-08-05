import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from _create_figs import create_m8_map, mark_sources

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
IMAGE_DIRECTORY = config["data_dir"]


# SET OUTPUT
output_folder = "results"
output_pdf = os.path.join(output_folder, "all_m8_maps.pdf")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
source_info = pd.read_csv("data/output/source_info.csv", index_col='Main')
source_info.index = source_info.index.str.casefold()

# SCRIPT
with PdfPages(output_pdf) as pdf:
    # get files
    files = glob.glob(f'{IMAGE_DIRECTORY}/*/*12co*.fits') + glob.glob(f'{IMAGE_DIRECTORY}/*/*spw39*.fits')

    for file in files:
        # set output name and output path
        target_name = os.path.dirname(file).split('/')[-1]

        try:
            target_info = source_info.loc[target_name]

            hdulist = fits.open(file)
            hdu = hdulist[0]

            # set center and size of cutout
            center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
            size = np.array([39, 39]) * u.arcsecond
            distance = target_info.iloc[0]['Dis']

            # create figure
            fig = create_m8_map(hdu, center, size, distance)
            fig.set_title(f"{target_name} 12CO M8")

            # add a marker at each source
            mark_sources(fig, target_info)

        except Exception as error:
            print(f"Error for {target_name}, {error}")

        # SAVE
        pdf.savefig(fig._figure, transparent=True)
        fig._figure.clf()
