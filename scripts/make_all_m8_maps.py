import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
from matplotlib.backends.backend_pdf import PdfPages
from create_figs import create_m8_map, mark_sources


### SCRIPT OPTIONS

# output
output_folder = "results"
output_pdf = os.path.join(output_folder, "all_m8_maps.pdf")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# this script assumes that the directory below contains a folder for each field
# and that within each field folder there is a 12CO image, which contains
# '12co' or 'spw39' in the filename
image_directory = "/Volumes/Alpha/Research/data/"


### SCRIPT

# read data
df = pd.read_csv('data/output/outflow_data.csv')

source_info = pd.read_csv("data/output/source_info.csv")
source_info['Main'] = source_info['Main'].apply(lambda x: str(x).casefold())
source_info.set_index('Main', inplace=True)


with PdfPages(output_pdf) as pdf:
    # get files
    files = glob.glob(f'{image_directory}/*/*12co*.fits') + glob.glob(f'{image_directory}/*/*spw39*.fits')

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
