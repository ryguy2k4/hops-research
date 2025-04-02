import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import aplpy
import glob
import os
import pandas as pd

from create_figs import create_m8_map
from create_figs import mark_sources

# script options
overwrite = False
output_folder = "results/m8_maps"

# read source info
source_info = pd.read_csv("data/output/source_info.csv")
source_info['Main'] = source_info['Main'].apply(lambda x: str(x).casefold())
source_info.set_index('Main', inplace=True)

# get files
files = glob.glob('/Volumes/Alpha/Research/data/*/*12co*.fits') + glob.glob('/Volumes/Alpha/Research/data/*/*spw39*.fits')

for file in files:

    # set output name and output path
    target_name = os.path.dirname(file).split('/')[-1]
    filename = f"{target_name}_m8.pdf"
    output_path = os.path.join(output_folder, filename)

    # verify output path exists and
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # skip already existing files if you don't want to overwrite them
    if (os.path.exists(output_path)) and (not overwrite):
        continue

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

        fig.savefig(output_path)

    except Exception as error:
        print(f"Error for {target_name}, {error}")