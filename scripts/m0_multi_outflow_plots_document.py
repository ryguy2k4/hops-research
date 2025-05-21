import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
from matplotlib.backends.backend_pdf import PdfPages

from create_figs import create_m0_map, mark_sources, plot_vector

# script options
output_folder = "results/m0_multi_outflow_plots"
output_pdf = os.path.join(output_folder, "all_m0_multi_outflow_plots.pdf")
# overwrite output files if they already exist
overwrite = False
# this script assumes that the directory below contains a folder for each field
# and that within each field folder there is a 12CO image, which contains
# '12co' or 'spw39' in the filename
image_directory = "/Volumes/Alpha/Research/data/"

# read data
df = pd.read_csv('data/output/outflow_data.csv')

source_info = pd.read_csv("data/output/source_info.csv")
source_info['Main'] = source_info['Main'].apply(lambda x: str(x).casefold())
source_info.set_index('Main', inplace=True)

# helper function to extract channel indices from the data table
def getIdx(listOf):
    ranges = []
    for string in listOf:
        if pd.isna(string):
            continue
        for part in string.split(', '):
            if '-' in part:
                start, end = map(int, part.split('-'))
                ranges.append(np.r_[start:end + 1])
            else:
                ranges.append(int(part))
    if not ranges:
        return np.r_[:]
    return np.r_[tuple(ranges)]

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

with PdfPages(output_pdf) as pdf:
    for i, field in df.groupby('field').agg('first').reset_index().iterrows():
        target_name = field['field']

        image_filename = (glob.glob(f'{image_directory}{target_name.casefold()}/*12co*.fits') +
                          glob.glob(f'{image_directory}{target_name.casefold()}/*spw39*.fits'))[0]
        hdulist = fits.open(image_filename)
        hdu = hdulist[0]

        center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
        size = np.array([39, 39]) * u.arcsecond
        distance = field['distance']

        channels = getIdx([field['red_channels'], field['blue_channels']])
        fig = create_m0_map(hdu, center, size, channels, 3, distance)
        # Set figure size to letter (8.5 x 11 inches)
        fig._figure.set_size_inches(8.5, 11)
        fig.set_title(f"{target_name} 12CO M0")

        center_origin = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]),
                                  np.mean([field['source_a_dec'], field['source_b_dec']])])
        separation_angle_north = field['binary_PA']
        plot_vector(fig, center_origin, separation_angle_north, color='white', length=0.005)
        plot_vector(fig, center_origin, separation_angle_north + 180, color='white', length=0.005)

        for j, source in df[df['field'] == target_name].reset_index(drop=True).iterrows():
            if source['outflow_source'] == 'both':
                outflow_origin = center_origin
            elif source['outflow_source'] == source['source_a']:
                outflow_origin = np.array([source['source_a_ra'], source['source_a_dec']])
            else:
                outflow_origin = np.array([source['source_b_ra'], source['source_b_dec']])

            outflow_angle_north = source['outflow_PA']
            plot_vector(fig, outflow_origin, outflow_angle_north, color='red', length=0.005)

            angle = np.abs(outflow_angle_north - separation_angle_north) % 180
            angle = np.min([angle, 180 - angle])

            source_label = source['outflow_source']
            if source_label == 'both':
                source_label = source['source_a'] + '+' + source['source_b']
            fig.ax.text(30, fig.ax.get_xlim()[1] - 50 - 40 * (j), f"{source_label} : {np.abs(angle):.2f}°")

        target_info = source_info.loc[target_name.casefold()]
        mark_sources(fig, target_info)

        # SAVE
        pdf.savefig(fig._figure, transparent=True)
        fig._figure.clf()