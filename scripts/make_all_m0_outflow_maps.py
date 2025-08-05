import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob
from matplotlib.backends.backend_pdf import PdfPages

from scripts._create_figs import create_m0_map, mark_sources, plot_vector, getIdx
from scripts._script_options import IMAGE_DIRECTORY


# SET OUTPUT
output_folder = "results"
output_pdf = os.path.join(output_folder, "m0_outflow_maps.pdf")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# READ DATA
outflow_data = pd.read_csv('data/output/outflow_data.csv')
source_info = pd.read_csv("../data/output/source_info.csv", index_col='Main')
source_info.index = source_info.index.str.casefold()

# SCRIPT
with PdfPages(output_pdf) as pdf:
    for i, field in outflow_data[~outflow_data['binary_PA'].isna()].groupby('field').agg('first').sort_values('source_a_ra').reset_index().iterrows():
        target_name = field['field']

        image_filename = (glob.glob(f'{IMAGE_DIRECTORY}{target_name.casefold()}/*12co*.fits') +
                          glob.glob(f'{IMAGE_DIRECTORY}{target_name.casefold()}/*spw39*.fits'))[0]
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

        for j, source in outflow_data[outflow_data['field'] == target_name].reset_index(drop=True).iterrows():
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