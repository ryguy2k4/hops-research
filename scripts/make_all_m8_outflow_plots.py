import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd

from create_figs import create_m8_map, mark_sources_2, plot_vector

# script options
overwrite = False
output_folder = "results/m8_outflow_plots"

# read data
df = pd.read_csv('data/output/outflow_data.csv')

def angle_between_vectors(v1, v2):
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    a = np.arccos(np.clip(cos_theta, -1.0, 1.0)) * (180 / np.pi)
    return a

# This loops through each source field with an angle measurement and creates
# a figure with the separation vector and outflow vector overlayed
for i, field in df.iterrows():

    # verify output path exists and
    # skip already existing files if you don't want to overwrite them
    target_name = field['field']
    filename = f"{target_name}_outflow.pdf"
    output_path = os.path.join(output_folder, filename)

    # verify output path exists and
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # skip already existing files if you don't want to overwrite them
    if (os.path.exists(output_path)) and (not overwrite):
        continue

    # open image
    hdulist = fits.open(f"/Volumes/Alpha/Research/data/{target_name.casefold()}/{target_name.casefold()}__s15__12co.fits")
    hdu = hdulist[0]

    # set center and size of cutout
    center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
    size = np.array([39, 39]) * u.arcsecond
    distance = field['distance']

    # create figure
    fig = create_m8_map(hdu, center, size, distance)
    fig.set_title(f"{target_name} 12CO M8")
    fig.show_colorscale(cmap='viridis', stretch='sqrt')

    # add a marker at each source with legend
    mark_sources_2(fig, field)

    ### VECTORS
    # define vector origin at the outflow source
    if field['outflow_source'] == 'both':
        outflow_origin = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    elif field['outflow_source'] == field['source_a']:
        outflow_origin = np.array([field['source_a_ra'], field['source_a_dec']])
    else:
        outflow_origin = np.array([field['source_b_ra'], field['source_b_dec']])

    # plot outflow angle
    outflow_angle_north = field['outflow_PA']
    plot_vector(fig, outflow_origin, outflow_angle_north, color='red', length=0.005)
    
    # plot separation angle
    separation_angle_north = field['binary_PA']
    # choose a separation vector that provides the smallest angle between vectors
    angle = np.abs(outflow_angle_north - separation_angle_north)
    if angle < 90:
        plot_vector(fig, outflow_origin, separation_angle_north, color='white', length=0.005)
    else:
        plot_vector(fig, outflow_origin, separation_angle_north + 180, color='white', length=0.005)
        angle = 180 - angle

    # display angle between outflow and separation in top left corner
    fig.ax.text(30,fig.ax.get_xlim()[1]-50, f"{np.abs(angle):.2f}°")

    # save image
    fig.savefig(output_path)