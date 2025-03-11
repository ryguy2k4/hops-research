import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd

from create_figs import create_m8_map

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

    # extract center coords
    center_ra = hdu.header['OBSRA']
    center_dec = hdu.header['OBSDEC']

    # set center and size of cutout
    center = SkyCoord(center_ra, center_dec, unit=u.degree)
    size = np.array([39, 39]) * u.arcsecond
    distance = field['distance']

    # create figure
    fig = create_m8_map(hdu, center, size, distance)
    fig.set_title(f"{target_name} 12CO M8")
    fig.show_colorscale(cmap='viridis', stretch='sqrt')

    # add a marker at each source with legend
    center_a = SkyCoord(field['source_a_ra'], field['source_a_dec'], unit=u.degree)
    center_b = SkyCoord(field['source_b_ra'], field['source_b_dec'], unit=u.degree)
    fig.show_markers(center_a.ra.deg, center_a.dec.deg, coords_frame='world', marker='x', s=25, c='black', linewidths=1, label=field['source_a'])
    fig.show_markers(center_b.ra.deg, center_b.dec.deg, coords_frame='world', marker='x', s=25, c='magenta', linewidths=1, label=field['source_b'])
    legend_handles = []
    legend_handles.append(mlines.Line2D([], [], color='black', marker='x', markersize=6, linestyle='None', label=field['source_a']))
    legend_handles.append(mlines.Line2D([], [], color='magenta', marker='x', markersize=6, linestyle='None', label=field['source_b']))
    fig.ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1,1.15))


    ### VECTORS
    def make_vector(origin, angle_north, length=0.005):
        angle_east_rad = np.radians(90 - angle_north)
        return origin[0] + length * np.cos(angle_east_rad), origin[1] + length * np.sin(angle_east_rad)

    # define vector origin at the outflow source
    if field['outflow_source'] == 'both':
        outflow_origin = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    elif field['outflow_source'] == field['source_a']:
        outflow_origin = np.array([field['source_a_ra'], field['source_a_dec']])
    else:
        outflow_origin = np.array([field['source_b_ra'], field['source_b_dec']])

    outflow_origin_pix = fig.world2pixel(outflow_origin[0], outflow_origin[1])

    # get outflow angle
    outflow_angle_north = field['outflow_angle']
    # create outflow vector
    outflow_tip = make_vector(outflow_origin, outflow_angle_north)
    outflow_tip_pix = fig.world2pixel(outflow_tip[0], outflow_tip[1])
    outflow_vector = np.array([outflow_tip_pix[0] - outflow_origin_pix[0], outflow_tip_pix[1] - outflow_origin_pix[1]])
    
    # get separation angle
    separation_angle_north = field['separation_angle']
    # choose a separation vector that provides the smallest angle between vectors
    angle = np.abs(outflow_angle_north - separation_angle_north)
    if angle < 90:
        separation_tip = make_vector(outflow_origin, separation_angle_north)
    else:
        separation_tip = make_vector(outflow_origin, separation_angle_north + 180)
        angle = 180 - angle
    separation_tip_pix = fig.world2pixel(separation_tip[0], separation_tip[1])
    separation_vector = np.array([separation_tip_pix[0] - outflow_origin_pix[0], separation_tip_pix[1] - outflow_origin_pix[1]])
    
    # plot separation vector
    fig.ax.quiver(outflow_origin_pix[0], outflow_origin_pix[1], separation_vector[0], separation_vector[1],
                angles='xy', scale_units='xy', scale=1, color='white', width=0.005)
    # plot outflow vector
    fig.ax.quiver(outflow_origin_pix[0], outflow_origin_pix[1], outflow_vector[0], outflow_vector[1],
                angles='xy', scale_units='xy', scale=1, color='red', width=0.005)

    # display angle between outflow and separation in top left corner
    fig.ax.text(30,fig.ax.get_xlim()[1]-50, f"{np.abs(angle):.2f}°")

    # save image
    fig.savefig(output_path)