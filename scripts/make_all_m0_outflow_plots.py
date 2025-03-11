import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd

from create_figs import create_m0_map

# script options
overwrite = False
output_folder = "results/m0_outflow_plots"

# read data
df = pd.read_csv('data/output/outflow_data.csv')

def angle_between_vectors(v1, v2):
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    a = np.arccos(np.clip(cos_theta, -1.0, 1.0)) * (180 / np.pi)
    return a

# helper function to extract channel indices from the data table
def getIdx(listOf):
    ranges = []
    for string in listOf:
        if pd.isna(string):
            continue
        for part in string.split(', '):
            if '-' in part:
                start, end = map(int, part.split('-'))
                ranges.append(np.r_[start:end+1])  # np.r_ includes the range
            else:
                ranges.append(int(part))
        for part in string.split(', '):
            if '-' in part:
                start, end = map(int, part.split('-'))
                ranges.append(np.r_[start:end+1])  # np.r_ includes the range
            else:
                ranges.append(int(part))

    # Flatten the result into a single NumPy array
    result = np.r_[tuple(ranges)]
    return result

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
    channels = getIdx([field['red_channels'], field['blue_channels']])
    fig = create_m0_map(hdu, center, size, channels, 3, distance)
    fig.set_title(f"{target_name} 12CO M0")
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

    ### calculate outflow vector
    angle_north = field['angle']
    angle_east_rad = np.radians(90 - angle_north)
    # get coordinate vectors
    # define outflow origin
    if field['outflow_source'] == 'both':
        outflow_origin = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    elif field['outflow_source'] == field['source_a']:
        outflow_origin = np.array([field['source_a_ra'], field['source_a_dec']])
    else:
        outflow_origin = np.array([field['source_b_ra'], field['source_b_dec']])
    outflow_tip_pix = fig.world2pixel(outflow_origin[0] + 0.005 * np.cos(angle_east_rad), outflow_origin[1] + 0.005 * np.sin(angle_east_rad))
    outflow_origin_pix = fig.world2pixel(outflow_origin[0], outflow_origin[1])
    # get outflow vector
    outflow_vector = np.array([outflow_tip_pix[0] - outflow_origin_pix[0], outflow_tip_pix[1] - outflow_origin_pix[1]])
    # plot outflow vector
    fig.ax.quiver(outflow_origin_pix[0], outflow_origin_pix[1], outflow_vector[0], outflow_vector[1],
                angles='xy', scale_units='xy', scale=1, color='red', width=0.005)


    ### calculate separation vector
    # get coordinate vectors
    star_a_pix = fig.world2pixel(center_a.ra.deg, center_a.dec.deg)
    star_b_pix = fig.world2pixel(center_b.ra.deg, center_b.dec.deg)
    # get separation vector
    separation_vector = np.array([star_a_pix[0] - star_b_pix[0], star_a_pix[1] - star_b_pix[1]])
    # normalize to same length as outflow vector
    separation_vector = np.linalg.norm(outflow_vector) * separation_vector / np.linalg.norm(separation_vector)
    # get smallest angle
    angle = angle_between_vectors(outflow_vector, separation_vector)
    if angle > 90:
        # recompute for the smaller angle
        separation_vector = np.array([star_b_pix[0] - star_a_pix[0], star_b_pix[1] - star_a_pix[1]])
        separation_vector = np.linalg.norm(outflow_vector) * separation_vector / np.linalg.norm(separation_vector)
        angle = angle_between_vectors(outflow_vector, separation_vector)
    # plot separation vector
    fig.ax.quiver(star_a_pix[0], star_a_pix[1], separation_vector[0], separation_vector[1],
                angles='xy', scale_units='xy', scale=1, color='white', width=0.005)


    # display angle between outflow and separation in top left corner
    fig.ax.text(30,fig.ax.get_xlim()[1]-50, f"{angle:.2f}°")

    # save image
    fig.savefig(output_path)