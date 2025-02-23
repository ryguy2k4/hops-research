import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd

from create_figs import create_m8_map

import warnings
warnings.filterwarnings("ignore")

def angle_between_vectors(v1, v2):
    """Calculate the angle (in degrees) between two vectors."""
    dot_product = np.dot(v1, v2)  # Compute dot product
    norm_v1 = np.linalg.norm(v1)  # Compute magnitude of v1
    norm_v2 = np.linalg.norm(v2)  # Compute magnitude of v2
    
    # Compute angle in radians and then convert to degrees
    angle_rad = np.arccos(dot_product / (norm_v1 * norm_v2))
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

# data
df = pd.read_csv('../data/outflow_data.csv')

# exclusions for semi-unclear outflows
exlude = ['HOPS-363', 'HOPS-213', 'HOPS-75', 'HOPS-384', 'HOPS-173', 'HOPS-77']
df = df[~df['field'].isin(exlude)]


# This loops through each source field with an angle measurement and creates
# a figure with the separation vector an outflow vector overlayed
overwrite = False
smallest_angles = []
for i, field in df.iterrows():

    # verify output path exists and
    # skip already existing files if you don't want to overwrite them
    target_name = field['field']
    output_path = f"../results/outflow_alignment2/{target_name}_outflow.png"
    if (os.path.exists(output_path)) and (not overwrite):
        continue

    # open image
    hdulist = fits.open(f"/Volumes/Alpha/Research/data/{target_name.casefold()}/{target_name.casefold()}__s15__12co.fits")
    hdu = hdulist[0]

    # crop field
    center_ra = hdu.header['OBSRA']
    center_dec = hdu.header['OBSDEC']
    center = SkyCoord(center_ra, center_dec, unit=u.degree)
    size = np.array([39, 39]) * u.arcsecond

    # create figure
    fig = create_m8_map(hdu, center, size, distance=field['distance'])
    fig.set_title(f"{target_name} 12CO M8")
    fig.show_colorscale(cmap='viridis', stretch='sqrt')

    # add a marker at each source with legend
    center_a = SkyCoord(field['source_a_ra'], field['source_a_dec'], unit=u.degree)
    center_b = SkyCoord(field['source_b_ra'], field['source_b_dec'], unit=u.degree)
    fig.show_markers(center_a.ra.deg, center_a.dec.deg, coords_frame='world', marker='x', s=25, c='black', linewidths=1, label=field['source_a'])
    fig.show_markers(center_b.ra.deg, center_b.dec.deg, coords_frame='world', marker='x', s=25, c='red', linewidths=1, label=field['source_b'])
    legend_handles = []
    legend_handles.append(mlines.Line2D([], [], color='black', marker='x', markersize=6, linestyle='None', label=field['source_a']))
    legend_handles.append(mlines.Line2D([], [], color='red', marker='x', markersize=6, linestyle='None', label=field['source_b']))
    fig.ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1,1.15))


    ### calculate outflow vector
    angle_north = field['angle']
    angle_east_rad = np.radians(90 - angle_north)
    outflow_tip_pix = fig.world2pixel(center_ra + 0.005 * np.cos(angle_east_rad), center_dec + 0.005 * np.sin(angle_east_rad))
    center_coords = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    center_pix = fig.world2pixel(center_coords[0], center_coords[1])
    outflow_vector = np.array([outflow_tip_pix[0] - center_pix[0], outflow_tip_pix[1] - center_pix[1]])
    # plot outflow vector
    fig.ax.quiver(center_pix[0], center_pix[1], outflow_vector[0], outflow_vector[1],
                angles='xy', scale_units='xy', scale=1, color='red', width=0.005)


    ### calculate separation vector
    star_a_pix = fig.world2pixel(center_a.ra.deg, center_a.dec.deg)
    star_b_pix = fig.world2pixel(center_b.ra.deg, center_b.dec.deg)
    separation_vector = np.array([star_a_pix[0] - star_b_pix[0], star_a_pix[1] - star_b_pix[1]])
    separation_vector = np.linalg.norm(outflow_vector) * separation_vector / np.linalg.norm(separation_vector)
    # get smallest angle
    angle = angle_between_vectors(outflow_vector, separation_vector)
    if angle > 90:
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
    smallest_angles.append({"field": target_name, "smallest_angle": angle})

# create histogram
df2 = pd.DataFrame(smallest_angles)
fig2, ax = plt.subplots()
ax.hist(df2['smallest_angle'], label=f"N = {len(df2)}")
ax.legend(loc='upper left')
ax.set_xlabel("smallest angle between binary separation and outflow")
ax.set_ylabel("count")
fig2.savefig("../results/outflow_alignment2/histogram.png")