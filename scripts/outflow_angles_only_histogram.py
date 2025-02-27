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
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.arccos(np.clip(cos_theta, -1.0, 1.0)) * (180 / np.pi)  # Convert to degrees

df = pd.read_csv('../data/outflow_data.csv')


# This loops through each source field with an angle measurement and creates
# a figure with the separation vector an outflow vector overlayed
smallest_angles = []
for i, field in df.iterrows():

    target_name = field['field']

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

    ### calculate outflow vector
    center_a = SkyCoord(field['source_a_ra'], field['source_a_dec'], unit=u.degree)
    center_b = SkyCoord(field['source_b_ra'], field['source_b_dec'], unit=u.degree)
    angle_north = field['angle']
    angle_east_rad = np.radians(90 - angle_north)
    outflow_tip_pix = fig.world2pixel(center_ra + 0.005 * np.cos(angle_east_rad), center_dec + 0.005 * np.sin(angle_east_rad))
    center_coords = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    center_pix = fig.world2pixel(center_coords[0], center_coords[1])
    outflow_vector = np.array([outflow_tip_pix[0] - center_pix[0], outflow_tip_pix[1] - center_pix[1]])

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

    smallest_angles.append({"field": target_name, "smallest_angle": angle})

# create histogram
df2 = pd.DataFrame(smallest_angles)
fig2, ax = plt.subplots()
ax.hist(df2['smallest_angle'], label=f"N = {len(df2)}")
ax.legend(loc='upper left')
ax.set_xlabel("smallest angle between binary separation and outflow")
ax.set_ylabel("count")
fig2.savefig("../results/histogram3.png")
df2.to_csv("../results/angles.csv", index=False)