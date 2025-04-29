import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
import pandas as pd
import glob

from create_figs import create_m0_map, mark_sources, plot_vector

# script options
overwrite = False
output_folder = "results/m0_multi_outflow_plots"

# read data
df = pd.read_csv('data/output/outflow_data.csv')

source_info = pd.read_csv("data/output/source_info.csv")
source_info['Main'] = source_info['Main'].apply(lambda x: str(x).casefold())
source_info.set_index('Main', inplace=True)

# helper function to extract channel indices from the data table
def getIdx(listOf):
    ranges = []
    
    for string in listOf:
        if pd.isna(string):  # Skip NaN values
            continue
        
        for part in string.split(', '):
            if '-' in part:  # Handle ranges
                start, end = map(int, part.split('-'))
                ranges.append(np.r_[start:end + 1])  # Use np.r_
            else:  # Handle single values
                ranges.append(int(part))  # Append as integer
    
    # Ensure `ranges` is not empty before passing to np.r_
    if not ranges:
        return np.r_[:]  # Returns an empty np.r_

    return np.r_[tuple(ranges)]  # Use `tuple(ranges)` to avoid errors

# This loops through each source field with an angle measurement and creates
# a figure with the separation vector and outflow vector overlayed
for i, field in df.groupby('field').agg('first').reset_index().iterrows():

    # verify output path exists and
    # skip already existing files if you don't want to overwrite them
    target_name = field['field']
    filename = f"{target_name}_outflow.png"
    output_path = os.path.join(output_folder, filename)

    # verify output path exists and
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # skip already existing files if you don't want to overwrite them
    if (os.path.exists(output_path)) and (not overwrite):
        continue

    # open image
    image_filename = (glob.glob(f'/Volumes/Alpha/Research/data/{target_name.casefold()}/*12co*.fits') + glob.glob(f'/Volumes/Alpha/Research/data/{target_name.casefold()}/*spw39*.fits'))[0]
    hdulist = fits.open(image_filename)
    hdu = hdulist[0]

    # set center and size of cutout
    center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
    size = np.array([39, 39]) * u.arcsecond
    distance = field['distance']

    # FIGURE
    channels = getIdx([field['red_channels'], field['blue_channels']])
    fig = create_m0_map(hdu, center, size, channels, 3, distance)
    fig.set_title(f"{target_name} 12CO M0")

    ### VECTORS
    # plot binary separation angle
    center_origin = np.array([np.mean([field['source_a_ra'], field['source_b_ra']]), np.mean([field['source_a_dec'], field['source_b_dec']])])
    separation_angle_north = field['binary_PA']
    # draw separation vector in both directions
    plot_vector(fig, center_origin, separation_angle_north, color='white', length=0.005)
    plot_vector(fig, center_origin, separation_angle_north + 180, color='white', length=0.005)
    
    # plot each outflow vector
    for j, source in df[df['field'] == target_name].reset_index(drop=True).iterrows():
        # define vector origin at the outflow source
        if source['outflow_source'] == 'both':
            outflow_origin = center_origin
        elif source['outflow_source'] == source['source_a']:
            outflow_origin = np.array([source['source_a_ra'], source['source_a_dec']])
        else:
            outflow_origin = np.array([source['source_b_ra'], source['source_b_dec']])

        # plot outflow vector
        outflow_angle_north = source['outflow_PA']
        plot_vector(fig, outflow_origin, outflow_angle_north, color='red', length=0.005)
        
        # calculate delta_PA
        angle = np.abs(outflow_angle_north - separation_angle_north) % 180
        angle = np.min([angle, 180 - angle])
            
        # display angle between outflow and separation in top left corner
        source_label = source['outflow_source']
        if source_label == 'both':
            source_label = source['source_a']+'+'+source['source_b']
        fig.ax.text(30,fig.ax.get_xlim()[1]-50-40*(j), f"{source_label} : {np.abs(angle):.2f}°")

    ### MARKERS
    # add a marker at each source with legend
    target_info = source_info.loc[target_name.casefold()]
    mark_sources(fig, target_info)

    fig.savefig(os.path.join(output_folder, f"{target_name}_outflow.png"), dpi=300, transparent=True)
    fig.savefig(os.path.join(output_folder, f"{target_name}_outflow.pdf"))