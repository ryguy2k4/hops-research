import warnings
warnings.filterwarnings("ignore")

import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import os
import pandas as pd
import glob

from create_figs import create_m0_map, create_m8_map, getIdx, plot_vector, mark_sources

### SCRIPT OPTIONS

# output
output_folder = "results/figs_for_paper"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# this script assumes that the directory below contains a folder for each field
# and that within each field folder there is a 12CO image, which contains
# '12co' or 'spw39' in the filename
image_directory = "/Volumes/Alpha/Research/data/"

### SCRIPT
# read data
outflow_data = pd.read_csv('data/output/outflow_data.csv')
outflow_data2 = outflow_data.groupby('field').first().reset_index()
by_field = pd.read_csv('data/output/data_by_field.csv')
source_info = pd.read_csv("data/output/source_info.csv")
source_info['Main'] = source_info['Main'].apply(lambda x: str(x).casefold())
source_info.set_index('Main', inplace=True)

# define fields in each figure
fig_1 = ['HOPS-32', 'HOPS-168', 'HOPS-281', 'Per-emb-17']
fig_2 = ['HOPS-290', 'HOPS-288']
no_outflows = ['HOPS-28', 'HOPS-163', 'HOPS-242', 'HOPS-255', 'HOPS-248', 'HOPS-357']
app_1 = by_field[~by_field['field'].isin(no_outflows)]['field'].values
app_2 = by_field[by_field['field'].isin(no_outflows)]['field'].values

def make_compound_plots(field_list, output_name, rows, cols, wspace=0.35, hspace=0.2):
    figure = plt.figure(figsize=(17/3*cols, 22/5*rows))

    # This loops through each source field with an angle measurement and creates
    # a figure with the separation vector and outflow vector overlayed

    data = by_field[by_field['field'].isin(field_list)]
    for i, field in data.sort_values(by='ra', ascending=True).reset_index().iterrows():

        # verify output path exists and
        # skip already existing files if you don't want to overwrite them
        target_name = field['field']

        # verify output path exists and
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        # open image
        image_filename = (glob.glob(f'/Volumes/Alpha/Research/data/{target_name.casefold()}/*12co*.fits') + glob.glob(f'/Volumes/Alpha/Research/data/{target_name.casefold()}/*spw39*.fits'))[0]
        hdulist = fits.open(image_filename)    
        hdu = hdulist[0]

        # set center and size of cutout
        center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
        size = np.array([39, 39]) * u.arcsecond
        distance = source_info.loc[target_name.casefold(), 'Dis'].iloc[0]

        if pd.isna(field['integrated_channels']):
            # m8 map
            # create figure
            fig = create_m8_map(hdu, center, size, distance, figure=figure, subplot=(rows, cols, i+1), multiimage=True)
            fig.set_title(f"{target_name} M8")   
        
        else:
            # m0 map
            # create figure
            outflow = outflow_data.loc[outflow_data['field'] == target_name].groupby('field').first().reset_index()
            channels = getIdx([outflow.at[0, 'red_channels'], outflow.at[0, 'blue_channels']])
            fig = create_m0_map(hdu, center, size, channels, 3, distance, figure, subplot=(rows, cols, i+1), multiimage=True)
            fig.set_title(f"{target_name} M0")

            ### VECTORS
            # plot binary separation angle
            center_origin = np.array([np.mean([outflow['source_a_ra'], outflow['source_b_ra']]), np.mean([outflow['source_a_dec'], outflow['source_b_dec']])])
            separation_angle_north = outflow['binary_PA']
            # draw separation vector in both directions
            plot_vector(fig, center_origin, separation_angle_north, color='white', length=0.005)
            plot_vector(fig, center_origin, separation_angle_north + 180, color='white', length=0.005)
            
            # plot each outflow vector
            for j, source in outflow_data[outflow_data['field'] == target_name].reset_index(drop=True).iterrows():
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

        # manage axis labels
        if i % cols != cols - 1:
            fig.colorbar.set_axis_label_text("")

        if i % cols != 0:
            fig.axis_labels.hide_y()

        if (rows > 1) & (i < rows*cols - cols):
            fig.axis_labels.hide_x()

        ### MARKERS
        # add a marker at each source with legend
        target_info = source_info.loc[target_name.casefold()]
        mark_sources(fig, target_info, use_short_label=True)

    plt.subplots_adjust(wspace=wspace, hspace=hspace)

    # save image
    figure.savefig(os.path.join(output_folder, f"{output_name}.pdf"), bbox_inches='tight')
    # figure.savefig(os.path.join(output_folder, f"{output_name}.png"), dpi=800, transparent=True, bbox_inches='tight')
    figure.show()

make_compound_plots(fig_1, "fig_1", rows=1, cols=4, wspace=0.4, hspace=0.1)
make_compound_plots(fig_2, "fig_2", rows=1, cols=4, wspace=0.4, hspace=0.1)
make_compound_plots(app_1[0:15], "appendix-1", rows=5, cols=3)
make_compound_plots(app_1[15:30], "appendix-2", rows=5, cols=3)
make_compound_plots(app_1[30:45], "appendix-3", rows=5, cols=3)
make_compound_plots(app_2, "appendix-4", rows=2, cols=3)