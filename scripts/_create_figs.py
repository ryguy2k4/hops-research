import astropy.io.fits as fits
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import aplpy
import astropy.wcs.wcs as wcs
import pandas as pd

"""

Creates a simple aplpy FITSFigure with a colorbar, scalebar, and beam

"""
def create_fig(img, distance=0, colorbar_label="Intensity (Jy/Beam)", figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False):
    # Create FITSFigure
    if not multiimage:
        figure.clear()

    fig = aplpy.FITSFigure(img, figure=figure, subplot=subplot)

    # Display image
    fig.show_colorscale(cmap='viridis', stretch='sqrt')

    # Title
    fig.set_title(img.header['OBJECT'])

    # Colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text(colorbar_label)
    fig.colorbar.set_location('right')
    fig.colorbar.set_pad(0.0)

    # Scalebar
    if distance == 0:
        fig.add_scalebar(1 * u.arcsecond)
        fig.scalebar.set_label('1"')
    else:
        fig.add_scalebar((500 / (distance)) * u.arcsecond)
        fig.scalebar.set_label('500 AU')
        fig.scalebar.set_font_size(8)

    fig.scalebar.set_linewidth(3)

    # Beam
    fig.add_beam()
    fig.beam.show()
    fig.beam.set_corner("bottom left")
    fig.beam.set_color("black")
    fig.beam.set_edgecolor("black")
    fig.beam.set_facecolor("white")

    return fig

"""

Creates a simple aplpy FITSFigure with a colorbar, scalebar, and beam
that has been cropped with parameters `center` and `size`

"""
def create_sub_fig(hdu, center, size, distance=0, colorbar_label="Intensity (Jy/Beam)", figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False):
    # get the wcs
    wcs_original = wcs.WCS(hdu)
    wcs_2d = wcs_original.celestial

    # make the cut
    cut = Cutout2D(hdu.data[0,0,:,:], center, size, wcs=wcs_2d)

    # create a header with the wcs
    # CAUTION: lots of header data not being preserved here
    new_header = cut.wcs.to_header()
    for key in ['BUNIT', 'OBJECT', 'TELESCOP', 'INSTRUME', 'BMAJ', 'BMIN', 'BPA']:
        if key in hdu.header:
            new_header[key] = hdu.header[key]

    # combine cut and wcs and an HDU
    hdu_cut = fits.PrimaryHDU(data=cut.data, header=new_header)
    
    return create_fig(hdu_cut, distance, colorbar_label="Intensity (Jy/Beam)", figure=figure, subplot=subplot, multiimage=multiimage)

"""

Takes image data and image header and returns a cropped image
packaged into a primary hdu, with a 2D WCS

"""
def cut_fig(data, header, center, size):
    # get the wcs
    wcs_original = wcs.WCS(header)
    wcs_2d = wcs_original.celestial

    # make the cut
    cut = Cutout2D(data, center, size, wcs=wcs_2d)

    # create a header with the wcs
    # CAUTION: lots of header data not being preserved here
    new_header = cut.wcs.to_header()
    for key in ['BUNIT', 'OBJECT', 'TELESCOP', 'INSTRUME', 'BMAJ', 'BMIN', 'BPA']:
        if key in header:
            new_header[key] = header[key]

    # combine cut and wcs and an HDU
    return fits.PrimaryHDU(data=cut.data, header=new_header)


"""

Takes an image and creates a moment 8 map

"""
def create_m8_map(hdu, center, size, distance=0, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False):
    # create map
    m8_map = np.max(hdu.data[0,:,:,:], axis=0)

    cut = cut_fig(m8_map, hdu.header, center, size)
    
    return create_fig(cut, distance, colorbar_label="Intensity (Jy/bm)", figure=figure, subplot=subplot, multiimage=multiimage)


"""

Takes an image and creates a moment 0 map
Included channels are specified with channel_idx and
sigma-clipping is specified with sigma; default is 3

"""
def create_m0_map(hdu, center, size, channel_idx, sigma=3, distance=0, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False):
    # sigma
    blank_channel = hdu.data[0,0,:,:]
    region_size = np.array([100, 100]) * u.pixel
    wcs_2d = wcs.WCS(hdu.header).celestial
    region = Cutout2D(blank_channel, center, region_size, wcs=wcs_2d)
    mean = np.mean(region.data)
    std = np.std(region.data)

    # create map
    channels = hdu.data[0,channel_idx,:,:]
    channels[channels < (mean + sigma*std)] = 0
    m0_map = np.sum(channels, axis=0)

    cut = cut_fig(m0_map, hdu.header, center, size)
    
    return create_fig(cut, distance, colorbar_label="Intensity (Jy/bm km/s)", figure=figure, subplot=subplot, multiimage=multiimage)


"""

Marks sources given rows from 'source_info.csv'

"""
def mark_sources(fig, source_rows, use_short_label=False):
    marker_colors = ['black', 'magenta', 'red', 'darkred', 'darkblue']
    legend_handles = []
    sources_to_mark = source_rows.reset_index()[0:4]
    for i, row in sources_to_mark.iterrows():
        center2 = SkyCoord(row['RA'], row['Dec'], unit=u.degree)
        fig.show_markers(center2.ra.deg, center2.dec.deg, coords_frame='world', marker='x', s=25, c=marker_colors[i], linewidths=1, label=row['Source'])

        # Create legend handle for this source (only if not already added)
        short_label = str(row['Source']).casefold().removeprefix(str(source_rows.index.tolist()[i]).casefold()+'-').upper()
        legend_handles.append(mlines.Line2D([], [], color=marker_colors[i], marker='x', markersize=4, linestyle='None', label=short_label if use_short_label else row['Source']))

    fig.ax.legend(handles=legend_handles, fontsize=6, loc='upper right', bbox_to_anchor=(0.98,1))


"""

Marks sources given a single row from 'outflow_data.csv'

"""
def mark_sources_2(fig, field):
    center_a = SkyCoord(field['source_a_ra'], field['source_a_dec'], unit=u.degree)
    center_b = SkyCoord(field['source_b_ra'], field['source_b_dec'], unit=u.degree)
    fig.show_markers(center_a.ra.deg, center_a.dec.deg, coords_frame='world', marker='x', s=25, c='black', linewidths=1, label=field['source_a'])
    fig.show_markers(center_b.ra.deg, center_b.dec.deg, coords_frame='world', marker='x', s=25, c='magenta', linewidths=1, label=field['source_b'])
    legend_handles = []
    legend_handles.append(mlines.Line2D([], [], color='black', marker='x', markersize=6, linestyle='None', label=field['source_a']))
    legend_handles.append(mlines.Line2D([], [], color='magenta', marker='x', markersize=6, linestyle='None', label=field['source_b']))
    fig.ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1,1.15))

def mark_sources_3(fig, field_rows):
    legend_handles = []
    for i, field in field_rows.iterrows():
        center_a = SkyCoord(field['source_a_ra'], field['source_a_dec'], unit=u.degree)
        fig.show_markers(center_a.ra.deg, center_a.dec.deg, coords_frame='world', marker='x', s=25, c='black', linewidths=1, label=field['source_a'])
        legend_handles.append(mlines.Line2D([], [], color='black', marker='x', markersize=6, linestyle='None', label=field['source_a']))
        if ~pd.isna(field['source_b']):
            center_b = SkyCoord(field['source_b_ra'], field['source_b_dec'], unit=u.degree)
            fig.show_markers(center_b.ra.deg, center_b.dec.deg, coords_frame='world', marker='x', s=25, c='magenta', linewidths=1, label=field['source_b'])
            legend_handles.append(mlines.Line2D([], [], color='magenta', marker='x', markersize=6, linestyle='None', label=field['source_b']))
    fig.ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1,1.15))

"""

Plots a vector on the figure, starting from `origin`, pointing in
the direction of `angle` with length `length` (default 0.005)

"""
def plot_vector(fig, origin, angle_north_deg, color, length=0.005):
    angle_east_rad = np.radians(90 - angle_north_deg)
    origin_pix = fig.world2pixel(origin[0], origin[1])
    tip = origin[0] + length * np.cos(angle_east_rad), origin[1] + length * np.sin(angle_east_rad)
    tip_pix = fig.world2pixel(tip[0], tip[1])
    outflow_vector = np.array([tip_pix[0] - origin_pix[0], tip_pix[1] - origin_pix[1]])
    fig.ax.quiver(origin_pix[0], origin_pix[1], outflow_vector[0], outflow_vector[1],
                angles='xy', scale_units='xy', scale=1, color=color, width=0.005)
    
def plot_dotted_vector(fig, origin, angle_north_deg, color, length=0.005):
    angle_east_rad = np.radians(90 - angle_north_deg)
    origin_pix = fig.world2pixel(origin[0], origin[1])
    tip = origin[0] + length * np.cos(angle_east_rad), origin[1] + length * np.sin(angle_east_rad)
    tip_pix = fig.world2pixel(tip[0], tip[1])
    outflow_vector = np.array([tip_pix[0] - origin_pix[0], tip_pix[1] - origin_pix[1]])
    fig.ax.plot(
        [origin_pix[0], origin_pix[0] + outflow_vector[0]], 
        [origin_pix[1], origin_pix[1] + outflow_vector[1]], 
        linestyle="dashed", color=color, linewidth=0.5  # Dotted line
    )

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