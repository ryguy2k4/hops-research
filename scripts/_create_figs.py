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
def create_fig(img, distance=0, scalebar_au=500, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False, use_offset_labels=False):
    # Create FITSFigure
    if not multiimage:
        figure.clear()

    if np.nanmax(img.data) < 0.1:
        # Jy → mJy
        units = "mJy/bm"
        img.data *= 1e3
    else:
        units = "Jy/bm"

    fig = aplpy.FITSFigure(img, figure=figure, subplot=subplot)

    # Display image
    fig.show_colorscale(cmap='inferno', stretch='sqrt')

    # Title
    fig.set_title(img.header['OBJECT'])

    # Colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text(f"Intensity ({units})")
    fig.colorbar.set_location('right')
    fig.colorbar.set_pad(0.0)

    # Scalebar
    if distance == 0:
        fig.add_scalebar(1 * u.arcsecond)
        fig.scalebar.set_label('1"')
    else:
        fig.add_scalebar((scalebar_au / (distance)) * u.arcsecond)
        fig.scalebar.set_label(f"{scalebar_au} au")
        fig.scalebar.set_font_size(8)

    fig.scalebar.set_linewidth(3)

    # Beam
    fig.add_beam()
    fig.beam.show()
    fig.beam.set_corner("bottom left")
    fig.beam.set_color("black")
    fig.beam.set_edgecolor("black")
    fig.beam.set_facecolor("white")

    # Offset Tick Labels
    if use_offset_labels:
        # hide aplpy ticks
        fig.tick_labels.hide()
        # enable plain Matplotlib ticks
        fig.ax.tick_params(axis='both', which='both', direction='out', length=5, labelsize=10)
        fig.ax.xaxis.set_visible(True)
        fig.ax.yaxis.set_visible(True)

        # create offset ticks
        tick_spacing = 10  # arcsec
        xticks_new = np.arange(-20, 21, tick_spacing)
        yticks_new = np.arange(-20, 21, tick_spacing)
        # convert offsets to pixel positions
        pixscale = np.mean(np.abs(fig._wcs.pixel_scale_matrix.diagonal())) * 3600
        x_center = fig._data.shape[1]/2
        y_center = fig._data.shape[0]/2
        xticks_pix_new = xticks_new / pixscale + x_center
        yticks_pix_new = yticks_new / pixscale + y_center

        # Apply new tick positions and labels
        fig.ax.set_xticks(xticks_pix_new)
        fig.ax.set_yticks(yticks_pix_new)
        fig.ax.set_xticklabels([f"{x:.0f}" for x in xticks_new])
        fig.ax.set_yticklabels([f"{y:.0f}" for y in yticks_new])
        fig.axis_labels.set_xtext("RA Offset (arcsec)")
        fig.axis_labels.set_ytext("Dec Offset (arcsec)")
        fig.axis_labels.set_xpad(3)
        fig.axis_labels.set_ypad(3)

    figure.show()
    return fig



"""

Takes image data and image header and returns a cropped image
packaged into a primary hdu, with a 2D WCS

"""
def cut_fig(data, header, center=None, size=None):
    # get the wcs
    wcs_original = wcs.WCS(header)
    wcs_2d = wcs_original.celestial

    # center and size
    if center == None:
        center = SkyCoord(header['OBSRA'], header['OBSDEC'], unit=u.degree)
    if size == None:
        size = np.array([39, 39]) * u.arcsecond

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


def create_cont_map(hdu, center=None, size=None, distance=0, scalebar_au=500, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False, use_offset_labels=False):

    # crop_center should determine the crop
    # field_center should determine the offset labels
    if center == None:
        field_center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
        crop_center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
    else:
        field_center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
        crop_center = center

    cut = cut_fig(hdu.data[0,0,:,:], hdu.header, crop_center, size)
    
    fig = create_fig(cut, distance=distance, scalebar_au=scalebar_au, figure=figure, subplot=subplot, multiimage=multiimage)
    fig.scalebar.set_color('white')

    if use_offset_labels:
        # need to do this separately because the scale is different here
        # hide aplpy ticks
        fig.tick_labels.hide()
        # enable plain Matplotlib ticks
        fig.ax.tick_params(axis='both', which='both', direction='out', length=5, labelsize=10)
        fig.ax.xaxis.set_visible(True)
        fig.ax.yaxis.set_visible(True)

        # create offset ticks
        tick_spacing = 0.5  # arcsec
        xticks_new = np.arange(-2, 3, tick_spacing)
        yticks_new = np.arange(-1.5, 2.5, tick_spacing)
        # convert offsets to pixel positions
        pixscale = np.mean(np.abs(fig._wcs.pixel_scale_matrix.diagonal())) * 3600
        x_center, y_center = fig.world2pixel(field_center.ra, field_center.dec)
        xticks_pix_new = xticks_new / pixscale + x_center
        yticks_pix_new = yticks_new / pixscale + y_center

        # Apply new tick positions and labels
        fig.ax.set_xticks(xticks_pix_new)
        fig.ax.set_yticks(yticks_pix_new)
        fig.ax.set_xticklabels([f"{x:.1f}" for x in xticks_new])
        fig.ax.set_yticklabels([f"{y:.1f}" for y in yticks_new])
        fig.axis_labels.set_xtext("RA Offset (arcsec)")
        fig.axis_labels.set_ytext("Dec Offset (arcsec)")
        fig.axis_labels.set_xpad(3)
        fig.axis_labels.set_ypad(3)

        # fix limits
        ny, nx = fig._data.shape
        fig.ax.set_xlim(-0.5, nx - 0.5)
        fig.ax.set_ylim(-0.5, ny - 0.5)

    return fig


"""

Takes an image and creates a moment 8 map

"""
def create_m8_map(hdu, center=None, size=None, distance=0, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False, use_offset_labels=False):
    # create map
    m8_map = np.max(hdu.data[0,:,:,:], axis=0)

    cut = cut_fig(m8_map, hdu.header, center, size)
    
    return create_fig(cut, distance=distance, figure=figure, subplot=subplot, multiimage=multiimage, use_offset_labels=use_offset_labels)


"""

Takes an image and creates a moment 0 map
Included channels are specified with channel_idx and
sigma-clipping is specified with sigma; default is 3

"""
def create_m0_map(hdu, channel_idx, center=None, size=None, sigma=3, distance=0, figure=plt.figure(figsize=(6,6)), subplot=(1,1,1), multiimage=False, use_offset_labels=False):
    # sigma
    blank_channel = hdu.data[0,0,:,:]
    region_size = np.array([100, 100]) * u.pixel
    wcs_2d = wcs.WCS(hdu.header).celestial
    if center == None:
        center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)
    region = Cutout2D(blank_channel, center, region_size, wcs=wcs_2d)
    mean = np.mean(region.data)
    std = np.std(region.data)

    # create map
    channels = hdu.data[0,channel_idx,:,:]
    channels[channels < (mean + sigma*std)] = 0
    m0_map = np.sum(channels, axis=0)

    cut = cut_fig(m0_map, hdu.header, center, size)

    if np.nanmax(cut.data) < 0.1:
        units = "mJy/bm km/s"
    else:
        units = "Jy/bm km/s"
    
    fig = create_fig(cut, distance=distance, figure=figure, subplot=subplot, multiimage=multiimage, use_offset_labels=use_offset_labels)

    fig.colorbar.set_axis_label_text(f"Intensity ({units})")
    return fig


"""

Marks sources given rows from 'source_info.csv'

"""
def mark_sources(fig, source_rows, use_short_label=False, fontsize=4):
    marker_colors = ["dodgerblue", "magenta", "chartreuse", "springgreen", "deepskyblue"]
    legend_handles = []
    sources_to_mark = source_rows.sort_values('Source').reset_index()[0:4]
    for i, row in sources_to_mark.iterrows():
        center2 = SkyCoord(row['RA'], row['Dec'], unit=u.degree)
        fig.show_markers(center2.ra.deg, center2.dec.deg, coords_frame='world', marker='x', s=50, c='k', linewidths=2, zorder=2, label=row['Source'])
        fig.show_markers(center2.ra.deg, center2.dec.deg, coords_frame='world', marker='x', s=50, c=marker_colors[i], linewidths=1, zorder=3, label=row['Source'])


        # Create legend handle for this source (only if not already added)
        short_label = str(row['Source']).casefold().removeprefix(str(source_rows.index.tolist()[i]).casefold()+'-').upper()
        legend_handles.append(mlines.Line2D([], [], color=marker_colors[i], marker='x', markersize=4, linestyle='None', label=short_label if use_short_label else row['Source']))

    fig.ax.legend(handles=legend_handles, fontsize=fontsize, loc='upper right', bbox_to_anchor=(0.98,1))


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
def plot_vector(fig, origin, angle_north_deg, color, length=0.0025):
    angle_east_rad = np.radians(90 - angle_north_deg)
    origin_pix = fig.world2pixel(origin[0], origin[1])
    tip = origin[0] + length * np.cos(angle_east_rad), origin[1] + length * np.sin(angle_east_rad)
    tip_pix = fig.world2pixel(tip[0], tip[1])
    outflow_vector = np.array([tip_pix[0] - origin_pix[0], tip_pix[1] - origin_pix[1]])
    fig.ax.quiver(origin_pix[0], origin_pix[1], outflow_vector[0], outflow_vector[1],
                angles='xy', scale_units='xy', scale=1, color=color, width=0.0075)
    
def plot_dotted_vector(fig, origin, angle_north_deg, color, length=0.0025):
    angle_east_rad = np.radians(90 - angle_north_deg)
    origin_pix = fig.world2pixel(origin[0], origin[1])
    tip = origin[0] + length * np.cos(angle_east_rad), origin[1] + length * np.sin(angle_east_rad)
    tip_pix = fig.world2pixel(tip[0], tip[1])
    outflow_vector = np.array([tip_pix[0] - origin_pix[0], tip_pix[1] - origin_pix[1]])
    fig.ax.plot(
        [origin_pix[0], origin_pix[0] + outflow_vector[0]], 
        [origin_pix[1], origin_pix[1] + outflow_vector[1]], 
        linestyle="dashed", color=color, linewidth=1  # Dotted line
    )

def plot_outflow_and_separation_vectors(fig, outflow_data, target_name, delta_pa_label=False):
    ### VECTORS
    # plot binary separation angle
    outflow = outflow_data.loc[outflow_data['field'] == target_name].groupby('field').first().reset_index()
    center_origin = np.array([np.mean([outflow['source_a_ra'], outflow['source_b_ra']]), np.mean([outflow['source_a_dec'], outflow['source_b_dec']])])
    separation_angle_north = outflow['binary_PA']
    # draw separation vector in both directions
    plot_vector(fig, center_origin, separation_angle_north, color='white')
    plot_vector(fig, center_origin, separation_angle_north + 180, color='white')
    
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
        plot_vector(fig, outflow_origin, outflow_angle_north, color='#00FFFF')

        # add delta PA label
        if delta_pa_label:
            angle = np.abs(outflow_angle_north - separation_angle_north) % 180
            angle = np.min([angle, 180 - angle])

            source_label = source['outflow_source']
            if source_label == 'both':
                source_label = source['source_a'] + '+' + source['source_b']
            fig.ax.text(30, fig.ax.get_xlim()[1] - 50 - 40 * (j), f"{source_label} : {np.abs(angle):.2f}°")

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