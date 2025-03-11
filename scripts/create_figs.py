import astropy.io.fits as fits
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import aplpy
import astropy.wcs.wcs as wcs

def create_fig(img, distance=0):
    # Create FITSFigure
    figure = plt.figure(figsize=(6,6))
    fig = aplpy.FITSFigure(img, figure=figure)

    # Display image
    fig.show_colorscale(cmap='viridis')

    # Title
    fig.set_title(img.header['OBJECT'])

    # Colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text("Intensity (Jy/Beam)")
    fig.colorbar.set_location('right')
    fig.colorbar.set_pad(0.0)

    # Scalebar
    if distance == 0:
        fig.add_scalebar(1 * u.arcsecond)
        fig.scalebar.set_label('1"')
    else:
        fig.add_scalebar((500 / (distance)) * u.arcsecond)
        fig.scalebar.set_label('500 AU')

    fig.scalebar.set_linewidth(3)

    # Beam
    fig.add_beam()
    fig.beam.show()
    fig.beam.set_corner("bottom left")
    fig.beam.set_color("black")
    fig.beam.set_edgecolor("black")
    fig.beam.set_facecolor("white")

    return fig

def create_sub_fig(hdu, center, size, distance=0):
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
    
    return create_fig(hdu_cut, distance)

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

def create_m8_map(hdu, center, size, distance=0):
    # create map
    m8_map = np.max(hdu.data[0,:,:,:], axis=0)

    cut = cut_fig(m8_map, hdu.header, center, size)
    
    return create_fig(cut, distance)

def create_m0_map(hdu, center, size, channel_idx, sigma=3, distance=0):
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
    
    return create_fig(cut, distance)

def mark_sources(fig, source_rows):
    marker_colors = ['black', 'red', 'magenta', 'darkred', 'darkblue']
    legend_handles = []
    sources_to_mark = source_rows.reset_index()[0:4]
    for i, row in sources_to_mark.iterrows():
        center2 = SkyCoord(row['RA'], row['Dec'], unit=u.degree)
        fig.show_markers(center2.ra.deg, center2.dec.deg, coords_frame='world', marker='x', s=25, c=marker_colors[i], linewidths=1, label=row['Source'])

        # Create legend handle for this source (only if not already added)
        legend_handles.append(mlines.Line2D([], [], color=marker_colors[i], marker='x', markersize=6, linestyle='None', label=row['Source']))
    fig.ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1.05,1+0.07*np.min(len(sources_to_mark))))