import astropy.io.fits as fits
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
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
        fig.add_scalebar((100 / (distance)) * u.arcsecond)
        fig.scalebar.set_label('100 AU')

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

def create_m8_map(hdu, center, size, distance=0):
    # create map
    m8_map = np.max(hdu.data[0,:,:,:], axis=0)

    # get the wcs
    wcs_original = wcs.WCS(hdu)
    wcs_2d = wcs_original.celestial

    # make the cut
    cut = Cutout2D(m8_map, center, size, wcs=wcs_2d)

    # create a header with the wcs
    # CAUTION: lots of header data not being preserved here
    new_header = cut.wcs.to_header()
    for key in ['BUNIT', 'OBJECT', 'TELESCOP', 'INSTRUME', 'BMAJ', 'BMIN', 'BPA']:
        if key in hdu.header:
            new_header[key] = hdu.header[key]

    # combine cut and wcs and an HDU
    hdu_cut = fits.PrimaryHDU(data=cut.data, header=new_header)
    
    return create_fig(hdu_cut, distance)

