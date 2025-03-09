import numpy as np
#import pyfits as pf
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import interpolation as inter
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io import fits as pf
from astropy.wcs import WCS
import math
#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()

# last updated: 22 nov 2019; fje

def Rotate_and_Sum(Filename,Rotated_Image, Rotation, Width, ra, dec):
    '''
    This function rotates an given image with an angle Rotation and then
    sums the pixels in y-direction in the range
    Ny/2-width/2:Ny/2+width/2, with Ny the length of the array in the y
    direction. It then returns an 1D array.
    Image = 2D array
    Rotation = Is defined positive from the x-axis to the y-axis
    (anti-clockwise), must be given in degrees.
    Width = Total width over wich there will be summed, must be given in
    pixels.
    '''

    Shape_Rotated_Image = np.shape(Rotated_Image)

    int1 = int(np.round(Shape_Rotated_Image[1]/2))
    int2 = int(np.round(Width/2))

    # Summed_Array = np.sum(Rotated_Image[(np.round(Shape_Rotated_Image[1]/2) - np.round(Width/2)):(np.round(Shape_Rotated_Image[1]/2) + np.round(Width/2))],0)
    Summed_Array = np.nansum(Rotated_Image[(int1 - int2):(int1 + int2)],0)

    temp = np.zeros([1,Shape_Rotated_Image[1]])
    temp[0,:] = Summed_Array

    #print("Summed_Array:\n{}".format(Summed_Array))
    return temp

def PV_diagram(ra,dec, File, Directory, rotation_deg, width, Object, Line,
               inclination, Velocity_Curve = False, mass = 0, mass_err = 0,
               v_source = 0, d_source = 0, Thindisk = False, Zoom = False,
               v_width = 5, arcsec_width = 20, Overlay_Contour = 'None',
               imagemin=-99.0, imagemax=-99.0, imagecut=None, along=True,
               contour_interval=10.0, savepdf=False, saveeps=False,
               saveps=True, savepng=True, suffix=""):
    '''
    This function creates an position velocity diagram for a given fits
    file.
    The rotation that has to be performed must be given in
    degrees. The width is given in pixels and is assumed to be around
    the middle. If Velocity_Curve is True then a velocity curve based on
    the given mass and mass error is produced. Both should be given in
    solar masses.

    - File = Name of the 3D datacube. [string]
    - Directory = Directory where the file is saved.  [string]
    - rotation_deg = The rotation in degrees that has to be performed to
                     align the axis over which there will be summed with
                     the y-axis. [float]
    - width = The total amount of pixels over which will be summed. [float]
    - Object = The name of the object, only used in title figure. [string]
    - Line = The molecular line of the data, only used in title figure. [string]
    - Velocity_Curve = If true, a keplerian rotation curve will be drawn
                       in the figure based on the parameters explained
                       below. [boolean]
    - Mass = The mass of the central object in solar masses, used for
             the velocity curve. [float]
    - Mass_Err = The error on the mass in solar masses, used for the
                 velocity curve. [float]
    - v_source = The source velocity in km s^-1, used for the velocity
                 curve and the zoom function.[float]
    - d_source = The distance of the source in parsec, used for the
                 velocity curve. [float]
    - Thindisk = Set to true if the data is from a thindisk model, so it
                 wont draw any contour lines based on the sigma's in the
                 image. [boolean]
    - Zoom = If true the function will zoom in. [boolean]
    - v_width = The total range in velocities (km s^-1) that will be shown
                when zoomed in. [float]
    - arcsec_width = The total range in arcsec that will be shown when
                     zoomed in. [float]
    - Overlay_Contour = A contour will be drawn over the figure of a
                        different dataset, which already must have a
                        PV-diagram made by this programm. If 'None', no
                        contour will be drawn. If you want to
                        add these contours, you have to line 191 to manually
                        add the directories and filenames of the data. [string]
    - imagecut = since scipy.ndimage.interpolate.rotate doesn't like an image
                 to have any NaNs, we must first cut the image a bit. this
                 should be a list, [x1cut, x2cut, y1cut, y2cut], where
                 x1,y1 are the first axe of the image and x2,y2 are the second.
                 in the data itself, the cut comes from Data[:,x1:x2,y1:y2]
    - along = whether we do all this business along the y-axis, which should
              be along the outflow, or if we do this perpendicular to the
              outflow. default is True.
    - suffix = some extra text to help distinguish the files. used only
               on the product filenames.
    '''

    font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 16}

    mpl.rc('font', **font)
    #---------------------------------------------------------------------------
    #Getting the data from the fits file.
    #---------------------------------------------------------------------------

    #getting the header.
    header = pf.getheader(Directory + File)

    #Getting the data.
    data = pf.getdata(Directory + File)

    # the datacubes coming from the thindisk model have 3 dimensions,
    # while the science datacubes have 4 dimension. So we have to account
    # for that.
    if len(np.shape(data)) == 4:
        Data = data[0,:,:,:]
    else:
        Data = data

    # hard-coding the zoom because scipy.ndimage.interpolate.rotate doesn't
    # play well with NaN's. so we zoom in to a value that cuts out most
    # all the NaN's and hopefully not too much of the stuff we want.
    if imagecut is not None:
      Data = Data[:, imagecut[0]:imagecut[1], imagecut[2]:imagecut[3]]

    # if the image has any NaNs, let's remove it.
    if np.isnan(Data).any():
        print("replacing NaNs with 0's")
        Data[np.isnan(Data)]    = 0
        undoNaNs                = True
    else:
        undoNaNs                = False

    #Determining the shape of the data.
    Shape_Data = np.shape(Data)

    #Determining the begin_v, delta_v and the amount of images N.
    N = header['NAXIS3']
    if header['CTYPE3']  == 'VELO-LSR':
        print('True')
        begin_v = header['LSTART']
        delta_v = header['LWIDTH']
    elif header['CTYPE3']  == 'FREQ':
        #Reading the data in frequencies. We have to bring this to velocities.
        begin_freq = header['CRVAL3']
        delta_freq = header['CDELT3']
        begin_pos  = header['CRPIX3'] - 1
        rest_freq  = header['RESTFRQ']

        #The speed of light is.
        c = const.c.value / 1000.0  #km s^-1

        # some temp variables to format the code better
        begin_v_numer  = rest_freq**2 - (begin_freq - delta_freq*begin_pos)**2
        begin_v_denom  = rest_freq**2 + (begin_freq - delta_freq*begin_pos)**2
        begin_v_numer1 = rest_freq**2 - (begin_freq-delta_freq*(begin_pos+1))**2
        begin_v_denom1 = rest_freq**2 + (begin_freq-delta_freq*(begin_pos+1))**2

        #Calculating the begin velocity.
        begin_v = c * begin_v_numer / begin_v_denom

        #Now we calculate the delta v
        begin_v_plus_one  = c * begin_v_numer1 / begin_v_denom1
        delta_v           = np.round(begin_v - begin_v_plus_one, 2)
        delta_v           = begin_v - begin_v_plus_one
        vgrid             = np.arange(begin_v,
                                      begin_v + delta_v*float(Shape_Data[0]),
                                      delta_v)
    else:
        print('I am not sure how to get the velocities.')
        raise KeyError

    PixelWidth_RA   = header['CDELT1']
    PixelWidth_DEC  = header['CDELT2']

    #The used units.
    Image_Units = header['BUNIT']

    #Creating the directory for saving the images, if it does not
    #already exists.
    save_directory = 'pv_diagrams/'
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    #---------------------------------------------------------------------------
    #Creating the array that contains the PV-diagram
    #---------------------------------------------------------------------------

    # if you rotated using ds9 or casa, and you got the outflow along the
    # y-axis, the that's not the real y-axis, but x-axis. so we'll fix that.
    if along:
        rotation_deg += 90
    rotation_rad = np.radians(rotation_deg)

    #Creating the array of the correct size.
    abs_rot = np.abs(rotation_rad)
    y_size  = int(np.round(np.abs(np.cos(abs_rot) * Shape_Data[1]) + \
                           np.abs(np.cos(np.pi/2 - abs_rot) * Shape_Data[2])))
    PV_Data = np.zeros([y_size, Shape_Data[0]])

    print("file loc:       {}".format(Directory+File))
    print("data shape:     {}".format(Data[:,:,:].shape))

    if imagecut is not None:
      imxcen = (imagecut[1] - imagecut[0])/2
      imycen = (imagecut[3] - imagecut[2])/2
    else:
      imxcen = Data[0,:,:].shape[0]/2
      imycen = Data[0,:,:].shape[0]/2

    w     = WCS(Directory+File)

    print("ra,dec:         {},{}".format(ra,dec))
    print("rotation:       {}".format(rotation_deg))

    xcen, ycen = w.sub(['longitude','latitude']).wcs_world2pix(ra, dec,0)
    if imagecut is not None:
      xcen = xcen - imagecut[0]
      ycen = ycen - imagecut[2]

    print("xcen, ycen:     {},{}".format(xcen, ycen))

    center_shift = [0, imycen-ycen, imxcen-xcen]

    print("cent shift:     {}".format(center_shift))

    '''

    if center_shift[0] >0:
       center_shift[0]=center_shift[0]*(-1.0)
    if center_shift[1] >0:
       center_shift[1]=center_shift[1]*(-1.0)
    print(center_shift)
    '''

    print("abs cent shift: {},{}".format(abs(center_shift[1]),
                                         abs(center_shift[2])))

    if (abs(center_shift[1]) > 0.5) or (abs(center_shift[2]) >0.5):
       print('Shifting Image to Center')
       Shifted_Image=inter.shift(Data, center_shift)

       print('Rotating Image')
       Rotated_Image = inter.rotate(Shifted_Image, rotation_deg, axes=(1,2),
                                    cval=np.nan)
       if undoNaNs:
           print("changing back the 0's to NaN's")
           Rotated_Image[Rotated_Image == 0] = np.nan
    else:
       print('_NOT_ Shifting Image')
       print('Rotating Image')
       Rotated_Image = inter.rotate(Data, rotation_deg, axes=(1,2), cval=np.nan)

    i = 0
    while i < Shape_Data[0]:
        PV_Data[:,i] = Rotate_and_Sum(Directory + File, Rotated_Image[i,:,:],
                                      rotation_deg, width, ra, dec)
        #print("std of current pv_data: {}".format(np.std(PV_Data)))
        i += 1

    print("Pv_Data shape: {}".format(PV_Data.shape))
    #---------------------------------------------------------------------------
    #Contour lines
    #---------------------------------------------------------------------------
    std_PV = np.nanstd(PV_Data[0:PV_Data.shape[0],0:Shape_Data[0]])
    print("std_PV: {}".format(std_PV))

    PV_Contour_Levels = np.arange(1, 50, 1)*contour_interval*std_PV
    print("contour levels: \n{}".format(PV_Contour_Levels))
    #PV_Contour_Levels = [1000.0,2000.0]

    #---------------------------------------------------------------------------
    #Plotting
    #---------------------------------------------------------------------------
    #First we create arrays for the correct ticks
    x_values = np.arange(begin_v, begin_v + delta_v*float(Shape_Data[0]), delta_v)

    #The total length in arcsec of the y axis in the new image.
    length_arcsec_new =  (np.abs(np.cos(np.abs(rotation_rad)))*np.abs(PixelWidth_RA)*Shape_Data[2]+np.abs(np.cos(np.pi/2.0-np.abs(rotation_rad)))*np.abs(PixelWidth_DEC)*Shape_Data[1])*3600
    y_values = np.arange(-length_arcsec_new/2.0, length_arcsec_new/2.0 + length_arcsec_new/10.0, length_arcsec_new/10.0)

    #Calculating the size of y pixel in the y direction in arcsec.
    pixel_size_y_arcsec = length_arcsec_new/y_size

    y_arcsec = np.arange(1, PV_Data.shape[0])*pixel_size_y_arcsec - \
                         length_arcsec_new/2.0


    ##################
    # start plotting #
    ##################

    print('Starting the plotting.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if imagemin !=-99.0 and imagemax !=-99.0:
       cax = ax.imshow(PV_Data, origin='lower', cmap='magma',
                       vmin=imagemin, vmax=imagemax, interpolation='nearest')
    else:
       cax = ax.imshow(PV_Data, origin='lower', cmap='magma',
                       interpolation='nearest')
       #cax=ax.pcolormesh(vgrid,y_arcsec,PV_data, cmap='magma')

    #The thindisk model does not contain noise, so no contourlines are needed.
    if not Thindisk:
        ax.contour(PV_Data, PV_Contour_Levels, colors='white',
                   linewidths=0.5)
        print("adding contours")

    def Make_Contour(Directory, File, Color, Line, Thindisk = False):
        '''
        This file opens the file given and plots the contours of this
        file.
        '''
        Data = pf.getdata(Directory + File)

        if Thindisk:
            #If we have thindisk model data, we have no noise. So we just
            #use percentages instead.
            Max = np.max(Data)
            Data_Contour_Levels = np.arange(0.2, 1, 0.2)*Max
        else:
            std_Data = np.std(Data[0:100,0:100])
            Data_Contour_Levels = np.arange(0, 30, 1)*std_Data
            Data_Contour_Levels[0] = 3*std_Data

        ax.contour(Data, Data_Contour_Levels, colors = Color)
        ax.axvline(-10000, color = Color, label = Line)
        #ax.legend(loc = 3)

    if Overlay_Contour == '13CO':
        Directory = '13CO_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-13CO_1.fits'
        Make_Contour(Directory, File_Name, 'yellow', '$^{13}$CO')
        ax.legend(loc = 3)
    if Overlay_Contour == 'C18O':
        Directory = 'C18O_2/PV-diagram/'
        File_Name = 'PV-Diagram_L1527_C18O_2.fits'
        Make_Contour(Directory, File_Name, 'white', 'C$^{18}$O')
        #ax.legend(loc = 3)
    if Overlay_Contour == 'Model':
        Directory =  'Optimization_Results/Round_3_d_414/Model_Best_38.3/PV-diagram/'
        File_Name = 'PV-Diagram_Optimization_Best_38.3_2_convolved_with_beam.fits'
        Make_Contour(Directory, File_Name, 'white', 'Model', Thindisk = True)
        #ax.legend(loc = 3)
    if Overlay_Contour == 'Both':
        Directory = 'C18O_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-C18O_1.fits'
        Make_Contour(Directory, File_Name, 'yellow', 'C$^{18}$O')

        Directory = '13CO_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-13CO_1.fits'
        Make_Contour(Directory, File_Name, 'red', '$^{13}$CO')

    plot_title = "{} {}".format(Object, Line)

    ax.set_title(plot_title)
    ax.set_xlabel('Velocity (km s$^{-1}$)')
    ax.set_ylabel('Offset ($^{\prime\prime}$)')#'Position [arcsec]')$\Delta$X

    offset = 0.0
    #If we have correct masses we can calculate the velocity curve.
    if Velocity_Curve:
        print(' the velocities curves.')
        #Calculating the extreme masses within the errors, do we can also plot
        #those.
        mass_min_err = mass - mass_err
        mass_plus_err = mass + mass_err

        #Creating an array containing the velocities in km s^-1.
        delta_v_plot=delta_v/10.0
        velocities = (np.arange(begin_v, begin_v + delta_v*float(Shape_Data[0]),
                                delta_v_plot) - v_source)
        print("printing the array containing the velocities")
        print("size of the array: {}".format(velocities.shape))
        print("velocities: {}".format(velocities))

        #This function returns for a given mass (solar masses), velocity (in
        #km s^-1) and distance to the source (in pc) the radius (in arcsec)
        #assuming Keplerian rotation.
        def Keplerian_Rotation(mass, velocity, Distance, inclination):
            radii_return =  np.sin(inclination*3.14/180.0)**2*const.G.value*mass*const.M_sun.value/(velocity*1000)/(velocity*1000)/(Distance*u.pc.to(u.m))*u.rad.to(u.arcsec)

            #All the positive radii.
            radii_positive = radii_return[velocity < 0]
            #We also have some negative radii, so thats why we have to do this.
            radii_negative = -1*radii_return[velocity > 0]
            return radii_positive, radii_negative

        #Calculate the radii.
        _tmp1, _tmp2 = Keplerian_Rotation(mass, velocities,
                                          d_source, inclination)
        radii_positive, radii_negative = _tmp1, _tmp2

        _tmp1, _tmp2 = Keplerian_Rotation(mass_min_err, velocities,
                                          d_source, inclination)
        radii_positive_min_err, radii_negative_min_err = _tmp1, _tmp2

        _tmp1, _tmp2 = Keplerian_Rotation(mass_plus_err, velocities,
                                          d_source, inclination)
        radii_positive_plus_err, radii_negative_plus_err = _tmp1, _tmp2

        #print(y_arcsec)
        #print(len(y_arcsec))
        #print(len(radii_positive))
        #print(len(radii_negative))
        #print(len(velocities))
        #print(radii_positive)
        #print(radii_negative)
        #Changing the radii to the correct pixel coordinates for correct
        #plotting. Plus bring the lines to the object.
        radii_positive_pixel_coor = radii_positive/pixel_size_y_arcsec + \
                                    (y_size - 1.0)/2.0 + \
                                    offset
        radii_negative_pixel_coor = radii_negative/pixel_size_y_arcsec + \
                                    (y_size - 1.0)/2.0  + \
                                    offset

        radii_positive_min_err_pixel_coor = radii_positive_min_err / \
                                            pixel_size_y_arcsec + \
                                            (y_size - 1.0)/2.0 + offset
        radii_negative_min_err_pixel_coor = radii_negative_min_err / \
                                            pixel_size_y_arcsec + \
                                            (y_size - 1.0)/2.0 + offset

        radii_positive_plus_err_pixel_coor = radii_positive_plus_err / \
                                             pixel_size_y_arcsec + \
                                             (y_size - 1.0)/2.0 + offset
        radii_negative_plus_err_pixel_coor = radii_negative_plus_err / \
                                             pixel_size_y_arcsec + \
                                             (y_size - 1.0)/2.0 + offset

        #print(radii_positive_pixel_coor)
        #print(radii_negative_pixel_coor)
        #Plotting the velocities
        #print(len(radii_positive_pixel_coor))
        #print(len(np.arange(0,len(radii_positive)/10.0, 1.0/10.0)))

        #print(len(radii_negative_pixel_coor))
        #print(len(np.arange(len(radii_positive)/10.0 , (len(velocities)-0.5)/10.0, 1.0/10.0)))

        length_negative_v       = len(radii_negative_pixel_coor)

        xaxis_negative_v        = np.arange(len(radii_positive)/10.0,
                                            (len(velocities))/10.0,
                                            1.0/10.0)
        length_xaxis_negative_v = len(xaxis_negative_v)

        if length_negative_v > length_xaxis_negative_v:
            xaxis_negative_v  = np.arange(len(radii_positive)/10.0,
                                          (len(velocities)+0.5)/10.0,
                                          1.0/10.0)
        elif length_negative_v < length_xaxis_negative_v:
            xaxis_negative_v  = np.arange(len(radii_positive)/10.0,
                                          (len(velocities)-0.5)/10.0,
                                          1.0/10.0)

        xaxis_pos_v         = np.arange(0,len(radii_positive)/10.0,
                                        1.0/10.0)
        length_pos_v        = len(radii_positive_pixel_coor)
        length_xaxis_pos_v  = len(np.arange(0,len(radii_positive)/10.0,
                                            1.0/10.0))
        if length_pos_v > length_xaxis_pos_v:
            xaxis_pos_v = np.arange(0,(len(radii_positive)+0.5)/10.0, 1.0/10.0)
        elif length_pos_v < length_xaxis_pos_v:
            xaxis_pos_v = np.arange(0,(len(radii_positive)-0.5)/10.0, 1.0/10.0)


        ax.plot(xaxis_pos_v, radii_positive_pixel_coor, color='white',
                label='Keplerian rotation', linestyle=':')
        ax.plot(xaxis_negative_v, radii_negative_pixel_coor, color='white',
                linestyle=':')

        ax.axhline(np.where(y_arcsec > 0)[0][0] - 0.5 + offset,
                   color='white', linestyle='--')

        hloc=float(np.where(velocities > 0)[0][0])/10.0
        ax.axvline(hloc -0.15 , color='white', linestyle='--',
                   label='v$_{source}$ = ' + str(v_source) + ' km s$^{-1}$')
        #ax.legend(loc = 3)

    #If zoom is true we zoom in else we just show the entire image.
    if Zoom:
        #First we determine at which pixel we have the source velocity.
        pix_v_source = float(np.abs((begin_v - v_source)/delta_v))
        print("pix with source vel: {}".format(pix_v_source))
        #Then we determine what half the width of the v slice must be.
        pix_v_shift = float(v_width/delta_v/2.0)
        print("half the pix width: {}".format(pix_v_shift))
        print("vel width: {}".format(v_width))
        print("delta v: {}".format(delta_v))
        #Now we determine the central pixel for the arcsec.
        pix_arcsec_central = float(y_size/2.0) - 1.0 + float(offset)
        pix_arcsec_shift = float(arcsec_width/pixel_size_y_arcsec/2.0)

        x = np.arange(pix_v_source - pix_v_shift,
                      pix_v_source + pix_v_shift + 2.0*pix_v_shift/10.0,
                      2.0*pix_v_shift/10.0)
        y = np.linspace(pix_arcsec_central - pix_arcsec_shift,
                        pix_arcsec_central + pix_arcsec_shift, 11.0)

        x_label = np.round(begin_v + delta_v*x, 1)
        y_label = np.linspace(-arcsec_width/2.0, arcsec_width/2.0, 11.0)




        from matplotlib.ticker import FormatStrFormatter
        #ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        y_label = np.array(["%.1f" % i for i in y_label])
        ax.set_xticks(x)
        ax.set_yticks(y)
        ax.set_xticklabels(x_label)
        ax.set_yticklabels(y_label)

        # doing the zooming by limiting the shown x and y coordinates.
        # move after the labels, because those functions automatically
        # adjust x and y limits.
        print(pix_v_source - pix_v_shift,
              pix_v_source + pix_v_shift,
              pix_v_source)
        print('!')
        ax.set_xlim([pix_v_source - pix_v_shift,
                     pix_v_source + pix_v_shift])
        ax.set_ylim([pix_arcsec_central - pix_arcsec_shift,
                     pix_arcsec_central + pix_arcsec_shift])
        print(x)
        ax.set_aspect(1.0*pix_v_shift/pix_arcsec_shift)

    else:
        x = np.arange(0,Shape_Data[0] + Shape_Data[0]/10, Shape_Data[0]/10)
        y = np.arange(0,y_size + np.round(y_size/10), np.round(y_size/10))
        x_label = np.round(x_values[0::Shape_Data[0]/10], 1)
        y_label = np.round(y_values,1)
        xticks=np.linspace(-4,1,11)
        ax.set_xticks(xticks)
        ax.set_yticks(y)
        #ax.set_xticklabels(x_label)
        ax.set_yticklabels(y_label)

        ax.set_xlim([0, Shape_Data[0]-1])
        ax.set_ylim([0,y_size-1])
        ax.set_aspect(Shape_Data[0]/y_size)
        #ax.set_aspect((delta_v*(Shape_Data[0] - 1) )/length_arcsec_new)


    #ax.grid(color='white', alpha=0.5, linestyle='solid')
    cbar = plt.colorbar(cax)
    cbar.set_label('Jy/Beam')

    #ax.set_xticks(color='white', exclude_overlapping=True)
    #ax.set_yticks(color='white', exclude_overlapping=True)

    fileprefix = "{}{}{}_{}".format(save_directory, "pv-diagram",
                                     File[:-5], suffix)
    if along:
      fileprefix = "{}_outflow".format(fileprefix)

    if savepdf:
        plt.savefig(fileprefix + '.pdf', bbox_inches='tight')
    if savepng:
        plt.savefig(fileprefix + '.png', bbox_inches='tight', dpi=400)
    if saveps:
        plt.savefig(fileprefix + '.ps',  bbox_inches='tight')
    if saveeps:
        plt.savefig(fileprefix + '.eps', bbox_inches='tight')

    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(PV_Data)
    HDR = HDU.header
    #Saving the correct header data.

    #Need to create header.
    HDU.writeto(fileprefix + '.fits', overwrite = True)

    #---------------------------------------------------------------------------
    #Saving data in .txt file
    #---------------------------------------------------------------------------
    np.savetxt(fileprefix + '.txt', PV_Data, delimiter = ' ')

if __name__ == '__main__':
    rotation = 0.0 #degrees

    Mass = 0.7#Solar masses
    Mass_err = 0#Solar masses
    v_source = 4.8#5#km s^-1
    distance_source = 300.#pc
    inclination = 70.

    #-------------------------------------------------------------------------------
    #C180
    #-------------------------------------------------------------------------------
    Object = 'BHR7'
    Line = 'C$^{18}$O'
    Line = '$^{13}$CO'
    width = 21#pixels
    File = 'BHR7-13CO.fits'
    Directory = 'EXT2/13CO/'


    PV_diagram(File, Directory, rotation, width, Object, Line, inclination, Velocity_Curve = True, mass = Mass, mass_err = Mass_err, v_source = v_source, d_source = distance_source, Thindisk = False, Zoom = True, v_width = 10, arcsec_width = 10,Overlay_Contour = 'None')

