from PV_diagram_2 import *
import numpy as np
import sys
from astropy.coordinates import SkyCoord


"""
i have all the .fits files in the same directory. they follow a naming
format as such:
hops-XXX__sYY__ZZZ.fits

where
XXX: the sources number (e.g., hops-70, hops-400, hops-75, etc.)
YY : the spectral window (e.g., 11, 13, 15, 17, etc..)
ZZZ: the molecular line (e.g., c18o, h2co_218.2, ch3oh, so, etc.)

take that into account when editing this file.
"""

#######################################
# check these numbers for each source #
#######################################

# coords
coords_Aa   = SkyCoord('05h35m18.337s', '-05d00m32.94s', frame='icrs')
coords_Ab   = SkyCoord('05h35m18.328s', '-05d00m33.18s', frame='icrs')
coords_b    = SkyCoord('05h35m18.270s', '-05d00m33.92s', frame='icrs')

# rotation (needed for the loc_cut ipython script)
rot_Aa      = 140.8
rot_Ab      = 140.8   # not given, assuming Aa rot.
rot_b       = 174.1

# needs rot, v_sour, and maybe zoom?
source      = "hops-92"     # used to get the filename in my naming convention.
directory   = "cubes/"      # where to place the final images
rotation    = 174.1         # johns; degrees; Aa:140.8, Ab:..., B:174.1
mass        = 0.5           # Solar masses
mass_err    = 0             # Solar masses
v_source    = 0.0           # km s^-1; the slice in the cube that is blank
source_dist = 400.0         # pc
inclination = 60            # degrees
ra          = coords_b.ra.deg
dec         = coords_b.dec.deg
suffix      = "objB__rot{}".format(rotation)


########################
# specific lines stuff #
########################

# used to grab the right cubes. edit as your naming convention changes
file_suffix = ['s11__c18o.fits',            's13__13co.fits',
               's15__12co.fits',            's17__h2co_218.2.fits',
               's19__h2co_218.475.fits',    's21__h2co_218.76.fits',
               's23__so.fits',              's25__ch3oh.fits']

# used for the image
lines       = ['C$^{18}$O (2-1)',
               '$^{13}$CO (2-1)',
               '$^{12}$CO (2-1)',
               'H$_2$CO (3$_{0,3}$-2$_{0,2}$)',
               'H$_2$CO (3$_{2,2}$-2$_{2,1}$)',
               'H$_2$CO (3$_{2,1}$-2$_{2,0}$)',
               'SO',
               'CH$_3$OH' ]

# total amount of pixels over which we'll be summing. (thickness)
widths      = [15, 15, 25, 15, 15, 15, 15, 15]
v_sources   = [10.75, 11.00, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0]

# cosmetics
vel_widths  = [15, 15, 15, 15, 15, 15, 15, 15]
arc_widths  = [6.0, 6.0, 6.0,  6.0, 6.0, 6.0, 6.0, 6.0] # image zoom
cont_intrs  = [5.0, 5.0, 10.0, 5.0, 5.0, 5.0, 5.0, 5.0]


########
# code #
########

if __name__ == '__main__':

  for i in range(len(lines)):

      filename    = "{}__{}".format(source, file_suffix[i])
      width       = widths[i]
      line        = lines[i]
      v_source    = v_sources[i]

      vel_width   = vel_widths[i]
      arc_width   = arc_widths[i]
      cont_intr   = cont_intrs[i]

      print("\n")
      print("doing: {}".format(filename))
      print("width: {}".format(width))
      print("line:  {}".format(line))
      print("vsour: {}".format(v_source))
      print("vel_w: {}".format(vel_width))
      print("arc_w: {}".format(arc_width))
      print("c_lvl: {}".format(cont_intr))

      if v_source == "no":
          print("skipping")
      else:
          # doing the cut along the outflow direction
          PV_diagram(ra, dec, filename, directory, rotation, width, source, line,
                     inclination, Velocity_Curve=True, mass=mass,
                     mass_err=mass_err, v_source=v_source,
                     d_source=source_dist, Thindisk=False, Zoom=True,
                     v_width=vel_width, arcsec_width=arc_width,
                     Overlay_Contour='None', imagemin=-99.0, imagemax=-99.0,
                     contour_interval=cont_intr, along=True, suffix=suffix)

          # doing the cut normal to the outflow direction.
          PV_diagram(ra, dec, filename, directory, rotation, width, source, line,
                     inclination, Velocity_Curve=True, mass=mass,
                     mass_err=mass_err, v_source=v_source,
                     d_source=source_dist, Thindisk=False, Zoom=True,
                     v_width=vel_width, arcsec_width=arc_width,
                     Overlay_Contour='None', imagemin=-99.0, imagemax=-99.0,
                     contour_interval=cont_intr, along=False, suffix=suffix)

  sys.exit()
