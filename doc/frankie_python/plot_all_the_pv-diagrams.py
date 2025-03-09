#!/usr/bin/python3

import os
import sys
import glob
from PIL import Image
import numpy as np
import astropy as ap
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.constants import c as c

#########
# usage #
########

# run it in the dir where it is to be used.


##################
# file selection #
##################

files = sorted(glob.glob("*"))


#############
# user info #
#############

# assuming you know what you're doing, so no error catching.

# shape of the subplots to make
plot_x = 4
plot_y = 4


########
# code #
########

# starting the plotting.
fig, axs = plt.subplots(plot_x, plot_y, sharex='col', sharey='row')

for ax,file in zip(axs.ravel(),files):

  print("opening {}".format(file))
  data = Image.open(file)

  ax.imshow(data)
  ax.axis('off')

#plt.tight_layout()
#plt.show()
plt.savefig("master_view.png", dpi=900)
sys.exit()
