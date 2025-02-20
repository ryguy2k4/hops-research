#!/usr/bin/python3

"""
open casa from the dir this file is in and run it from in
there. it'll make a new directory called "bnr" and place
the blue.fits and red.fits in there. you give it a range
of slices for each side. does them all in one go but doesnt
have to.
"""


import os
import sys
import glob
import numpy as np
import astropy as ap
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.constants import c as c

################
# user options #
################

"""
s#	spw		molecule
11	25		c18o
13	27		13co
15	37		12co
17	29		h2co
19	31		h2co
21	33		h2co
23	35		so
25	31		ch3oh
27	39		??
"""

# which molecules you want to do this for. see the table above.
#slist  = [11, 13, 15, 17, 19, 21, 23, 25, 27]
slist = [11, 13, 15, 17, 19, 21, 23]

# will overwrite any previous blue and red images made.
overwrite = True

# if you want to mess with the chan selection, there's a
# function down below for that.

##########
##########
## code ##
##########
##########

#######
# def #
#######

def molist(s):
  if s == 11:
    return "c18o"
  elif s == 13:
    return "13co"
  elif s == 15:
    return "12co"
  elif s == 17:
    return "h2co_218.2"
  elif s == 19:
    return "h2co_218.475"
  elif s == 21:
    return "h2co_218.76"
  elif s == 23:
    return "so"
  elif s == 25:
    return "ch3oh"
  elif s == 27:
    return ""

def bNr(s):
  if s == 11:
    return ("65~68", "72~75")
  elif s == 13:
    return ("63~67", "73~77")
  elif s == 15:
    return ("96~101", "105~113")
  elif s == 17:
    return ("64~68", "72~75")
  elif s == 19:
    return ("68~70", "73~75")
  elif s == 21:
    return ("68~70", "72~74")
  elif s == 23:
    return ("68~70", "72~74")
  else:
    return (None, None)


################
# dir checking #
################

if os.path.isdir("bnr"):
  print("dir 'bnr' exists")
else:
  os.system("mkdir bnr")
  print("made dir 'bnr'")


######################################
# making the blue and red fits files #
######################################

for i,s in enumerate(slist):

  imagename = glob.glob("data_reduction/oussid.s{}_*cube*.image".format(s))

  #print("using file: {}".format(imagename[0]))

  bluchan, redchan = bNr(s)

  ########
  # blue #
  ########

  if bluchan != None:
    newimage  = "bnr/blu__s{}__{}__{}".format(s, molist(s), bluchan)
    fitsimage = "{}.fits".format(newimage)

    if overwrite:
      try:
        os.system("rm -rf {}".format(newimage))
        print("removed old {} (if it existed)".format(newimage))
      except:
        pass

    try:
      immoments(imagename=imagename[0], moments=0, chans=bluchan,
                outfile=newimage)
      print("made the {} image".format(molist(s)))

      exportfits(imagename=newimage, fitsimage=fitsimage, overwrite=overwrite)
      print("made the {} fits file".format(molist(s)))


    except:
      print("couldn't make fits file. try overwrite=True")

  #######
  # red #
  #######

  if redchan != None:
    newimage  = "bnr/red__s{}__{}__{}".format(s, molist(s), redchan)
    fitsimage = "{}.fits".format(newimage)

    if overwrite:
      try:
        os.system("rm -rf {}".format(newimage))
        print("removed old {} (if it existed)".format(newimage))
      except:
        pass

    try:
      immoments(imagename=imagename[0], moments=0, chans=redchan,
                outfile=newimage)
      print("made the {} image".format(molist(s)))

      exportfits(imagename=newimage, fitsimage=fitsimage, overwrite=overwrite)
      print("made the {} fits file".format(molist(s)))

    except:
      print("couldn't make image. try overwrite=True or something.")

  ########
  # none #
  ########

  if bluchan != None or redchan != None:
    print("finished: {}".format(molist(s)))
  else:
    print("skipping: {}".format(molist(s)))

print("done")
