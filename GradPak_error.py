#! /usr/bin/python
#
##################################################
#
#
# This script is used to produce error images by propagating a raw error
# (.sig) file through the standard GradPak pipeline. The raw error image is
# typically something produced by the IRAF task rawimerr, which is part of
# Matt Bershady's ifupkg.
#
# The calling syntax is easy:
#
# > python GradPak_error.py RAW_IMAGE.sig.fits
#
# The result of which will be a file called RAW_IMAGE.me_rf_lin.fits.
#
# The reduction steps that are performed, along with the suffix produced by
# that step are:
#
# 1. Aperture extraction and flat-field calibration using DOHYDRA (.me.fits)
# 2. Dispersion correction using DISPCOR (.me_lin.fits)
# 3. Flux calibration using CALIBRATE (.me_rf_lin.fits)
#
# Note that sky subtraction is not performed (see below).
#
# For this script to work properly all the parameters for DOHYDRA and
# CALIBRATE must be set (in IRAF) exactly as they were when the raw data
# images were reduced. In other words, the values of, e.g., dohydra.flat,
# dohydra.apidtab, or calibrate.sensiti must be the same as when your data
# were reduced. The typical usage is to reduce your GradPak data through flux
# calibration and then simply run rawimerr followed by this script to produce
# a final error image in two easy steps.
#
# Assumptions and Simplifications:
#
# During the error propagation a few simplifying assumptions are made:
#
# 1. The flat error is negligible compared to the object error. This
#    assumption is strengthened by the fact that during flat-field corrections
#    the error term from the flat itself depends on the inverse of the square
#    of the flat signal and can therefore be thought of as a second order
#    correction. This assumption is consistent the ifupkg.mkmes task.
#
# 2. The dispersion solution is essentially linear. This allows a much simpler
#    treating of the error propagation resulting from the DISPCOR task and
#    should be valid for most GradPak spectra. This assuption makes this
#    script less accurate than ifupkg.dispcor_err, but the difference should
#    be negligable.
#
# 3. The error is not propagated through the sky subtraction step. This is
#    akin to assuming the combined sky fibers have negligable error compared
#    to the object fibers and is probably the most egregious assumption made
#    by this scipt. It will probably be fixed sometime in the future.
#
# History:
#      v1 - A. Eigenbrot Jan. 2015
#
##################################################

import glob
import sys
import os

import numpy as np
import pyfits

#Load the IRAF packages we'll need
try:
    current_dir = os.getcwd()
    if os.getlogin() == 'Arthur':
            os.chdir('/Users/Arthur/Ureka/iraf/local')
    from pyraf import iraf
    os.chdir(current_dir)
    iraf.imred(_doprint=0)
    iraf.hydra(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.onedspec(_doprint=0)
except Exception as e:
    print "Failure: could not find pyraf/iraf"
    sys.exit(1)

glob = glob.glob

if os.getlogin() == 'Arthur':
    APIDTABLE = '/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf'
else:
    APIDTABLE = '/usr/users/eigenbrot/research/Pak/gradpak_sizes.iraf'


def pow_image(inputname, outputname, power):
    '''
    Take an input image, raise all of its pixesl to the specified power and
    write the output image.
    '''
    print 'raising {} to power of {}'.format(inputname,power)
    
    h = pyfits.open(inputname)[0]
    
    pyfits.PrimaryHDU(h.data**power,h.header).writeto(outputname)
    return

def create_tmps(errname, flatname):
    '''
    Set up some temporary, squared version of the raw error and flat
    files. These are then used as inputs to dohydra.
    '''
    sq_errname = 'tmp_sq{}'.format(errname)
    sq_flatname = 'tmp_sq{}'.format(flatname)

    pow_image(errname,sq_errname,2.)
    pow_image(flatname,sq_flatname,2.)
    
    return sq_errname, sq_flatname

def find_msname(rawname):
    '''
    Given a raw file (that was input to dohydra), find the name that dohydra
    gave its processed result. This is used to get to correct naming
    convention for the master flat as the specific type of IRAF scrunching can
    vary from system to system.
    '''
    flist = glob('{}*.ms.fits'.format(rawname.split('.fits')[0]))
    
    if len(flist) == 0:
        print "WARNING: dohydra'd flat corresponding to {} not found".\
            format(rawname)
        print "dying..."
        sys.exit(1)
    elif len(flist) > 1:
        print "WARNING: I don't know which master flat to use:"
        for f in flist:
            print '\t{}'.format(f)
        print "dying..."
        sys.exit(1)
    else:
        return flist[0]

def dohydra_err(errname):
    '''
    Extract apertures, perform flat correction, and apply wavelength solution
    header values to the raw (.sig) error file.
    
    Errors are propagated assuming that the raw flat errors are minimal. This
    causes the full error term,

    dMS' = \sqrt( (dF*MS/F**2)**2 + (dMS/F)**2 ),

    to reduce to

    dMS' = dMS/F,

    where MS is the extracted, .ms file, F is the flat and d represents errors
    on that particular image. The error on a particular aperture in the .ms
    file is simply

    dMS = \sqrt( \sum_i(dI**2) ),

    where dI is the raw error (.sig) file and the sum_i is over all columns
    used in a fiber.

    Given the above expressions the final error (dMS') is achieved by simply
    passing dI**2 and F**2 into dohydra and taking the square root of the
    result.
    '''
    #Save the OG flat so we can set it back once we're done with the squared
    #flat
    pd = iraf.dohydra.getParDict() normal_flat = pd['flat'].get()

    msflat = find_msname(normal_flat)
    sq_errname, sq_flatname = create_tmps(errname,msflat)

    iraf.dohydra(sq_errname,
                 flat='tmp_sq{}'.format(normal_flat),
                 readnoise=3.9,
                 gain=0.438,
                 fibers=109,
                 width=5,
                 minsep=1,
                 maxsep=10,
                 apidtable=APIDTABLE,
                 scatter=False,
                 clean=False,
                 dispcor=True,
                 savearc=False,
                 skyalign=False,
                 skysubt=False,
                 skyedit=False,
                 savesky=False,
                 splot=False,
                 redo=False,
                 update=False,
                 batch=False,
                 listonl=False)
                 
    dohydra_output = '{}.ms.fits'.format(sq_errname.split('.fits')[0])
    final_image = dohydra_output.replace('tmp_sq','').\
        replace('sig.ms.fits','me.fits')
    pow_image(dohydra_output,final_image,0.5)

    #Set the flat back to whatever it was before
    pd['flat'].set(normal_flat)
    iraf.dohydra.saveParList()

    return final_image

def dispcor_err(msfile):
    '''
    Read wavelength solution information from the FITS header and resample an
    image to a linear wavelength scale.

    The error is propagated by simply passing a squared .ms.fits file to
    DISPCOR and taking the square root of the result. This method makes the
    following assumptions:

    1. The wavelength solution in the header is close enough to linear that
       spline3 interpolation used by DISPCOR essentially becomes a linear
       interpolation. This means the error on each output pixel is just the
       quadrature sum of the errors of the input pixels that went into that
       output pixel (divided by the number of pixels).

    2. Any effects of fractional pixels is minimal.
    '''
    tmpname = 'tmp_sq{}'.format(msfile)
    outputname =  msfile.replace('me.fits','me_lin.fits')
    tmpoutput = 'tmp_sq{}'.format(outputname)
    pow_image(msfile,tmpname,2.)

    iraf.dispcor(tmpname,
                 tmpoutput,
                 linearize=True,
                 samedisp=True)

    pow_image(tmpoutput,outputname,0.5)

    return outputname

def calibrate_err(mefile):
    '''
    Apply a sensitivity function to a linearized multispectrum file.
    
    The error is propagated simply by dividing it by the sensitivity
    function. This makes the assumption that there is no error in the
    sensitivity function itself.
    '''
    outputname = mefile.replace('me_','me_rf_')
    iraf.calibrate(mefile,
                   outputname,
                   extinct=False,
                   flux=True,
                   ignoreap=True,
                   sensiti='sens',
                   fnu=False)

    return outputname

def main():
    '''
    Take the first argument passed to this program and run it through aperture
    extraction, flat-fielding, linearization, and flux calibration.
    '''
    errimage = sys.argv[1]
    print 'Running dohydra on {}'.format(errimage)
    hydra = dohydra_err(errimage)
    
    print 'Linearizing {}'.format(hydra)
    cor = dispcor_err(hydra)
    
    print 'Calibrating {}'.format(cor)
    final = calibrate_err(cor)
    print 'Created {}'.format(final)

    print 'Cleaning intermediates'
    os.system('rm tmp_sq*')
    return 0

if __name__ == '__main__':
    sys.exit(main())
