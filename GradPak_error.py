#! /usr/bin/python

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
except Exception as e:
    print "Failure: could not find pyraf/iraf"
    sys.exit(1)

glob = glob.glob

if os.getlogin() == 'Arthur':
    APIDTABLE = '/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf'
else:
    APIDTABLE = '/usr/users/eigenbrot/research/Pak/gradpak_sizes.iraf'


def pow_image(inputname, outputname, power):
    
    print 'raising {} to power of {}'.format(inputname,power)
    
    h = pyfits.open(inputname)[0]
    
    pyfits.PrimaryHDU(h.data**power,h.header).writeto(outputname)
    return

def create_tmps(errname, flatname):

    sq_errname = 'tmp_sq{}'.format(errname)
    sq_flatname = 'tmp_sq{}'.format(flatname)

    pow_image(errname,sq_errname,2.)
    pow_image(flatname,sq_flatname,2.)
    
    return sq_errname, sq_flatname

def find_msname(rawname):

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

    pd = iraf.dohydra.getParDict()
    normal_flat = pd['flat'].get()

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

    pd['flat'].set(normal_flat)
    iraf.dohydra.saveParList()

    return final_image

def dispcor(msfile):

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

def main():

    errimage = sys.argv[1]
    print 'Running dohydra on {}'.format(errimage)
    hydra = dohydra_err(errimage)
    
    print 'Linearizing {}'.format(hydra)
    final = dispcor(hydra)
    
    print 'Created {}'.format(final)

    return 0

if __name__ == '__main__':
    sys.exit(main())
