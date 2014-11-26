#! /usr/bin/python
#
##################################################
#
# This script sets up multiple flats so they can be used as a single
# input to dohydra. You will want to do this because an exposure that
# puts enough signal in the 200 micron fibers will push the larger
# fibers into the non-linear regime.
#
# The calling syntax is:
# 
# >$ python GradPak_flatfu.py Flat1 Flat2... Flatn pivot1 pivot2... pivotn-1
#
# Where the pivots define the aperture at which to cut the flat,
# inclusive. For example, if the call is:
#
# >$ python GradPak_flatfu.py Flat1.fits Flat2.fits 70
#
# Then the resulting master flat will have apertures 1 - 70 from
# Flat1.fits and 71 - 109 from Flat2.fits.
#
# Flat1 is the flat used for the aperture extraction of the entire
# run, so it should probably be the one with the longest exposure. 
#
# All the flats are scaled by their exposure times and then cut up as
# specified by the user. The final, master flat lives only as a .ms
# multispec file. For this reason it is imperative that you set the
# "flat" parameter to EXACTLY dFlat_master.fits. Nothing else will
# work. Even when you do this dohydra will complain that it can't find
# dFlat_master.fits. This warning is OK; dohydra will still do
# everything you want it to.
#
# History:
#      v1 - A. Eigenbrot Nov. 2014
#
####################################################


import os
import sys
import pyfits
import numpy as np

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

if os.getlogin() == 'Arthur':
    APIDTABLE = '/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf',
else:
    APIDTABLE = '/usr/users/eigenbrot/research/Pak/gradpak_sizes.iraf'

def scale_images(hdulist):
    '''
    Take in a list of fits HDUS and scale the data in all of them to
    the exposure time of the first HDU
    '''
    exptimes = [h.header['EXPTIME'] for h in hdulist]
    scales = exptimes[0]/np.array(exptimes)

    print 'Scaling flats...'
    for h, scale in zip(hdulist, scales):
        print '\t{} {}'.format(h.header['OBJECT'],scale)
        h.data *= scale

    return hdulist

def setup_files(flatlist):
    '''
    Take in a list of names of flat fields, scale the data by exposure
    time, and then write them out to appropriatly named files.
    Also create a dummy spectrum that will be "reduced" by dohydra.
    '''
    hdulist = [pyfits.open(i)[0] for i in flatlist]
    hdulist = scale_images(hdulist)

    scalednames = []
    print 'Writing flats...'
    for h, flat in zip(hdulist, flatlist):
        newname = '{}_scale.fits'.format(flat.split('.fits')[0])
        print '\t{}'.format(newname)
        h.writeto(newname,clobber=True)
        scalednames.append(newname)
        
    make_tmp(hdulist[0])

    return scalednames

def make_tmp(hdu):
    '''
    Make the dummy spectrum to be reduced by dohydra. This spectrum is
    all ones so that the affect of the flat scaling can be seen.
    '''
    tmpdata = np.ones(hdu.data.shape)
    tmphdu = pyfits.PrimaryHDU(tmpdata,hdu.header)
    tmphdu.header['OBJECT'] = 'FlatFu tmp'

    print 'Writing tmp spectrum'
    tmphdu.writeto('ffTmp.fits',clobber=True)

    return

def initial_run(scalednames,traceflat,fitflat=True):
    '''
    Use dohydra to extract .ms multispec files from the input flat
    fields. This is done so that we can stitch them back together in
    aperture space rather than pixel space. To save time we don't
    bother with any sort of wavelength solution or other garbage.
    '''
    print 'Doing initial flat aperture extraction...'
    idtable = os.path.basename(iraf.dohydra.apidtable)
    print '\tusing idtable {}'.format(idtable)
    outputnames = []
    for flat in scalednames:
        print '\trunning {}'.format(flat)
        try:
            iraf.dohydra('ffTmp.fits',
                         apref=traceflat,
                         flat=flat,
                         readnoise=3.9,
                         gain=0.438,
                         fibers=109,
                         width=5,
                         minsep=1,
                         maxsep=10,
                         apidtable=APIDTABLE,
                         scatter=False,
                         fitflat=fitflat,
                         clean=False,
                         dispcor=False,
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
        except iraf.IrafError:
            '''For some reason, if you run this from the shell
            (non-interactively) dohydra will always fail the first
            time when trying to copy the database file. If/when this
            happens we'll just return false so main() knows just to
            try it again.
            '''
            print 'Fucked up, trying again'
            return False
        outputnames.append('{}{}.ms.fits'.format(flat.split('.fits')[0],idtable))
        os.remove('ffTmp.ms.fits')

    return outputnames

def stitch_flats(outputnames,pivots):
    '''
    Take a list of names of multispec flat files produced by dohydra
    and a list of pivot apertures and stitch together a master flat.

    The pivot values are inclusive, so if pivots = [72] then the
    master flat will contain apertures 1 - 72 from flat #1 and 73 -
    109 from flat #2
    '''
    pivots = [0] + pivots + [109]

    tmpfiles = []
    print 'Extracting flat apertures...'
    for i, flat in enumerate(outputnames):
        print '\ttaking {} from {} to {}'.format(flat,pivots[i]+1,pivots[i+1])
        name = 'tmp{}'.format(flat)
        iraf.scopy(flat,name,
                   apertur='{}-{}'.format(pivots[i]+1,pivots[i+1]),
                   w1='INDEF',
                   w2='INDEF',
                   format='multispec',
                   verbose=False)
        tmpfiles.append(name)

    print 'Stitching together master flat'
    iraf.scombine(','.join(tmpfiles),'dFlat_master{}.ms.fits'.format(os.path.basename(iraf.dohydra.apidtable)),
                  apertur='',
                  group='apertures',
                  first=True,
                  w1='INDEF',
                  w2='INDEF',
                  dw='INDEF',
                  nw='INDEF',
                  log=False,
                  scale='none',
                  zero='none',
                  weight='none',
                  logfile='flatfu.log')

    for tmp in tmpfiles:
        os.remove(tmp)

    '''I would really like this to work in no-interactive mode, but
    for some reason pyraf doesn't actually write it's params to the
    uparam file. Oh well.
    '''
    iraf.dohydra.flat = 'dFlat_master.fits'

    return

def parse_input(inputlist):

    flat_list = []
    pivot_list = []
    fitflat = True
    traceflat = False

    for token in inputlist:
        if token == '-nf':
            fitflat = False
        elif token[0:2] == '-t':
            traceflat = token[2:]
            flat_list.append(token[2:])
        elif '.fits' in token:
            flat_list.append(token)
        else:
            try:
                pivot_list.append(int(token))
            except ValueError:
                print 'Could not parse option {}'.format(token)
    if not traceflat:
        traceflat = flat_list[0]
        
    return flat_list, pivot_list, traceflat, fitflat

def main():
    '''
    Parse the user inputs, check that there are the right number of
    pivot points, and run through the steps in the script.
    '''
    flat_list, pivot_list, traceflat, fitflat = parse_input(sys.argv[1:])

    print 'Flat list is {}'.format(flat_list)
    print 'Pivots are {}'.format(pivot_list)
    print 'Tracing using {}'.format(traceflat)
    
    if len(pivot_list) != len(flat_list) - 1:
        print "There must be one less pivot than flats"
        return 1

    '''Run the script'''
    sl = setup_files(flat_list)
    msl = initial_run(sl, traceflat, fitflat)
    if not msl:
        '''Here is where we catch IRAF being bad'''
        msl = initial_run(sl, traceflat, fitflat)
    stitch_flats(msl,pivot_list)

    return 0
            
if __name__ == '__main__':

    try:
        sys.exit(main())
    except Exception as e:
        print "The request was made but it was not good"
        sys.exit(1)
