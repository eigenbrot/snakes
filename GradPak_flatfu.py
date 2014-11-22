#! /usr/bin/python

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

def scale_images(hdulist):

    exptimes = [h.header['EXPTIME'] for h in hdulist]
    scales = exptimes[0]/np.array(exptimes)

    print 'Scaling flats...'
    for h, scale in zip(hdulist, scales):
        print '\t{} {}'.format(h.header['OBJECT'],scale)
        h.data *= scale

    return hdulist

def setup_files(flatlist):

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

    tmpdata = np.ones(hdu.data.shape)
    tmphdu = pyfits.PrimaryHDU(tmpdata,hdu.header)
    tmphdu.header['OBJECT'] = 'FlatFu tmp'

    print 'Writing tmp spectrum'
    tmphdu.writeto('ffTmp.fits',clobber=True)

    return

def initial_run(scalednames):

    print 'Doing initial flat aperture extraction...'
    idtable = os.path.basename(iraf.dohydra.apidtable)
    print '\tusing idtable {}'.format(idtable)
    outputnames = []
    for flat in scalednames:
        print '\trunning {}'.format(flat)
        try:
            iraf.dohydra('ffTmp.fits',
                         apref=scalednames[0],
                         flat=flat,
                         readnoise=3.9,
                         gain=0.438,
                         fibers=109,
                         width=5,
                         minsep=1,
                         maxsep=10,
                         apidtable='/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf',
                         scatter=False,
                         fitflat=False,
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
            print 'Fucked up, trying again'
            return False
        outputnames.append('{}{}.ms.fits'.format(flat.split('.fits')[0],idtable))
        os.remove('ffTmp.ms.fits')

    return outputnames

def stitch_flats(outputnames,pivots):

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
                  logfile='STDOUT')

    for tmp in tmpfiles:
        os.remove(tmp)

    iraf.dohydra.flat = 'dFlat_master.fits'

    return

def main():

    flat_list = []
    pivot_list = []

    for token in sys.argv[1:]:
        try:
            pivot_list.append(int(token))
        except ValueError:
            flat_list.append(token)

    print 'Flat list is {}'.format(flat_list)
    print 'Pivots are {}'.format(pivot_list)

    if len(pivot_list) != len(flat_list) - 1:
        print "There must be one less pivot than flats"
        return 1

    sl = setup_files(flat_list)
    msl = initial_run(sl)
    if not msl:
        msl = initial_run(sl)
    stitch_flats(msl,pivot_list)

    return 0
            
if __name__ == '__main__':

    main()
