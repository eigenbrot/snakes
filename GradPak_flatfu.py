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
# > python GradPak_flatfu.py Flat1 Flat2... Flatn pivot1 pivot2... pivotn-1
#
# Where the pivots define the aperture at which to cut the flat,
# inclusive. For example, if the call is:
#
# > python GradPak_flatfu.py Flat1.fits Flat2.fits 70
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
#      v1.1 - A. Eigenbrot Dec. 2014
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
    APIDTABLE = '/Users/Arthur/Documents/School/MetaPak/gradpak_sizes.iraf'
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
        h.header['EXPTIME'] *= scale

    return hdulist

def scale_spectra(imagelist):

    fiber1 = 44
    fiber2 = 62

    hdulist = [pyfits.open(image)[0] for image in imagelist]
    means = [np.mean(h.data[fiber1 - 1:fiber2 - 1,:]) for h in hdulist]
    print means
    scales = means[0]/np.array(means)

    print 'Scaling extracted flats...'
    outputnames = ['{}_scale.ms.fits'.format(image.split('.ms.fits')[0]) for image in imagelist]
    for h, scale, name in zip(hdulist, scales, outputnames):
        print '\t{} {}'.format(name, scale)
        h.data *= scale
        h.writeto(name,clobber=True)

    return outputnames

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
        print '\t{}'.format(newname), np.mean(h.data[:,3:70]/hdulist[0].data[:,3:70])
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

def initial_run(scalednames,traceflat,throughput=''):
    '''
    Use dohydra to extract .ms multispec files from the input flat
    fields. This is done so that we can stitch them back together in
    aperture space rather than pixel space. To save time we don't
    bother with any sort of wavelength solution or other garbage.
    '''
    
    oldlog = iraf.hydra.logfile
    print 'Doing initial flat aperture extraction...'
    idtable = os.path.basename(iraf.dohydra.apidtable)
    print '\tusing idtable {}'.format(idtable)
    outputnames = []
    outputscales = []
    for flat in scalednames:
        print '\trunning {}'.format(flat)
        iraf.hydra.logfile = '{}.log'.format(flat)
        try:
            iraf.dohydra('ffTmp.fits',
                         apref=traceflat,
                         flat=flat,
                         through=throughput,
                         readnoise=3.9,
                         gain=0.438,
                         fibers=109,
                         width=5,
                         minsep=1,
                         maxsep=10,
                         apidtable=APIDTABLE,
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
                         listonl=False,
                         Stdout=1)
        except iraf.IrafError:
            '''For some reason, if you run this from the shell
            (non-interactively) dohydra will always fail the first
            time when trying to copy the database file. If/when this
            happens we'll just return false so main() knows just to
            try it again.
            '''
            print 'Fucked up, trying again'
            return False, False
        f = open('{}.log'.format(flat),'r')
        o = f.readlines()
        for dump in o:
            if 'bscale' in dump:
                scale = dump.split()[-1]
                outputscales.append(float(scale))
            if 'Create the normalized response' in dump:
                outname = '{}.fits'.format(dump.split()[-1])
                if outname not in outputnames:
                    outputnames.append('{}.fits'.format(dump.split()[-1]))
        os.remove('ffTmp.ms.fits')
        f.close()

    print "Extracted flats:"
    for out,scale in zip(outputnames,outputscales):
        print '\t{}  {}'.format(out,scale)
    iraf.hydra.logfile = oldlog
    return outputnames, outputscales

def stitch_flats(outputnames,pivots,outstring):
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


    mastername = 'dFlat_master{}.ms.fits'.format(outstring)
    print 'Stitching together master flat {}'.format(mastername)    
        
    iraf.scombine(','.join(tmpfiles),mastername,
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
    return mastername

def get_scrunch(flatname, msname):
    '''
    Take in two strings and return the difference between them, modulo
    the difference between .ms.fits and .fits. 

    This is used to figure out what scrunching scheme IRAF used to
    name the aperture-extracted flats. This function exists because
    IRAF has different behavior when it comes to scrunching on
    different systems.
    '''
    basefname = flatname.split('.fits')[0]
    basemsname = msname.split('.ms.fits')[0]
    
    return basemsname.replace(basefname,'')

def mean_scale(mslist,scalelist):
    '''
    Take in a list of fits file names and scale each file by the corresponding
    value in the scalelist. The input file is overwritten.

    When constructing the aperture-extracted flat that will be applied
    to all data apertures IRAF's last step is to normalize the entire
    multispec file to a unity mean. Because we will be stitching
    together a few different flats we need to undo this scaling so the
    relative strengths of the individual flats is maintained. This
    assumes that these flats were already properly scaled (see
    scale_images()).
    '''
    print 'Undoing IRAF flat scaling...'
    mshdulist = [pyfits.open(i)[0] for i in mslist]
    for h, name, scale in zip(mshdulist, mslist, scalelist):
        print '\t{:} scale: {:5.4f}'.format(name,scale)
        h.data *= scale
        h.writeto(name,clobber=True)

    return

def fit_flat(mastername):

    print 'Dividing out average spectrum'
    h = pyfits.open(mastername)[0]
    h.data /= np.mean(h.data,axis=0)
    h.writeto(mastername,clobber=True)

    return

def normalize(mastername):
    '''
    Take the name of a fits file and normalize it so the mean over all pixels is 1. The input file is overwritten.

    This function exists to keep reduction as close as possible to
    IRAF's own routines. In msresp1d the last step is to normalize the
    flat so the mean is one. In mean_scale() we undid this scaling so
    that the flats could be stitched together. This function reapplys
    this normalization to the stitched together flat.
    '''
    print 'Renormalizing final flat'
    h = pyfits.open(mastername)[0]
    h.data /= np.mean(h.data)
    h.writeto(mastername,clobber=True)

    return

def parse_input(inputlist):
    '''
    Take the list of arguments passed to the terminal when this
    program was called and parse them into useful arguments and
    options.
    '''
    flat_list = []
    pivot_list = []
    fitflat = True
    traceflat = False
    throughput = ''

    for i, token in enumerate(inputlist):
        if token == '-nf':
            fitflat = False
        elif token == '-t':
            traceflat = inputlist[i+1]
            flat_list.append(inputlist[i+1])
        elif token == '-r':
            throughput = inputlist[i+1]
            del inputlist[i+1]
        elif '.fits' in token:
            if token not in flat_list:
                flat_list.append(token)
        else:
            try:
                pivot_list.append(int(token))
            except ValueError:
                print 'Could not parse option {}'.format(token)
    if not traceflat:
        traceflat = flat_list[0]
        
    return flat_list, pivot_list, traceflat, throughput, fitflat

def print_help():
    '''
    Print usage help.
    '''
    print """Calling syntax:
    
> GradPak_flatfu.py flat1.fits flat2.fits ... flatn.fits pivot

    Options:
            -t   Use the next flat as the trace flat. Default is 
                 the first flat given.

            -r   Use the next flat as the throughput (sky) image.
                 This flat is not used to make the master flat.

            -nf  Do not fit a average spectral function to the flats.
                 Equivalent to fitflat- in IRAF.

    Common Example:
             > GradPak_flatfu.py dFlat_4s.fits -t dFlat_1s.fits -r sFlat.fits -nf 43

             This will create a master flat consisting of apertures 1
             - 43 from dFlat_4s.fits and apertures 44 - 109 from
             dFlat_1s.fits. No functions will be fit to the flats but
             a fiber-by-fiber scaling will be applied based on the
             total counts in each aperture extracted from sFlat.fits
    """
    return

def main():
    '''
    Parse the user inputs, check that there are the right number of
    pivot points, and run through the steps in the script.
    '''
    flat_list, pivot_list, traceflat, throughput, fitflat = parse_input(sys.argv[1:])

    print 'Flat list is {}'.format(flat_list)
    print 'Pivots are {}'.format(pivot_list)
    print 'Tracing using {}'.format(traceflat)
    print 'Throughput file is {}'.format(throughput)
    
    if len(pivot_list) != len(flat_list) - 1:
        print "There must be one less pivot than flats"
        return 1

    '''Run the script'''
    # sl = setup_files(flat_list)
    make_tmp(pyfits.open(flat_list[0])[0])
    msl, scales = initial_run(flat_list, traceflat, throughput)
    if not msl:
        '''Here is where we catch IRAF being bad'''
        msl, scales = initial_run(sl, traceflat, throughput)
    outstring = get_scrunch(flat_list[0],msl[0])
    mean_scale(msl,scales)
    msl = scale_spectra(msl)
    master = stitch_flats(msl,pivot_list,outstring)
    if fitflat:
        fit_flat(master)
    normalize(master)

    '''Finally, set some IRAF parameters to make it easier to run dohydra with
    the master flat'''
    pd = iraf.dohydra.getParDict()
    pd['apref'].set(traceflat)
    pd['flat'].set('dFlat_master.fits')
    pd['through'].set(throughput)
    iraf.dohydra.saveParList()

    return 0
            
if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "The request was made but it was not good"
        sys.exit(1)
    elif sys.argv[1] == '-h':
        sys.exit(print_help())
    # try:
    #     sys.exit(main())
    # except Exception as e:
    #     print "The request was made but it was not good"
    #     sys.exit(1)
    sys.exit(main())
