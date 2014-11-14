#!/usr/bin/python
#
##################################################
#
# This is a script used to do sky subtraction on multispec files produced by
# the WIYN Bench Spectrograph and the GradPak IFU. To make use of this script
# the spectra need to be reduced usig the gradpak_sizes.iraf aperture
# table. This table defines different beam numbers for each fiber size such
# that the beam number is fibersize(in microns)/100. The sky fibers are
# assigned beam numbers as fibersize/100*11.
#
# History:
#     v1 - A. Eigenbrot Nov. 2014
#
###################################################

import sys
import os

#Load the IRAF packages we'll need
try:
    from pyraf import iraf
    iraf.imred(_doprint=0)
    iraf.hydra(_doprint=0)
except Exception as e:
    print "Failure: could not find pyraf/iraf"
    sys.exit(1)

def skysub(imagename,fibersize):
    '''
    Take a multispec file with GradPak beam numbers and do a sky subtraction
    of a single fiber size. The fiber size is an integer that is the actual
    size of the fibers (in microns) divided by 100. So the 200 micron fibers
    have fibersize=2.

    The result is a .ms file that contains sky subtracted spectra for only
    that particular fiber size. The name of this file is returned.
    '''
    skyname = '{:}.ms_s{:n}_lin.fits'.\
        format(imagename.split('.ms_lin.fits')[0],fibersize)
    print "Generating {} micron sky sub image {}".format(fibersize*100,skyname)

    iraf.skysub(imagename,
                output=skyname,
                objbeam='{},{}'.format(fibersize,fibersize*11),
                skybeam='{}'.format(fibersize*11),
                skyedit='yes',
                combine='average',
                reject='avsigclip',
                scale='none',
                savesky='yes',
                logfile='spool.txt')

    return skyname

def recombine(imagelist,outputimage):
    '''
    Take a list of single-fiber-size, sky-subtracted .ms files and put them
    back into a single multispec file with all fibers.

    This uses the scombine task with group set to apertures. Because the
    aperture information is preserved when doing sky subtraction all this does
    is concatenate the individual .ms files into a single file.
    '''
    print 'Combining appertures from:'
    for image in imagelist: print '\t'+image

    iraf.scombine(','.join(imagelist),outputimage,
                  apertures='',
                  group='apertures',
                  combine='average', #These are irrelivant b/c
                  reject='avsigclip',# we aren't actually combining anything
                  first='yes', #Very important
                  w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',
                  log='no',
                  gain=0.438,
                  rdnoise=3.9,
                  logfile='spool.txt')

    return

def cleanup(imagelist):
    '''
    Delete all the intermediate, single-fiber-size, sky-subtracted images
    '''
    print 'cleaning intermediates:'
    for image in imagelist:
        print '\t'+image
        os.remove(image)

    return
                  
def main():
    '''
    Take the input file and do a sky subtraction for each of the GradPak fiber
    sizes. Then combine everything back into the output file name.
    '''

    if len(sys.argv) != 3:
        #We don't want to overwhelm the user with verbosity
        return "The request was made, but it was not good"

    imagename = sys.argv[1]
    outputname = sys.argv[2]

    #We do this check now because scombine doesn't have a clobber option
    if os.path.exists(outputname):
        print "Warning, you are about to overwrite an existing {}".\
            format(outputname)
        clobber = raw_input('Continue? (Y/n): ')
        if clobber.lower() == 'n':
            sys.exit()
        else:
            os.remove(outputname)
    
    #Fuck this line is cool
    skylist = [skysub(imagename,f) for f in range(2,7)]
    recombine(skylist,outputname)
    cleanup(skylist)

    return 0

if __name__ == '__main__':
    sys.exit(main())
