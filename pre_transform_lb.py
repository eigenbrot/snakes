import numpy as np
import glob
import pyfits
import os
import ADEUtils as ADE
import matplotlib.pyplot as plt
import scipy.misc.pilutil as pl
from pyraf import iraf
from pyraf.iraf import stsdas
from datetime import datetime
import time

debug=1
fignum=0
frd_sav=0
_ = np.seterr(invalid='ignore')


def fat_angel(findstr, num_ap, EXTEN=0, OUTPUT=0,
              PIXELPITCH=1, FIBERSIZE=500):
    '''
    Description:
        Fat_angel is deisgned to process a large number of data from
        the laser bench. It takes an input string, findstr, and computes
        the width and radius of the ring in each data image. Surface brightness
        profiles are found using Annulize in the ADEUtils package with
        the number of annuli specfied in the num_ap input. Output can be
        directed to a file for plotting bliss.
    
    Inputs:
        findstr- Str
                 A string containing the names of the data files to use.
                 Wildcards are allowed so you can chose a single file or
                 as many as want. Files are assumed to be FITS files.
        num_ap - Int
                 The number of annuli to use when constructing the surface
                 brightness profile.
        EXTEN  - Int
                 The FITS extension where the primary data is stored. As of
                 right now there is no way to specify different extensions
                 for different files.
        OUTPUT - Str
                 Name of output file. Output contains angle, ring radius, and
                 ring width in column form
        PIXELPITCH - Float
                 Size of the camera pixels in micrometers. This is only used to 
                 get the correct scale on the debug plots.
    
    Output:
        Output is a tuple of Numpy vectors containing angle, ring radius, and
        ring width.
    
    Example:
        Assume you have a bunch of data called X_red.fits where X is some
        data iterator.

        >>lb.fat_angel('*_red.fits',150,EXTEN=1,OUTPUT='my_data.dat')
    '''


    file_list = glob.glob(findstr)
    numfiles = len(file_list)

    'initialize some data arrays'
    widths = np.zeros(numfiles)
    radii = np.zeros(numfiles)
    angles = np.zeros(numfiles)
    r1_vec = np.zeros(numfiles)
    r2_vec = np.zeros(numfiles)
    frd_widths = np.zeros(numfiles)
    xcent = np.zeros(numfiles)
    ycent = np.zeros(numfiles)

    t1 = time.time()
    for i in range(numfiles):
        print(file_list[i])

        ADE.mediclean(file_list[i],'clean'+file_list[i],exten=EXTEN)

        huds = pyfits.open('clean'+file_list[i])
        data = huds[EXTEN].data
        
        'get the angle from the FITS header. MUCH easier in python than IDL!'
        angles[i] = huds[EXTEN].header['ANGLE']

        '''poor man's background subtraction'''
        mode = ADE.mode(data)[0][0]
        if debug: print 'Mode = '+str( mode)

 
        'Annulize!'
        if debug:
            fig0 = plt.figure(0)
            plt.clf()
            t0 = time.time()
            (r_vec,fluxes,center) = ADE.annulize(data,num_ap,NOREAD=1)#,MODE=mode)
            print "Annulize took "+str(time.time() - t0)+" seconds"
            plt.plot(r_vec*PIXELPITCH,fluxes)
            fig0.show()
        else:
            (r_vec,fluxes,center) = ADE.annulize(data,num_ap,NOREAD=1)#,MODE=mode) 

        '''Compute the cdf from the pdf (fluxes). Working with the CDF
        allows us to assume that the annulus is gaussian-ish without 
        having to worry about the particulars'''
        (rm_idx,r1,r2) = find_peak(r_vec,fluxes)
        
        rm = r_vec[rm_idx]
 #       r1 = r_vec[r1_idx]
  #      r2 = r_vec[r2_idx]
        
        'Now deconvolve and find the width of the FRD smearing kernel'
#        (frd_width,frd) = decon(r_vec,fluxes,rm_idx,r2,r1,FIBERSIZE,PIXELPITCH)

#        global frd_sav
#        frd_sav = frd
        
        if debug:
            fig1 = plt.figure(fignum)
            plt.clf()
            sp0 = fig1.add_subplot(222)
            sp0.plot(r_vec*PIXELPITCH, 
                     np.cumsum(fluxes)/np.max(np.cumsum(fluxes)))
            sp0.axvline(x=rm*PIXELPITCH,ls='--',lw=0.3)
            sp0.axvline(x=r1*PIXELPITCH,ls='--',lw=0.3)
            sp0.axvline(x=r2*PIXELPITCH,ls='--',lw=0.3)
            sp0.set_xlabel("Radius (um)")
            sp0.set_ylabel("% of total counts")
            sp0.set_title("Normalized CDF")

#            plt.figure(0)
            sp1 = fig1.add_subplot(221)
            sp1.plot(r_vec*PIXELPITCH,fluxes)
            sp1.axvline(x=rm*PIXELPITCH,ls='--',lw=0.3)
            sp1.axvline(x=r1*PIXELPITCH,ls='--',lw=0.3)
            sp1.axvline(x=r2*PIXELPITCH,ls='--',lw=0.3)
            sp1.set_xlabel("Radius (um)")
            sp1.set_ylabel("Counts")
            sp1.set_title("Ring Profile")

#            sp2 = fig1.add_subplot(224)
#            sp2.plot(r_vec*PIXELPITCH, frd)
#            sp2.set_xlabel("Radius (um)")
#            sp2.set_ylabel("??")
#            sp2.set_title("FRD kernel")

            plt.suptitle("Angle = "+str(angles[i])+" degrees\n"+file_list[i])
            fig1.show()

            print "Center: "+str(center)

            if numfiles > 1: raw_input("press enter to continue...\n")


        widths[i] = r2 - r1
        radii[i] = rm
        r1_vec[i] = r1
        r2_vec[i] = r2
#        frd_widths[i] = frd_width
        xcent[i] = center[0]
        ycent[i] = center[1]

    
    print "Total annulize time was "+str(time.time()-t1)+" seconds"
    'We sort the data just make the output a little more readable'
    sort_idx = np.argsort(angles)

    widths *= PIXELPITCH
    radii *= PIXELPITCH
    r1_vec *= PIXELPITCH
    r2_vec *= PIXELPITCH
    frd_widths *= PIXELPITCH
    
    if OUTPUT:
        f = open(OUTPUT, 'w')
        f.write('#angle       radius     width     r1        r2          frd width        center\n')
        for i in range(angles.shape[0]):
            np.array([angles[sort_idx][i],
                      radii[sort_idx][i],
                      widths[sort_idx][i],
                      r1_vec[sort_idx][i],
                      r2_vec[sort_idx][i],
                      frd_widths[sort_idx][i],
                      xcent[sort_idx][i],
                      ycent[sort_idx][i]]).\
                      tofile(f, sep='   ',format='%3.4f')
            f.write('\n')
        
    return (angles[sort_idx],
            radii[sort_idx],
            widths[sort_idx],
            r1_vec[sort_idx],
            r2_vec[sort_idx],
            frd_widths[sort_idx],
            xcent[sort_idx],
            ycent[sort_idx])

def sodom(findstr, exten=0):

    angles = {}
    file_list = glob.glob(findstr)

    for image in file_list:
        angle = pyfits.open(image)[exten].header['ANGLE']
        if angle not in angles.keys():
            angles[angle] = []

        angles[angle].append(image)

    iraf.imcombine.combine = 'average'
    iraf.imcombin.reject = 'none'

    for theta in angles.keys():
        if angles[theta][0].find('../') > -1: s = 3
        else: s = 0
        name = 'final_'+angles[theta][0][s:angles[theta][0].rfind('_')]+'.fits'
        
        iraf.imcombine(','.join([s+'['+str(exten)+']' for s in angles[theta]]),name)
        ADE.mediclean(name,name)

    return

def decon(r, f, idx, r1, r2, fiber_size, pixel_pitch):

#    x1 = np.where(r >= r1)[0][0]
#    x2 = np.where(r >= r2)[0][0]

    'in pixel units'
    radius = (fiber_size/pixel_pitch)/(2)#*(r2-r1)/(x2-x1))

#    idx = np.where(r >= peak)[0][0]
    idx += chi(r,f,radius, idx)

    size = f.shape[0]

    _ = np.seterr(invalid='ignore')
    
    theta = np.arccos(np.abs(r - r[idx])/radius)
    pn = radius*np.nan_to_num(np.sin(theta))

    zpn = pn[np.where(pn != 0)[0]]
#    zpn = pn
#    zpn = np.zeros(f.size)
#    zpn[0:f.size] = pn[np.where(pn != 0)[0][0]]

    ff = np.copy(f)/np.sum(f[np.where(f > 0)[0]])
    pnf = np.copy(zpn)/np.sum(zpn)


    if debug:
        fig1 = plt.figure(fignum)
        plt.clf()

        sp3 = fig1.add_subplot(223)

        sp3.plot(r*pixel_pitch,ff,label='Ring profile')
        sp3.plot(r*pixel_pitch,pn/pn.sum(),label='Ideal beam profile')
        sp3.legend()
        sp3.set_title("Area normalized data and ideal beam profile")
        sp3.set_xlabel("Radius (um)")
        sp3.set_ylabel("Normalized Counts")

    'so we can parallelize'
    temp_str = 'temp_'+str(datetime.now().microsecond)
    idl_name = temp_str+'_idl.fits'

    pyfits.PrimaryHDU(ff).writeto(temp_str+'_ff.fits')
    pyfits.PrimaryHDU(pnf).writeto(temp_str+'_pn.fits')

    stsdas()
    stsdas.analysis()
    stsdas.analysis.restore()
    stsdas.analysis.restore.lucy\
        (temp_str+'_ff.fits',temp_str+'_pn.fits',temp_str+'_frd.fits',1,0)
    
    command = "\"MWRFITS, MRDFITS('"+temp_str+"_frd.fits'),'"+idl_name+"'\""

    os.system('idl -quiet -e '+command)

    frd = pyfits.open(idl_name)[0].data

    global frd_sav
    frd_save = frd

    (peak, r1, r2) = find_peak(r,frd)

    os.system('rm ./'+temp_str+'*.fits')

    return ((r2 - r1),frd)


def rapture(datalist,output):
    '''
    Calling:
        rapture(datalist, output)
    
    Description:
        Rapture is designed to bring together multiple data ranges
        from the laser bench and modify them to all be on the same scale.
        It then outputs one master file (bringing all the angels together)
        that has all the data for one laser over multiple data ranges.

    Inputs:
        datalist - Str
            A pointer to the data files to combine. Wildcards are
            allowed and indeed necessary because you need to have at
            least two data files for Rapture to be of any use. The 
            files are expected to be output files from fat_angel.
        output   - Str
            Filename where output is directed. Format is columns with angle,
            ring radius, and ring width.

    Output:
        None. Output is directed to the file specified by output.

    Example:
        First create the data files to be used
        
        >>lb.fat_angel('7.0*.fits',100,OUTPUT='7.0.dat')
        >>lb.fat_angel('7.2*.fits',100,OUTPUT='7.2.dat')
        
        Now bring the rapture

        >>lb.rapture('7.*.dat','my_rapture.cat')
    '''

    file_list=glob.glob(datalist)
    
    '''We have to sort the file names because only consecutive data ranges
     will have overlapping data. This means that the whole program is
     highly dependent on you not fucking up the naming'''
    file_list.sort()
        
    f = open(output, 'w')
    f.write('# angle   radius     width     frd width\n')
    
    '''Read in the first file and write it to disk. All other data will be
    shifted to have the same scale as the first file'''
    data1 = np.loadtxt(file_list[0],skiprows=1)
    print(file_list[0])
    for j in range(data1[:,0].shape[0]):
        np.append(data1[j][0:3],data1[j][5]).tofile(f,sep='   ', format='%3.4f')
        f.write('\n')


    '''Now go through the other files and shift them. Shifts are made
     assuming a linear pixel scale where a pixel Y on data2 maps to pixel
     X on data1 by Y = mX + b'''
    for i in range(len(file_list)-1):
        print(file_list[i+1])
        data2 = np.loadtxt(file_list[i+1],skiprows=1)

#        oldData = np.copy(data2[:,1:3])
        
        m1 = get_slope(data1, data2, 1)
        data2[:,1] *= m1
        print(m1)
        
        for j in range(data2[:,0].shape[0]):
            np.append(data2[j][0:3],data2[j][5]).\
                tofile(f, sep='   ', format='%3.4f')
            f.write('\n')
        
        '''As data2 has been modified we are essentially extending the scale 
         of the first data through the other data files one by one.'''
        data1 = data2
    return

def zion(findstr, catfile, outstr, EXTEN=0):
    
    file_list = glob.glob(findstr)
    catData = np.loadtxt(catfile, skiprows=1)

    sort_idx = np.argsort(catData[:,0])
    catData = catData[sort_idx]

    for i in range(len(file_list)):
        print(file_list[i])
        hdu = pyfits.open(file_list[i])
        image = hdu[EXTEN].data
        imAngle = hdu[EXTEN].header['ANGLE']

        idx = np.where(catData[:,0] == imAngle)[0]
        radius = np.average(catData[:,3][idx])
        width = np.average(catData[:,4][idx])

        cent = ADE.centroid(image)
        dims = image.shape

        imRad = radius + width

        cutIm = image[cent[0] - imRad:cent[0] + imRad,\
                          cent[1] - imRad:cent[1] + imRad]

        bigIm = pl.imresize(cutIm, (1028,1028))

        frameNumber = idx[0]
 
        name = str(frameNumber).zfill(3)+'_'+str(imAngle)+outstr+'.fits'

        if len(glob.glob(name)) == 0: pyfits.PrimaryHDU(bigIm).writeto(name)

    return

def find_peak(r,f):

    CDF = np.cumsum(f)
    CDF /= np.max(CDF)

    '''First find the peak (0.5 in the CDF) and then find the
    limits of the fwhm'''
    rp_idx = np.where(CDF >= 0.50)[0][0]
    r1_idx = np.where(f >= 0.5*f[rp_idx])[0][0]
    r2_idx = np.where(f >= 0.5*f[rp_idx])[0][-1]

    '''Now intertpolate to get a better width estimate'''
    m1 = (f[r1_idx] - f[r1_idx - 1])/(r[r1_idx] - r[r1_idx - 1])
    r1 = r[r1_idx - 1] + (0.5*f[rp_idx] - f[r1_idx - 1])/m1
    m2 = (f[r2_idx + 1] - f[r2_idx])/(r[r2_idx + 1] - r[r2_idx])
    r2 = r[r2_idx] + (0.5*f[rp_idx] - f[r2_idx])/m2

    return (rp_idx, r1, r2)

def get_slope(data1, data2, dim):
    '''
    get_slope takes in two data arrays that contain data with different
    slopes and returns m such that Y = mX + b where Y represents points
    in data1 and X represents points in data2
    '''
    

    (GT,LT) = format_data(data1, data2, dim)
    vec1GT = GT[0]
    vec2GT = GT[1]
    
    vec1LT = LT[0]
    vec2LT = LT[1]
    
    mGT = np.average(vec1GT/vec2GT)
    mLT = np.average(vec1LT/vec2LT)

    return (mGT + mLT)/2

def get_shift(data1, data2, dim):
    '''
    get_shit takes in two sets of data THAT ARE ASSUMED TO HAVE THE SAME
    SLOPE and returns the offset between them.
    '''
    
    (GT,LT) = format_data(data1, data2, dim)
    vec1GT = GT[0]
    vec2GT = GT[1]
    
    vec1LT = LT[0]
    vec2LT = LT[1]

    GT_diffs = vec1GT - vec2GT
    LT_diffs = vec1LT - vec2LT

    return (np.average(GT_diffs) + np.average(LT_diffs))/2

def format_data(data1, data2, dim):
    '''
    format_data is used to take data from the laser bench and return specially
    formatted data for use by get_slope and get_shift. This program assumes
    that the [:,1] dimension of the data corresponds to ring radius and that
    the data contains both negative and positive angles.
    '''
    
    (OneInTwo, TwoInOne) = ADE.multi_where(data1[:,0], data2[:,0])

    data1Overlap = data1[OneInTwo]
    data2Overlap = data2[TwoInOne]

    data1GT = data1Overlap[:,dim][np.where(data1Overlap[:,0] >= 0)]
    data1LT = data1Overlap[:,dim][np.where(data1Overlap[:,0] < 0)]
    
    data2GT = data2Overlap[:,dim][np.where(data2Overlap[:,0] >= 0)]
    data2LT = data2Overlap[:,dim][np.where(data2Overlap[:,0] < 0)]

    data1GT.sort()
    data1LT.sort()
    
    data2GT.sort()
    data2LT.sort()

    return (np.array([data1GT, data2GT]), np.array([data1LT,data2LT]))

def chi(r, f, radius, idx):

    fn = np.copy(f)/np.sum(f)
    size = np.float64(r.shape)

    numtest = 1000
    best = 99999999.99
    goodshift = -999
    for i in range(numtest):
        i -= numtest/2
        shift = i/10.
        testcent = min(idx + shift, r.size - 1)
        
        theta = np.arccos(np.abs(r - r[testcent])/radius)
        pn = radius*np.nan_to_num(np.sin(theta))
        
        pn /= np.sum(pn)
        
        nonz = np.where(pn > 0)[0]
        
        chisq = np.sum(((fn[nonz] - pn[nonz])**2)/pn[nonz])
        
        if chisq < best: 
            best = chisq
            goodshift = shift
            
    return goodshift
