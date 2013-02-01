import numpy as np
import glob
import pyfits
import os
import ADEUtils as ADE
import matplotlib
import scipy.misc.pilutil as pl
import time
import scipy.optimize as spo
from multiprocessing import Queue, Process
import sys
import time
plt = matplotlib.pyplot


debug=0
fignum=0
frd_sav=0
_ = np.seterr(invalid='ignore')


def fat_angel(findstr, num_ap, phi, dp, f, exten=0, 
              output=0, pixelpitch=1, fibersize=500):
    '''
    Description:
        fat_angel is deisgned to process a large number of data from
        the laser bench. It takes an input string, findstr, and computes
        the width and radius of the ring in each data image. The images
        are first cleaned using a median subtraction algorithm, which
        assumes that your S/N is very good. Surface brightness
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
        phi    - Float
                 The angle between the screen normal and the detector normal
                 in RADIANS.
        dp     - Float
                 The distance, in millimeters, from the center of the screen
                 to the front glass of the camera lens.
        f      - Float
                 The nominal focal length of the camera lens in millimeters.
                 Read this number off of the lens body.
        exten  - Int
                 The FITS extension where the primary data is stored. As of
                 right now there is no way to specify different extensions
                 for different files.
        output - Str
                 Name of output file. Output contains angle, ring radius, and
                 ring width in column form
        pixelpitch - Float
                 Size of the camera pixels in micrometers.
    
    Output:
        Output is a tuple of Numpy vectors that each contain the folling info:
             Field:       Description:
               0          input angle
               1          ring radius (mm)
               2          ring width (mm)
               3          inner ring radius (mm)
               4          outer ring radius (mm)
               5,6        the x and y coordinates of the ring center (pixels)

    
    Example:
        Assume you have a bunch of data called X_red.fits where X is some
        data iterator.

        >>lb.fat_angel('*_red.fits',150,0.319,1084,26,exten=1,output='my_data.dat')
    '''


    file_list = glob.glob(findstr)
    numfiles = len(file_list)

    'initialize some data arrays'
    widths = np.zeros(numfiles)
    radii = np.zeros(numfiles)
    angles = np.zeros(numfiles)
    r1_vec = np.zeros(numfiles)
    r2_vec = np.zeros(numfiles)
    xcent = np.zeros(numfiles)
    ycent = np.zeros(numfiles)

    t1 = time.time()
    for i in range(numfiles):
        print(file_list[i])


        huds = pyfits.open(file_list[i])
        data = np.float32(huds[exten].data)

        'get the angle from the FITS header. MUCH easier in python than IDL!'
        angles[i] = huds[exten].header['ANGLE']

        '''get rid of the background noise with mediclean. This algorithm
        is very good if you have bodacious S/N'''
        data = ADE.mediclean(data)

        '''given the nominal focal length, f, and object distance, dp, we
        can approximate the image distance, sp, using the thin lens eq.'''
        sp = (1/float(f) - 1/float(dp))**-1

        '''find the center of the image by minimizing the reported ring width'''
        center = cent_test(data,phi,dp,sp,pixelpitch)
#        center = ADE.centroid(data)

        '''t_dist holds the transformation from the detector space to screen space'''
        t_dist = metatron(data,center,phi,dp,sp,pixelpitch)

        'Annulize!'
        if debug:
            t0 = time.time()
            (r_vec,fluxes) = ADE.annulize(data,num_ap,distances=t_dist)
            print "Annulize took "+str(time.time() - t0)+" seconds"
        else:
            (r_vec,fluxes) = ADE.annulize(data,num_ap,distances=t_dist)

        '''find_peak uses the CDF of the fluxes to find the peak
        and interpolation to find the limits of the FWHM. Working with 
        the CDF allows us to assume that the annulus is gaussian-ish 
        without having to worry about the particulars'''
        (rm_idx,r1,r2) = find_peak(r_vec,fluxes)
        
        rm = r_vec[rm_idx]
        
        if debug:
            '''plot a bunch of stuff'''
            fig1 = plt.figure(fignum)
            plt.clf()
            sp0 = fig1.add_subplot(212)
            sp0.plot(r_vec, 
                     np.cumsum(fluxes)/np.max(np.cumsum(fluxes)))
            sp0.axvline(x=rm,ls='--',lw=0.3)
            sp0.axvline(x=r1,ls='--',lw=0.3)
            sp0.axvline(x=r2,ls='--',lw=0.3)
            sp0.set_xlabel("Radius (mm)")
            sp0.set_ylabel("% of total counts")
            sp0.set_title("Normalized CDF")

            sp1 = fig1.add_subplot(211)
            sp1.plot(r_vec,fluxes)
            sp1.axvline(x=rm,ls='--',lw=0.3)
            sp1.axvline(x=r1,ls='--',lw=0.3)
            sp1.axvline(x=r2,ls='--',lw=0.3)
            sp1.set_xlabel("Radius (mm)")
            sp1.set_ylabel("Counts")
            sp1.set_title("Ring Profile")

            plt.suptitle("Angle = "+str(angles[i])+" degrees\n"+file_list[i])
            fig1.show()

            print "Center: "+str(center)

            if numfiles > 1: raw_input("press enter to continue...\n")


        widths[i] = r2 - r1
        radii[i] = rm
        r1_vec[i] = r1
        r2_vec[i] = r2
        xcent[i] = center[0]
        ycent[i] = center[1]

    
    print "Total annulize time was "+str(time.time()-t1)+" seconds"
    
    'We sort the data by angle to make the output a little more readable'
    sort_idx = np.argsort(angles)

    if output:
        f = open(output, 'w')
        f.write('#angle       radius     width     r1        r2       center\n')
        for i in range(angles.shape[0]):
            np.array([angles[sort_idx][i],
                      radii[sort_idx][i],
                      widths[sort_idx][i],
                      r1_vec[sort_idx][i],
                      r2_vec[sort_idx][i],
                      xcent[sort_idx][i],
                      ycent[sort_idx][i]]).\
                      tofile(f, sep='   ',format='%3.4f')
            f.write('\n')
        
    return (angles[sort_idx],
            radii[sort_idx],
            widths[sort_idx],
            r1_vec[sort_idx],
            r2_vec[sort_idx],
            xcent[sort_idx],
            ycent[sort_idx])

def metatron(data, center,  phi, dp, sp, pixelsize):
    '''
    Description:
        metatron uses the relevant physical parameters from the laser bench
        setup and computes a transform that goes from the detector plane
        (what we measure) to the screen plane (what we want to measure).
        It is a helper function only intended to be called by fat_angel.

    Inputs:
        data      - ndarray
                    The data to be transformed. This is only used to make sure the
                    transform has the right dimensions.
        center    - Tuple
                    The center of the ring in pixels
        phi       - Float
                    The angle between the detector normal and the screen normal
                    in radians.
        dp        - Float
                    The distance from the center of the screen to the front glass
                    of the camera lens, in millimeters.
        sp        - Float
                    The image distance given the object distance and lens focal
                    length, in millimeters.
        pixelsize - Float
                    The size of the detector pixels, in microns.

    Output:
        The output is a ndarray that contains the real physical distance of
        each pixel from the center of the ring at the SCREEN PLANE. In other words,
        all the pixels of the same value correspond to a perfect circle on the screen.
    '''

    '''first get the distances in the detector plane'''
    r = ADE.dist_gen(data,center)*pixelsize*1e-3
    alpha = ADE.angle_gen(data,center)

    '''now transform. See my notebook pages 101-103 for an explanation'''
    den = sp**2*(np.cos(2*alpha) - 2*np.cos(alpha)**2*np.cos(2*phi))\
        + 4*r**2 - 3*sp**2
    bignum = 2*dp*r*(\
        (sp-r)*(sp+r)\
            *(3 + np.cos(2*phi) - 2*np.cos(2*alpha)*np.sin(phi)**2)\
            )**0.5
    t_dist = (4*dp*r**2*np.cos(alpha)*np.sin(phi) - bignum)/den
    
    return t_dist

def find_peak(r,f):
    '''
    Description:
        find_peak finds the location of the peak value and the limits of the FWHM
        of the data in f. It assumes that the general shape of f is something
        gaussian-ish (that is, a single, isolated peak with low background).

    Inputs:
        f - ndarray
            The data you want to find peaks and limits of. Should be 1D.
        r - ndarray
            The abcissa for above data. This is used to provide a better
            estimation of the FWHM through interpolation.

    Outputs:
        rp_idx - Int
                 The index that correspond to the peak value of f
        r1     - Float
                 The actual value of the lower limit of the FWHM
        r2     - Float
                 The actual value of the upper limit of the FWHM
    '''
    
    'compute the CDF and normalize it'
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

    '''we do this check becuase otherwise r1 would be set to r[-1] which is just
    about as wrong as you can get'''
    if r1_idx == 0: r1 = r[0]

    return (rp_idx, r1, r2)

def cent_test(data,phi,dp,sp,pixelsize):
    '''
    Description:
        cent_test really trys hard to find the center of a ring image that
        has been disorted because it is on a screen. It uses your favorite
        non-derivative minimization algorithm to find where the FWHM of the
        annulized ring is smallest. The logic here is that as we go off center
        the FWHM gets larger because each TRUE ring annulus is sampled at a
        few differet radii. The reason we can't use a trasform produced by 
        metatron is that the transform itself is heavily dependant on knowing 
        the correct center.

    Inputs:
        data - ndarray
               The ring data
        phi  - Float
               The angle between the detector normal and the screen normal,
               in radians.
        dp   - Float
               The distance from the screen center to the front glass of the
               camera lens, in millimeters.
        sp   - Float
               The image distance of the lense, in millimeters.

    Output:
        A tuple containing the x and y coordinates of the best fit center,
        in pixel units.
    '''
    
    'a centroid is still a good starting point'
    cent = ADE.centroid(data)

    if debug: print "initial center: "+str(cent[0])+','+str(cent[1])

    '''here is where the actual minimization happens, the values for xtol 
    and ftol were chosen with speed in mind (also, we don't really need 
    precision to better than 0.1 pixels).'''
    lstsq = spo.fmin_powell(func,np.array([cent[0],cent[1]])\
                                ,args=(data,phi,dp,sp,pixelsize)\
                                ,xtol=0.01,ftol=0.01)

    return lstsq
    
def func(p,data,phi,dp,sp,pixelsize):
    '''func is the function to be minimized by get_cent. It returns the value
    of the FWHM'''

    '''seems a little dumb to essentially have to run the entire reduction 
    routine every iteration, but that's just life I guess.'''
    t_dist = metatron(data,(p[0],p[1]),phi,dp,sp,pixelsize)

    (r,f) = ADE.annulize(data,90,distances=t_dist)

    (_,r1,r2) = find_peak(r,f)
    diff = r2 - r1

    return diff

def grid_plot(full_dir):
    
    file_list = glob.glob(full_dir+'/*.dat')

    fig = plt.figure(2)
    plt.clf()

    font = {'family' : 'sans-serif',
            'weight' : 'normal',
            'size' : '5'}
    matplotlib.rc('font',**font)

    print 'There are '+str(len(file_list))+' profiles to plot.'
    r = raw_input('How many rows?: ')
    c = raw_input('and how many columns?: ')
    while int(r)*int(c) < len(file_list):
        print 'Fool! '+str(r)+' times '+str(c)+' is less than '\
            +str(len(file_list))+'. Try again.'
        r = raw_input('How many rows?: ')
        c = raw_input('and how many columns?: ')

    pp = 0
    for full in file_list:

        f = open(full,'r')
        a,p = np.loadtxt(f,unpack=True)
        

        sp = fig.add_subplot(r,c,pp+1)
        sp.plot(a,np.log10(p))
        sp.set_ylabel('Log(Counts)')
        sp.set_xlabel('Output angle [deg]')
        sp.set_xlim(0,30)
        sp.set_ylim(0,6)
        sp.text(10,5,full[full.find('/')+1:])
#        sp2 = sp.twiny()
#        sp2.plot(theta,output[i,1]/1e4,alpha=0.0)
#        if np.sign(angles[i]) < 0:
#            a,b = sp2.get_xlim()
#            print (a,b)
#            sp2.set_xlim(b-abs(abs(b)-abs(a)),a-abs(abs(b)-abs(a)))
#        sp2.axvline(x=np.arctan(moms[i,2]/305)*180/np.pi,ls='--',lw=0.3)
        pp += 1

    plt.suptitle(time.asctime(time.localtime()))
    fig.show()

    scratch = raw_input("any key press wipes the plots and exits")
    
    return

if __name__ == '__main__': 
    print "hello"
    main()
