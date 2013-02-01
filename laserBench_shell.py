#! /usr/bin/env python

import numpy as np
import pyfits
import ADEUtils as ADE
import scipy.optimize as spo
from multiprocessing import Queue, Process, Lock
import sys
from datetime import datetime
import time
import ConfigParser
import os

version = 1.2
full_output = False
full_dir = './'

def fat_angel(name, num_ap, phi, dp, f, fsg, pixelpitch, exten, q, lock):
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

    print name

    huds = pyfits.open(name)
    data = np.float32(huds[exten].data)

    'get the angle from the FITS header. MUCH easier in python than IDL!'
    angle = huds[exten].header['ANGLE']

    '''get rid of the background noise with mediclean. This algorithm
    is very good if you have bodacious S/N'''
    data = ADE.mediclean(data)

    '''given the nominal focal length, f, and object distance, dp, we
    can approximate the image distance, sp, using the thin lens eq.'''
    sp = (1/float(f) - 1/float(dp))**-1

    '''find the center of the image by minimizing the reported ring width'''
    center = cent_test(data,phi,dp,sp,pixelpitch) #slow and right
#    center = ADE.centroid(data) #fast and wrong

    '''t_dist holds the transformation from the detector space to 
    screen space'''
    t_dist = metatron(data,center,phi,dp,sp,pixelpitch)

    'Annulize!'
    (r_vec,fluxes,errors) = ADE.annulize(data,num_ap,distances=t_dist)

    '''find_peak uses the CDF of the fluxes to find the peak
    and interpolation to find the limits of the FWHM. Working with 
    the CDF allows us to assume that the annulus is gaussian-ish 
    without having to worry about the particulars'''
    (rp,r1,r2) = find_peak(r_vec,fluxes)
        
    ap = np.arctan(rp/fsg)*180/np.pi
    a1 = np.arctan(r1/fsg)*180/np.pi
    a2 = np.arctan(r2/fsg)*180/np.pi
        
    awidth = a2 - a1
    rwidth = r2 - r1
    xcent = center[0]
    ycent = center[1]

    if full_output:
        'write power profile info for use by theModule'
        if not os.path.exists(full_dir): os.makedirs(full_dir)
        angs = np.arctan(r_vec/fsg)*180/np.pi
        write_full(angs,fluxes,angle,ap)

    q.put((angle,ap,rp,awidth,rwidth,r1,r2,xcent,ycent))
    
    return

def metatron(data, center,  phi, dp, sp, pixelsize):
    '''
    Description:
        metatron uses the relevant physical parameters from the laser bench
        setup and computes a transform that goes from the detector plane
        (what we measure) to the screen plane (what we want to measure).
        It is a helper function only intended to be called by fat_angel.

    Inputs:
        data      - ndarray
                    The data to be transformed. This is only used to make 
                    sure the transform has the right dimensions.
        center    - Tuple
                    The center of the ring in pixels
        phi       - Float
                    The angle between the detector normal and the screen normal
                    in radians.
        dp        - Float
                    The distance from the center of the screen to the front 
                    glass of the camera lens, in millimeters.
        sp        - Float
                    The image distance given the object distance and lens focal
                    length, in millimeters.
        pixelsize - Float
                    The size of the detector pixels, in microns.

    Output:
        The output is a ndarray that contains the real physical distance of
        each pixel from the center of the ring at the SCREEN PLANE. In other 
        words, all the pixels of the same value correspond to a perfect 
        circle on the screen.
    '''

    '''first get the distances in the detector plane'''
    r = ADE.dist_gen(data,center)*pixelsize
    alpha = ADE.angle_gen(data,center)

    '''now transform. See my notebook pages 101-103 for an explanation'''
#    den = sp**2*(np.cos(2*alpha) - 2*np.cos(alpha)**2*np.cos(2*phi))\
#        + 4*r**2 - 3*sp**2
#    bignum = 2*dp*r*(\
#        (sp-r)*(sp+r)\
#            *(3 + np.cos(2*phi) - 2*np.cos(2*alpha)*np.sin(phi)**2)\
#            )**0.5
#    t_dist = (4*dp*r**2*np.cos(alpha)*np.sin(phi) - bignum)/den

    den = sp*(np.cos(alpha)**2*np.cos(phi)**2 + np.sin(alpha)**2)**0.5 \
        + r*np.cos(alpha)*np.sin(phi)

    t_dist = dp*r/den
    
    return t_dist

def find_peak(r,f):
    '''
    Description:
        find_peak finds the location of the peak value and the limits of the
        FWHM of the data in f. It assumes that the general shape of f is 
        something gaussian-ish (that is, a single, isolated peak with low 
        background).

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
    rp_idx = np.interp(0.5,CDF,np.arange(CDF.size))
    rp = np.interp(rp_idx,np.arange(r.size),r)
    peak_val = np.interp(rp,r,f)


    r1_idx = np.where(f >= 0.5*peak_val)[0][0]
    r2_idx = np.where(f >= 0.5*peak_val)[0][-1]

    '''Now intertpolate to get a better width estimate'''
    m1 = (f[r1_idx] - f[r1_idx - 1])/(r[r1_idx] - r[r1_idx - 1])
    r1 = r[r1_idx - 1] + (0.5*peak_val - f[r1_idx - 1])/m1
    m2 = (f[r2_idx + 1] - f[r2_idx])/(r[r2_idx + 1] - r[r2_idx])
    r2 = r[r2_idx] + (0.5*peak_val - f[r2_idx])/m2

    '''we do this check becuase otherwise r1 would be set to r[-1] which 
    is just about as wrong as you can get'''
    if r1_idx == 0: 
        r1 = r[0]

    return (rp, r1, r2)

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

    (r,f,_) = ADE.annulize(data,90,distances=t_dist)

    (_,r1,r2) = find_peak(r,f)
    diff = r2 - r1

    return diff

def write_output(q,out_name,out_info):
    '''write the final data file. We have this as separate function because
    we need to wait for all the threads to finish before we write the output.
    '''
    angles = np.array([])
    aradii = np.array([])
    awidths = np.array([])
    radii = np.array([])
    widths = np.array([])
    r1_vec = np.array([])
    r2_vec = np.array([])
    xcent = np.array([])
    ycent = np.array([])
    
    'now fill them in'
    while q.qsize() != 0:
        pickel = q.get()
        angles = np.append(angles,pickel[0])
        aradii = np.append(aradii,pickel[1])
        radii = np.append(radii,pickel[2])
        awidths = np.append(awidths,pickel[3])
        widths = np.append(widths,pickel[4])
        r1_vec = np.append(r1_vec,pickel[5])
        r2_vec = np.append(r2_vec,pickel[6])
        xcent = np.append(xcent,pickel[7])
        ycent = np.append(ycent,pickel[8])
        
    sort_idx = np.argsort(angles)
    
    (sgd,fsg,foc,phi,num_ap,ini_name) = out_info

    if sgd % 1 == 0.0: sgd = int(sgd)
    if fsg % 1 == 0.0: fsg = int(fsg)
    if foc % 1 == 0.0: foc = int(foc)

    f = open(out_name, 'w')
    f.write('# Generated by fat_angel v.'+str(version)+'\n'
            +'# Output written on '+datetime.now().isoformat(' ')+'\n'
            +'# Input INI file: '+ini_name+'\n'
            +'# Screen-glass distance: '+str(sgd)+' mm\n'
            +'# Fibre-screen distance: '+str(fsg)+' mm\n'
            +'# Nominal lens focal length: '+str(foc)+' mm\n'
            +'# Fibre-detector angle: '+str(phi)+' rad = '
            +str(phi*180/np.pi)+' deg\n'
            +'# Number of sampling annuli: '+str(num_ap)+'\n'
            +'#\n'
            +'# a    = laser input angle (deg)\n'
            +'# r_a  = ring radius (deg)\n'
            +'# r    = ring radius on screen (mm)\n'
            +'# w_a  = ring width (deg)\n'
            +'# w    = ring width on screen (mm)\n'
            +'# r1   = inner ring radius (mm)\n'
            +'# r2   = outer ring radius (mm)\n'
            +'# c_x  = x center of ring on detector (pix)\n'
            +'# c_y  = y center of ring on detector (pix)\n'
            +'#\n'
            +str('# {0:<10}{1:<11}{2:<11}{3:<11}{4:<11}{5:<11}'
                 +'{6:<11}{7:<11}{8:<11}')\
                .format('a','r_a','r','w_a','w','r1','r2','c_x','c_y')+'\n'
            +str('# {0:<10}{1:<11}{2:<11}{3:<11}{4:<11}{5:<11}'
                 +'{6:<11}{7:<11}{8:<11}').format(1,2,3,4,5,6,7,8,9)+'\n')


    for i in range(angles.size):
        np.array([angles[sort_idx][i],
                  aradii[sort_idx][i],
                  radii[sort_idx][i],
                  awidths[sort_idx][i],
                  widths[sort_idx][i],
                  r1_vec[sort_idx][i],
                  r2_vec[sort_idx][i],
                  xcent[sort_idx][i],
                  ycent[sort_idx][i]]).\
                  tofile(f, sep='  ',format='% -9.4f')
        f.write('\n')
    f.close()
    return

def write_full(ang,p,in_angle,out_angle):
    name = full_dir+'/full_'+str(in_angle)+'.dat'
    f = open(name, 'w')
    f.write('# Written on '+datetime.now().isoformat(' ')+'\n'
            +'# Input angle: '+str(in_angle)+'\n'
            +'# Output angle: '+str(out_angle)+'\n'
            +'#\n'
            +'# {0:<8}{1:<11}\n'.format('Angle','Counts')
            +'# {0:<8}{1:<11}\n'.format('[deg]','[ADU]'))
    for i in range(ang.size):
        np.array([ang[i],p[i]]).tofile(f, sep='  ',format='%7.4f')
        f.write('\n')
    f.close()
    return

def main():
    '''reads in the .ini file and sets up all the processes that parallelizes
    reduction'''
    
    options = ConfigParser.ConfigParser()
    options.read(sys.argv[1])

    sgd = options.getfloat('Data','screen_glass_distance')
    fsg = options.getfloat('Data','fiber_screen_distance')
    f = options.getfloat('Data','lens_focal_length')
    phi = options.getfloat('Data','fiber_detector_angle')
    pixelsize = options.getfloat('Data','pixel_size')
    fits_exten = options.getint('Data','fits_exten')
    num_ap = options.getint('Options','num_an')
    out_name = options.get('Options','output_name')
    max_proc = options.getint('Options','max_proc')
    global full_output
    global full_dir
    full_output = options.getboolean('Options','full_output')
    full_dir = options.get('Options','full_directory')

    print "max_proc = "+str(max_proc)

    file_list = sys.argv[2:]

    jobs = []

    '''I don't think we need this lock anymore'''
    lock = Lock()
    resultq = Queue()

    for name in file_list:
        args = (name,num_ap,phi,sgd,f,fsg,pixelsize,fits_exten,resultq,lock)
        jobs.append(Process(target=fat_angel,args=args))

    '''now jobs has a bunch of processes in it ready to be run'''

    '''this loop ensures that we only have max_proc jobs running at once'''
    i = 0
    while i < len(jobs):
        while [x.is_alive() for x in jobs].count(True) < max_proc:
            jobs[i].start()
            xi += 1
            if i >= len(jobs): break
        '''without this sleep statement the program would check how many
        jobs are alive as fast as it could, which would take 100% of a
        processor'''
        time.sleep(0.01)

    '''this loop makes sure all the jobs are done before continuing'''
    while any([x.is_alive() for x in jobs]):
        time.sleep(0.01)
    
    output_info = (sgd, fsg, f, phi, num_ap, sys.argv[1])
    write_output(resultq,out_name,output_info)

    '''just a cleanup step. Not strictly necessary'''
    for p in jobs: p.join()

    return

if __name__ == '__main__': 
    sys.exit(main())
