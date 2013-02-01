#! /usr/bin/env python

import numpy as np
from pyraf import iraf
import pyfits
import ADEUtils as ADE
import os
import sys
import glob
from datetime import datetime
import ConfigParser
import matplotlib.pyplot as plt

exten=0

def FReD(direct_image, fiber_image, num_ap, OUTPUT=0, FL=50, FR=4.2):
    '''
    Description:
        FReD is primary reduction tool for the FRD Bench. It takes in two
        images, one for the direct beam and one for the fiber beam, and
        produces outputs suitable for FRD analysis. The main tasks that
        FReD performs are:
            1) Producing curves of growth that show the enclosed energy
            as a function of radius
            
            2) Correcting both beam's data for effects that make the
            direct beam non-ideal
            
            3) Computing effective f-ratios for each beam. The effective
            f-ratio is the f-ratio an ideal lens would have if its enclosed
            energy was EE(r) at r.

     Input:
         direct_image - str
             The name of the FITS file containing the direct beam data
         fiber_image  - stf
             The name of the FITS file containing the fiber beam data
         num_ap - Int
             The number of apertures to use in creating the curve of growth
         exten  - Int
             The FITS extension where the primary data is storred. This
             must be the same for both the direct and fiber beams
         OUTPUT - str
             The name of the output file. FReD produce no output besides
             this file, so if you don't set it then nothing will happen.
         FL     - Float
             The focal length of the lens used (L2)
         FR     - Float
             The focal ratio (f-number) of the lens used

     Output:
         The only output produced by FReD is the output file specified by the
         OUTPUT keyword. This file is a bunch of vectors suitable for plotting
         by an external package, like supermongo. See the output file header
         for more information on the file's contents.
     
     Version History:
         1.0 - 1.12.2011
         1.1 - 1.13.2011 - Changed from appetize to annulize_sb followed
                           by a cumsum.
         1.2 - 1.17.2011 - Changed to a corrected version of annulize (not
                           annulize_sb). The correction accounts for
                           incomplete data at large radii without having
                           to use surface brightness.
         1.3 - 4.28.2011 - Corrected a slight error in the fiber correction
                           calculation. Instead of f_r_c = f_r - abs(x - y)
                           we have f_r_c = f_r - (y - x). This shouldn't
                           make too much of a difference b/c usually y > x.
      
     '''


    version = 1.3

    # Shit's gonna get real real later on with variable names
    # just remember that d_ is for variables relating to the direct beam
    # and f_ is for variables relating to the fiber beam.
    #
    # Comments preceded by '#' relate to variable naming conventions
    
    (d_rvec, d_sb) = ADE.annulize(direct_image,num_ap,EXTEN=exten)
    (f_rvec, f_sb) = ADE.annulize(fiber_image,num_ap,EXTEN=exten)

    '''Turn pixels into an physical length. The SBIG STL-1001E has 24 micron
    square pixels'''
    d_rvec *= 0.024 #-> mm
    f_rvec *= 0.024

    d_flux = np.cumsum(d_sb)
    f_flux = np.cumsum(f_sb)

    '''Now we normalize the fluxes so we are talking about EE, the enclosed
    energy'''
    
    print np.max(d_flux)/np.max(f_flux)
    
    d_EE = d_flux/np.max(d_flux)
    f_EE = f_flux/np.max(f_flux)

    '''Now we need to use the difference between an ideal beam and the direct
    beam (which should be ideal) to create a correction that we can apply
    to both the direct and fiber data'''
    
    # _correction will correspond to corrections that will be applied to
    # get corrected values and _c will correspond to values that have 
    # been corrected
    
    f_r_correction = np.zeros(f_EE.size,dtype=float)
    d_r_c = np.zeros(d_EE.size,dtype=float)

    j=0
    for k in range(f_r_correction.size):
        # Naming conventions here match what is in my notebook on pages
        # 47 through 51
        
        '''First the direct beam'''
        d_r_c[k] = (d_EE[k]*(FL/(2*FR))**2)**0.5

        '''Now the fiber beam'''
        f_r_i = (f_EE[k]*(FL/(2*FR))**2)**0.5
        
        '''find the closest d_EE value that is less than f_EE[k]'''
        while d_EE[j] < f_EE[k]: j += 1
 
        '''interpolate'''
        m = (d_EE[j] - d_EE[j-1])/(d_rvec[j] - d_rvec[j-1])
        r_d = (f_EE[k] - d_EE[j-1])/m + d_rvec[j-1]

        '''f_dr2 is f_dr**2'''
        f_dr2 = r_d**2 - f_r_i**2
        f_r_correction[k] = f_dr2

        '''We do this to fix some weirdness that might happen at really large
        radii. It's more of a visual appeal thing than anything else'''
        if ((f_rvec[k]**2 - f_r_correction[k])**0.5) <\
                ((f_rvec[k-1]**2 - f_r_correction[k-1])**0.5):
            
            f_r_correction[k] = f_r_correction[k-1]

    '''Actually perform the correction on the fiber data'''
    f_r_c = ((f_rvec**2 - (f_r_correction)))**0.5

    #############

    '''For the various plots we want to make we need to have the f-numbers
    which is pretty easy b/c it only depends on radius'''
    # N stands for f-number, a la my notes
    d_N = FL/(2*d_rvec)
    f_N = FL/(2*f_rvec)
    d_N_c = FL/(2*d_r_c)
    f_N_c = FL/(2*f_r_c)

    '''We also need the EFFECTIVE f-number, which is the what the f-number
    of the lens in an ideal system would be if the enclosed energy was 
    EE[k] at r[k]'''
    #_e is for an effective quantity
    d_N_e = d_N*(d_EE)**0.5
    f_N_e = f_N*(f_EE)**0.5
    d_N_e_c = d_N_c*(d_EE)**0.5
    f_N_e_c = f_N_c*(f_EE)**0.5

    ##################
    '''Grab some information about the data for the header'''
    try:
        filt = pyfits.open(direct_image)[exten].header['FILTER']
    except KeyError:
        filt = 'NA'
    try:
        polish = pyfits.open(direct_image)[exten].header['TELESCOP']
    except KeyError:
        polish = 'NA'

    if OUTPUT:
        f = open(OUTPUT,'w')
        f.write('# Generated by FReD v.'+str(version)+'\n'
                +'# Output writen on: '+datetime.now().isoformat(' ')+'\n'
                +'# Input file (direct beam): '+direct_image+'\n'
                +'# Input file (fiber beam): '+fiber_image+'\n'
                +'# Focal length and beam speed: '+str(FL)+'mm '+str(FR)+'\n'
                +'# Filter: '+filt+'\n'
                +'# Polish: '+polish+'\n'
                +'#\n'
                +'# d_r     = aperture radius of direct beam (mm)\n'
                +'# f_r     = aperture radius of fiber beam (mm)\n'
                +'# d_N     = f-ratio of direct beam\n'
                +'# f_N     = f-ratio of fiber beam\n'
                +'# d_r_c   = corrected direct beam radius (mm)\n'
                +'# f_r_c   = corrected fiber beam radius (mm)\n'
                +'# d_N_c   = corrected direct f-ratio\n'
                +'# f_N_c   = corrected fiber f-ratio\n'
                +'# d_EE    = normalized enclosed energy of direct beam\n'
                +'# f_EE    = normalized enclosed energy of fiber beam\n'
                +'# d_N_e   = effective f-ratio for direct beam\n'
                +'# f_N_e   = effective f-ratio for fiber beam\n'
                +'# d_N_e_c = effective f-ratio for corrected direct beam\n'
                +'# f_N_e_c = effective f-ratio for corrected fiber beam\n'
                +'#\n'
                +'#     d_r        f_r        d_N        f_N      d_r_c'
                +'      f_r_c      d_N_c      f_N_c       d_EE       f_EE'
                +'      d_N_e      f_N_e    d_N_e_c    f_N_e_c\n'
                +'#       1          2          3          4          5'
                +'          6          7          8          9         10'
                +'         11         12         13         14\n')
        
        for i in range(min(d_rvec.size,f_rvec.size)):
            np.array([d_rvec[i],
                      f_rvec[i],
                      d_N[i],
                      f_N[i],
                      d_r_c[i],
                      f_r_c[i],
                      d_N_c[i],
                      f_N_c[i],
                      d_EE[i],
                      f_EE[i],
                      d_N_e[i],
                      f_N_e[i],
                      d_N_e_c[i],
                      f_N_e_c[i]]).tofile(f,sep='  ',format='%9.3E')
            f.write('\n')
        f.close()
    return

def friends(findstr, num_ap, outstr, d_suffix='d', f_suffix='f',
            FL=50, FR=4.2, GOGO=False):
    '''
    Description:
        friends is the batch processing front end to FReD. It takes
        in a string describing the reduced data images for multiple
        filters and runs each set of image (direct and fiber) through
        FReD for each filter. Each filter's output is directed to a
        file whos suffix is specified by the user.

    Inputs:
        findstr - str
            The string that specifies what all the filenames have in common.
            Standard wildcards are allowed and indeed necessary.
        num_ap  - Int
            The number of annuli used by FReD in comuting the curve of growth.
        outstr  - str
            The suffix appended to file names when creating the output file.
            See Output section below for usage.
        d_suffix- str
            The suffix of the input files that differentiates data related
            to the direct beam.
        f_suffix- str
            The suffix of the input files that differentiates data related
            to the fiber beam.
        exten   - Int
            The FITS extension where the primary data array is found
        FL      - Float
            The focal length of the lens used on the bench (L2)
        FR      - Float
            The focal ratio (f-number) of the lens used

    Output:
        For each filter friends writes a standard FReD output file. The name
        of this file is the name of the data files, minus either the
        d_ or f_ suffix and appended by the outstr.

    Example:
        Assume we have data for two filter (R and I) across 4 files:
        
        SR_IF42d.fits  SR_RF42d.fits
        SR_IF42f.fits  SR_RF42f.fits

        then running:

        >>friends('SR_?F42*',150,'.dat')

        will result in the following output files:

        SR_IF42.dat  SR_RF42.dat

    '''

    if GOGO: file_list = findstr
    else: file_list = glob.glob(findstr)

    while len(file_list) > 0:
        
        if file_list[0].find(d_suffix+'.fits') > 0:
            name = file_list[0][:file_list[0].find(d_suffix+'.fits')]
        else:
            name = file_list[0][:file_list[0].find(f_suffix+'.fits')]
            
        direct_name = name+d_suffix+'.fits'
        fiber_name = name+f_suffix+'.fits'

        print '\nDirect file: '+direct_name
        print 'Fiber file:  '+fiber_name

        FReD(direct_name,fiber_name,num_ap,FL=FL,FR=FR,
             OUTPUT=name+outstr)
        
        file_list.remove(direct_name)
        file_list.remove(fiber_name)

    return

def stew(raw_directory, num_ap, prefilter='_',EXTEN=0,CLEANUP=True, 
         GOGO=False, FL=50, FR=4.2):
    '''
    Description:
        stew is the highest level reduction program for the FRD Bench. In
        its simplist form it takes in raw bench data and performs the
        standard IRAF reduction that produces two images per filter
        (direct and fiber). It can also take these reduced images and
        run them through FReD to produce data suitable for plotting and 
        publication. It is THE FRD reduction program.

    Inputs:
        raw_directory - str
            the path to the directory where the raw bench files are stored
            stew assumes that the dark frames have DARK somewhere in their
            file names (they do by default if you are using CCDSoft).
        num_ap - Int
            The number of apertures to use when constructing the curve of growth
        pre_filter - str
            the character or sequence of characters that comes right before
            the filter indicator in the data file names. stew assumes that
            the filter indicators are all one character long
        exten - Int
            the fits extension where the primary data live
        CLEANUP - Bool
            Set to False if you want to keep all the intermediate IRAF images,
            i.e. the combined darks and the dark subtracted data image
        GOGO - Bool
            Set to True if you want stew to continue the analysis after the
            final, reduced images are produced. This essential triggers a call
            to FReD through friends.

    Example:
        In the shell:
        
        : pwd
        /home/users/astro/data/reduced
        :ls ../Raw
        DARK_8sec.00001.FIT
        ...
        SR_vF42d.00001.FIT
        SR_vF42f.00001.FIT
        ...
        
        In python
        >>> stew('../Raw',150,GOGO=True)
        
        Shell:
        :ls
        SR_vF42d.fits
        SR_vF42f.fits
        SR_vF42.dat
        ...

        '''
    
    data_dir = config.get('Data','raw_directory')
    dark_dir = config.get('Data','dark_directory')
    usehead = config.getboolean('Data','use_header_info')
    fr = config.getfloat('Data','f-ratio')
    fl = config.getfloat('Data','focal_length')
    gogo = config.getboolean('Options','gogo')
    clean = config.getboolean('Options','cleanup')
    prefilter = config.get('Data','prefilter')
    exten = config.getint('Data','fits_exten')
    data_string = config.get('Data','data_string')
    num_ap = config.get('Options','num_ap')
        
    data_dict = get_data(data_dir, data_string)
    dark_dict = get_darks(dark_dir)

    print "Finding data..."
    needed_darks = data_dict.keys()
    for k in dark_dict.keys(): needed_darks.remove(k)

    if len(needed_darks) != 0:
        print "Generating darks..."
        dark_dict = gen_darks(dark_dict, needed_darks, dark_dir)
        
    print "\nSubtracting darks..."
    rm_dark(data_dict, dark_dict)
    
    print "Sorting data..."
    focal_dict = sort_data(data_dir)

    print "Combining images..."
    name_dict = combine(focal_dict, data_string)

    if clean:
        print "Cleaning up intermediate steps..."
        os.system('rm ./*_ds.fits')

    print name_dict

    if gogo:
        print "Starting analysis..."
        for f_ratio in name_dict.keys():
            friends(name_dict[f_ratio], num_ap, '.dat', FR=f_ratio, GOGO=True)
        os.system('cp ~/snakes/plot_FRD.sm .')
    
    print "\nAll done!"
 
    return

def get_data(raw_dir, data_string):

    data_list = glob.glob(raw_dir+'/'+data_string)
    exp_times = {}

    for data in data_list:

        exptime = pyfits.open(data)[exten].header['EXPTIME']

        if exptime not in exp_times.keys():
            exp_times[exptime] = []
        
        exp_times[exptime].append(data)

    return exp_times

def get_darks(dark_dir):

    dark_list = glob.glob(dark_dir+'/Combined*DARK.fits')
    dark_dict = {}
    
    for dark in dark_list:
        exptime = pyfits.open(dark)[exten].header['EXPTIME']
        dark_dict[exptime] = dark
    
    return dark_dict

def gen_darks(dark_dict, needs, dark_dir):
    
    dark_list = glob.glob(dark_dir+'/DARK*.FIT')
    dark_dict = {}

    for dark in dark_list:
        exptime = pyfits.open(dark)[exten].header['EXPTIME']

        if exptime not in dark_dict.keys():
            dark_dict[exptime] = []

        dark_dict[exptime].append(dark)

    iraf.imcombine.combine = 'median'
    iraf.imcombine.reject = 'none'

    for exp in needs:
        iraf.imcombine(','.join(dark_dict[exp]),\
                           dark_dir+'/Combined_'+str(exp)+'_DARK.fits')
        dark_dict[exp] = dark_dir+'/Combined_'+str(exp)+'_DARK.fits'

    return dark_dict

def rm_dark(data_dict, dark_dict):
    
    for exp in data_dict.keys():
        inputlist = data_dict[exp]
        
        if len(inputlist) > 10:
            while len(inputlist) > 0:
                inputstring = ','.join(inputlist[:10])
                outputstring = ','.join([s[s.rfind('/')+1:s.find('.FIT')]+\
                                             '_ds.fits'\
                                             for s in inputlist[:10]])
                iraf.imarith(inputstring,'-',\
                                 dark_dict[exp],outputstring)
                inputlist = inputlist[10:]
                
        else:
            iraf.imarith(','.join(inputlist),\
                             '-',\
                             dark_dict[exp],\
                             ','.join([s[s.rfind('/')+1:s.find('.FIT')]+\
                                           '_ds.fits' for s in inputlist]))
    return

def sort_data(data_dir):

    ds_list = glob.glob('*_ds.fits')
    focal = {}
    
    for ds in ds_list:
        focal_length = pyfits.open(ds)[exten].header['FOCALLEN']
        diameter = float(pyfits.open(ds)[exten].header['APTDIA'])
        focal_ratio = focal_length/diameter
        filt = pyfits.open(ds)[exten].header['FILTER']
        ftype = pyfits.open(ds)[exten].header['OBSERVER']

        if focal_ratio not in focal.keys():
            focal[focal_ratio] = {}
        if filt not in focal[focal_ratio].keys():
            focal[focal_ratio][filt] = {}
        if ftype not in focal[focal_ratio][filt].keys():
            focal[focal_ratio][filt][ftype] = []
        
        focal[focal_ratio][filt][ftype].append(ds)


    return focal

def combine(fdict, data_string):

    iraf.imcombine.combine = 'average'
    iraf.imcombine.reject = 'avsigclip'
    iraf.imcombine.lsigma = 3

    name_dict = {}
    if data_string[data_string.find('*')-1] == '_':
        fill = ''
    else:
        fill = '_'

    for f_ratio in fdict.keys():
        name_dict[f_ratio] = []
        
        for filt in fdict[f_ratio].keys():
            name = data_string[:data_string.find('*')]+fill+filt+'F'+\
                str(int(f_ratio*10))
            for ftype in fdict[f_ratio][filt].keys():
                iraf.imcombine(','.join(fdict[f_ratio][filt][ftype]),\
                                   name+ftype[0]+'.fits')
                
                name_dict[f_ratio].append(name+ftype[0]+'.fits')

    return name_dict

def main():
    
    options = ConfigParser.ConfigParser()
    options.read(sys.argv[1])
    
    global exten
    exten = options.getint('Data','fits_exten')

    stew(options)
    
    return
    
if __name__ == '__main__':
    sys.exit(main())

def clean_head(string, keyword):

    dlist = glob.glob(string)

    for data in dlist:
        dd = pyfits.open(data)[0]
        dd.header.update('APTDIA',keyword)
        dd.writeto(data+'.clean')
    return
        

def tput(findstr,xlim,ylim,exten=0,offset=0):
    
    file_list = glob.glob(findstr)

    time_list = []
    mean_list = []
    std_list = []
    
    for image in file_list:
        

        HDU = pyfits.open(image)[exten]
        data = np.float64(HDU.data)
        to = HDU.header['TIME-OBS']

        time = float(to[:to.find(':')])*3600\
            +float(to[to.find(':')+1:to.rfind(':')])*60\
            +float(to[to.rfind(':')+1:])

        mean = np.mean(data[ylim[0]-1:ylim[1]-1,xlim[0]-1:xlim[1]-1])
        std = data[ylim[0]-1:ylim[1]-1,xlim[0]-1:xlim[1]-1].std()

        print image +'    '+str(std)

        time_list.append(time)
        mean_list.append(mean)
        std_list.append(std)


    time_vec = np.array(time_list)
    mean_vec = np.array(mean_list)
    std_vec = np.array(std_list)

    time_vec -= min(time_vec)
    time_vec += offset
#    time_vec /= 60
    p = np.polyfit(time_vec,mean_vec,3)
    fit = p[0]*time_vec**3 + p[1]*time_vec**2 + p[2]*time_vec + p[3]
#    fit = p[0]*time_vec + p[1]

    fig = plt.figure(0)
    fig.clf()
    ax = fig.add_subplot(211)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Counts')
    ax.plot(time_vec,mean_vec,'.k')
    ax.plot(time_vec,fit,'-b')
#    ax.errorbar(time_vec,mean_vec,fmt='.k',yerr=std_vec)

    ax2 = fig.add_subplot(212)
    ax2.plot(time_vec,fit-mean_vec,'.k')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Residuals')
    ax2.set_title('$\sigma$= '+str(np.std(mean_vec - fit)))

    fig.show()

    return (time_vec,mean_vec)
