#! /usr/bin/env python

import numpy as np
from pyraf import iraf
import pyfits
try:
    import ADEUtils as ADE
    SLOWAN = False
except ImportError:
    print "WARNING: Could not load ADEUtils, falling back on slow version of annulize"
    SLOWAN = True
import os
import sys
import glob
from datetime import datetime
import ConfigParser
import matplotlib.pyplot as plt
import pickle

debug = False

def FReD(direct_image, fiber_image, num_ap, pot, filt, dir_cut,\
             EXTEN=0, OUTPUT=0, FL=50, FR=4.2):
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
         EXTEN  - Int
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
         1.3 - 10.2011   - Added FRD metrics
     '''


    version = 1.3
    # Shit's gonna get real real later on with variable names
    # just remember that d_ is for variables relating to the direct beam
    # and f_ is for variables relating to the fiber beam.
    #
    # Comments preceded by '#' relate to variable naming conventions
    
    direct_HDU = pyfits.open(direct_image)[EXTEN]
    fiber_HDU = pyfits.open(fiber_image)[EXTEN]
    
    d_exptime = direct_HDU.header['EXPTIME']
    f_exptime = fiber_HDU.header['EXPTIME']    
    
    direct = direct_HDU.data
    fiber = fiber_HDU.data
    
    direct[np.where(direct < dir_cut)] = 0.0
#    fiber[np.where(fiber < 300)] = 0.0
    
    if pot:
        direct_time_str = direct_HDU.header['TIME-OBS']
        fiber_time_str = fiber_HDU.header['TIME-OBS']

        direct_time = np.float(direct_time_str[6:])\
            + np.float(direct_time_str[3:5])*60.\
            + np.float(direct_time_str[0:2])*3600.

        fiber_time = np.float(fiber_time_str[6:])\
            + np.float(fiber_time_str[3:5])*60.\
            + np.float(fiber_time_str[0:2])*3600.
        
        fcorrect = pot.get_correction(fiber_time,filt)
        dcorrect = pot.get_correction(direct_time,filt)
        print '   Direct throughput correction is '+str(dcorrect)
        print '   Fiber throughput correction is '+str(fcorrect)
        
        direct *= dcorrect
        fiber *= fcorrect

    if debug: print '    Annulizing...'
    d_rvec, d_sb, d_sberr = annulize(direct,num_ap)
    f_rvec, f_sb, f_sberr = annulize(fiber,num_ap)

    if debug:
        plt.clf()
        fig = plt.figure(0)
        ax1 = fig.add_subplot(221)
        ax1.plot(d_rvec,d_sb)
        ax1.plot(f_rvec,f_sb)
        fig.show()
        raw_input('     hit enter or something')

    '''Turn pixels into an physical length. The SBIG STL-1001E has 24 micron
    square pixels'''
    d_rvec *= 0.024 #-> mm
    f_rvec *= 0.024

    if debug: print '    Cumsumming...'
    d_flux = np.cumsum(d_sb)
    f_flux = np.cumsum(f_sb)
    d_ferr = (np.cumsum(d_sberr**2))**0.5
    f_ferr = (np.cumsum(f_sberr**2))**0.5

    if debug:
        ax2 = fig.add_subplot(222)
        ax2.plot(d_rvec, d_flux)
        ax2.plot(f_rvec, f_flux)
        fig.show()
        raw_input('     Cumsum')
        

    '''Now we normalize the fluxes so we are talking about EE, the enclosed
    energy'''
    if debug: print '    EEing...'
    d_max = np.max(d_flux)
    f_max = np.max(f_flux)
    d_EE = d_flux/d_max
    f_EE = f_flux/f_max
#    d_fmaxerr = d_ferr[np.where(d_flux == d_max)[0]]
#    f_fmaxerr = f_ferr[np.where(f_flux == f_max)[0]]
#    d_EEerr = ((d_ferr/d_max)**2 + (d_fmaxerr*d_flux/(d_max**2))**2)**0.5
#    f_EEerr = ((f_ferr/f_max)**2 + (f_fmaxerr*f_flux/(f_max**2))**2)**0.5

    '''Now we need to use the difference between an ideal beam and the direct
    beam (which should be ideal) to create a correction that we can apply
    to both the direct and fiber data'''
    
    # _correction will correspond to corrections that will be applied to
    # get corrected values and _c will correspond to values that have 
    # been corrected
    
    f_r_correction = np.zeros(f_EE.size,dtype=float)
    d_r_c = np.zeros(d_EE.size,dtype=float)

    if debug: print '    Correcting...'
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

#        if debug: print (f_rvec[k]**2 - f_r_correction[k])**0.5, f_rvec[k]**2, f_r_correction[k]

        '''We do this to fix some weirdness that might happen at really large
        radii. It's more of a visual appeal thing than anything else'''
        if (np.abs(f_rvec[k]**2 - f_r_correction[k]))**0.5 <\
                (np.abs(f_rvec[k-1]**2 - f_r_correction[k-1]))**0.5 \
                or (f_rvec[k]**2 - f_r_correction[k]) < 0:
#            if debug: print ' here'
            f_r_correction[k] = f_r_correction[k-1]

    '''Actually perform the correction on the fiber data'''
    f_r_c = (f_rvec**2 - f_r_correction)**0.5

    d_rerr = (np.abs(d_r_c - d_rvec))**0.5
    f_rerr = (np.abs(f_r_c - f_rvec))**0.5
    
    if debug: print (np.abs(f_r_c - f_rvec))**0.5

    ################

    if debug: print '    Computing...'
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

    '''Before we go any further we will compute some metrics of throughput.
    For this we'll first need to put fluxes in units of counts/time so that
    different exposure times do not affect the result'''
    
    d_ADUrate = d_flux/d_exptime
    f_ADUrate = f_flux/f_exptime

    '''We will compute the relative ADU/time at various percentages of the
    total'''
    
    tput_100 = float(f_ADUrate.max() / d_ADUrate.max())

    rf5 = FL/10.
    rf4 = FL/8.
    if rf5 > d_r_c.max():
        d_ADUf5 = d_ADUrate.max()
    else: d_ADUf5 = np.interp(rf5,d_r_c,d_ADUrate)
        
    if rf4 > d_r_c.max():
        d_ADUf4 = d_ADUrate.max()
    else: d_ADUf4 = np.interp(rf4,d_r_c,d_ADUrate)
    
    tput_f5 = np.interp(rf5,f_r_c,f_ADUrate)/d_ADUf5
    tput_f4 = np.interp(rf4,f_r_c,f_ADUrate)/d_ADUf4
    tput_f5b = np.interp(rf5,f_r_c,f_EE)
    tput_f4b = np.interp(rf4,f_r_c,f_EE)
    
    d_r_90 = np.interp(0.9, d_EE, d_r_c)
    tput_90 = np.interp(d_r_90, f_r_c, f_ADUrate) / (d_ADUrate.max()*0.9)
    
    d_r_80 = np.interp(0.8, d_EE, d_r_c)
    tput_80 = np.interp(d_r_80, f_r_c, f_ADUrate) / (d_ADUrate.max()*0.8)

    r_ideal = FL/(2*FR)
    r_ideal_test = d_r_c[np.where(f_ADUrate == f_ADUrate.max())[0]]
    metric = np.interp(r_ideal, f_r_c, f_EE)
    metric80 = FL/(2*np.interp(0.8,f_EE,f_r_c))
    metric90 = FL/(2*np.interp(0.9,f_EE,f_r_c))
    
    if debug:
        ax3 = fig.add_subplot(223)
        ax3.errorbar(f_r_c,f_ADUrate,yerr=f_EEerr/f_exptime,xerr=f_rerr,fmt='g')
        ax3.errorbar(d_r_c,d_ADUrate,yerr=d_EEerr/d_exptime,xerr=d_rerr,fmt='b')
        ax3.axvline(ls='--',x=rf5)
        ax3.axvline(x=d_r_90)
        ax3.axhline(ls='--',y=np.interp(rf5,f_r_c,f_ADUrate))
        ax3.axhline(ls='--',y=np.interp(rf5,d_r_c,d_ADUrate))
        fig.show()
        print "    Fiber value at f/5 ("+str(rf5)+"):  "+str(np.interp(rf5,f_r_c,f_ADUrate))
        print "    Direct value at f/5 ("+str(rf5)+"): "+str(np.interp(rf5,d_r_c,d_ADUrate))
        print "    Full direct is:  "+str(float(d_ADUrate.max()))
        print "    Max d_r_c value is:  "+str(d_r_c.max())
        print "    r_ideal is:  "+str(r_ideal)
        print "    r_max is:    "+str(r_ideal_test)
        raw_input("    ADUrate")

    ##################
    '''Grab some information about the data for the header'''
    try:
        filt = pyfits.open(direct_image)[EXTEN].header['FILTER']
    except KeyError:
        filt = 'NA'
    try:
        fiber = pyfits.open(direct_image)[EXTEN].header['OBSERVER']
    except KeyError:
        fiber = 'NA'
    try:
        polish = pyfits.open(direct_image)[EXTEN].header['TELESCOP']
    except KeyError:
        polish = 'NA'

    if OUTPUT:
        f = open(OUTPUT,'w')
        f.write('# Generated by FReD v.'+str(version)+'\n'
                +'# Output writen on: '+datetime.now().isoformat(' ')+'\n'
                +'# Input file (direct beam): '+direct_image+'\n'
                +'# Input file (fiber beam): '+fiber_image+'\n'
                +'# Number of apertures: '+str(num_ap)+'\n'
                +'# Focal length and beam speed: '+str(FL)+'mm '+str(FR)+'\n'
                +'# Filter: '+filt+'\n'
                +'# Fiber: '+fiber+'\n'
                +'# Polish: '+polish+'\n'
                +'# Throughput: '+str(tput_100)+'\n'
                +'# Direct image cutoff: '+str(dir_cut)+'\n'
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
        
    return (metric,metric90,metric80,tput_100,tput_90,tput_80,
            tput_f5,tput_f4,tput_f5b,tput_f4b)

def soba(nood,num_ap,dir_cut,exten,pot,mfile):

    '''
    Description:
        soba takes a reduced Noodle object and runs it through FReD for
        analysis. The Noodle must have been reduced or the soba will not
        taste good! All inputs are handled by main as soba is not intended
        to be run alone.

    Inputs:
        nood   - The reduced Noodle object
        num_ap - The number of apertures to use when reducing the data
        exten  - The FITS extension that holds the data

    Output:
        soba runs FReD on every filter/f-ratio combination found in nood.
        As a result, the output is a bunch of .dat files produced by FReD.

    '''
    
    print '\nAnalyzing data...'
    f = open(mfile,'wb')
    f.write('# Output writen on: '+datetime.now().isoformat(' ')+'\n'
            +'# Reduction directory: '+os.getcwd()+'\n'
            +'# INI file: '+sys.argv[1]+'\n'
            +'#\n'
            +'# {:4}= '.format('FR')+'f-ratio\n'
            +'# {:4}= '.format('filt')+'filter\n'
            +'# {:4}= '.format('m')+'fiber EE at f/6.3\n'
            +'# {:4}= '.format('N90')+'fiber f/# at EE90\n'
            +'# {:4}= '.format('N80')+'fiber f/# at EE80\n'
            +'# {:4}= '.format('tput')+'total throughput\n'
            +'# {:4}= '.format('tput90')+'throughput at direct beam EE90\n'
            +'# {:4}= '.format('tput80')+'throughput at direct beam EE80\n'
            +'# {:4}= '.format('tput5')+'throughput at f/5\n'
            +'# {:4}= '.format('tput4')+'throughput at f/4\n'
            +'# {:4}= '.format('EE5')+'fiber EE at f/5\n'
            +'# {:4}= '.format('EE4')+'fiber EE at f/4\n'
            +'#\n'
            +str('#{0:>8}{1:>9}{2:>9}{3:>9}{4:>9}{5:>9}{6:>9}'
                 +'{7:>9}{8:>9}{9:>9}{10:>9}{11:>9}\n')\
                .format('FR','filt','m','N90','N80','tput','tput90','tput80',\
                            'tput5','tput4','EE4','EE5')
            +str('#{0:>8}{1:>9}{2:>9}{3:>9}{4:>9}{5:>9}{6:>9}'
                 +'{7:>9}{8:>9}{9:>9}{10:>9}{11:>9}\n')\
                .format(1,2,3,4,5,6,7,8,9,10,11,12))

    for f_ratio in nood.keys():

        print 'F'+str(f_ratio)
        focal_length = nood[f_ratio]['focal_length']
        print 'Focal length = '+str(focal_length)+'mm'
    
        for filt in nood[f_ratio]['data'].keys():
            direct_name = nood[f_ratio]['data'][filt]['direct']['final']
            fiber_name = nood[f_ratio]['data'][filt]['fiber']['final']
            name = fiber_name[:fiber_name.rfind('f.')]

            print ' Filter: '+str(filt)
            print '  Direct image is '+direct_name
            print '  Fiber image is  '+fiber_name

            metric = FReD(direct_name,fiber_name,num_ap,pot,filt,dir_cut,EXTEN=exten,\
                     FL=focal_length,FR=f_ratio,OUTPUT=name+'.dat')
            
            f.write('{0:9.1f}{1:>9}'.format(f_ratio,filt))
            for m in metric: f.write('{:9.4f}'.format(m))
            f.write('\n')

            '''now generate the FRD plots using supermongo'''
            sm = open('tmp_'+name+'.sm','wb')
            sm.write('verbose 0\n'
                     +'macro read plot_FRD.sm\n'
                     +'device postencap_color '+name+'.ps\n'
                     +'plot_FRD '+name+'.dat\ndevice nodevice\n')
            sm.close()
            os.system('sm < tmp_'+name+'.sm')
            os.system('convert -density 200 '+name+'.ps -quality 92 '+name+'.jpg')
            os.system('rm tmp_'+name+'.sm')

    f.close()
    return

def main():
    
    options = ConfigParser.ConfigParser()
    options.read(sys.argv[1])
    resume = options.getboolean('Options','resume')

    if resume:
        noodfile = options.get('Options','noodle_file')
        print 'Loading previous data run from '+noodfile
        N = pickle.load(open(noodfile,'rb'))
    else:
        N = Noodle(options)
        N.build_run()
        noodsave = options.get('Options','noodle_save')
        print 'Saving data run to '+noodsave
        pickle.dump(N,open(noodsave,'wb'))

    T = thePot(options)

    num_ap = options.getint('Options','num_ap')
    dir_cut = options.getint('Options','direct_cutoff')
    gogo = options.getboolean('Options','gogo')
    mfile = options.get('Options','metric_file')
    hname = options.get('Options','html_name')
    html = options.getboolean('Options','html_go')
    global debug
    debug = options.getboolean('Options','debug')

    if gogo: 
        os.system('cp ~/snakes/plot_FRD.sm .')
        soba(N.ratios,num_ap,dir_cut,N.exten,T,mfile)
#        soba(N,num_ap,dir_cut,0,T,mfile)
        if html: webit(hname)

    return
    
class Noodle:
    ''' Noodle is a data structure that keeps track of all the different
    data files associated with the full-beam bench. It also can reduce
    these files to produce images ready for analysis by FReD.
    '''


    def __init__(self, config):
        '''The data is intiallized with a .ini file that is passed
        to Noodle throuh main(). All initialization does is read in the
        relevant parameters and create the empty data structures.
        '''
        self.fiber_dir = config.get('Data','fiber_directory')
        self.direct_dir = config.get('Data','direct_directory')
        self.dark_dir = config.get('Data','dark_directory')
        self.exten = config.getint('Data','fits_exten')
        self.clean = config.getboolean('Options','cleanup')
        self.darkcombine = config.get('Options','darkcombine')
        self.datacombine = config.get('Options','datacombine')
        self.datareject = config.get('Options','datareject')
        self.rejectpar = config.get('Options','rejectpar')
        
        self.darks = {}
        self.fiber = {}
        self.direct = {}
        self.ratios = {}

    def build_run(self):
        '''build_run is the highest level method for Noodle and really the
        only one you would need to call directly. It assembles the raw data
        and reduces it. It also calls a lot of helper functions.
        '''
        print '\nBuilding data run...'
        self.get_darks()
        print ' dark-subtracting fiber images...'
        self.reduce_data(self.fiber_dir,self.fiber)
        print ' dark-subtracting direct images...'
        self.reduce_data(self.direct_dir,self.direct)
        self.clean_up(self.ratios)
        print ' combining images...'
        self.combine()
        if self.clean: os.system('rm *_ds.fits')

    def reduce_data(self,directory,dictionary):
        
        self.fill_dict(directory,dictionary)
        self.sub_darks(dictionary)
        self.add_to_ratios(dictionary)

    def fill_dict(self,directory,dictionary):
        data_list = glob.glob(directory+'/*FIT')
        
        for data in data_list:
            exptime = pyfits.open(data)[self.exten].header['EXPTIME']
            if exptime not in dictionary.keys():
                dictionary[exptime] = {'raw':[],'ds':[]}
                
            dictionary[exptime]['raw'].append(data)
            
    def get_direct(self):
        direct_list = glob.glob(self.direct_dir+'/*')

        for data in direct_list:
            
            head = pyfits.open(data)[self.exten].header
            focal_length = head['FOCALLEN']
            diameter = float(head['APTDIA'])
            focal_ratio = focal_length/diameter
            filt = head['FILTER']
            exptime = head['EXPTIME']
            
            if focal_ratio not in self.direct.keys():
                self.direct[focal_ratio] = {}
            if filt not in self.direct[focal_ratio].keys():
                self.direct[focal_ratio][filt] = {'raw':{},'combined':None}
            if data.find('Combined') >= 0:
                self.direct[focal_ratio][filt]['combined'] = data
            elif data.find('Combined') == -1:
                if exptime not in self.direct[focal_ratio][filt]['raw'].keys():
                    self.direct[focal_ratio][filt]['raw'][exptime] =\
                        {'raw':[],'ds':[]}
                self.direct[focal_ratio][filt]['raw'][exptime]['raw']\
                    .append(data)

    def gen_direct(f_ratio,filt):
        iraf.imcombine.combine = self.datacombine
        iraf.imcombine.reject = self.datareject
        iraf.imcombine.lsigma = self.rejectpar

        self.sub_darks(self.direct[f_ratio][filt]['raw'])
        name = self.direct_dir+\
            '/Combined_'+filt+'F'+str(int(focal_ratio*10))+'d.fits'
        if len(self.direct[f_ratio][filt]['raw'].keys()) > 1:
            print 'WARNING, DIRECT IMAGES MAY HAVE MULTIPLE EXPS PER FILT'
        exp = self.direct[f_ratio][filt]['raw'].keys()[0]
        iraf.imcombine(','.join(self.direct[f_ratio][filt]['raw'][exp]['ds']),\
                           name)
        self.direct[f_ratio][filt]['combined'] = name

    def get_darks(self):
        print ' finding darks...'
        dark_list = glob.glob(self.dark_dir+'/*')
    
        for dark in dark_list:
            exptime = pyfits.open(dark)[self.exten].header['EXPTIME']
            if exptime not in self.darks.keys():
                self.darks[exptime] = {'raw':[],'combined':None}

            if dark.find('Combined') >= 0:
                self.darks[exptime]['combined'] = dark
            else:
                self.darks[exptime]['raw'].append(dark)

    def gen_dark(self,exptime):
            iraf.imcombine.reject = 'none'
            iraf.imcombine.combine = self.darkcombine

            name = self.dark_dir+'/Combined_'+str(exptime)+'_DARK.fits'
            iraf.imcombine(','.join(self.darks[exptime]['raw']),name)
            self.darks[exptime]['combined'] = name

    def sub_darks(self,data):
        for exp in data.keys():
            inputlist = data[exp]['raw']
            outputlist = [s[s.rfind('/')+1:s.find('.FIT')]+\
                              '_ds.fits'\
                              for s in inputlist]
            
            if self.darks[exp]['combined'] == None: self.gen_dark(exp)
            if len(inputlist) > 10:
                while len(inputlist) > 0:
                    inputstring = ','.join(inputlist[:10])
                    outputstring = ','.join(outputlist[:10])
                    iraf.imarith(inputstring,'-',\
                                     self.darks[exp]['combined'],\
                                     outputstring)
                    data[exp]['ds'] += outputlist[:10]
                    inputlist = inputlist[10:]
                    outputlist = outputlist[10:]
                
            else:
                iraf.imarith(','.join(inputlist),\
                                 '-',\
                                 self.darks[exp]['combined'],\
                                 ','.join([s[s.rfind('/')+1:s.find('.FIT')]+\
                                               '_ds.fits' for s in inputlist]))
                data[exp]['ds'] += outputlist
    
    def add_to_ratios(self, data):
        for exp in data.keys():
            for ds in data[exp]['ds']:
                head = pyfits.open(ds)[self.exten].header
                focal_length = head['FOCALLEN']
                diameter = float(head['APTDIA'])
                focal_ratio = round(focal_length/diameter,1)
                filt = head['FILTER']
                ftype = head['OBSERVER']
            
                if focal_ratio not in self.ratios.keys():
                    self.ratios[focal_ratio] =\
                        {'data':{},'focal_length':focal_length}
                if filt not in self.ratios[focal_ratio]['data'].keys():
                    self.ratios[focal_ratio]['data'][filt] = {}
                if ftype not in self.ratios[focal_ratio]['data'][filt].keys():
                    self.ratios[focal_ratio]['data'][filt][ftype] =\
                        {'raw':[],'final':None}
                
                self.ratios[focal_ratio]['data'][filt][ftype]['raw'].append(ds)

    def clean_up(self,data):
        for fratio in data.keys():
            for filt in data[fratio]['data'].keys():
                if 'fiber' not in data[fratio]['data'][filt].keys():
                    del data[fratio]['data'][filt]

    def combine(self):
        iraf.imcombine.combine = self.datacombine
        iraf.imcombine.reject = self.datareject
        iraf.imcombine.lsigma = self.rejectpar

        for f_ratio in self.ratios.keys():
            for filt in self.ratios[f_ratio]['data'].keys():
                for ftype in self.ratios[f_ratio]['data'][filt].keys():
                    name = self.ratios[f_ratio]['data'][filt][ftype]['raw'][0]
                    name = name[:name.rfind('.0')]+'.fits'
                    iraf.imcombine(\
                        ','.join(self.ratios\
                                     [f_ratio]['data'][filt][ftype]['raw']),\
                            name)
                    self.ratios[f_ratio]['data'][filt][ftype]['final'] = name

class thePot:

    def __init__(self, config):

        self.functions = {}
        self.ref_levels = {}

        t_file = config.get('Data','Tput_file')
        order = config.getint('Data','Tput_fit_order')

        self.times, self.levels = np.loadtxt(t_file,usecols=(0,1),unpack=True,
                                   converters={0: self.format_time})
        self.filters = np.loadtxt(t_file,usecols=(2,),dtype=np.str)

        for filt in np.unique(self.filters):
            fidx = np.where(self.filters == filt)
            firstidx = np.where(self.times[fidx] == self.times[fidx].min())
            self.functions[filt] = \
                np.polyfit(self.times[fidx],self.levels[fidx],order)
            self.ref_levels[filt] = \
                self.make_func(self.times[fidx],self.functions[filt])

    def get_correction(self, time, filt):

        try: level = self.make_func(time,self.functions[filt])
        except KeyError:
            print 'Oops! For some reason I am trying to get a correction'+\
                'for a filter that isn\'t in my database.'
            return
        return (level/self.ref_levels[filt])[0]
    
    def format_time(self,string):
        
        h = np.float(string[0:2])*3600.
        m = np.float(string[3:5])*60.
        s = np.float(string[6:])
        
        return h+m+s

    def make_func(self,x,p):
        
        out = 0.
        for i in range(p.size):
            out += p[i]*x**(p.size-1-i)
            
        return out
    
    def test(self,filt):
        fig = plt.figure(0)
        plt.clf()
        ax1 = fig.add_subplot(211)
        
        idx = np.where(self.filters == filt)
        t = np.linspace(self.times[idx].min(),self.times[idx].max())
        l = self.make_func(t,self.functions[filt])
        
        ax1.plot(self.times[idx],self.levels[idx],'.',t,l,'-')
        
        ax2 = fig.add_subplot(212)
        residuals = self.levels[idx]\
                     - self.make_func(self.times[idx],self.functions[filt])
        ax2.plot(self.times[idx],residuals,'.')

        ax2.set_title('Residuals\n$\sigma$='+str(np.std(residuals)))

        fig.show()
        return

def ABABA(filter_list,prestr):
    Nood = {6.3: {'focal_length': 50, 'data':{} } }

    Avec = np.array([[0,0,0]])

    for filt in filter_list:
        
        A_data = np.empty((3,1024,1024))
        A_times = np.empty((3,2))
        for i in [1,2,3]:
            name = glob.glob('*_'+filt+'*d_A'+str(i)+'.fits')[0]
            A_HDU = pyfits.open(name)
            A_timestr = A_HDU[0].header['TIME-OBS']
            A_time = float(A_timestr[:A_timestr.find(':')])*3600 +\
                float(A_timestr[A_timestr.find(':')+1:A_timestr.rfind(':')])*\
                60+float(A_timestr[A_timestr.rfind(':')+1:])
            A_time = A_time + (A_HDU[0].header['EXPTIME'] + 5.5)*\
                A_HDU[0].header['NCOMBINE']/2
            idx = np.where(A_HDU[0].data > 5000)
            A_data[i-1] = A_HDU[0].data
            A_times[i-1] = np.array([A_time, np.mean(A_HDU[0].data[idx])])
            Avec = np.append(Avec,[[A_time, np.mean(A_HDU[0].data[idx]),\
                                        np.std(A_HDU[0].data[idx])]],axis=0)

            mask = np.zeros(A_data[i-1].shape)
            mask[idx] = 1.0
            pyfits.PrimaryHDU(mask).writeto('mask_'+name)
            
        B_data = np.empty((2,1024,1024))
        B_times = np.empty((2))
        for j in [1,2]:
            name = glob.glob('*_'+filt+'*f_B'+str(j)+'.fits')[0]
            B_HDU = pyfits.open(name)
            B_timestr = B_HDU[0].header['TIME-OBS']
            B_time = float(B_timestr[:B_timestr.find(':')])*3600 +\
                float(B_timestr[B_timestr.find(':')+1:B_timestr.rfind(':')])*\
                60+float(B_timestr[B_timestr.rfind(':')+1:])
            B_time = B_time + (B_HDU[0].header['EXPTIME'] + 5.5)*\
                B_HDU[0].header['NCOMBINE']/2
            B_data[j-1] = B_HDU[0].data
            B_times[j-1] = B_time

        B1_corr = np.interp(B_times[0],np.array([A_times[0,0],A_times[1,0]]),\
                                np.array([A_times[0,1],A_times[1,1]]))/\
                                A_times[0,1]
        B_data[0] /= B1_corr
        
        B2_corr = np.interp(B_times[1],np.array([A_times[1,0],A_times[2,0]]),\
                                np.array([A_times[2,1],A_times[2,1]]))/\
                                A_times[0,1]
        B_data[1] /= B2_corr

        
        A_data[1] /= A_times[1,1] / A_times[0,1]
        A_data[2] /= A_times[2,1] / A_times[0,1]

        fits = pyfits.PrimaryHDU(np.mean(B_data,axis=0))
        fits.header.update('EXPTIME',B_HDU[0].header['EXPTIME'])
        fits.header.update('FILTER',filt)
        dits = pyfits.PrimaryHDU(np.mean(A_data,axis=0))
        dits.header.update('EXPTIME',A_HDU[0].header['EXPTIME'])
        dits.header.update('FILTER',filt)

        fits.writeto(prestr+'ABABA_'+filt+'f.fits')
        dits.writeto(prestr+'ABABA_'+filt+'d.fits')

        Nood[6.3]['data'][filt] = {\
            'direct': {'final': prestr+'ABABA_'+filt+'d.fits'},\
                'fiber': {'final': prestr+'ABABA_'+filt+'f.fits'}}

    pickle.dump(Nood,open('ABABA.pkl','wb'))
    
    Avec = Avec[1:]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('ADU')
    ax.errorbar(Avec[:3,0],Avec[:3,1],yerr=Avec[:3,2],label=filter_list[0])
    ax.errorbar(Avec[3:6,0],Avec[3:6,1],yerr=Avec[3:6,2],label=filter_list[1])
    ax.errorbar(Avec[6:,0],Avec[6:,1],yerr=Avec[6:,2],label=filter_list[2])
    ax.legend(loc=0)
    ax.set_ylim(0,35000)
    ax.set_xlim(Avec[0,0]-100,Avec[-1,0]+100)
    ax.set_title(prestr+'\n'+datetime.now().isoformat(' '))
    plt.savefig('ABABA_correction.pdf')

    return Avec

def webit(name):
    
    path = '/d/www/eigenbrot/polish/'+name
    g = glob.glob('*_face_*.fits')
    idd = g[0][g[0].rfind('_')+1:g[0].rfind('.')]
    os.system('mkdir '+path)
    os.system('cp *.dat '+path)
    os.system('cp *.ps '+path)
    os.system('cp *.jpg '+path)
    os.system('cp *_face_'+idd+'.fits '+path)
    os.system('convert '+path+'/Output_face_'+idd+'.fits '
              +path+'/Output_face_'+name+'.jpg')
    os.system('convert '+path+'/Input_face_'+idd+'.fits '
              +path+'/Input_face_'+name+'.jpg')
    os.system('rename _'+idd+' _'+name+' '+path+'/*face_*.fits')
    ht = open(path+'/index.html','wb')
    ht.write('<HTML><TITLE>FRD Bench data output '+name+'</TITLE>\n'
             +'<body bgcolor="EEEEEE" link="#333FF" vlink="#5500BB">\n'
             +'<h1>Data products for '+name+'\n'
             +'<h2>Web summary generated on '+datetime.now().isoformat(' ')
             +'</h2></h1>\n<hr noshade>\n'
             +'<h3>Direct images</h3>\n'
             +'<table cellpadding=10><tr><td>\n'
             +'<a href=Input_face_'+name+'.fits><img src=Input_face_'
             +name+'.jpg height=150 width=200></a><br>\n'
             +'Input face (click for FITS)\n</td><td>\n'
             +'<a href=Output_face_'+name+'.fits><img src=Output_face_'
             +name+'.jpg height=150 width=200></a><br>\n'
             +'Output face (click for FITS)\n</td></tr></table><br>\n\n'
             +'<h3>Data files</h3>\n')

    dat_list = glob.glob(path+'/*.dat')
    for dat in dat_list:
        link = dat[dat.rfind('/')+1:]
        ht.write('<a href='+link+' type="text/plain">'+link+'</a><br>\n')

    ht.write('<h3>FRD plots</h3>Click for postscript<br><br>\n')
    plt_list = [i[i.rfind('/')+1:i.rfind('.')] for i in glob.glob(path+'/*.ps')]
    for plot in plt_list:
        ht.write(plot+'<br>\n')
        ht.write('<a href='+plot+'.ps>ps</a> <a href='+plot+'.jpg>jpg</a><br>\n')
        ht.write('<img src='+plot+'.jpg height=300 width=232><br><br>\n')
    ht.write('</body></HTML>')
    ht.close()
    
    old_m = open('/d/www/eigenbrot/polish/index.html')
    lines = old_m.readlines()
    old_m.close()

    lines[3] = '<h2>Last updated '+datetime.now().isoformat(' ')+'</h2></h1>\n'
    m = open('/d/www/eigenbrot/polish/index.html','wb')
    m.writelines([l for l in lines[:-1]])
    m.write('<a href='+name+'>'+name+'</a><br>\n</body></HTML>')
    m.close()


def annulize(data, num_an, distances=np.array([0])):
    '''
    Description:
        annulize takes in a FITS image and computes the total power 
        contained within progressivly larger annuli. The number of 
        annuli used is set by the user and, unlike annulize_sb, it 
        is the width of the annuli that remains constant, not the area.

    Inputs:
        data      - ndarray
                    The data to be annulized
        num_an    - Int
                    The number of annuli to use
        distances - ndarray
                    This is expected to be a transformed distance array

    Output:
        r_vec - ndarray
                vector of radii where each entry is the radius of the
                middle of the corresponding annulus.
        fluxes- ndarray
                each entry is the flux contained within that annulus.
        errors- ndarray
                The standard deviation of the pixels in each annulus.
    '''

    if debug: print'annulizing...'
    if debug: print"data type is "+str(data.dtype)

    dims = data.shape

    '''check to see if we got some transformed distances and generate
    a basic distance array if we didn't'''
    if not(distances.any()):
        center = centroid(data)
        if debug: print "center at "+str(center)
        distances = dist_gen(data, center)


    rlimit = distances.max()
    rstep = rlimit/num_an

    'initialize future output arrays. I think that maybe with python'
    'these can be left empty, but w/e'
    fluxes = np.zeros(num_an, dtype=np.float32)
    r_vec = np.zeros(num_an, dtype=np.float32)
    outarray = np.zeros(dims, dtype=np.float32)
    area_array = np.zeros(dims, dtype=np.float32) + 1
    errors = np.zeros(num_an, dtype=np.float32)
    r1 = 0.0
    r2 = rstep
    
    for i in range(num_an):
        idx = np.where((distances <= r2) & (distances > r1))
        
        '''The correction is designed to correct annuli that are not
        entirely contained in the data array, but I still can't decide
        if it screws up the data or not'''
        correction = math.pi*(r2**2 - r1**2)/np.sum(area_array[idx])
        fluxes[i] = np.sum(data[idx])#*correction
        errors[i] = np.std(data[idx])

        if debug == 3: print(i,r1,r2,fluxes[i],correction) 
        
        '''this is only used during debug to show where the annuli
        lie on the image'''
        outarray[idx] = np.mod(i,2) + 1

        r_mid = (r1 + r2)*0.5
        r_vec[i] = r_mid
        r1 = r2
        r2 = r1 + rstep

    if debug == 2:
        pyfits.PrimaryHDU(data+outarray*data.max()/2)\
            .writeto('an_show.fits',clobber=True)
                    
    return(r_vec,fluxes,errors)

if not SLOWAN:
    annulize = ADE.fast_annulize

if __name__ == '__main__':
    sys.exit(main())
