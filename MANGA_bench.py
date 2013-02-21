#! /usr/bin/env python

print 'Loading module...'
import numpy as np
import os
import matplotlib
if os.popen('echo $DISPLAY').readline() == 'localhost:10.0\n': 
    print 'Deactivating display...'
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import AxesGrid as AG
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib.collections import PatchCollection
print 'syncing mesh...'
from pyraf import iraf
import pyfits
print 'calculating splines...'
import ADEUtils as ADE
import time
print 'initializing goodness...'
import sys
import glob
from datetime import datetime
import ConfigParser
import pickle
import MANGAmap as mmp
print 'load complete!'

debug = False

def FReD(direct_image, fiber_image, num_ap, pot, filt, dir_cut,\
             EXTEN=0, OUTPUT=0, FL=50, FR=4.2, FP='99,99'):
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
        # direct_time_str = direct_HDU.header['TIME-OBS']
        # fiber_time_str = fiber_HDU.header['TIME-OBS']

        # direct_time = np.float(direct_time_str[6:])\
        #     + np.float(direct_time_str[3:5])*60.\
        #     + np.float(direct_time_str[0:2])*3600.

        # fiber_time = np.float(fiber_time_str[6:])\
        #     + np.float(fiber_time_str[3:5])*60.\
        #     + np.float(fiber_time_str[0:2])*3600.
        
        # fcorrect = pot.get_correction(fiber_time,filt)
        # dcorrect = pot.get_correction(direct_time,filt)
        direct_start_time = direct_HDU.header['STARTIME']
        direct_end_time = direct_HDU.header['ENDTIME']
        fiber_start_time = fiber_HDU.header['STARTIME']
        fiber_end_time = fiber_HDU.header['ENDTIME']
        
#        filt = 'V'

        pot.set_ref(filt,direct_start_time,direct_end_time)

        dcorrect = pot.get_correction2(direct_start_time,direct_end_time,filt)
        fcorrect = pot.get_correction2(fiber_start_time,fiber_end_time,filt)
        print '   Direct throughput correction is '+str(dcorrect)
        print '   Fiber throughput correction is '+str(fcorrect)
        
        # direct *= dcorrect
        # fiber *= fcorrect
        
        if debug: 
            pot.plot(filt,[direct_start_time,direct_end_time,fiber_start_time,fiber_end_time])
            raw_input('    Tput plots...')


    if debug: print '    Annulizing...'
    d_rvec, d_sb, d_sberr = ADE.fast_annulize(direct,num_ap)
    f_rvec, f_sb, f_sberr = ADE.fast_annulize(fiber,num_ap)

    if debug:
        plt.clf()
        fig = plt.figure(0)
        ax1 = fig.add_subplot(221)
        ax1.plot(d_rvec,d_sb)
        ax1.plot(f_rvec,f_sb)
        ax1.set_xlabel('Radius [px]')
        ax1.set_ylabel('Counts [ADU]')
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
        ax2.set_xlabel('Radius [mm]')
        ax2.set_ylabel('Cumulative flux [ADU]')
        fig.show()
        raw_input('     Cumsum')
        

    '''Now we normalize the fluxes so we are talking about EE, the enclosed
    energy'''
    if debug: print '    EEing...'
    d_max = np.max(d_flux)
    f_max = np.max(f_flux)
    d_EE = d_flux/d_max
    f_EE = f_flux/f_max
    '''we'll use f_EE_nd for plotting later'''
    f_EE_nd = f_flux/d_max
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
    
#    if debug: print (np.abs(f_r_c - f_rvec))**0.5

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
    
    sloanf = np.interp(rf4,f_rvec,f_ADUrate)
    sloand = np.interp(rf4,d_rvec,d_ADUrate)
    sloan_m = sloanf/sloand
    
    r_ideal = FL/(2*FR)
    r_ideal_test = d_r_c[np.where(f_ADUrate == f_ADUrate.max())[0]]
    metric = np.interp(r_ideal, f_r_c, f_EE)
    metric80 = FL/(2*np.interp(0.8,f_EE,f_r_c))
    metric90 = FL/(2*np.interp(0.9,f_EE,f_r_c))
    
    if debug:
        try:
            direct_V = pot.get_voltage(direct_start_time,direct_end_time,filt)
            fiber_V = pot.get_voltage(fiber_start_time,fiber_end_time,filt)
        except KeyError:
            direct_V,drawV = (1.0,1.0)
            fiber_V,frawV = (1.0,1.0)
        ax3 = fig.add_subplot(223)
        ax3.plot(f_rvec,f_ADUrate/fiber_V,'g')
        ax3.plot(d_rvec,d_ADUrate/direct_V,'b')
        ax3.set_xlabel('Radius [mm]')
        ax3.set_ylabel('Count rate/voltage [ADU/(V s)]')
        ax3.axvline(ls=':',x=rf4,color='b')
        ax3.axvline(ls='--',x=rf4,color='g')
        ax3.axhline(ls='--',y=sloanf/fiber_V,color='g')
        ax3.axhline(ls=':',y=sloand/direct_V,color='b')
        infostr = "\tFiber ratio at f/4 ({:}mm): {:8.4E}".format(rf4,sloanf/fiber_V)+\
            "\n\tDirect ratio at f/4 ({}mm): {:8.4E}".format(rf4,sloand/direct_V)+\
            "\n\tFiber counts/s at f/4: {:8.4E}".format(sloanf)+\
            "\n\tDirect counts/s at f/4: {:8.4E}".format(sloand)+\
            "\n\tFiber voltages are : {:4.4f}".format(fiber_V)+\
            "\n\tDirect voltages are : {:4.4f}".format(direct_V)+\
            "\n\tSloan metric is {}".format((sloanf*direct_V)/(sloand*fiber_V))+\
            "\n\tFull direct is:  {:8.4E}".format(float(d_ADUrate.max()))+\
            "\n\tMax d_r_c value is: {:4.2f}".format(d_r_c.max())+\
            "\n\tr_ideal is: {:4.2f}".format(r_ideal)+\
            "\n\tr_max is: {}".format(r_ideal_test)
        ax3.text(1.1,0.4,infostr,transform=ax3.transAxes,ha='left',va='center')
        print infostr
        fig.suptitle('{} - {}'.format(os.popen('pwd').readlines()[0],datetime.now().isoformat(' ')))
        fig.show()
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
                +'# Fiber position: '+FP+'\n'
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
        
    return ((metric90,metric80,tput_100,sloan_m,
            tput_f5,tput_f4,tput_f5b,tput_f4b),(d_N,d_N_c,f_N,f_N_c,d_EE,f_EE,f_EE_nd))

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
            +'# {:10}= '.format('Fiber_pos')+'fiber input position\n'
            +'# {:10}= '.format('Out_pos')+'fiber output position\n'
            +'# {:10}= '.format('filt')+'filter\n'
            +'# {:10}= '.format('N90')+'fiber f/# at EE90\n'
            +'# {:10}= '.format('N80')+'fiber f/# at EE80\n'
            +'# {:10}= '.format('tput')+'total throughput\n'
            +'# {:10}= '.format('sloan')+'fiber within f/4 / direct within f/4 (uncorrected)\n'
            +'# {:10}= '.format('tput5')+'throughput at f/5\n'
            +'# {:10}= '.format('tput4')+'throughput at f/4\n'
            +'# {:10}= '.format('EE5')+'fiber EE at f/5\n'
            +'# {:10}= '.format('EE4')+'fiber EE at f/4\n'
            +'#\n'
            +str('#{0:>10}{1:>10}{2:>9}{3:>9}{4:>9}{5:>9}{6:>9}{7:>9}'
                 +'{8:>9}{9:>9}{10:>9}\n')\
                .format('Fiber_pos','Out_pos','filt','N90','N80','tput','sloan',\
                            'tput5','tput4','EE5','EE4')
            +str('#{0:>10}{1:>10}{2:>9}{3:>9}{4:>9}{5:>9}{6:>9}{7:>9}'
                 +'{8:>9}{9:>9}{10:>9}\n')\
                .format(1,2,3,4,5,6,7,8,9,10,11))

    hdulist = []

    for fiber_pos in nood.keys():

        print 'Fiber ('+str(fiber_pos)+')'
        try: outpos = nood[fiber_pos]['outpos']
        except KeyError: outpos = 'NA'
        print 'Output position is ('+outpos+')'
        focal_length = nood[fiber_pos]['focal_length']
        focal_ratio = nood[fiber_pos]['focal_ratio']
        print 'Camera focal length = '+str(focal_length)+' mm'
        print 'Fiber fed at f/{}'.format(int(focal_ratio))
    
        for filt in nood[fiber_pos]['data'].keys():
            direct_name = nood[fiber_pos]['data'][filt]['direct']['final']
            fiber_name = nood[fiber_pos]['data'][filt]['fiber']['final']
            name = fiber_name[:fiber_name.rfind('f.')]

            print ' Filter: '+str(filt)
            print '  Direct image is '+direct_name
            print '  Fiber image is  '+fiber_name

            (metric,plot_data) = FReD(direct_name,fiber_name,num_ap,pot,filt,\
                                          dir_cut,EXTEN=exten,FL=focal_length,\
                                          FR=focal_ratio,OUTPUT=name+'.dat',FP=fiber_pos)
            
            hdu = pyfits.ImageHDU(np.array(plot_data))
            hdu.header.update('FIBERPOS',fiber_pos)
            hdu.header.update('OUTPOS',outpos)
            hdu.header.update('SLOAN',metric[3])
            hdulist.append(hdu)
            f.write('{0:>11}{1:>10}{2:>9}'.format(fiber_pos,outpos,filt))
            for m in metric: f.write('{:9.4f}'.format(m))
            f.write('\n')

            '''now generate the FRD plots using supermongo'''
            sm = open('tmp_'+name+'.sm','wb')
            sm.write('verbose 0\n'
                     +'macro read manga_FRD.sm\n'
                     +'device postencap_color '+name+'.ps\n'
                     +'manga_FRD "'+name+'.dat"\ndevice nodevice\n')
            sm.close()
            os.system('sm < tmp_'+name+'.sm')
            os.system('convert -density 200 '+name+'.ps -quality 92 '+name+'.jpg')
            os.system('rm tmp_'+name+'.sm')

    try:
        pyfits.HDUList([pyfits.PrimaryHDU(None)]+hdulist).writeto('plotdata.fits')
    except IOError:
        if raw_input("Overwrite file 'plotdata.fits'?: (Y/n)").lower() in ['y','']:
            os.system('rm plotdata.fits')
            pyfits.HDUList([pyfits.PrimaryHDU(None)]+hdulist).writeto('plotdata.fits')
            # plot_helper('plotdata.fits','allfibers.pdf','')
            # os.system('convert -density 200 allfibers.pdf -quality 92 allfibers.jpg')
    f.close()
    return

def main():
    
    options = ConfigParser.ConfigParser()
    options.read(sys.argv[1])
    print "Reading settings from "+sys.argv[1]
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

    num_ap = options.getint('Options','num_ap')
    dir_cut = options.getint('Options','direct_cutoff')
    gogo = options.getboolean('Options','gogo')
    mfile = options.get('Options','metric_file')
    hname = options.get('Options','html_name')
    html = options.getboolean('Options','html_go')
    global debug
    debug = options.getboolean('Options','debug')

    if gogo: 
        doT = options.get('Data','Tput_file')
        if doT.lower() == 'false': T = False
        else: T = thePot(doT)
        os.system('cp /d/monk/eigenbrot/MANGA/manga_FRD.sm .')
        soba(N.ratios,num_ap,dir_cut,N.exten,T,mfile)
        ring_helper(mfile)
        plot_helper('plotdata.fits',
                    'allfibers.pdf',
                    hname+'\n'+datetime.now().isoformat(' '))
        os.system('convert -density 200 allfibers.pdf -quality 92 allfibers.jpg')

#        soba(N,num_ap,dir_cut,0,T,mfile) #use if running an ABABA run

        if html: 
            print 'Producing web product...'
            webit(hname)
        print 'Reduction completed'

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
        self.mapping_scheme = config.get('Options','inoutmap')
        
        self.darks = {}
        self.fiber = {}
        self.direct = {}
        '''I'll keep calling the top-level distinction a ratio because it will
        avoid any confusion that might arise by calling it fibers'''
        self.ratios = {}
        self.mapdb = self.genmap(self.mapping_scheme)

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
        self.fill_dict(self.direct_dir,self.direct)
        self.sub_darks(self.direct)
        self.direct_to_ratios(self.direct)
        self.clean_up(self.ratios)
        print ' combining images...'
        self.combine()
        self.uphead(self.mapdb)
        if self.clean: os.system('rm *_ds.fits')

    def genmap(self,scheme):
        '''All this does is return the correct input-output mapping. It doesn't
        actually compute anything'''
        
        try:
            return mmp.master[scheme]
        except KeyError:
            print "Warning: mapping scheme not in database"
            return False

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
            fiber_pos = head['ORIGIN']
            filt = head['FILTER']
            exptime = head['EXPTIME']
            
            if fiber_pos not in self.direct.keys():
                self.direct[fiber_pos] = {}
            if filt not in self.direct[fiber_pos].keys():
                self.direct[fiber_pos][filt] = {'raw':{},'combined':None}
            if data.find('Combined') >= 0:
                self.direct[fiber_pos][filt]['combined'] = data
            elif data.find('Combined') == -1:
                if exptime not in self.direct[fiber_pos][filt]['raw'].keys():
                    self.direct[fiber_pos][filt]['raw'][exptime] =\
                        {'raw':[],'ds':[]}
                self.direct[fiber_pos][filt]['raw'][exptime]['raw']\
                    .append(data)

    def get_darks(self):
        print ' finding darks...'
        dark_list = glob.glob(self.dark_dir+'/*.[Ff][Ii][Tt]*')
    
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
                    iraf.nhedit(inputstring,'FIBERPOS','"(ORIGIN)"','Fiber input position',add=True,addonly=True)
                    iraf.imarith(inputstring,'-',\
                                     self.darks[exp]['combined'],\
                                     outputstring)
                    data[exp]['ds'] += outputlist[:10]
                    inputlist = inputlist[10:]
                    outputlist = outputlist[10:]
                
            else:
                iraf.nhedit(','.join(inputlist),'FIBERPOS','"(ORIGIN)"','Fiber input position',add=True,addonly=True)
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
                fiber_pos = head['FIBERPOS']
                L2_focal_length = head['APTAREA']
                diameter = float(head['APTDIA'])
                focal_ratio = round(L2_focal_length/diameter,1)
                filt = head['FILTER']
                ftype = head['OBSERVER']
                timestr = head['TIME-OBS']
                obstime = np.float(timestr[6:])\
                    + np.float(timestr[3:5])*60.\
                    + np.float(timestr[0:2])*3600.
                
                if fiber_pos not in self.ratios.keys():
                    self.ratios[fiber_pos] =\
                        {'data':{},'focal_ratio':focal_ratio}
                if filt not in self.ratios[fiber_pos]['data'].keys():
                    self.ratios[fiber_pos]['data'][filt] = {}
                if ftype not in self.ratios[fiber_pos]['data'][filt].keys():
                    self.ratios[fiber_pos]['data'][filt][ftype] =\
                        {'raw':[],'obstimes':[],'exptime':exp,'final':None}
                
                self.ratios[fiber_pos]['data'][filt][ftype]['raw'].append(ds)
                self.ratios[fiber_pos]['data'][filt][ftype]['obstimes'].append(obstime)

    def direct_to_ratios(self, data):
        for exp in data.keys():
            for ds in data[exp]['ds']:
                head = pyfits.open(ds)[self.exten].header
                L3_focal_length = head['FOCALLEN']
                filt = head['FILTER']
                ftype = head['OBSERVER']
                timestr = head['TIME-OBS']
                obstime = np.float(timestr[6:])\
                    + np.float(timestr[3:5])*60.\
                    + np.float(timestr[0:2])*3600.

                for pos in self.ratios.keys():
                    self.ratios[pos]['focal_length'] = L3_focal_length
                    try: 
                        self.ratios[pos]['data'][filt][ftype]['raw'].append(ds)
                        self.ratios[pos]['data'][filt][ftype]['obstimes'].append(obstime)
                    except KeyError:
                        self.ratios[pos]['data'][filt][ftype] =\
                            {'raw':[ds],'obstimes':[obstime],'exptime':exp,'final':None}

    def clean_up(self,data):
        final_direct = []
        for fratio in data.keys():
            for filt in data[fratio]['data'].keys():
                if 'fiber' not in data[fratio]['data'][filt].keys():
                    del data[fratio]['data'][filt]

    def combine(self):
        iraf.imcombine.combine = self.datacombine
        iraf.imcombine.reject = self.datareject
        iraf.imcombine.lsigma = self.rejectpar

        for fiber_pos in self.ratios.keys():
            for filt in self.ratios[fiber_pos]['data'].keys():
                for ftype in self.ratios[fiber_pos]['data'][filt].keys():
                    name = self.ratios[fiber_pos]['data'][filt][ftype]['raw'][0]
                    name = name[:name.rfind('.0')]+'.fits'
                    mintime = min(self.ratios[fiber_pos]['data'][filt][ftype]['obstimes'])
                    maxtime = max(self.ratios[fiber_pos]['data'][filt][ftype]['obstimes'])
                    maxtime += float(self.ratios[fiber_pos]['data'][filt][ftype]['exptime'])

                    iraf.imcombine(\
                        ','.join(self.ratios\
                                     [fiber_pos]['data'][filt][ftype]['raw']),\
                            name)
                    self.ratios[fiber_pos]['data'][filt][ftype]['final'] = name
                    iraf.nhedit(name,'STARTIME',mintime,'Start time of first combined image',addonly=True)
                    iraf.nhedit(name,'ENDTIME',maxtime,'End time of last combined image',addonly=True)

    def uphead(self,mapping):
        
        if mapping == False: 
            print "Warning: Input to output mapping not loaded correctly.\n"+\
                "FITS headers will NOT be updated"
            return
        
        for fiber_pos in self.ratios.keys():
            for filt in self.ratios[fiber_pos]['data'].keys():
                outpos = mapping[fiber_pos]
                iraf.nhedit(self.ratios[fiber_pos]['data'][filt]['fiber']['final'],\
                                "OUTPOS",outpos,'Fiber output position',\
                                add=True,addonly=False,after='FIBERPOS')
                self.ratios[fiber_pos]['outpos'] = outpos

class thePot:

    def __init__(self, t_file):

        self.ref_levels = {}

##        t_file = config.get('Data','Tput_file')
 ##       self.order = config.getint('Data','Tput_fit_order')

        try: self.times, self.levels, self.errs = np.loadtxt(t_file,usecols=(0,1,2),
                                                             unpack=True,
                                                             converters={0: self.format_time})
        except ValueError:
            self.times, self.levels = np.loadtxt(t_file,usecols=(0,1),unpack=True,
                                                 converters={0: self.format_time})

        self.filters = np.loadtxt(t_file,usecols=(-1,),dtype=np.str)


    def set_ref(self, filt, time1, time2):

        if time2 > self.times.max() or time1 < self.times.min():
            self.ref_levels[filt] = {'level': 1.0, 'times': [time1,time2]}
            return

        tidx = np.where((self.times >= time1) & (self.times <= time2))
        avg_level = np.mean(self.levels[tidx])
        self.ref_levels[filt] = {'level': avg_level, 'times': [time1,time2]}
#        self.plot(filt,times=[time1,time2])

    
    def get_correction2(self, time1, time2, filt):

        if time2 > self.times.max() or time1 < self.times.min():
            print "Asking for time outside of range"
            return 1.0
        
        tidx = np.where((self.times >= time1) & (self.times <= time2))
        avg_level = np.mean(self.levels[tidx])
        return self.ref_levels[filt]['level']/avg_level

    def get_voltage(self,time1,time2,filt,std=False):
        
        idx = np.where((self.times >= time1) & (self.times <= time2))
        rawlevel = np.mean(self.levels[idx])

        if std:
            return rawlevel, np.std(self.levels[idx])
        
        else: return rawlevel
                
    def format_time(self,string):
        
        h = np.float(string[0:2])*3600.
        m = np.float(string[3:5])*60.
        s = np.float(string[6:])
        
        return h+m+s
    
    def plot(self,filt,times=False):
        fig = plt.figure(0)
        plt.clf()
        ax1 = fig.add_subplot(111)
        
        idx = np.where(self.filters == filt)
        
        ax1.plot(self.times[idx],self.levels[idx],'.')#,t,l,'-')
        try: ax1.plot(self.ref_levels[filt]['times'],\
                          [self.ref_levels[filt]['level'] for i in range(2)],linestyle='',marker='s')
        except KeyError: pass
        ax1.set_ylabel('$V$')
        ax1.set_xlabel('time [s]')
        fig.suptitle(os.popen('pwd').readlines()[0]+datetime.now().isoformat(' '))
        
        if times:
            for time in times:
                ax1.axvline(ls=':',x=time,color='r')
        fig.show()
        return

def webit(name):
    
    path = '/d/www/eigenbrot/MANGA/'+name
    os.system('mkdir '+path)
    os.system('cp *.dat '+path)
    os.system('cp *.txt '+path)
    os.system('cp *.ps '+path)
    os.system('cp *.jpg '+path)
    os.system('cp *.pdf '+path)
    ht = open(path+'/index.html','wb')
    ht.write('<HTML><TITLE>FRD Bench data output '+name+'</TITLE>\n'
             +'<body bgcolor="EEEEEE" link="#333FF" vlink="#5500BB">\n'
             +'<h1>Data products for '+name+'\n'
             +'<h2>Web summary generated on '+datetime.now().isoformat(' ')
             +'</h2></h1>\n<hr noshade>\n')

    dat_list = glob.glob(path+'/*.txt')+glob.glob(path+'/*.dat')
    for dat in dat_list:
        link = dat[dat.rfind('/')+1:]
        ht.write('<a href='+link+' type="text/plain">'+link+'</a><br>\n')

    ht.write('<h3>FRD plots</h3>Click images for jpg\n')
    ht.write('<h4>All-fiber grid</h4>\n'
             +'<a href=allfibers.pdf>pdf</a> <a href=allfibers.jpg>jpg</a><br>\n'
             +'<a href=allfibers.jpg><img src=allfibers.jpg height=300 width=300></a>\n')
    plt_list = [i[i.rfind('/')+1:i.rfind('.')] for i in glob.glob(path+'/*.ps')]
    ht.write('<h4>Individual fibers</h4>\n<table>\n<tr>')
    i = 1
    for plot in plt_list:
        ht.write('<td>\n'+plot+':<br>\n')
        ht.write('<a href='+plot+'.ps>ps</a> <a href='+plot+'.jpg>jpg</a><br>\n')
        ht.write('<a href='+plot+'.jpg><img src='+plot+'.jpg height=300 width=232></a>\n</td>\n')
        if i % 5 == 0: ht.write('</tr>\n<tr>')
        i += 1
    ht.write('</tr></table>\n</body>\n</HTML>')
    ht.close()
    
    old_m = open('/d/www/eigenbrot/MANGA/index.html')
    lines = old_m.readlines()
    old_m.close()

    lines[3] = '<h2>Last updated '+datetime.now().isoformat(' ')+'</h2></h1>\n'
    m = open('/d/www/eigenbrot/MANGA/index.html','wb')
    m.writelines([l for l in lines[:-1]])
    m.write('<a href='+name+'>'+name+'</a><br>\n</body></HTML>')
    m.close()

def plot_helper(datafile,filename,title):
    
    pp = PDF(filename)

    rc('text', usetex=False)
    rc('font', family='serif')
    rc('font', size=9.0)
    rc('axes', linewidth=0.4)
    
    fiber_plot_mapping = dict([('0,-2',11),
                              ('1,-2',10),
                              ('2,-2',9),
                              ('-1,-1',12),
                              ('0,-1',3),
                              ('1,-1',2),
                              ('2,-1',8),
                              ('-2,0',13),
                              ('-1,0',4),
                              ('0,0',1),
                              ('1,0',7),
                              ('2,0',19),
                              ('-2,1',14),
                              ('-1,1',5),
                              ('0,1',6),
                              ('1,1',18),
                              ('-2,2',15),
                              ('-1,2',16),
                              ('0,2',17)])

    hdus = pyfits.open(datafile)[1:]
    
    numpages = len(hdus)/20
    if len(hdus) % 20 != 0: numpages += 1

    hi = 0
    for p in range(numpages):
        fig = plt.figure()
        grid = AG(fig,111,
                  nrows_ncols = (4,5),
                  axes_pad = 0.0,
                  label_mode = 'L',
                  aspect = False,
                  share_all = False)
        
        gg = 0
        for ax in grid:
            
            try: h = hdus[hi]
            except IndexError: continue
            plot_data = h.data
            sloan_metric = h.header['SLOAN']
            fiber_pos = h.header['OUTPOS']
            if fiber_pos == 'NA':
                fiber_pos = h.header['FIBERPOS']
                posstring = '(fiber input position)'
            else:
                posstring = '(fiber output position)'
            try: plot_pos = fiber_plot_mapping[fiber_pos] - 1
            except KeyError: plot_pos = gg
#            ax = grid[plot_pos]
            ax.plot(plot_data[1],plot_data[4],'r',lw=0.4)
            ax.plot(plot_data[3],plot_data[5],'b',lw=0.4)
#            ax.plot(plot_data[2],plot_data[5],'b:',lw=0.4)
#            ax.plot(plot_data[0],plot_data[4],'r:',lw=0.4)
            ax.plot(plot_data[3],plot_data[6],linestyle='-.',color='b',lw=0.4)
            ax.text(9,0.7,'({0})\n{1:4.3f}'.format(fiber_pos,sloan_metric))
            ax.set_xlim(2,30)
            ax.set_ylim(0,1.1)
            ax.set_xscale('log')
            ax.set_xticks([5,10,20])
            ax.set_xticklabels(['$5$','$10$','$20$'])
            hi += 1
            gg += 1

        if p == range(numpages)[-1]:
            grid[-1].set_xlim(2,30)
            grid[-1].set_ylim(0,1.1)
            grid[-1].set_xscale('log')
            grid[-1].set_xticks([])
            grid[-1].text(2.5,0.7,posstring+'\nthroughput value',fontsize=8)
            grid[-1].hlines([0.5,0.4,0.3],3,5,colors=['b','b','r'],linestyles=['solid','dashdot','solid'])
            grid[-1].text(6,0.5,'Fiber normalized',fontsize=5.0)
            grid[-1].text(6,0.4,'Fiber relative',fontsize=5.0)
            grid[-1].text(6,0.3,'Direct normalized',fontsize=5.0)
            
        fig.text(0.5,0.04,'$f$-ratio',fontsize=12.0)
        fig.text(0.06,0.65,'Normalized Encircled Energy',rotation='vertical',fontsize=11.0)
        fig.suptitle(title)
    
        pp.savefig(fig)
    
    pp.close()

    return

def ring_helper(metric_file):


    inpos, outpos = np.loadtxt(metric_file,usecols=(0,1),unpack=True,dtype=str)
    sloan = np.loadtxt(metric_file,usecols=(6,),unpack=True)
    if type(inpos) != np.ndarray: inpos = np.array([inpos])
    if type(outpos) != np.ndarray: outpos = np.array([outpos])
    if sloan.shape == (): sloan = np.array([sloan])

    if outpos[0].find(',') == -1: out = False
    else: out = True

    dinring = {}
    inringarr = np.array([])
    if out:
        doutring = {}
        outringarr = np.array([])
        ringavg = np.array([])

    for i in range(inpos.size):
        inring = findring(inpos[i])
        inringarr = np.append(inringarr,inring)
        try: dinring[inring].append(sloan[i])
        except KeyError: dinring[inring] = [sloan[i]]

        if out:
            outring = findring(outpos[i])
            outringarr = np.append(outringarr,outring)
            ringavg = np.append(ringavg,(inring + outring)/2.)
            try: doutring[outring].append(sloan[i])
            except KeyError: doutring[outring] = [sloan[i]]

    outname = metric_file.split('.')[0]+'_summary.txt'
    f = open(outname,'w')
    f.write('# Summary written on: '+datetime.now().isoformat(' ')+'\n'
            +'# Metric file: '+metric_file+'\n'
            +'#\n')

    for iring in dinring.keys():
        f.write('# Input Ring {:n} mean sloan: {:5.4f}\n'.format(iring,np.mean(dinring[iring])))

    if out:
        f.write('#\n')
        for oring in doutring.keys():
            f.write('# Output Ring {:n} mean sloan: {:5.4f}\n'.format(oring,np.mean(doutring[oring])))

    f.write('#\n# {:10}= '.format('Fiber_pos')+'fiber input position\n'
            +'# {:10}= '.format('Out_pos')+'fiber output position\n'
            +'# {:10}= '.format('InRing')+'Input Ring\n')
    if out: 
        f.write('# {:10}= '.format('OutRing')+'Output Ring\n'
                + '# {:10}= '.format('RingAvg')+'Average of Input and Output rings\n')
    f.write('# {:10}= '.format('sloan')+'fiber within f/4 / direct within f/4 (uncorrected)\n'
            +'#\n')
    if out: 
        f.write(str('#{:>10}{:>10}{:>9}{:>9}{:>9}{:>9}\n'\
                        .format('Fiber_pos','Out_pos','InRing','OutRing','RingAvg','sloan'))
                +str('#{:>10}{:>10}{:>9}{:>9}{:>9}{:>9}\n'.format(*range(6))))
    else:
        f.write(str('#{:>10}{:>10}{:>9}{:>9}\n'\
                        .format('Fiber_pos','Out_pos','InRing','sloan'))
                +str('#{:>10}{:>10}{:>9}{:>9}\n'.format(*range(4))))


    if out:
        sidx = np.argsort(ringavg)
        for k in range(inringarr.size):
            f.write('{:>11}{:>10}{:>9n}{:>9n}{:9.1f}{:9.4f}\n'\
                        .format(inpos[sidx][k],outpos[sidx][k],inringarr[sidx][k],
                                outringarr[sidx][k],ringavg[sidx][k],sloan[sidx][k]))
    else:
        sidx = np.argsort(inringarr)
        for k in range(inringarr.size):
            f.write('{:>11}{:>10}{:>9n}{:9.4f}\n'\
                        .format(inpos[sidx][k],outpos[sidx][k],inringarr[sidx][k],sloan[sidx][k]))
    return

def findring(pos):
    xpos, ypos = map(int,pos.split(','))
    
    return int(np.ceil(((xpos**2. + ypos**2.)/2.)**0.5))
                 
def hex_plot_helper(datafile,outputfile,title,numfibers):
    
    if numfibers == 127:
        rank = 6
    elif numfibers == 61:
        rank = 4
    elif numfibers == 19:
        rank = 2

    low = 0.90
    high = 0.96

    in_ax = plot_hex(rank)
    out_ax = plot_hex(rank)
    in_plist = []
    out_plist = []
    in_tlist = []
    out_tlist = []

    hdus = pyfits.open(datafile)[1:]
    for h in hdus:
        fiberpos = h.header['FIBERPOS']
        outpos = h.header['OUTPOS']
        sloan = h.header['SLOAN']
        in_ax = shade_circle(in_ax,fiberpos,sloan,in_plist,in_tlist)
        out_ax = shade_circle(out_ax,outpos,sloan,out_plist,out_tlist)

    pp = PDF(outputfile)
    in_ax.figure.suptitle(title+' - INPUT\n'+datetime.now().isoformat(' '))
    in_p = PatchCollection(in_plist,cmap=plt.get_cmap('RdYlGn'),
                           norm=matplotlib.colors.Normalize(vmin=low,vmax=high),
                           alpha=0.3)
    in_p.set_array(np.array(in_tlist))
    in_ax.add_collection(in_p)
    in_ax.figure.colorbar(in_p)

    out_ax.figure.suptitle(title+' - OUTPUT\n'+datetime.now().isoformat(' '))
    out_p = PatchCollection(out_plist,cmap=plt.get_cmap('RdYlGn'),
                            norm=matplotlib.colors.Normalize(vmin=low,
                                                             vmax=high),
                            alpha=0.3)
    out_p.set_array(np.array(out_tlist))
    out_ax.add_collection(out_p)
    out_ax.figure.colorbar(out_p)

    pp.savefig(in_ax.figure)
    pp.savefig(out_ax.figure)
    pp.close()
    return

def plot_hex(rank):

    offset = (rank % 2)*0.5
    row_lengths = np.concatenate((
        np.arange(rank) + rank + 1,
        [rank*2 + 1],
        (np.arange(rank) + rank + 1)[::-1]))
    y_offsets = np.arange(rank*2 + 1)*2 - rank*2
    
    fig = plt.figure()
    ax = fig.add_subplot(111,aspect='equal')
    ax.set_axis_off()
    theta = np.arange(0,2*np.pi,0.01)

    for (row,y) in zip(row_lengths,y_offsets):
         for x in (np.arange(row) - row/2 + offset)*2:
            ax.plot(np.cos(theta)+x,np.sin(theta)+y,'k')

         offset = abs(offset - 0.5)

    return ax

def shade_circle(ax,coords,tput,plist,tlist):

    x_coord, y_coord = [int(i) for i in coords.split(',')[0:2]]
    
    x_coord = x_coord*2 + y_coord
    y_coord *= 2

    circ = Circle((x_coord,y_coord),radius=1)
    plist.append(circ)
    tlist.append(tput)
    ax.text(x_coord,y_coord,'{:4.3f}'.format(tput),
            fontsize=9,ha='center',va='center')

    return ax


if __name__ == '__main__':
    sys.exit(main())
