import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import ADEUtils as ADE
from MANGA_bench import Noodle as N
from MANGA_bench import thePot
from ConfigParser import ConfigParser
import os

def gen_resp():

#    level_list = [4,5,6,7,8,9,10,11,12]
    level_list = [1,2,3,4,5,6,7]
   
    PD_vals = np.array([])
    fits_vals = np.array([])

    PD_err = np.array([])
    fits_err = np.array([])

    for level in level_list:

        PD_signal = np.loadtxt('PD_{}.txt'.format(level),usecols=(1,),unpack=True)
        fits_raw = pyfits.open('PD_{}_FINAL.fits'.format(level))[0].data

        center = ADE.centroid(fits_raw)
        dist = ADE.dist_gen(fits_raw,center) * 0.024
        idx = np.where(dist < 6.75)
        fits_data = fits_raw[idx]

        PD_vals = np.append(PD_vals,np.mean(PD_signal))
        PD_err = np.append(PD_err,np.std(PD_signal))

        fits_vals = np.append(fits_vals,np.mean(fits_data))
        fits_err = np.append(fits_err,np.std(fits_data))

    level_vec = np.array(level_list)/12.
    fig = plt.figure()
    ax1 = fig.add_subplot(211)

    PD_plot = PD_vals/PD_vals[-1]
    PD_err_plot = PD_err/PD_vals[-1]

    fits_plot = fits_vals/fits_vals[-1]
    fits_err_plot = fits_err/fits_vals[-1]

#    ax1.errorbar(PD_vals,fits_vals,xerr=PD_err,yerr=fits_err,linestyle='x')
    ax1.plot(fits_vals,PD_vals,'.')
    # ax1.axhline(y=95693.7799/2,linestyle=':')
    # ax1.text(1.,95000,'Full Well')

    ax1.set_ylabel('PD Voltage')
    ax1.set_xlabel('CCD Mean Counts')
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')

    ax2 = fig.add_subplot(212)
    ax2.plot(fits_vals,fits_vals/PD_vals,'.')
    ax2.set_xlabel('CCD Mean Counts')
    ax2.set_ylabel('CCD/PD')
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')

    fig.show()

    return level_vec, fits_vals, fits_err, PD_vals, PD_err

def get_stats():

    level_list = [4,5,6,7,8,9,10,11,12]
    
    fits_vals = np.array([])
    fits_err = np.array([])
    
    for level in level_list:
        
        fits_list = glob.glob('PD_{}_*_ds.fits'.format(level))
        
        imarray = np.array([])
        
        for image in fits_list:
            print image
            data = pyfits.open(image)[0].data
            
            center = ADE.centroid(data)
            dist = ADE.dist_gen(data,center) * 0.024
            idx = np.where(dist < 6.75)
            fits_data = data[idx]

            imarray = np.append(imarray,np.mean(fits_data))
            
        fits_vals = np.append(fits_vals,np.mean(imarray))
        fits_err = np.append(fits_err,np.std(imarray))

    return fits_vals, fits_err

def bias_resp():

    levels = np.array([])
    errs = np.array([])
    times = np.array([])

    file_list = glob.glob('BIAS*.FIT')

    for image in file_list:
        print image,
        hdu = pyfits.open(image)[0]
        data = hdu.data
        timestr = hdu.header['TIME-OBS']
        print timestr
        time = np.float(timestr[6:]) + np.float(timestr[3:5])*60. +\
            np.float(timestr[0:2])*3600.
        
        times = np.append(times,time)
        levels = np.append(levels,np.mean(data))
        errs = np.append(errs,np.std(data))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(times,levels,yerr=errs,linestyle='.')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Counts [ADU]')
    ax.set_title('Mean mean: {:4.3f}\nMean std: {:4.3f}'.format(np.mean(levels),np.std(levels)))
    fig.show()

    return times, levels, errs


def jump_test(inifile):
    ''' to be run in /d/monk/eigenbrot/MANGA/20121015 '''

    options = ConfigParser()
    options.read(inifile)

    finald = []
    
    for d in ['d1','d2','d3','d4']:
        nood = N(options)
        nood.get_darks()
        nood.fill_dict(d,nood.direct)
        nood.sub_darks(nood.direct)
        nood.ratios = {'direct': {'data': {'V':{}}}}
        nood.direct_to_ratios(nood.direct)
        nood.combine()
        finald.append(nood.ratios['direct']['data']['V']['direct']['final'])
        os.system('rm *_ds.fits')
    
    print "reduction finished"

    countarr = np.array([])
    voltarr = np.array([])

    pot = thePot(options)

    for image in finald:
        
        hdu = pyfits.open(image)[0]
        fits_raw = hdu.data

        center = ADE.centroid(fits_raw)
        dist = ADE.dist_gen(fits_raw,center) * 0.024
        idx = np.where(dist < 6.75)
        fits_data = fits_raw[idx]
        countarr = np.append(countarr,np.mean(fits_data))
        
        stime = hdu.header['STARTIME']
        etime = hdu.header['ENDTIME']
        voltage = pot.get_voltage(stime,etime,'V')
        voltarr = np.append(voltarr,voltage)

    print "\nImage  :{:>15}{:>15}{:>15}{:>15}\nCounts :{:15.3E}{:15.3E}{:15.3E}{:15.3E}\nVoltage:{:15.2f}{:15.2f}{:15.2f}{:15.2f}".format(*finald+countarr.tolist()+voltarr.tolist())
    print "-"*(4*15+8)
    print "Ratio  :{:15.3E}{:15.3E}{:15.3E}{:15.3E}\n".format(*(countarr/voltarr).tolist())

    return countarr, voltarr, countarr/voltarr
