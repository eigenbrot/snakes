# originally from /d/monk/eigenbrot/SALT/2011-3-UW_002/20120215/flat_test

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import bottleneck as bn
from datetime import datetime

def make_ratios():

    '''read in the data'''
    raw_data = pyfits.open('../reduced/tiESO_z0_MgI.fits')[0].data
    old_data = pyfits.open('to45.fits')[0].data
    MB_data = pyfits.open('tn45.fits')[0].data

    # raw_data = pyfits.open('tsky_sub.fits')[0].data
    # old_data = pyfits.open('tosky_sub.fits')[0].data
    # MB_data = pyfits.open('tnsky_sub.fits')[0].data

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)

    for data, name in zip([raw_data,old_data,MB_data], 
                          ['No flat','Nightly flat','Reconstructed flat']):
        
        x, line1, line2 = get_lines(data,20)

        ratio = line1/line2

        zeroed_ratio = ratio/ bn.nanmean(ratio[30:])
        
        ax.plot(x,ratio,label=name)
        ax1.plot(x,zeroed_ratio,label=name)

    ax.set_xlabel('$\mathrm{Slit\ position\ [px]}$')
    ax.set_ylabel('$f(4861\mathrm{nm})/f(5198\mathrm{nm})$')
    ax.legend(loc=0)
#    ax.set_ylim(0,2.5)
#    ax.set_title(datetime.now().isoformat(' '))

    ax1.set_xlabel('$\mathrm{Slit\ position\ [px]}$')
    ax1.set_ylabel('$f(4861\mathrm{nm})/f(5198\mathrm{nm})$')
    ax1.legend(loc=0,numpoints=1)
    # ax1.set_ylim(0.5,1.5)
    # ax1.set_title('Sky-subtracted object')
    # ax1.text(250,1.4,'Plot 2',fontsize=17)
    ax1.set_ylim(0.9,1.1)
    ax1.set_title('Object frame only')
    ax1.text(250,1.07,'Plot 1',fontsize=17)

#    fig.show()
    fig1.show()
    return

def get_lines(data,smooth_amount):
    
    line1 = np.mean(data[:,770:790],axis=1)
    bkgrd1 = np.mean(data[:,750:770] + data[:,791:811],axis=1)
    line2 = np.mean(data[:,2435:2453],axis=1)
    bkgrd2 = np.mean(data[:,2417:2435] + data[:,2454:2472],axis=1)

    # line1 -= bn.move_median(bkgrd1,smooth_amount)
    # line2 -= bn.move_median(bkgrd2,smooth_amount)

    # line1 -= bkgrd1
    # line2 -= bkgrd2

    line1_smooth = bn.move_median(line1,smooth_amount)
    line2_smooth = bn.move_median(line2,smooth_amount)
    
    x = np.arange(line1_smooth.size)

    return x, line1_smooth, line2_smooth
