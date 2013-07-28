import numpy as np
#import ADEUtils as ADE
import ADESALT as sa
#import gc
#import Salty2 as salty
import moment_tests as mt
#import bottleneck as bn
from datetime import datetime
import matplotlib
plt = matplotlib.pyplot
rc = matplotlib.rc
import pyfits
from glob import glob
from matplotlib.backends.backend_pdf import PdfPages as PDF

'''This file is for making grid plots of data and various morphologies.
'''


def grid_plot(msfiles, flips, sims, simleg):

    font = {'size':4}
    rc('font',**font)

    for msfile, flip, simfiles in zip(msfiles,flips,sims):
        print msfile
        height = int(msfile.split('_')[1][1:])
        if height == 5:
            height = 0.5
        name = msfile.split('.ms.fits')[0] + '_binplot.pdf'
        simname = '_'.join(simfiles[0].split('_')[0:2]) + '_binplot.pdf'
        pp = PDF(name)
        simpp = PDF(simname)
        tempr, _, _ = sa.openslay(msfile.split('.ms')[0] + '.slay.fits',
                                  flip=flip)        
        print tempr
        tmphead = pyfits.open(msfile)[0].header
        ap = 1
        fig0 = plt.figure(figsize=(11,8.5))
        fig1 = plt.figure(figsize=(11,8.5))
        for radius in np.sort(tempr):
            print '\t{}'.format(radius)
            tmprstr = tmphead['APNUM{}'.format(ap)].split(' ')
            tmpwidth = int(tmprstr[3]) - int(tmprstr[2])
            tmpwidth *= 0.118*8*34.1e3/206265
            
            #############
            ax0 = fig0.add_subplot(3,4,ap)
            ax0.set_title('{}\n{}'.\
                              format(msfile,datetime.now().isoformat(' ')),
                          fontsize=4)
            ax0.set_xlabel('Velocity [km/s]')
            ax0.set_ylabel('ADU/s')
            ax0.text(0.1,0.95,'z ~ {:}$h_z$\nr = {:5.3f}kpc\ndr = {:5.3f}kpc'\
                         .format(height,radius,tmpwidth),
                     transform=ax0.transAxes,
                     ha='left',va='top')
            ax0.set_xlim(-600,600)
            sa.plot_line(msfile,radius,ax=ax0,plot=False,flip=flip,velo=True,
                         baseline=1,linewidth=0.4)

            #############
            ax1 = fig1.add_subplot(3,4,ap)
            ax1.set_title('{}\n{}'.\
                              format(simname.split('_binplot')[0]
                                     ,datetime.now().isoformat(' ')),
                          fontsize=4)
            ax1.text(0.1,0.65,'z = {:}$h_z$\nr = {:5.3f}kpc\ndr = {:5.3f}kpc'\
                         .format(height,radius,tmpwidth),
                     transform=ax1.transAxes,
                     ha='left',va='top')
            ax1.set_xlabel('Velocity [km/s]')
            ax1.set_ylabel('Normalized flux')
            ax1.set_xlim(-600,600)
            
            '''The negative 1 is there because all my sims are created
            backwards, by convention
            '''
            for simfile in simfiles:
                variable = pyfits.open(simfile)[0].header[simleg.upper()]
                mt.do_line(simfile,-1*radius,1,ax=ax1,plot=False,
                           label=str(variable),rwidth=tmpwidth,linewidth=0.8)
            ax1.legend(loc=0,title=simleg)
            ap += 1

        pp.savefig(fig0)
        pp.close()
        simpp.savefig(fig1)
        simpp.close()

    return

def linear_flare_script():

    mslist = glob('*_bin??.ms.fits')
    simlist1 = glob('sim*083.fits')
    simlist2 = glob('sim*290.fits')
    simlist3 = glob('sim*025.fits')
    
    mslist.sort()
    simlist1.sort()
    simlist2.sort()
    simlist3.sort()
    
    print mslist, simlist1, simlist2, simlist3

    grid_plot(mslist,[True,False,False,False],
              zip(simlist1, simlist2, simlist3),
              'h_zR')

    return

def ring_script():
    
    mslist = glob('*_bin??.ms.fits')
    simlist1 = glob('sim*w010.fits')
    simlist2 = glob('sim*w085.fits')
    
    mslist.sort()
    simlist1.sort()
    simlist2.sort()

    print mslist, simlist1, simlist2
    
    grid_plot(mslist,[True,False,False,False],
              zip(simlist1,simlist2),
              'r_w')

    return

def ring_script_radius():

    mslist = glob('*_bin??.ms.fits')
    simlist1 = glob('sim*R0100.fits')
    simlist2 = glob('sim*R2100.fits')
    
    mslist.sort()
    simlist1.sort()
    simlist2.sort()

    print mslist, simlist1, simlist2
    
    grid_plot(mslist,[True,False,False,False],
              zip(simlist1,simlist2),
              'r_R')

    return
