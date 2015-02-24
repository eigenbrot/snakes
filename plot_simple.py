import numpy as np
import GradPak_plot as GPP
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

rc('text', usetex=False)
rc('font', family='serif')
rc('font', size=10.0)
rc('axes', linewidth=1.0)
rc('lines', linewidth=0.4)
rc('ps', usedistiller='Xpdf')
rc('xtick', labelsize=10.0)
rc('ytick', labelsize=10.0)

AGES = np.array([0.0050119, 0.0251188, 0.101518, 0.28611901, 0.64054298,
                 0.9047920, 1.434, 2.5, 5.0, 10.0])

def plot_idl(inputfile):

    data = np.loadtxt(inputfile)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fibers = data[:,0]
    for i in np.arange(data.shape[1] - 1) + 1:
        
        ax.plot(fibers,data[:,i],'.',label='{} Myr'.format(AGES[i-1]))
        
    ax.legend(loc=0,numpoints=1)
    
    fig.show()

    return data
    
def plot_age(inputfile):

    data = np.loadtxt(inputfile)
    fibers = data[:,0]
    ageidx = np.argmax(data[:,1:],axis=1)
    
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('Fiber number')
    ax.set_ylabel('Age of highest light fraction population [Myr]')
    ax.set_yticks(range(10))
    ax.set_yticklabels(AGES)
    ax.set_ylim(-0.5,11)
    ax.plot(fibers,ageidx,'.')
    
    ax.figure.show()

    return data

def plot_age_hist(inputfile, outputfile, exclude=[]):

    data = np.loadtxt(inputfile)
    fibers = data[:,0]

    agelabels = ['{:4.0e}'.format(i) if np.log10(i) < -1 
                 else '{:3.1f}'.format(i)
                 for i in AGES]
    print agelabels

    pp = PDF(outputfile)
    for i in range(fibers.size):
        print i
        ax = plt.figure().add_subplot(111)
        ax.bar(np.arange(AGES.size),data[i,1:-2],align='center',width=1.0)
        ax.set_ylabel('Light fraction')
        ax.set_xlim(-1,AGES.size)
        ax.set_xticks(np.arange(AGES.size))
        ax.set_xticklabels(agelabels)
        ax.set_xlabel('Age [Gyr]')
        MLWA = data[i,-2]
        ax.set_title('Fiber {}\nMLWA = {:4.3f} Gyr'.format(i+1,MLWA))
        pp.savefig(ax.figure)
        plt.close(ax.figure)
    
    pp.close()

    return

def plot_maps(inputfile, outputfile, eps=False, exclude=[], sky=False,
              labelfibers = True, MLWA = True):

    if MLWA:
        data = np.loadtxt(inputfile,usecols=(12,),unpack=True)
        label = 'Mean Light Weighted Age [Gyr]'
        minval = np.nanmin(data)
        maxval = np.nanmax(data)
    else:
        data = np.loadtxt(inputfile,usecols=(11,),unpack=True)
        label = 'Mean Mass Weighted Age [Gyr]'
        minval = np.nanmin(data)#7
        maxval = np.nanmax(data)#10

    map_ax = GPP.plot_img(data,
                          fitsfile=\
                          '/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits',
                          imrot=67.0
                          pa=295.787,
                          center = [35.6034125,42.32349444],
                          clabel=label,
                          method='cubic',
                          cmap='gnuplot2',
                          exclude=exclude,
                          sky=sky,
                          minval=minval,
                          maxval=maxval)
    
    fiber_ax = GPP.plot(data,
                        fitsfile='/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits',
                        imrot=67.0,
                        pa=295.787,
                        center = [35.6034125,42.32349444],
                        clabel=label,
                        cmap='gnuplot2',
                        labelfibers=labelfibers,
                        exclude=exclude,
                        sky=sky,
                        minval=minval,
                        maxval=maxval)

    if eps:
        fiber_name = outputfile+'_fibers.eps'
        map_name = outputfile+'_map.eps'

        fiber_ax.figure.savefig(fiber_name,format='eps')
        map_ax.figure.savefig(map_name,format='eps')
    else:
        pp = PDF(outputfile)
        pp.savefig(fiber_ax.figure)
        pp.savefig(map_ax.figure)
        pp.close()

    # plt.close(fiber_ax.figure)
    # plt.close(map_ax.figure)
    
    return fiber_ax

def plot_heights(inputfile, outputfile, title=''):

    
    MMWA, MLWA, SNR = np.loadtxt(inputfile,usecols=(11,12,14),unpack=True)
    
    pp = PDF(outputfile)
    
    ax = GPP.plot_rows(MMWA, weights=SNR, 
                       ylabel='Age [Gyr]',label='MMWA', kpc_scale = 0.0485)
    GPP.plot_rows(MLWA, weights=SNR, label='MLWA', ax = ax,
                  kpc_scale = 0.0485)
    ax.legend(loc=0,scatterpoints=1,numpoints=1,frameon=False)
    ax.set_xlim(-0.5,3)
    ax.set_title(title)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)

    return
