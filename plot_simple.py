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
                          imrot=67.0,
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

def all_maps(output,MLWA=True,labelfibers=False):

    pp = PDF(output)
    centers = [
        [35.603413,42.323494],
        [35.615154,42.3426],
        [35.636325,42.37935],
        [35.625617,42.35935],
        [35.591642,42.304247],
        [35.648196,42.401217]]

    ax = None
    for i in range(6):
        print i
        inputfile = 'P{}.dat'.format(i+1)
        if MLWA:
            data = np.loadtxt(inputfile,usecols=(12,),unpack=True)
            label = 'Mean Light Weighted Age [Gyr]'
            minval = 0
            maxval = 10#np.nanmax(data)
        else:
            data = np.loadtxt(inputfile,usecols=(11,),unpack=True)
            label = 'Mean Mass Weighted Age [Gyr]'
            minval = 5#np.nanmin(data)#7
            maxval = 10#np.nanmax(data)#10

        ax = GPP.plot(data,
                      ax=ax,
                      figsize=(8,4),
                      fitsfile=\
                      '/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits',
                      centpos=[35.63689,42.34633],
                      imrot=67,#64.213,
                      pa=295.787,
                      center = centers[i],
                      clabel=label,
                      cmap='gnuplot2',
                      labelfibers=labelfibers,
                      exclude=[],
                      sky=False,
                      minval=minval,
                      maxval=maxval)
        
    ax.set_xlim(160,755)
    ax.set_ylim(350,600)
    pp.savefig(ax.figure)
    pp.close()
    return ax

def all_heights(output,err=True):

    datname = output.split('.pdf')[0]+'.dat'
    f = open(datname,'w')
    f.write('#{:5}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}\n#\n'.\
            format('','height [kpc]','age [Gyr]',
                   'age eom','age stderr','AV','AV eom','AV stderr'))

    pp = PDF(output)
    symblist = ["o","^","v","s","*","x"]
    rlist = [4.4,0.8,-6.4,-2.3,8.1,-10.0]

    ax = None
    AVax = None
    for i in range(6):
        print i+1
        inputfile = 'P{}.dat'.format(i+1)
        MMWA, MLWA, TAUV, SNR = np.loadtxt(inputfile,usecols=(11,12,13,14),
                                           unpack=True)
        
        ax, tmpz, tmpage, tmperr, tmpstd =  GPP.plot_rows(MLWA, 
                                                          weights=SNR,
                                                          ax=ax,
                                                          label='{}'.\
                                                          format(rlist[i]),
                                                          fullout=True,
                                                          kpc_scale = 0.0485,
                                                          err=err, 
                                                          marker=symblist[i], 
                                                          linestyle='', 
                                                          color='k')
        AVax, _, tmpAV, tmpAVerr, tmpAVstd = GPP.plot_rows(TAUV*1.086,
                                                           weights=SNR,
                                                           ax=AVax,
                                                           label='{}'.\
                                                           format(rlist[i]),
                                                           fullout=True,
                                                           kpc_scale=0.0458,
                                                           err=err,
                                                           marker=symblist[i],
                                                           linestyle='',
                                                           color='k')

        f.write('\n# P{} '.format(i+1)+92*'#'+'\n# r ~ {} kpc\n'.\
                format(rlist[i]))
        for i in range(tmpz.size):
            f.write('{:6}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}\n'.\
                    format('',tmpz[i],tmpage[i],tmperr[i],
                           tmpstd[i],tmpAV[i],tmpAVerr[i],tmpAVstd[i]))

        try:
            z = np.vstack((z,tmpz))
            age = np.vstack((age,tmpage))
            err = np.vstack((err,tmperr))
            std = np.vstack((std,tmpstd))
            AV = np.vstack((AV,tmpAV))
            AVstd = np.vstack((AVstd,tmpAVstd))
            AVerr = np.vstack((AVerr,tmpAVerr))
        except UnboundLocalError:
            z = tmpz
            age = tmpage
            err = tmperr
            std = tmpstd
            AV = tmpAV
            AVstd = tmpAVstd
            AVerr = tmpAVerr

    bigz = np.mean(z,axis=0)
    bigage = np.mean(age,axis=0)
    bigerr = np.sqrt(
        np.sum(err*(age - bigage)**2,axis=0)/
        ((err.shape[0] - 1.)/(err.shape[0]) * np.sum(err,axis=0)))
    bigAV = np.mean(AV,axis=0)
    bigAVerr = np.sqrt(
        np.sum(AVerr*(AV - bigAV)**2,axis=0)/
        ((AVerr.shape[0] - 1.)/(AVerr.shape[0]) * np.sum(AVerr,axis=0)))

    ax.plot(bigz, bigage)
    ax.fill_between(bigz,bigage-bigerr,bigage+bigerr,alpha=0.1)
    ax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    ax.set_xlim(-0.1,2.6)
    ax.set_ylabel('Light-weighted age [Gyr]')

    AVax.plot(bigz,bigAV)
    AVax.fill_between(bigz,bigAV-bigAVerr,bigAV+bigAVerr,alpha=0.1)
    AVax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    AVax.set_xlim(-0.1,2.6)
    AVax.set_ylabel(r'$A_V$')

    pp.savefig(ax.figure)
    pp.savefig(AVax.figure)
    plt.close(ax.figure)

    pp.close()
    f.close()

    return
