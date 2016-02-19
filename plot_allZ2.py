import numpy as np
import bottleneck as bn
import GradPak_plot as GPP
import scipy.ndimage as spnd
import scipy.interpolate as spi
#import model_A as mA
import matplotlib.pyplot as plt
import pyfits
#import pywcs
import time
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib import rc
from matplotlib import colors as mplcolors
from matplotlib.ticker import ScalarFormatter
import glob
import time
glob = glob.glob

rc('text', usetex=False)
rc('font', family='serif')
rc('font', size=10.0)
rc('axes', linewidth=1.0)
rc('lines', linewidth=0.4)
rc('patch', linewidth=0.1)
rc('ps', usedistiller='Xpdf')
rc('xtick', labelsize=10.0)
rc('ytick', labelsize=10.0)

def all_heights(output, inputprefix='NGC_891', err=True, binned=True):

    datname = output.split('.pdf')[0]+'.dat'
    f = open(datname,'w')
    f.write('#{:5}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}{:>13}\n#\n'.\
            format('','height [kpc]','age [Gyr]',
                   'age eom','age stderr',
                   'AV','AV eom','AV stderr',
                   'Z/Z_sol','Z eom', 'Z stderr'))

    pp = PDF(output)
    symblist = ["o","^","v","s","*","x"]
    colorlist = ['b','c','g','y','m','r']
    kpc_scale = 0.0485
    plist = [6,3,4,2,1,5]
    zcorrlist = [-0.01, -0.113, -0.225, -0.07, -0.142, -0.004]
    rlist = [-10.0, -6.4, -2.3, 0.8, 4.4, 8.1]

    ageax = None
    AVax = None
    Zax = None

    for i in range(6):
        
        inputfile = glob('{}*P{}*allz2*.dat'.format(inputprefix,plist[i]))[0]
        print inputfile
        if binned:
            fitsname = glob('{}*P{}*.ms*fits'.format(inputprefix,plist[i]))[0]
            print fitsname
            binhead = pyfits.open(fitsname)[0].header
        else:
            binhead = False
            
        MLWA, MLWZ, TAUV, SNR = np.loadtxt(inputfile, usecols=(62,64,66,67),
                                           unpack=True)

        MLWZ = np.log10(MLWZ)
        ageax, tmpz, tmpage, tmperr, tmpstd = GPP.plot_rows(MLWA, 
                                                            binheader=binhead,
                                                            weights=SNR,
                                                            ax=ageax,
                                                            label='{}'.\
                                                            format(rlist[i]),
                                                            fullout=True,
                                                            kpc_scale=kpc_scale,
                                                            zcorr=zcorrlist[i],
                                                            err=err, 
                                                            marker=symblist[i], 
                                                            linestyle='',
                                                            color=colorlist[i])

        AVax, _, tmpAV, tmpAVerr, tmpAVstd = GPP.plot_rows(TAUV*1.086,
                                                           binheader = binhead,
                                                           weights=SNR,
                                                           ax=AVax,
                                                           label='{}'.\
                                                           format(rlist[i]),
                                                           fullout=True,
                                                           kpc_scale=kpc_scale,
                                                           zcorr=zcorrlist[i],
                                                           err=err,
                                                           marker=symblist[i],
                                                           linestyle='',
                                                           color=colorlist[i])

        Zax, _, tmpZ, tmpZerr, tmpZstd = GPP.plot_rows(MLWZ, 
                                                       binheader = binhead,
                                                       weights=SNR,
                                                       ax=Zax,
                                                       label='{}'.\
                                                       format(rlist[i]),
                                                       fullout=True,
                                                       kpc_scale=kpc_scale,
                                                       zcorr=zcorrlist[i],
                                                       err=err,
                                                       marker=symblist[i],
                                                       linestyle='',
                                                       color=colorlist[i])
        
        
        f.write('\n# P{} '.format(plist[i])+92*'#'+'\n# r ~ {} kpc\n'.\
                format(rlist[i]))
        for i in range(tmpz.size):
            f.write('{:6}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}\n'.\
                    format('',tmpz[i],tmpage[i],tmperr[i],tmpstd[i],
                           tmpAV[i],tmpAVerr[i],tmpAVstd[i],
                           tmpZ[i],tmpZerr[i],tmpZstd[i]))
        
        try:
            z = np.vstack((z,tmpz))
            age = np.vstack((age,tmpage))
            err = np.vstack((err,tmperr))
            std = np.vstack((std,tmpstd))
            AV = np.vstack((AV,tmpAV))
            AVstd = np.vstack((AVstd,tmpAVstd))
            AVerr = np.vstack((AVerr,tmpAVerr))
            Z = np.vstack((Z,tmpZ))
            Zerr = np.vstack((Zerr,tmpZerr))
            Zstd = np.vstack((Zstd,tmpZstd))
        except UnboundLocalError:
            z = tmpz
            age = tmpage
            err = tmperr
            std = tmpstd
            AV = tmpAV
            AVstd = tmpAVstd
            AVerr = tmpAVerr
            Z = tmpZ
            Zerr = tmpZerr
            Zstd = tmpZstd

    bigz = np.nanmean(z,axis=0)
    bigage = np.nanmean(age,axis=0)
    bigerr = np.sqrt(
        np.nansum(err*(age - bigage)**2,axis=0)/
        ((err.shape[0] - 1.)/(err.shape[0]) * np.nansum(err,axis=0)))
    bigAV = np.nanmean(AV,axis=0)
    bigAVerr = np.sqrt(
        np.nansum(AVerr*(AV - bigAV)**2,axis=0)/
        ((AVerr.shape[0] - 1.)/(AVerr.shape[0]) * np.nansum(AVerr,axis=0)))
    bigZ = np.nanmean(Z,axis=0)
    bigZerr = np.sqrt(
        np.nansum(Zerr*(Z - bigZ)**2,axis=0)/
        ((Zerr.shape[0] - 1.)/(Zerr.shape[0]) * np.nansum(Zerr,axis=0)))

    with open('{}_means.dat'.format(datname.split('.dat')[0]),'w') as fm:
        fm.write(str('#{:>9}'+'{:>10}'*6+'\n').format('height',
                                                      'MLWA',
                                                      'MLWA err',
                                                      'Av',
                                                      'Av err',
                                                      'Z',
                                                      'Z err'))

        for i in range(bigz.size):
              fm.write('{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}{:10.4f}\n'.\
                       format(bigz[i],bigage[i],bigerr[i],bigAV[i],bigAVerr[i],bigZ[i],bigZerr[i]))

    ageax.plot(bigz, bigage)
    ageax.fill_between(bigz,bigage-bigerr,bigage+bigerr,alpha=0.1)
    ageax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    ageax.set_xlim(-0.6,2.6)
    ageax.set_ylim(0,13)
    ageax.set_ylabel('Light-weighted age [Gyr]')
    ageax.set_title('Generated on {}'.format(time.asctime()))

    AVax.plot(bigz,bigAV)
    AVax.fill_between(bigz,bigAV-bigAVerr,bigAV+bigAVerr,alpha=0.1)
    AVax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    AVax.set_xlim(-0.6,2.6)
    AVax.set_ylim(0,6)
    AVax.set_ylabel(r'$A_V$')

    Zax.plot(bigz,bigZ)
    Zax.fill_between(bigz, bigZ-bigZerr, bigZ+bigZerr,alpha=0.1)
    Zax.legend(loc=1,title='radius [kpc]',scatterpoints=1,
               numpoints=1,frameon=False)
    Zax.set_xlim(-0.6,2.6)
#    Zax.set_ylim(-1.5,3.0)
    Zax.set_ylabel(r'Log( MLWZ [$Z_{\odot}$] )')

    pp.savefig(ageax.figure)
    pp.savefig(AVax.figure)
    pp.savefig(Zax.figure)
    plt.close(ageax.figure)
    plt.close(AVax.figure)
    plt.close(Zax.figure)

    pp.close()
    f.close()

    return

def simple_plot(inputsuffix='allz2.dat', label='Mean Light Weighted Age [Gyr]', 
                col=62, order=5, ylims=None, exclude=[[],[],[],[],[],[]]):

    zz = np.array([])
    dd = np.array([])

    axlist = []

    bigax = plt.figure().add_subplot(111)
    bigax.set_xlabel('|Height [kpc]|')
    bigax.set_ylabel(label)
    
    for i in range(6):
        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('|Height [kpc]|')
        ax.set_ylabel(label)
        ax.set_title('{}\nP{}'.format(time.asctime(),i+1))

        dat = glob('*P{}*{}'.format(i+1, inputsuffix))[0]
        print dat
        loc = glob('*P{}*locations.dat'.format(i+1))[0]
        print loc
    
        td = np.loadtxt(dat, usecols=(col,), unpack=True)
        r, tz = np.loadtxt(loc, usecols=(4,5), unpack=True)
        
        exarr = np.array(exclude[i])-1 #becuase aps are 1-indexed
        td = np.delete(td,exarr)
        t = np.delete(r,exarr)
        tz = np.delete(tz,exarr)

        z = np.array([])
        d = np.array([])
        e = np.array([])
        while tz.size > 0:
            zi = tz[0]
            idx = np.where(np.abs(tz - zi) < 0.05)
            d = np.r_[d,np.mean(td[idx])]
            e = np.r_[e,np.std(td[idx])]
            z = np.r_[z,np.abs(zi)]
            tz = np.delete(tz, idx)
            td = np.delete(td, idx)

        gidx = d == d
        d = d[gidx]
        z = z[gidx]
        e = e[gidx]
        
        sidx = np.argsort(z)
        mean = bn.move_mean(d[sidx],order)
        std = bn.move_std(d[sidx],order)
        spl = spi.UnivariateSpline(z[sidx],d[sidx])
        mean = spl(z[sidx])
        # mean = np.convolve(d[sidx],np.ones(order)/order,'same')
        # std = np.sqrt(np.convolve((d - mean)**2,np.ones(order)/order,'same'))

        l = bigax.plot(z[sidx],mean, label='P{}'.format(i+1))[0]
        bigax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=l.get_color())
        
        ax.plot(z[sidx],mean,color=l.get_color())
        ax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=l.get_color())

        ax.errorbar(z, d, yerr=e, fmt='.', color=l.get_color())
        ax.set_xlim(-0.1,2.6)
        
        if ylims is not None:
            ax.set_ylim(*ylims)
        
        axlist.append(ax)

    bigax.legend(loc=0, numpoints=1, scatterpoints=1)

    bigax.set_title(time.asctime())
    bigax.set_xlim(-0.1,2.6)

    if ylims is not None:
        bigax.set_ylim(*ylims)

    axlist.append(bigax)

    return axlist

def simple_batch(suffix, order=5, exclude=[[],[],[],[],[],[]]):

    clist = [62,61,66,63]
    llist = ['Mean Light Weighted Age [Gyr]',
             'Mean Mass Weighted Age [Gyr]',
             r'$\tau_V$',
             'Mean Light Weighted Metallicity [Z$_{\odot}$]']
    yllist = [[0,11],[0,11],[-1,6],[0,3]]
    sllist = ['MLWA','MMWA','TauV','MLWZ']

    for c, l, sl, yl in zip(clist, llist, sllist, yllist):
        pp = PDF('{}_{}.pdf'.format(sl,suffix))
        for x in simple_plot(col=c,label=l,ylims=yl,exclude=exclude,order=order):
            pp.savefig(x.figure)

        pp.close()

        # pp.savefig(simple_plot(col=62,label='Mean Light Weighted Age [Gyr]',
        #                        ylims=[0,11],order=order,exclude=exclude))
        # pp.savefig(simple_plot(col=61,label='Mean Mass Weighted Age [Gyr]',
        #                        ylims=[0,11],order=order,exclude=exclude))
        # pp.savefig(simple_plot(col=66,label=r'$\tau_V$',
        #                        ylims=[-1,6],order=order,exclude=exclude))
        # pp.savefig(simple_plot(col=63,label='Mean Light Weighted Metallicity [Z$_{\odot}$]',
        #                        ylims=[0,3],order=order,exclude=exclude))

    plt.close('all')

    return
