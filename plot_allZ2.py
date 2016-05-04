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

exclude = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]
e2 = [[5, 7, 9, 30, 34], [1, 2, 3, 11, 32, 33, 34, 35, 37], 
      [1, 3, 4, 50, 53, 54, 55, 59], [1, 2, 15], [1, 2, 3, 5, 24, 27, 28, 29,5], [35, 36, 37, 38]]

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
                col=62, order=5, ylims=None, labelr=False, bigpoints=False,
                exclude=[[],[],[],[],[],[]]):

    zz = np.array([])
    dd = np.array([])

    axlist = []

    bigax = plt.figure().add_subplot(111)
    bigax.set_xlabel('|Height [kpc]|')
    bigax.set_ylabel(label)
    
    plist = [6,3,4,2,1,5]
    #color_list = ['blue','turquoise','chartreuse','yellow','tomato','red']
    color_list = ['blue','seagreen','sienna','sienna','seagreen','blue']
    style_list = ['-','-','-','--','--','--']

    for i in range(6):
        pointing = plist[i]
        color = color_list[i]
        style = style_list[i]

        dat = glob('*P{}*{}'.format(pointing, inputsuffix))[0]
        print dat
        loc = glob('*P{}*locations.dat'.format(pointing))[0]
        print loc
        print 'Excluding: ', exclude[pointing-1]
    
        td = np.loadtxt(dat, usecols=(col,), unpack=True)
        r, tz = np.loadtxt(loc, usecols=(4,5), unpack=True)
        avgr = np.mean(r)

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('|Height [kpc]|')
        ax.set_ylabel(label)
        if labelr:
            ax.set_title('{:4.0f} kpc'.format(avgr))
            linelabel = '{:4.0f} kpc'.format(avgr)
        else:
            ax.set_title('{}\nP{}'.format(time.asctime(),pointing))
            linelabel = 'P{}'.format(pointing)

        exarr = np.array(exclude[pointing-1])-1 #becuase aps are 1-indexed
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
        dp = np.r_[d[sidx][order::-1],d[sidx]]
        zp = np.r_[z[sidx][order::-1],z[sidx]]
        mean = bn.move_mean(dp,order)[order+1:]
        std = bn.move_std(dp,order)[order+1:]
        spl = spi.UnivariateSpline(z[sidx],d[sidx])
        mean = spl(z[sidx])
        # mean = np.convolve(d[sidx],np.ones(order)/order,'same')
        # std = np.sqrt(np.convolve((d - mean)**2,np.ones(order)/order,'same'))

        bigax.plot(z[sidx],mean, label=linelabel, color=color, ls=style)
        bigax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=color)
        if bigpoints:
            bigax.errorbar(z, d, yerr=e, fmt='.', color=color, alpha=0.6, capsize=0)
        
        ax.plot(z[sidx],mean,color=color, ls=style)
        ax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=color)

        ax.errorbar(z, d, yerr=e, fmt='.', color=color)
        ax.set_xlim(-0.1,2.6)
        
        if ylims is not None:
            ax.set_ylim(*ylims)
        
        axlist.append(ax)

    bigax.legend(loc=0, numpoints=1, scatterpoints=1)

    bigax.set_title(time.asctime())
    bigax.set_xlim(-0.1,2.6)

    if ylims is not None:
        bigax.set_ylim(*ylims)

    axlist = [bigax] + axlist

    return axlist

def plot_heights_with_err(inputsuffix,label=r'$\tau_{\mathrm{V,Balm}}$',
                          col=1, errcol=2, order=5, bigorder=60, s=None,
                          ylims=None, labelr=False, bigpoints=False,
                          plotfit=True, exclude=[[],[],[],[],[],[]]):

    zz = np.array([])
    dd = np.array([])
    ee = np.array([])
    axlist = []

    bigax = plt.figure().add_subplot(111)
    bigax.set_xlabel('|Height [kpc]|')
    bigax.set_ylabel(label)
    
    plist = [6,3,4,2,1,5]
    color_list = ['blue','seagreen','sienna','orange','yellowgreen','darkturquoise']
    style_list = ['-','-','-','--','--','--']

    for i in range(6):
        pointing = plist[i]
        color = color_list[i]
        style = style_list[i]

        dat = glob('*P{}*{}'.format(pointing, inputsuffix))[0]
        print dat
        loc = glob('*P{}*locations.dat'.format(pointing))[0]
        print loc
        print 'Excluding: ', exclude[pointing-1]
    
        if errcol is not None:
            data, err = np.loadtxt(dat, usecols=(col,errcol), unpack=True)
        else:
            data = np.loadtxt(dat, usecols=(col,), unpack=True)
            err = np.ones(data.size)*0.01

        r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
        avgr = np.mean(r)

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('|Height [kpc]|')
        ax.set_ylabel(label)
        if labelr:
            ax.set_title('{:4.0f} kpc'.format(avgr))
            linelabel = '{:4.0f} kpc'.format(avgr)
        else:
            ax.set_title('{}\nP{}'.format(time.asctime(),pointing))
            linelabel = 'P{}'.format(pointing)

        exarr = np.array(exclude[pointing-1])-1 #becuase aps are 1-indexed
        data = np.delete(data,exarr)
        err = np.delete(err,exarr)
        r = np.delete(r,exarr)
        z = np.delete(z,exarr)

        gidx = data == data
        data = data[gidx]
        z = z[gidx]
        err = err[gidx]
        
        zz = np.r_[zz,z]
        dd = np.r_[dd,data]
        ee = np.r_[ee,err]

        sidx = np.argsort(z)
        data_pad = np.r_[data[sidx][order::-1],data[sidx]]
        z_pad = np.r_[z[sidx][order::-1],z[sidx]]
        # mean = bn.move_mean(data_pad,order)[order+1:]
        std = bn.move_std(data_pad,order)[order+1:]
        spl = spi.UnivariateSpline(z[sidx],data[sidx])
        mean = spl(z[sidx])
        # mean = np.convolve(d[sidx],np.ones(order)/order,'same')
        # std = np.sqrt(np.convolve((d - mean)**2,np.ones(order)/order,'same'))

        bigax.errorbar(z, data, yerr=err, fmt='.', label=linelabel, color=color, capsize=0)
        
        # ax.plot(z[sidx],mean,color=color, ls=style)
        # ax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=color)

        ax.errorbar(z, data, yerr=err, fmt='.', color=color, capsize=0)
        ax.set_xlim(-0.1,2.6)
        
        if ylims is not None:
            ax.set_ylim(*ylims)
        
        axlist.append(ax)
        
    plot_title = time.asctime()
    if plotfit:
        sidx = np.argsort(zz)
        big_data_pad = np.r_[dd[sidx][bigorder::-1],dd[sidx]]
        big_z_pad = np.r_[zz[sidx][bigorder::-1],zz[sidx]]
        big_e_pad = np.r_[ee[sidx][bigorder::-1],ee[sidx]]
        big_sum = bn.move_sum(big_data_pad/big_e_pad,bigorder)[bigorder+1:]
        big_weight = bn.move_sum(1./big_e_pad,bigorder)[bigorder+1:]
        big_mean = big_sum/big_weight

        # std = bn.move_std(data_pad,order)[order+1:]
        # big_spl = spi.UnivariateSpline(zz[sidx],dd[sidx],w = 1./ee[sidx]**2, k=k, s=s)
        # big_mean = big_spl(zz[sidx])
        # big_pc = np.polyfit(zz[sidx], dd[sidx], polydeg, w=1./ee[sidx]**2)
        # big_poly = np.poly1d(big_pc)
        # big_mean = big_poly(zz[sidx])
        
        p = np.poly1d(np.polyfit(zz[sidx],big_mean,1))
        print p.coeffs
        
        bigax.plot(zz[sidx],big_mean,'-k',lw=2)
        bigax.plot(zz[sidx],p(zz[sidx]),':k',lw=1.5)
        plot_title += '\n'+label+'$={:4.2f}z{:+4.2f}$'.format(p.coeffs[0],p.coeffs[1])

    bigax.set_title(plot_title)
    bigax.legend(loc=0, numpoints=1, scatterpoints=1)
    bigax.set_xlim(-0.1,2.6)

    print zz.size

    if ylims is not None:
        bigax.set_ylim(*ylims)

    axlist = [bigax] + axlist
    
    return axlist

def height_plot_across_folders(folder_list, inputsuffix='allz2.dat', 
                               label='Mean Light Weighted Age [Gyr]', 
                               col=6, order=5, ylims=None, bigpoints=False,
                               binz=True, combine_all=False, plot_std=False,
                               exclude=[[],[],[],[],[],[]]):

    axlist = []
    
    plist = [6,3,4,2,1,5]
    #color_list = ['blue','turquoise','chartreuse','yellow','tomato','red']
    color_list = ['blue','seagreen','darkorange','crimson','dimgray','mediumorchid','lightblue']
    style_list = ['-','-','-','-','-','-','-']

    if not isinstance(col,list):
        col = [col] * len(folder_list)

    for i in range(6):                
        pointing = plist[i]

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('|Height [kpc]|')
        ax.set_ylabel(label)
        ax.set_title('{}\nP{}'.format(time.asctime(),pointing))

        for f, folder in enumerate(folder_list):
            color = color_list[f]
            style = style_list[f]
            
            dat = glob('{}/*P{}*{}'.format(folder, pointing, inputsuffix))[0]
            print dat
            loc = glob('{}/*P{}*locations.dat'.format(folder, pointing))[0]
            print loc
            print 'Excluding: ', exclude[pointing-1]
    
            td = np.loadtxt(dat, usecols=(col[f],), unpack=True)
            r, tz = np.loadtxt(loc, usecols=(4,5), unpack=True)

            exarr = np.array(exclude[pointing-1])-1 #becuase aps are 1-indexed
            td = np.delete(td,exarr)
            r = np.delete(r,exarr)
            tz = np.delete(tz,exarr)

            alpha=1.0
            if combine_all and f == 0:
                bigD = np.zeros(td.size)
                alpha=0.3
            
            if binz:
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
            else:
                z = tz
                d = td
                e = np.zeros(tz.size)

            if combine_all:
                bigD = np.vstack((bigD,d))
                bigz = z

            gidx = d == d
            d = d[gidx]
            z = z[gidx]
            e = e[gidx]

            sidx = np.argsort(z)
            dp = np.r_[d[sidx][order::-1],d[sidx]]
            zp = np.r_[z[sidx][order::-1],z[sidx]]
            mean = bn.move_mean(dp,order)[order+1:]
            std = bn.move_std(dp,order)[order+1:]
            spl = spi.UnivariateSpline(z[sidx],d[sidx])
            mean = spl(z[sidx])
            # mean = np.convolve(d[sidx],np.ones(order)/order,'same')
            # std = np.sqrt(np.convolve((d - mean)**2,np.ones(order)/order,'same'))
        
            ax.plot(z[sidx],mean,color=color, ls=style, label=folder, alpha=alpha)
            ax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=color)

            ax.errorbar(z, d, yerr=e, fmt='.', color=color,alpha=alpha,capsize=0)

        ax.set_xlim(-0.1,2.6)
        
        if ylims is not None:
            ax.set_ylim(*ylims)
        ax.legend(loc=0,numpoints=1)

        if combine_all:
            sidx = np.argsort(bigz)
            bigD = bigD[1:]
            bigMean = bn.nanmean(bigD,axis=0)
            bigStd = bn.nanstd(bigD,axis=0)
            bigspl = spi.UnivariateSpline(bigz[sidx],bigMean[sidx])
            bigFit = bigspl(bigz[sidx])
            
            ax.plot(bigz[sidx], bigFit, 'k-', lw=2)
            ax.errorbar(bigz, bigMean, yerr=bigStd, fmt='.', color='k',capsize=0)

        axlist.append(ax)
    
        if combine_all and plot_std:
            ax2 = plt.figure().add_subplot(111)
            ax2.set_xlabel('|Height [kpc]|')
            ax2.set_ylabel('$\delta$'+label)
            ax2.set_title(ax.get_title())
            ax2.plot(bigz, bigStd, 'k')
            axlist.append(ax2)

    return axlist

def histogram_across_folders(folder_list, inputsuffix='allz2.dat', 
                             label='Mean Light Weighted Age [Gyr]', 
                             col=62, bins=10, ylims=None, xlims=None,
                             delta=False,exclude=[[],[],[],[],[],[]]):
    
    axlist = []
    
    plist = [6,3,4,2,1,5]
    #color_list = ['blue','turquoise','chartreuse','yellow','tomato','red']
    color_list = ['blue','seagreen','darkorange','crimson','dimgray','mediumorchid','maroon']
    style_list = ['-','-','-','-','-','-','-']

    if not isinstance(col,list):
        col = [col] * len(folder_list)
    
    bigD = dict(zip(folder_list,[np.array([]) for i in range(len(folder_list))]))

    for i in range(6):                
        pointing = plist[i]

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel(label)
        ax.set_ylabel('N')
        ax.set_title('{}\nP{}'.format(time.asctime(),pointing))

        for f, folder in enumerate(folder_list):
            color = color_list[f]
            style = style_list[f]

            dat = glob('{}/*P{}*{}'.format(folder, pointing, inputsuffix))[0]
            print dat
            print 'Excluding: ', exclude[pointing-1]
    
            d = np.loadtxt(dat, usecols=(col[f],), unpack=True)
            
            exarr = np.array(exclude[pointing-1])-1 #becuase aps are 1-indexed
            d = np.delete(d,exarr)
                            
            gidx = d == d
            d = d[gidx]

            if delta:
                d -= np.min(d)

            ax.hist(d, bins=bins, histtype='stepfilled', color=color, label=folder,alpha=0.3)
            
            bigD[folder] = np.r_[bigD[folder], d]

        if xlims is not None:
            ax.set_xlim(*xlims)
        
        if ylims is not None:
            ax.set_ylim(*ylims)
        ax.legend(loc=0,numpoints=1)
        axlist.append(ax)

    bigax = plt.figure().add_subplot(111)
    bigax.set_xlabel(label)
    bigax.set_ylabel('N')
    bigax.set_title(time.asctime())
        
    for f, folder in enumerate(folder_list):
        color = color_list[f]
        style = style_list[f]
        
        bigax.hist(bigD[folder], bins=bins, histtype='stepfilled', color=color, label=folder,alpha=0.3)
        
    bigax.legend(loc=0,numpoints=1)
    if xlims is not None:
        bigax.set_xlim(*xlims)
        
    if ylims is not None:
        bigax.set_ylim(*ylims)
        
    axlist = [bigax] + axlist

    return axlist


def simple_batch(prefix, order=5, exclude=[[],[],[],[],[],[]]):

    clist = [62,61,66,63]
    llist = ['Mean Light Weighted Age [Gyr]',
             'Mean Mass Weighted Age [Gyr]',
             r'$\tau_V$',
             'Mean Light Weighted Metallicity [Z$_{\odot}$]']
    yllist = [[0,11],[0,11],[-1,6],[0,3]]
    sllist = ['MLWA','MMWA','TauV','MLWZ']

    for c, l, sl, yl in zip(clist, llist, sllist, yllist):
        pp = PDF('{}_{}_heights.pdf'.format(prefix,sl))
        for x in simple_plot(col=c,label=l,ylims=yl,exclude=exclude,order=order):
            pp.savefig(x.figure)

        pp.close()

    plt.close('all')

    return

def dfk_batch(prefix, order=5, exclude=[[],[],[],[],[],[]], 
              offset=0, labelr=False, bigpoints=False):

    clist = [6,5,10,7,8]
    llist = ['Mean Light Weighted Age [Gyr]',
             'Mean Mass Weighted Age [Gyr]',
             r'$\tau_V$',
             'Mean Mass Weighted Metallicity [Z$_{\odot}$]',
             'Mean Light Weighted Metallicity [Z$_{\odot}$]']
    yllist = [[0,11],[0,11],[-1,6],[0,1.3],[0,1.3]]
    sllist = ['MLWA','MMWA','TauV','MMWZ','MLWZ']

    for c, l, sl, yl in zip(clist, llist, sllist, yllist):
        pp = PDF('{}_{}_heights.pdf'.format(prefix,sl))
        for x in simple_plot(col=c+offset,label=l,ylims=yl,
                             exclude=exclude,order=order,
                             labelr=labelr,bigpoints=bigpoints):
            pp.savefig(x.figure)

        pp.close()

    plt.close('all')

    return
