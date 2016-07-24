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
                          col=1, errcol=2, lowhigh=False, order=5, bigorder=60, 
                          s=None, ylims=None, labelr=False, bigpoints=False,
                          plotfit=True, exclude=[[],[],[],[],[],[]]):

    zz = np.array([])
    dd = np.array([])
    if lowhigh:
        ee = np.array([[],[]])
    else:
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
            if lowhigh:
                data, Lerr, Herr = np.loadtxt(dat, usecols=(col,errcol,errcol+1), unpack=True)
                err = np.vstack((Lerr,Herr))
            else:
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
        r = np.delete(r,exarr)
        z = np.delete(z,exarr)

        gidx = data == data
        data = data[gidx]
        z = z[gidx]
        if lowhigh:
            err = np.delete(err,exarr,axis=1)
            err = err[:,gidx]
            ee = np.hstack((ee,err))
        else:
            err = np.delete(err,exarr)
            err = err[gidx]
            ee = np.r_[ee,err]
        
        zz = np.r_[zz,z]
        dd = np.r_[dd,data]
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
                               col=6, errcol=None, lowhigh=False, 
                               order=5, ylims=None, bigpoints=False,
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
    
            if errcol == None:
                td = np.loadtxt(dat, usecols=(col[f],), unpack=True)
            else:
                if lowhigh:
                    td, low, high = np.loadtxt(dat, usecols=(col[f],errcol,errcol+1), unpack=True)
                    te = np.vstack((low,high))
                else:
                    td, te = np.loadtxt(dat, usecols=(col[f],errcol), unpack=True)                
            r, tz = np.loadtxt(loc, usecols=(4,5), unpack=True)

            exarr = np.array(exclude[pointing-1])-1 #becuase aps are 1-indexed
            td = np.delete(td,exarr)
            r = np.delete(r,exarr)
            tz = np.delete(tz,exarr)
            if errcol != None:
                if lowhigh:
                    te = np.delete(te,exarr,axis=1)
                else:
                    te = np.delete(te,exarr)

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
                if errcol == None:
                    e = np.zeros(tz.size)
                else:
                    e = te

            if combine_all:
                bigD = np.vstack((bigD,d))
                bigz = z

            gidx = d == d
            d = d[gidx]
            z = z[gidx]
            if lowhigh:
                e = e[:,gidx]
            else:
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
        
            # ax.plot(z[sidx],mean,color=color, ls=style, label=folder, alpha=alpha)
            # ax.fill_between(z[sidx],mean-std,mean+std, alpha=0.1, color=color)

            # print d.shape, np.sum(e,axis=0).shape
            # d = d/np.sum(e,axis=0)
            # e = np.diff(e,axis=0)[0]
            # print e.shape

            ax.errorbar(z, d, yerr=e, fmt='.', color=color,alpha=alpha,capsize=0, label=folder)

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
    color_list = ['blue','turquoise','seagreen','darkorange','crimson','dimgray','mediumorchid','maroon']
    style_list = ['-','-','-','-','-','-','-','-']

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

def coef_height_plot(field_name, output, err_name=None, ylim=[0,12],
                     exclude=exclude, suffix='coef', oplotdir=None):

        plist = [6,3,4,2,1,5]
        color_list = ['blue','seagreen','sienna','orange','yellowgreen','darkturquoise']
        # style_list = ['-','-','-','--','--','--']
        style_list = ['.'] * 6

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'$|z| \mathrm{[kpc]}$')
        ax.set_xlim(0,2.6)
        ax.set_ylim(*ylim)
        ax.set_ylabel(field_name)
        ax.set_title(time.asctime())
        
        for p, color, ex in zip(plist, color_list, exclude):
            
            loc = 'NGC_891_P{}_bin30_locations.dat'.format(p)
            coef = 'NGC_891_P{}_bin30_allz2.{}.fits'.format(p,suffix)
            z = np.loadtxt(loc,usecols=(5,),unpack=True)
            z = np.abs(z)
            
            data = pyfits.open(coef)[1].data
            values = data[field_name]

            exarr = np.array(ex) - 1
            z = np.delete(z, exarr)
            values = np.delete(values, exarr)

            if oplotdir is not None:
                oplotdata = pyfits.open('{}/{}'.format(oplotdir,coef))[1].data
                oplot = oplotdata[field_name]
                oplot = np.delete(oplot, exarr)

            if err_name is None:
                ax.plot(z, values, '.', color=color, label='P{}'.format(p))
            else:
                err = data[err_name]
                err = np.delete(err, exarr)
                ax.errorbar(z, values, yerr=err, fmt='.', color=color, capsize=0, label='P{}'.format(p))

            if oplotdir is not None:
                ax.scatter(z, oplot, marker='s', color=color, alpha=0.4)

        ax.legend(loc=0,scatterpoints=1,numpoints=1)

        pp = PDF(output)
        pp.savefig(fig)
        pp.close()

        return

def coef_covar(field1, field2, output, plotmono=False, field3=False):

    plist = [6,3,4,2,1,5]
    color_list = ['blue','seagreen','sienna','orange','yellowgreen','darkturquoise']
    style_list = ['.'] * 6
    
    for p, color in zip(plist, color_list):
        N = 1
        pp = PDF('{}_P{}.pdf'.format(output,p))
        
        if plotmono:
            monobase = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/no_weight'
            monod = {}
            for Z in ['0.2Z','0.4Z','1Z','2.5Z']:
                monod[Z] = pyfits.open('{}/{}/NGC_891_P{}_bin30_allz2.coef.fits'.\
                                       format(monobase,Z,p))[1].data

        while True:
            try:
                coef = 'MCdir/NGC_891_P{}_bin30_allz2.MC{:03n}.fits'.format(p,N)
                data = pyfits.open(coef)[1].data
            except IOError:
                break

            print coef
            values1 = data[field1]
            values2 = data[field2]
            if field3:
                c = data[field3]
            else:
                c = 'b'

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel(field1)
            ax.set_ylabel(field2)
            ax.set_title('P{}.{}\n{}'.format(p,N,time.asctime()))
            ax.scatter(values1, values2, c=c, alpha=1, linewidth=0, cmap=plt.cm.gnuplot)

            if plotmono:
                for Z in ['0.2Z','0.4Z','1Z','2.5Z']:
                    ax.scatter(monod[Z][field1][N-1],
                               monod[Z][field2][N-1], marker='s', color='r',
                               linewidths=0)

            pp.savefig(fig)
            plt.close(fig)
            N += 1

        pp.close()            
                
    return

def coef_MCcovar_3panel(field1, field2, field3, output,
                        label1 = False, label2 = False, label3 = False,
                        lim1 = [], lim2 = [], lim3 = [], colorbar=True,
                        basedir='.', zlims=[0,0.4], plist=[1,2,3,4,5,6]):

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)

    for p in plist:
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p)
        z = np.loadtxt(loc, usecols=(5,), unpack=True)
        z = np.abs(z)
        idx = np.where((z >= zlims[0]) & (z < zlims[1]))[0]
        N = 1
        while True:
            try:
                coef = '{}/MCdir/NGC_891_P{:}_bin30_allz2.MC{:03n}.fits'.format(basedir,p,N)
                c = pyfits.open(coef)[1].data
            except IOError:
                break

            if N-1 in idx:
                print coef, z[N-1]
                p1 = c[field1] - np.mean(c[field1])
                p2 = c[field2] - np.mean(c[field2])
                p3 = c[field3] - np.mean(c[field3])
                ax1.scatter(p1,p2,alpha=0.4,linewidth=0,c=[z[N-1]]*c[field2].size,vmin=0,vmax=2.6,cmap=plt.cm.gnuplot)
                ax2.scatter(p1,p3,alpha=0.4,linewidth=0,c=[z[N-1]]*c[field2].size,vmin=0,vmax=2.6,cmap=plt.cm.gnuplot)
                scat = ax3.scatter(p2,p3,alpha=0.4,linewidth=0,c=[z[N-1]]*c[field2].size,vmin=0,vmax=2.6,cmap=plt.cm.gnuplot)
            
            N += 1 

    ax1.set_xticklabels([])
    ax1.set_xlim(*lim1)
    ax1.set_ylim(*lim2)
    ax1.set_ylabel(label2 or '{0} - <{0}>'.format(field2))

    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*lim3)
    ax2.set_xlabel(label1 or '{0} - <{0}>'.format(field1))
    ax2.set_ylabel(label3 or '{0} - <{0}>'.format(field3))

    ax3.set_ylim(*ax2.get_ylim())
    ax3.set_xlim(*lim2)
    ax3.set_yticklabels([])
    ax3.set_xlabel(label2 or '{0} - <{0}>'.format(field2))
    
    fig.subplots_adjust(hspace=0.001,wspace=0.001)

    if colorbar:
        cb = fig.colorbar(scat, ax=[ax1,ax2,ax3])
        cb.set_alpha(1)
        cb.set_label('$|z| \mathrm{[kpc]}$')
        cb.draw_all()

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def coef_MCcovar_contour(field1, field2, field3, output,
                         label1 = False, label2 = False, label3 = False,
                         lim1 = [], lim2 = [], lim3 = [], bins=30,
                         basedir='.', zlims=[0,0.4], plist=[1,2,3,4,5,6]):

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)

    data1_list = []
    data2_list = []
    data3_list = []

    for p in plist:
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p)
        z = np.loadtxt(loc, usecols=(5,), unpack=True)
        z = np.abs(z)
        idx = np.where((z >= zlims[0]) & (z < zlims[1]))[0]
        N = 1
        while True:
            try:
                coef = '{}/MCdir/NGC_891_P{:}_bin30_allz2.MC{:03n}.fits'.format(basedir,p,N)
                c = pyfits.open(coef)[1].data
            except IOError:
                break

            if N-1 in idx:
                print coef, z[N-1]
                p1 = c[field1] - np.mean(c[field1])
                p2 = c[field2] - np.mean(c[field2])
                p3 = c[field3] - np.mean(c[field3])
                data1_list.append(p1)
                data2_list.append(p2)
                data3_list.append(p3)
            
            N += 1 

    data1 = np.hstack(data1_list)
    data2 = np.hstack(data2_list)
    data3 = np.hstack(data3_list)

    print data1.shape

    H1, xe1, ye1 = np.histogram2d(data1, data2, bins=bins)
    x1 = 0.5*(xe1[:-1] + xe1[1:])
    y1 = 0.5*(ye1[:-1] + ye1[1:])

    H2, xe2, ye2 = np.histogram2d(data1, data3, bins=bins)
    x2 = 0.5*(xe2[:-1] + xe2[1:])
    y2 = 0.5*(ye2[:-1] + ye2[1:])

    H3, xe3, ye3 = np.histogram2d(data2, data3, bins=bins)
    x3 = 0.5*(xe3[:-1] + xe3[1:])
    y3 = 0.5*(ye3[:-1] + ye3[1:])

    H1 /= np.max(H1)
    H2 /= np.max(H2)
    H3 /= np.max(H3)    

    levels = [0.005, 0.1, 0.4, 0.7, 0.9]

    ax1.contour(x1,y1,H1.T,levels,colors='k')
    ax2.contour(x2,y2,H2.T,levels,colors='k')
    ax3.contour(x3,y3,H3.T,levels,colors='k')

    ax1.set_xticklabels([])
    ax1.set_xlim(*lim1)
    ax1.set_ylim(*lim2)
    ax1.set_ylabel(label2 or '{0} - <{0}>'.format(field2))

    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*lim3)
    ax2.set_xlabel(label1 or '{0} - <{0}>'.format(field1))
    ax2.set_ylabel(label3 or '{0} - <{0}>'.format(field3))

    ax3.set_ylim(*ax2.get_ylim())
    ax3.set_xlim(*lim2)
    ax3.set_yticklabels([])
    ax3.set_xlabel(label2 or '{0} - <{0}>'.format(field2))
    
    fig.subplots_adjust(hspace=0.001,wspace=0.001)

    if not output:
        return fig
    
    else:
        pp = PDF(output)
        pp.savefig(fig)
        pp.close()
        plt.close(fig)

    return H1

def err_histogram(output, basedir='.',bins=10, field='MLWA', err='dMLWA', suffix='coef',
                  label=r'$\delta\tau_{L,\mathrm{fit}}/\tau_L$',exclude=exclude):

    ratio_list = []

    for p in range(6):
        coef = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,p+1,suffix)
        print coef
        c = pyfits.open(coef)[1].data
        tmp = c[err]
        if field == 'TAUV':
            tmp *= 1.086
        else:
            tmp /= c[field]            
        tmp = np.delete(tmp,np.array(exclude[p]) - 1)
        ratio_list.append(tmp)

    ratio = np.hstack(ratio_list)
    ratio = ratio[ratio == ratio]
    ratio = ratio[np.where(ratio < 0.8)[0]]
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(label)
    ax.set_ylabel('N')
    ax.hist(ratio, bins=bins, histtype='step', color='k')

    pp = PDF(output)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)

    return

def light_weight_cuts(output, basedir='.', exclude=exclude, 
                      rcuts=[3,8], zcuts=[0.4,1], numZ=6, numAge=4):

    fig = plt.figure()
    lax = fig.add_subplot(111, label='bigax')
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_ylabel(r'$Z/Z_\odot$')
    lax.tick_params(axis='both',pad=20,length=0)

    bigz = [0] + zcuts + [2.6]
    bigr = [0] + rcuts + [11]
    
    i = 1
    for z in range(len(zcuts) + 1):
        zc = [bigz[-z-2], bigz[-z-1]]
        for r in range(len(rcuts) + 1):
            rc = [bigr[r], bigr[r+1]]
            print zc, rc
            bigD = np.zeros((numZ,numAge))
            for p in range(6):
                coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
                print coeffile
                coefs = pyfits.open(coeffile)[1].data
                loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
                print loc
                r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
                r = np.abs(r)
                z = np.abs(z)
                
                exarr = np.array(exclude[p]) - 1
                r = np.delete(r, exarr)
                z = np.delete(z, exarr)
                lw = np.delete(coefs['LIGHT_WEIGHT'], exarr, axis=0)
                
                lw /= np.sum(lw,axis=1)[:,None]

                idx = np.where((z >= zc[0]) & (z < zc[1])
                               & (r >= rc[0]) & (r < rc[1]))[0]
                lw = lw[idx,:]
                print lw.shape
                lw = np.reshape(lw,(idx.size,numZ,numAge))
                print lw.shape
                print bigD.shape
                bigD += np.sum(lw,axis=0)

            ax = fig.add_subplot(len(zcuts)+1,len(rcuts)+1,i)
            ax.imshow(bigD,origin='lower',cmap='Blues',interpolation='none')
            ax.set_xticks(range(numAge))
            ax.set_yticks(range(numZ))
            ax.set_xticklabels([])
            ax.set_yticklabels([])

            if i > (len(zcuts)) * (len(rcuts) + 1):
                ax.set_xticklabels(['Y','I1','I2','O'])
            if i % (len(rcuts)+1) == 1:
                ax.set_yticklabels([0.005,0.02,0.2,0.4,1,2.5])

            if i % (len(rcuts) + 1) == 0:
                ax.text(1.15,0.5,'${}\leq |z| <{}$ kpc'.format(*zc),rotation=90,
                        ha='center',va='center',transform=ax.transAxes)

            if i <= len(rcuts) + 1:
                ax.text(0.5,1.1,'${}\leq |r| <{}$ kpc'.format(*rc),
                        ha='center',va='center',transform=ax.transAxes)

            i += 1

    fig.subplots_adjust(hspace=0.00001,wspace=0.0001)
    
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
            
    return
