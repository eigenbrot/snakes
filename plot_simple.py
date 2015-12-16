import numpy as np
import GradPak_plot as GPP
import model_A as mA
import matplotlib.pyplot as plt
import pyfits
import pywcs
import time
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib import rc
from matplotlib import colors as mplcolors
from matplotlib.ticker import ScalarFormatter
import glob
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
    try:
        fibers = data[:,0]
    except IndexError:
        fibers = np.array([data[0]])
        data = np.array([data])

    form_age = 12.0 #Gyr
    logt = np.log10(np.r_[1e-99,AGES,form_age])
    tdiff = np.diff(logt)
    borders = 10**(logt[1:] - tdiff/2.)
    borders[0] = 1e-99
    borders[-1] = form_age

    pp = PDF(outputfile)
    for i in range(fibers.size):
        print i
        print np.log10(data[i,1:11])
        ax = plt.figure().add_subplot(111)
        ax.plot(AGES,np.log10(data[i,1:11]),'.g')
        ax.set_ylabel(r'Log($\int\psi (t)dt$ )')
        for b in borders:
            ax.axvline(b,alpha=0.2,ls='--')
        # ax.set_xlim(-1,AGES.size)
        # ax.set_xticks(np.arange(AGES.size))
        # ax.set_xticklabels(agelabels)
        ax.set_xlabel('Lookback time [Gyr]')
        ax.set_xlim(13,-1)
        ymin, ymax = ax.get_ylim()
#        ax.set_ylim(5,11)
        MMWA = data[i,11]
        ax.set_title('Fiber {}\nMMWA = {:4.3f} Gyr'.format(i+1,MMWA))
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

def all_maps(output,col=12,inputprefix='NGC_891',inputsuffix='fit.dat',labelfibers=False,
             label='Mean Light Weighted Age [Gyr]', log=False, plotbins=False,
             minval = None, maxval = None, exclude = None, binned=False, cmap='gnuplot2'):

    pp = PDF(output)
    centers = [
        [35.603413,42.323494],
        [35.615154,42.3426],
        [35.636325,42.37935],
        [35.625617,42.35935],
        [35.591642,42.304247],
        [35.648196,42.401217]]

    ax = None
    if exclude is None:
        exclude = [[],[],[],[],[],[]]
    figlist = []
    for i in range(6):
        print i
        try:
            inputfile = glob('{}*P{}*{}'.format(inputprefix,i+1,inputsuffix))[0]
        except IndexError:
            print 'Could not find P{} data, skipping'.format(i+1)
            continue

        print inputfile

        if binned:
            fitsname = glob('{}*P{}*ms*fits'.format(inputprefix,i+1))[0]
            print fitsname
            binhead = pyfits.open(fitsname)[0].header
        else:
            binhead = None
        
        data = np.loadtxt(inputfile,usecols=(col,),unpack=True)
        if col==12:
            label = 'Mean Light Weighted Age [Gyr]'
            minval = 0
            maxval = 13#np.nanmax(data)
        elif col==11:
            label = 'Mean Mass Weighted Age [Gyr]'
            minval = 0#np.nanmin(data)#7
            maxval = 10#np.nanmax(data)#10
        elif col == 13:
            label = '$A_V$'
            data *= 1.086
            minval = 0
            maxval = 6
        elif col==19:
            data = np.log10(data)
            label = r'Log( Metallicity [$Z_{\odot}$] )'
            # cmap = mplcolors.ListedColormap(['black','blue','green','yellow','red','white'])
            # norm = mplcolors.BoundaryNorm([0,0.005,0.02,0.2,0.4,1,2.5], cmap.N)
            minval = -4
            maxval = 1
        elif col==65:
            data *= 1.086
        
        if log:
            data = np.log10(data)
            if i == 0:
                label = 'Log( {} )'.format(label)
        print 'excluding', exclude[i]
        ax = GPP.plot(data,
                      ax=ax,
                      binheader=binhead,
                      plotbins=plotbins,
                      figsize=(8,4),
                      fitsfile=\
                      '/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits',
                      wcsax=False,
                      imrot=67,#64.213,
                      pa=295.787,
                      center = centers[i],
                      clabel=label,
                      cmap=cmap,
                      labelfibers=labelfibers,
                      exclude=exclude[i],
                      sky=False,
                      minval=minval,
                      maxval=maxval)

        #Do it again, for the single-pointing figure
        sax = GPP.plot(data,
                       ax=None,
                       binheader=binhead,
                       plotbins=plotbins,
                       figsize=(8,4),
                       fitsfile=\
                       '/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits',
                       wcsax=False,
                       imrot=67,#64.213,
                       pa=295.787,
                       center = centers[i],
                       clabel=label,
                       cmap=cmap,
                       labelfibers=labelfibers,
                       exclude=exclude[i],
                       sky=False,
                       minval=minval,
                       maxval=maxval)
        sax.text(0.5, 0.8, 'P{}'.format(i+1), color='r', fontsize=20,
                 transform=sax.transAxes, ha='center', va='center')
        figlist.append(sax.figure)

    ax.set_xlim(160,755)
    ax.set_ylim(350,600)
    header = pyfits.open('/d/monk/eigenbrot/WIYN/14B-0456/NGC_891.fits')[0].header
    hw = pywcs.WCS(header)
    centpx = hw.wcs_sky2pix([[35.63689,42.34633]],0)[0]
    xticks = np.arange(-300,301,100)/(header['CDELT1']*3600.) + centpx[0]
    yticks = np.arange(-100,201,100)/(header['CDELT2']*3600.) + centpx[1]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    xlabs = ['{:3.0f}'.format((centpx[0]-i)*header['CDELT1']*3600.) for i in ax.get_xticks()]
    ylabs = ['{:3.0f}'.format((i - centpx[1])*header['CDELT2']*3600.) for i in ax.get_yticks()]
    ax.set_xticklabels(xlabs)
    ax.set_yticklabels(ylabs)
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    pp.savefig(ax.figure)
    for f in figlist:
        pp.savefig(f)
    pp.close()
    plt.close('all')
    return ax

def all_heights(output,inputprefix='NGC_891',err=True,binned=False,reg=True):

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
    rlist = [-10.0, -6.4, -2.3, 0.8, 4.4, 8.1]
#    rlist = [4.4,0.8,-6.4,-2.3,8.1,-10.0]

    ax = None
    AVax = None
    metalax = None
    if reg:
        psi0ax = None
        tausfax = None
        invtausfax = None
        tformax = None
        massax = None
    for i in range(6):

        inputfile = glob('{}*P{}*fit.dat'.format(inputprefix,plist[i]))[0]
        print inputfile
        if binned:
            fitsname = glob('{}*P{}*ms*fits'.format(inputprefix,plist[i]))[0]
            print fitsname
            binhead = pyfits.open(fitsname)[0].header
        else:
            binhead = False
        if reg:
            MMWA, MLWA, TAUV, SNR, Z, PSI0, TAUSF, TFORM = np.loadtxt(inputfile,
                                                                      usecols=(11,12,13,14,17,18,19,20),
                                                                      unpack=True)

            mdata = np.loadtxt(inputfile)
            MASS = np.sum(mdata[:,1:11],axis=1)
        else:
            MMWA, MLWA, TAUV, SNR, Z = np.loadtxt(inputfile,
                                                  usecols=(11,12,13,14,19),
                                                  unpack=True)

        Z = np.log10(Z)
        ax, tmpz, tmpage, tmperr, tmpstd =  GPP.plot_rows(MLWA, 
                                                          binheader = binhead,
                                                          weights=SNR,
                                                          ax=ax,
                                                          label='{}'.\
                                                          format(rlist[i]),
                                                          fullout=True,
                                                          kpc_scale = kpc_scale,
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
                                                           err=err,
                                                           marker=symblist[i],
                                                           linestyle='',
                                                           color=colorlist[i])

        metalax, _, tmpmetal, tmpmetalerr, tmpmetalstd = GPP.plot_rows(Z, 
                                                                       binheader = binhead,
                                                                       weights=SNR,
                                                                       ax=metalax,
                                                                       label='{}'.\
                                                                       format(rlist[i]),
                                                                       fullout=True,
                                                                       kpc_scale=kpc_scale,
                                                                       err=err,
                                                                       marker=symblist[i],
                                                                       linestyle='',
                                                                       color=colorlist[i])

        if reg:
            psi0ax, _, tmppsi0, tmppsi0err, tmppsi0std = GPP.plot_rows(np.log10(PSI0), 
                                                                       binheader = binhead,
                                                                       weights=SNR,
                                                                       ax=psi0ax,
                                                                       label='{}'.\
                                                                       format(rlist[i]),
                                                                       fullout=True,
                                                                       kpc_scale=kpc_scale,
                                                                       err=err,
                                                                       marker=symblist[i],
                                                                       linestyle='',
                                                                       color=colorlist[i])
            
            tausfax, _, tmptausf, tmptausferr, tmptausfstd = GPP.plot_rows(TAUSF, 
                                                                           binheader = binhead,
                                                                           weights=SNR,
                                                                           ax=tausfax,
                                                                           label='{}'.\
                                                                           format(rlist[i]),
                                                                           fullout=True,
                                                                           kpc_scale=kpc_scale,
                                                                           err=err,
                                                                           marker=symblist[i],
                                                                           linestyle='',
                                                                           color=colorlist[i])
            
            invtausfax, _, tmpinvtausf, tmpinvtausferr, tmpinvtausfstd = GPP.plot_rows(1./TAUSF, 
                                                                                       binheader = binhead,
                                                                                       weights=SNR,
                                                                                       ax=invtausfax,
                                                                                       label='{}'.\
                                                                                       format(rlist[i]),
                                                                                       fullout=True,
                                                                                       kpc_scale=kpc_scale,
                                                                                       err=err,
                                                                                       marker=symblist[i],
                                                                                       linestyle='',
                                                                                       color=colorlist[i])
            
            tformax, _, tmptform, tmptformerr, tmptformstd = GPP.plot_rows(TFORM, 
                                                                           binheader = binhead,
                                                                           weights=SNR,
                                                                           ax=tformax,
                                                                           label='{}'.\
                                                                           format(rlist[i]),
                                                                           fullout=True,
                                                                           kpc_scale=kpc_scale,
                                                                           err=err,
                                                                           marker=symblist[i],
                                                                           linestyle='',
                                                                           color=colorlist[i])
            
            massax, _, tmpmass, tmpmasserr, tmpmassstd = GPP.plot_rows(np.log10(MASS), 
                                                                       binheader = binhead,
                                                                       weights=SNR,
                                                                       ax=massax,
                                                                       label='{}'.\
                                                                       format(rlist[i]),
                                                                       fullout=True,
                                                                       kpc_scale=kpc_scale,
                                                                       err=err,
                                                                       marker=symblist[i],
                                                                       linestyle='',
                                                                       color=colorlist[i])


        f.write('\n# P{} '.format(i+1)+92*'#'+'\n# r ~ {} kpc\n'.\
                format(rlist[i]))
        for i in range(tmpz.size):
            f.write('{:6}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}{:13.3f}\n'.\
                    format('',tmpz[i],tmpage[i],tmperr[i],tmpstd[i],
                           tmpAV[i],tmpAVerr[i],tmpAVstd[i],
                           tmpmetal[i],tmpmetalerr[i],tmpmetalstd[i]))

        try:
            z = np.vstack((z,tmpz))
            age = np.vstack((age,tmpage))
            err = np.vstack((err,tmperr))
            std = np.vstack((std,tmpstd))
            AV = np.vstack((AV,tmpAV))
            AVstd = np.vstack((AVstd,tmpAVstd))
            AVerr = np.vstack((AVerr,tmpAVerr))
            metal = np.vstack((metal,tmpmetal))
            metalerr = np.vstack((metalerr,tmpmetalerr))
            metalstd = np.vstack((metalstd,tmpmetalstd))
            if reg:
                psi0 = np.vstack((psi0,tmppsi0))
                psi0err = np.vstack((psi0err,tmppsi0err))
                psi0std = np.vstack((psi0std,tmppsi0std))
                tausf = np.vstack((tausf,tmptausf))
                tausferr = np.vstack((tausferr,tmptausferr))
                tausfstd = np.vstack((tausfstd,tmptausfstd))
                invtausf = np.vstack((invtausf,tmpinvtausf))
                invtausferr = np.vstack((invtausferr,tmpinvtausferr))
                invtausfstd = np.vstack((invtausfstd,tmpinvtausfstd))
                tform = np.vstack((tform,tmptform))
                tformerr = np.vstack((tformerr,tmptformerr))
                tformstd = np.vstack((tformstd,tmptformstd))
                mass = np.vstack((mass,tmpmass))
                masserr = np.vstack((masserr,tmpmasserr))
                massstd = np.vstack((massstd,tmpmassstd))
        except UnboundLocalError:
            z = tmpz
            age = tmpage
            err = tmperr
            std = tmpstd
            AV = tmpAV
            AVstd = tmpAVstd
            AVerr = tmpAVerr
            metal = tmpmetal
            metalerr = tmpmetalerr
            metalstd = tmpmetalstd
            if reg:
                psi0 = tmppsi0
                psi0err = tmppsi0err
                psi0std = tmppsi0std
                tausf = tmptausf
                tausferr = tmptausferr
                tausfstd = tmptausfstd            
                invtausf = tmpinvtausf
                invtausferr = tmpinvtausferr
                invtausfstd = tmpinvtausfstd            
                tform = tmptform
                tformerr = tmptformerr
                tformstd = tmptformstd
                mass = tmpmass
                masserr = tmpmasserr
                massstd = tmpmassstd

    bigz = np.nanmean(z,axis=0)
    bigage = np.nanmean(age,axis=0)
    bigerr = np.sqrt(
        np.nansum(err*(age - bigage)**2,axis=0)/
        ((err.shape[0] - 1.)/(err.shape[0]) * np.nansum(err,axis=0)))
    bigAV = np.nanmean(AV,axis=0)
    bigAVerr = np.sqrt(
        np.nansum(AVerr*(AV - bigAV)**2,axis=0)/
        ((AVerr.shape[0] - 1.)/(AVerr.shape[0]) * np.nansum(AVerr,axis=0)))
    bigmetal = np.nanmean(metal,axis=0)
    bigmetalerr = np.sqrt(
        np.nansum(metalerr*(metal - bigmetal)**2,axis=0)/
        ((metalerr.shape[0] - 1.)/(metalerr.shape[0]) * np.nansum(metalerr,axis=0)))
    
    if reg:
        bigpsi0 = np.nanmean(psi0,axis=0)
        bigpsi0err = np.sqrt(np.abs(
            np.nansum(psi0err*(psi0 - bigpsi0)**2,axis=0)/
            ((psi0err.shape[0] - 1.)/(psi0err.shape[0]) * np.nansum(psi0err,axis=0))))
        bigtausf = np.nanmean(tausf,axis=0)
        bigtausferr = np.sqrt(
            np.nansum(tausferr*(tausf - bigtausf)**2,axis=0)/
            ((tausferr.shape[0] - 1.)/(tausferr.shape[0]) * np.nansum(tausferr,axis=0)))
        biginvtausf = np.nanmean(invtausf,axis=0)
        biginvtausferr = np.sqrt(
            np.nansum(invtausferr*(invtausf - biginvtausf)**2,axis=0)/
            ((invtausferr.shape[0] - 1.)/(invtausferr.shape[0]) * np.nansum(invtausferr,axis=0)))
        bigtform = np.nanmean(tform,axis=0)
        bigtformerr = np.sqrt(
            np.nansum(tformerr*(tform - bigtform)**2,axis=0)/
            ((tformerr.shape[0] - 1.)/(tformerr.shape[0]) * np.nansum(tformerr,axis=0)))
        bigmass = np.nanmean(mass,axis=0)
        bigmasserr = np.sqrt(
            np.nansum(masserr*(mass - bigmass)**2,axis=0)/
            ((masserr.shape[0] - 1.)/(masserr.shape[0]) * np.nansum(masserr,axis=0)))

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
                       format(bigz[i],bigage[i],bigerr[i],bigAV[i],bigAVerr[i],bigmetal[i],bigmetalerr[i]))

    ax.plot(bigz, bigage)
    ax.fill_between(bigz,bigage-bigerr,bigage+bigerr,alpha=0.1)
    ax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    ax.set_xlim(-0.1,2.6)
    ax.set_ylim(-2,10)
    ax.set_ylabel('Light-weighted age [Gyr]')
    ax.set_title('Generated on {}'.format(time.asctime()))

    AVax.plot(bigz,bigAV)
    AVax.fill_between(bigz,bigAV-bigAVerr,bigAV+bigAVerr,alpha=0.1)
    AVax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    AVax.set_xlim(-0.1,2.6)
    AVax.set_ylim(0,6)
    AVax.set_ylabel(r'$A_V$')
    # plotz = np.linspace(0,2.5,20)
    # for t in range(len(rlist)):
    #     modelAv = mA.A_vec(np.abs(rlist[t]),plotz,zd=1.0,tau0=0.85,hd=7.68)
    #     AVax.plot(plotz,modelAv,color=colorlist[t])

    metalax.plot(bigz,bigmetal)
    metalax.fill_between(bigz, bigmetal-bigmetalerr, bigmetal+bigmetalerr,alpha=0.1)
    metalax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
    metalax.set_xlim(-0.1,2.6)
#    metalax.set_ylim(-1.5,3.0)
    metalax.set_ylabel(r'Log( $Z/Z_{\odot}$ )')

    if reg:
        psi0ax.plot(bigz,bigpsi0)
        psi0ax.fill_between(bigz, bigpsi0-bigpsi0err, bigpsi0+bigpsi0err,alpha=0.1)
        psi0ax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
        psi0ax.set_xlim(-0.1,2.6)
        psi0ax.set_ylabel(r'Log($\psi_0$)')
        
        tausfax.plot(bigz,bigtausf)
        tausfax.fill_between(bigz, bigtausf-bigtausferr, bigtausf+bigtausferr,alpha=0.1)
        tausfax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
        tausfax.set_ylim(-1,5)
        tausfax.set_xlim(-0.1,2.6)
        tausfax.set_ylabel(r'$\tau_{\mathrm{sf}}$')
        
        invtausfax.plot(bigz,biginvtausf)
        invtausfax.fill_between(bigz, biginvtausf-biginvtausferr, biginvtausf+biginvtausferr,alpha=0.1)
        invtausfax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
        invtausfax.set_ylim(-1,5)
        invtausfax.set_xlim(-0.1,2.6)
        invtausfax.set_ylabel(r'$1/\tau_{\mathrm{sf}}$')
        
        tformax.plot(bigz,bigtform)
        tformax.fill_between(bigz, bigtform-bigtformerr, bigtform+bigtformerr,alpha=0.1)
        tformax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
        tformax.set_xlim(-0.1,2.6)
        tformax.set_ylabel(r'$t_{\mathrm{form}}$')
        
        massax.plot(bigz,bigmass)
        massax.fill_between(bigz, bigmass-bigmasserr, bigmass+bigmasserr,alpha=0.1)
        massax.legend(loc=1,title='radius [kpc]',scatterpoints=1,numpoints=1,frameon=False)
        massax.set_xlim(-0.1,2.6)
        massax.set_ylabel(r'Log($M_{\mathrm{total}}$)')

    pp.savefig(ax.figure)
    pp.savefig(AVax.figure)
    pp.savefig(metalax.figure)
    plt.close(ax.figure)
    plt.close(AVax.figure)
    plt.close(metalax.figure)
    
    if reg:
        pp.savefig(psi0ax.figure)
        pp.savefig(tausfax.figure)
        pp.savefig(tformax.figure)
        pp.savefig(invtausfax.figure)
        pp.savefig(massax.figure)
        plt.close(psi0ax.figure)
        plt.close(tausfax.figure)
        plt.close(invtausfax.figure)
        plt.close(tformax.figure)
        plt.close(massax.figure)
    
    pp.close()
    f.close()
 
    return

def radial_gradient(input_file, output, col=1):

    f = open(input_file,'r')
    lines = f.readlines()

    radii = np.array([])
    values = np.array([])
    errs = np.array([])
    
    tmp = np.array([])
    for line in lines:
        
        if line == '\n':
            continue
        if line[0] == '#':
            if len(line) < 3:
                continue
            elif line[2] == 'r':
                print tmp
                values = np.append(values, np.mean(tmp))
                errs = np.append(errs, np.std(tmp))
                tmp = np.array([])
                radii = np.append(radii, float(line.split()[3]))
                continue
            else:
                continue

        cols = line.split()
        Av = float(cols[4])
        if Av <= 1.0:
            tmp = np.append(tmp, float(cols[col]))

    print tmp
    values = np.append(values, np.mean(tmp))
    errs = np.append(errs, np.std(tmp))
    f.close()
    
    values = values[1:]
    errs = errs[1:]
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('~r [kpc]')
    ax.set_ylabel('MLWA [Gyr]')
    ax.set_xlim(-11, 11)
    ax.errorbar(radii, values, yerr=errs, marker='.', ms=15, ls='')

    return radii, values, ax

def age_metal():

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])

    prefix = '/d/monk/eigenbrot/WIYN/14B-0456/bin/metal_dat/NGC_891_P2_bin30'

    age = np.zeros(fraclist.size)
    metal = np.zeros(fraclist.size)
    chi = np.zeros(fraclist.size)

    for i, z in enumerate(fraclist):
        
        name = '{}_Z{:04}.dat'.format(prefix,int(z*1000))        
        data = np.loadtxt(name)
        age[i] = data[5,12]
        chi[i] = data[5,15]
        metal[i] = z

    return age, metal, chi
