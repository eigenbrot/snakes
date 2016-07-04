from glob import glob
import os
import re
import numpy as np
import pyfits
import tau_model as tm
from pyraf import iraf
import time
import prep_balmer as pb
import scipy.stats as ss
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()
iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

deftlst = [-1,-3,-9,13,5,3,1]
deftmlwa = [0.287,0.691,1.31,2.82,4.71,7.17,11.2]
excl = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]
ma11_fraclist = np.array([0.0001, 0.001, 0.01, 0.02, 0.04])/0.02
bc03_fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])

def make_galaxies(tausf_list = deftlst, ma11 = False):

    if ma11:
        fraclist = ma11_fraclist
        modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/MA11_models/ma11_cha_{}.fits'.format(i) for i in ['0005Z','005Z','05Z','1Z','2Z']]
    else:
        fraclist = bc03_fraclist
        modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    
    for z in range(len(modellist)):
        #get the length of the wavelength vector
        tmp = tm.make_galaxy('tmp',SSPs=modellist[z],makeplot=False,writeoutput=False)
        nwave = tmp['wave'].size
        output = np.zeros((len(tausf_list), nwave))
        outhdu = pyfits.PrimaryHDU()
        for i, t in enumerate(tausf_list):
            if ma11:
                galname = 'ma11_gal_t{}'.format(t)
            else:
                galname = 'bc03_gal_t{}'.format(t)
            gal = tm.make_galaxy(galname,tau_sf = t,
                                 SSPs=modellist[z],makeplot=True,
                                 writeoutput=False,vdisp=0.0)
            output[i,:] = gal['flux']/np.mean(gal['flux'])
            outhdu.header.update('TSF{:02n}'.format(i+1),t)

        outhdu.data = output
        outhdu.header.update('Z',fraclist[z])
        outhdu.header.update('CTYPE1','LINEAR')
        outhdu.header.update('CRPIX1',1)
        outhdu.header.update('CRVAL1',gal['wave'].min())
        outhdu.header.update('CDELT1',np.mean(np.diff(gal['wave'])))
        outhdu.header.update('CTYPE2','LINEAR')
        outhdu.header.update('CRPIX2',1)
        outhdu.header.update('CRVAL2',len(tausf_list))
        outhdu.header.update('CDELT2',1)
        if ma11:
            outhdu.writeto('MA11_Z{:04n}_tau.fits'.format(fraclist[z]*1000),clobber=True)
        else:
            outhdu.writeto('BC03_Z{:04n}_tau.fits'.format(fraclist[z]*1000),clobber=True)
        
    return

def data_prep(datafile, velocity, output, isdata=True, emcorr=False):
    
    wavemin=3800.
    wavemax=6800.

    vel = np.loadtxt(velocity,usecols=(1,),unpack=True)

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    header = hdu.header
    
    if isdata:
        cdelt = header['CDELT1']
        crpix = header['CRPIX1']
        crval = header['CRVAL1']

        if emcorr:
            #WARNING: The emission correction done below only corrects Ha and
            # HB, to get all the balmer lines check out
            # prep_balmer.do_all(*args,balmer=True)

            print "Correcting for Ha and HB emission"
            base = os.path.basename(datafile)
            dirn = os.path.dirname(datafile)
            pointing = int(re.search('_P([1-6])_',base).groups()[0])
            fitfile = os.path.join(dirn,'{}_allz2.fit.fits'.format(base.split('.')[0]))
            contsub = os.path.join(dirn,'{}_contsub.ms.fits'.format(base.split('.')[0]))
            pb.prep_spectra(datafile, fitfile, contsub, velocity)
            pb.do_fitprof(contsub, pointing)
            emline = pyfits.open(os.path.join(dirn,'P{}_HB_fits.fits'.format(pointing)))[0].data
    else:
        #Hacky hack for fit files. These values should not really change, so I think it's OK
        cdelt = 2.1
        crpix = 1
        crval = 3800.
        header.update('CDELT1',cdelt)
        header.update('CRPIX',crpix)

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval
    idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
    
    wave = wave[idx]
    data = data[:,idx]

    shift = np.vstack([np.interp(wave,wave*(1 - vel[i]/3e5),data[i,:]) for i in range(data.shape[0])])

    if isdata and emcorr:
#        emline = emline[:,idx]
        print shift.shape
        print emline.shape
        shift -= emline/1e17

    header.update('CRVAL1', 3800.)
    pyfits.PrimaryHDU(shift,header).writeto(output,clobber=True)

    return

def prep_all_data(emcorr=False):

    if not emcorr:
        bs = '_balmsub'
    else:
        bs = ''

    for i in range(6):
        
        output = 'NGC_891_P{}_bin30.msoz.fits'.format(i+1)
        data_prep('NGC_891_P{}_bin30{}.mso.fits'.format(i+1,bs),
                  'NGC_891_P{}_bin30_velocities.dat'.format(i+1),output,emcorr=emcorr)

    return

def prep_all_fits():
    
    for i in range(6):
        
        output = 'NGC_891_P{}_bin30_allz2.fitz.fits'.format(i+1)
        data_prep('NGC_891_P{}_bin30_allz2.fit.fits'.format(i+1),
                  'NGC_891_P{}_bin30_velocities.dat'.format(i+1),output,isdata=False)

    return

def run_sbands(findstr, 
               bands='/d/monk/eigenbrot/WIYN/14B-0456/anal/LICK.bands', 
               clobber=True):
    
    inputlist = glob(findstr)
    
    for data in inputlist:
        output = data.split('.fits')[0]+'.bands.dat'
        print '{} -> {}'.format(data,output)
        if clobber and os.path.exists(output):
            os.system('rm {}'.format(output))
        iraf.sbands(data,output,bands,
                    normali=True,
                    mag=False,
                    verbose=True,
                    _save=0)
        
    return

def combine_sbands(output,numaps=len(deftlst),ma11 = False):

    numbands = 8
    if ma11:
        fraclist = ma11_fraclist
    else:
        fraclist = bc03_fraclist
    fraclist = np.sort(fraclist)

    results = np.zeros((numaps, fraclist.size, numbands))
    for f, frac in enumerate(fraclist):
        
        if ma11:
            fracfile = 'MA11_Z{:04n}_tau.bands.dat'.format(frac*1000)
        else:
            fracfile = 'BC03_Z{:04n}_tau.bands.dat'.format(frac*1000)
        data = np.loadtxt(fracfile,usecols=(0,1,5,6),
                          dtype={'names':('aps','bands','index','eqwidth'),
                                 'formats':('S50','S11','f4','f4')})

        aps = np.unique(data['aps'])
        sar = [int(s.split(',')[-1].split(']')[0]) for s in aps]
        sidx = np.argsort(sar)
        print aps[sidx]
        for a, ap in enumerate(aps[sidx]):
            idx = np.where(data['aps'] == ap)
            results[a,f,:] = eat_index(data['eqwidth'][idx])

    
    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','aperture')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','HB')
    outhdu.header.update('i1','HdA')
    outhdu.header.update('i2','HgA')
    outhdu.header.update('i3','HdF')
    outhdu.header.update('i4','HgF')
    outhdu.header.update('i5','FeAvg')
    outhdu.header.update('i6','MgFe')
    outhdu.writeto(output,clobber=True)

    return

def quick_eat(datafile, numbands=8):

    data = np.loadtxt(datafile,usecols=(0,1,5,6),
                      dtype={'names':('aps','bands','index','eqwidth'),
                             'formats':('S80','S11','f4','f4')},
                      converters={5: eat_indef, 6: eat_indef})
    aps = np.unique(data['aps'])
    sar = [int(s.split(',')[-1].split(']')[0]) for s in aps]
    sidx = np.argsort(sar)
    numaps = aps.size
    results = np.zeros((numaps,numbands))
    for a, ap in enumerate(aps[sidx]):
        idx = np.where(data['aps'] == ap)
        results[a,:] = eat_index(data['eqwidth'][idx])
    
    return results

def eat_indef(s):

    try:
        res = np.float(s)
    except ValueError:
        res = np.nan

    return res

def bc03_compare(databands, bc03bands, outputprefix):

    pp = PDF('{}.pdf'.format(outputprefix))
    f = open('{}.dat'.format(outputprefix), 'w')
    f.write('#{:>7}{:>10}{:>10}\n'.format('Apnum','Age [Gyr]','Z/Z_sol'))

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    
    agelist = np.array([  5.01186000e+06,   2.51188000e+07,   1.01518000e+08,
                          2.86119008e+08,   6.40542976e+08,   9.04792000e+08,
                          1.43400000e+09,   2.50000000e+09,   5.00000000e+09,
                          1.00000000e+10])/1e9
    
    bc03data = pyfits.open(bc03bands)[0].data
    numages, numZ, numindex = bc03data.shape

    galdata = np.loadtxt(databands,usecols=(0,1,5,6),
                         dtype={'names':('aps','bands','index','eqwidth'),
                         'formats':('S50','S11','f4','f4')})

    aps = np.unique(galdata['aps'])
    sar = [int(s.split('(')[-1].split(')')[0]) for s in aps]
    sidx = np.argsort(sar)

    outdata = np.zeros((aps.size,2))

    for a, ap in enumerate(aps[sidx]):
        idx = np.where(galdata['aps'] == ap)
        bands = eat_index(galdata['eqwidth'][idx])
        compdata = np.tile(bands,(numages,numZ,1))
        chigrid = np.sum((compdata - bc03data)**2/bc03data,axis=2)
        bestage, bestZ = np.unravel_index(np.argmin(chigrid),chigrid.shape)
        outdata[a,0] = agelist[bestage]
        outdata[a,1] = fraclist[bestZ]
        
        f.write('{:8n}{:10.2f}{:10.0e}\n'.format(a+1,agelist[bestage],fraclist[bestZ]))

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        I = ax.imshow(chigrid,cmap=plt.cm.gnuplot,interpolation='nearest',
                      origin='lower')
        ax.set_title('Aperture {:}\nZ = {:5.2e} Z/Z$_{{\odot}}$\nAge = {:5.3f} Gyr'.format(a,outdata[a,1],outdata[a,0]))
        ax.set_yticks(np.arange(agelist.size))
        ax.set_yticklabels(['{:3.3f}'.format(a) for a in agelist])
        ax.set_xlabel('Z/Z$_{\odot}$')
        ax.set_xticks(np.arange(fraclist.size))
        ax.set_xticklabels(fraclist)
        ax.set_ylabel('SSP Age [Gyr]')

        ax.axvline(x=bestZ,color='white')
        ax.axhline(y=bestage,color='white')
        pp.savefig(fig)
        plt.close(fig)

    f.close()
    pp.close()
    return outdata

def eat_index(index):
    #Assumes order is:
    #0 HdA
    #1 HgA
    #2 HdF
    #3 HgF
    #4 Mgb
    #5 Fe5270
    #6 Fe5335
    #7 HB
    #
    # Output is:
    # 0 HB
    # 1 HdA
    # 2 HgA
    # 3 HdF
    # 4 HgF
    # 5 <Fe>
    # 6 <MgFe>
    # 7 Mgb

    FeAvg = 0.5*(index[5] + index[6])
    MgFe = np.sqrt(index[4]*(0.72*index[5] + 0.28*index[6]))
    
    return np.array([index[7], index[0], index[1], 
                     index[2], index[3], 
                     FeAvg, MgFe, index[4]])

def plot_model_grid(model_data_file, ax, band1, band2, alpha=1,
                    ma11 = False, plotlabels=True, isochrones=False):

    if ma11:
        fraclist = ma11_fraclist
    else:
        fraclist = bc03_fraclist
    fraclist = np.sort(fraclist)
    
    tausf_list = np.array(deftlst)
    mlwa_list = np.array(deftmlwa)

    modeldata = pyfits.open(model_data_file)[0].data
    numtau, numZ, numindex = modeldata.shape

    if isochrones:
        colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d']
        for t in range(numtau)[::-1]:
            ax.plot(modeldata[t,:,band1],
                    modeldata[t,:,band2],
                    '-',alpha=alpha,zorder=0,color=colors[t],lw=1.4)
            if ax.get_subplotspec().get_geometry()[2] == 2 and t == numtau - 1:
                ax.text(modeldata[0,0,band1],
                        modeldata[0,0,band2],
                        '{:6.3f} Z/Z$_{{\odot}}$'.format(fraclist[0]),fontsize=10,ha='left',va='center')
                ax.text(modeldata[-1,-1,band1],
                        modeldata[-1,-1,band2],
                        '{:4.1f} Z/Z$_{{\odot}}$'.format(fraclist[-1]),fontsize=10,ha='left',va='center')

    else:
        for t in range(numtau)[::-1]:
            ax.plot(modeldata[t,:,band1],
                    modeldata[t,:,band2],
                    '-k',alpha=alpha,zorder=0)
        
            if ax.get_subplotspec().get_geometry()[2] == 2 and plotlabels:
                if t < 2:
                    ax.text(modeldata[t,-1,band1],
                            modeldata[t,-1,band2]+0.25*(2-t),
                            '{:4.1f} Gyr'.format(mlwa_list[t]),fontsize=8,ha='left')
                else:
                    ax.text(modeldata[t,-1,band1],
                            modeldata[t,-1,band2],
                            '{:4.1f} Gyr'.format(mlwa_list[t]),fontsize=8,ha='left')

        for z in range(numZ):
            ax.plot(modeldata[:,z,band1],
                    modeldata[:,z,band2],
                    ':k',alpha=alpha,zorder=0)

            if ax.get_subplotspec().get_geometry()[2] == 2 and plotlabels:
                ax.text(modeldata[-1,z,band1],
                        modeldata[-1,z,band2],
                        '{:4.2f} Z/Z$_{{\odot}}$'.format(fraclist[z]),fontsize=8,ha='right',va='top')

    return

def plot_yanny_on_grid(parfile, ax, band1, band2):

    par = yanny(parfile,np=True)
    
    scat = ax.scatter(par['APINFO'][band1], par['APINFO'][band2],
                      c=np.abs(par['APINFO']['z']), s=40, linewidths=0,
                      alpha=0.7, cmap=plt.cm.gnuplot2)
    
    return scat

def plot_quick_on_grid(datafile, ax, band1, band2, exclude=[], basedir='.', plot_r=False, spy=False,
                       err=True, nocolor=False, zcut=[-99,99], rcut=[-99,99], size=40, marker='o', alpha=0.8):
    
    if spy:
        res = np.loadtxt(datafile)
    else:
        res = quick_eat(datafile)

    pointing = int(re.search('_P([1-6])_',datafile).groups()[0])
    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc,usecols=(4,5),unpack=True)
    fullr = r
    r = np.abs(r)
    z = np.abs(z)

    exar = np.array(exclude) - 1
    res = np.delete(res,exar,axis=0)
    z = np.delete(z,exar)
    r = np.delete(r,exar)
    if plot_r:
        d = np.abs(r)
        vmx = 10
    else:
        d = np.abs(z)
        vmx = 2.5

    idx = np.where((z >= zcut[0]) & (z < zcut[1])
                   & (r >= rcut[0]) & (r < rcut[1]))[0]
    res = res[idx,:]
    d = d[idx]
    fullr = fullr[idx]
    posidx = np.where(fullr >= 0)[0]
    negidx = np.where(fullr < 0)[0]
    
    if nocolor:
        d = 'k'

    scat = ax.scatter(res[posidx,band1], res[posidx,band2], s=size, linewidths=0,
                      marker=marker, vmin=-0.1, vmax=vmx, zorder=100,
                      c=d, alpha=alpha, cmap=plt.cm.gnuplot2)
    scat = ax.scatter(res[negidx,band1], res[negidx,band2], s=size, zorder=100, linewidths=1.2,
                      marker=marker, vmin=-0.1, vmax=vmx, facecolors='w',
                      alpha=alpha, cmap=plt.cm.gnuplot2)
    if spy and err:
        ax.errorbar(res[:,band1], res[:,band2], xerr=res[:,band1+1], yerr=res[:,band2+1],
                    fmt='none',capsize=0, ecolor='gray',elinewidth=1)

    return scat

def prep_contour_data(datafile, band1, band2, exclude=[], zcut=[-99,99]):

    res = quick_eat(datafile)
    pointing = int(re.search('_P([1-6])_',datafile).groups()[0])
    loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)
    z = np.loadtxt(loc,usecols=(5,),unpack=True)

    res = np.delete(res, exclude, axis=0)
    z = np.delete(z, exclude, axis=0)
    zidx = np.where((z >= zcut[0]) & (z < zcut[1]))[0]

    b1 = res[zidx,band1]
    b2 = res[zidx,band2]
    idx = np.where((b1 == b1) & (b2 == b2))
    b1 = b1[idx]
    b2 = b2[idx]

    return b1, b2

def plot_contour_on_grid(ax, data1, data2, color='k'):

    H, xbins, ybins = np.histogram2d(data1,data2)
    H /= np.max(H)
    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    ax.contour(H.T,[0.9,0.63,0.5,0.25],extent=extent,colors=color)
    #ax.contour(H.T,extent=extent,colors=color)

    return

def plot_index_grid(model_data_file,data_file,output,exclude=[],ma11=False,bestfits=None):
    
    fig = plt.figure(figsize=(12,11))
    
    bandlist = ['Hb','HdA','HgA','HdF','HgF','Fe','MgFe']

    ab = ['<Fe>','<MgFe>']
    o = [r'$H\beta$',r'$H_{\delta,A}$','$H_{\gamma,A}$']

    axes = []
    for p in range(6):
        ax = fig.add_subplot(3,2,p+1)
        axes.append(ax)
        plot_model_grid(model_data_file,ax,
                        5 + (p % 2),
                        p/2,ma11=ma11)
        if bestfits is not None:
            _ = plot_quick_on_grid(bestfits, ax,
                                   5 + (p % 2),
                                   p/2,
                                   exclude=exclude,
                                   marker='s',size=20)

        scat = plot_quick_on_grid(data_file, ax,
                                  5 + (p % 2),
                                  p/2,
                                  exclude=exclude,
                                  marker='o',size=40)
        ax.set_xlabel(ab[p%2])
        ax.set_ylabel(o[p/2])
        # ax.set_xlim(0.87,1.05)
        # ax.set_ylim(0.728,1.197)
        # ax.set_xlim(-0.060,0.021)
        # ax.set_ylim(-0.138,0.078)
        if p < 4:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        if p % 2 == 1:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    fig.subplots_adjust(hspace=0.0001,wspace=0.0001)
    
    cb = fig.colorbar(scat,ax=axes)
    cb.set_label('|z| [kpc]')
    fig.suptitle('{}\nGenerated on {}'.format(data_file,time.asctime()))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
    return 

def plot_all_pointing_grid(output, plotdata=True, plotfits=False, 
                           ma11=False, exclude=excl, contour=False,
                           zcut=[-99,99]):

    if ma11:
        model_file = 'MA11_data.fits'
        outpre = 'ma11'
    else:
        model_file = 'BC03_data.fits'
        outpre = 'bc03'

    fig = plt.figure(figsize=(12,11))
    
    bandlist = ['Hb','HdA','HgA','HdF','HgF','Fe','MgFe']

    ab = ['<Fe>','<MgFe>']
    o = [r'$H\beta$',r'$H_{\delta,A}$','$H_{\gamma,A}$']

    axes = []
    
    for p in range(6):
        ax = fig.add_subplot(3,2,p+1)
        axes.append(ax)
        plot_model_grid(model_file,ax,
                        5 + (p % 2),
                        p/2,ma11=ma11)
        if contour:
            data1 = np.array([])
            data2 = np.array([])
            fits1 = np.array([])
            fits2 = np.array([])

        for pointing in range(6):

            if plotfits:
                bestfits = 'NGC_891_P{}_bin30_allz2.fitz.bands.dat'.format(pointing+1)
                if contour:
                    f1, f2 = prep_contour_data(bestfits,
                                               5 + (p%2),
                                               p/2, zcut=zcut,
                                               exclude=exclude[pointing])
                    fits1 = np.r_[fits1,f1]
                    fits2 = np.r_[fits2,f2]

                else:
                    scat = plot_quick_on_grid(bestfits, ax,
                                              5 + (p % 2),
                                              p/2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='s',size=20,alpha=0.7)
            
            if plotdata:
                data_file = 'NGC_891_P{}_bin30.msoz.bands.dat'.format(pointing+1)
                if contour:
                    d1, d2 = prep_contour_data(data_file,
                                               5 + (p%2),
                                               p/2, zcut=zcut,
                                               exclude=exclude[pointing])
                    data1 = np.r_[data1, d1]
                    data2 = np.r_[data2, d2]
                    
                    scat = plot_quick_on_grid(data_file, ax,
                                              5 + (p % 2),
                                              p/2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='o',size=40,alpha=0.2)


                else:
                    scat = plot_quick_on_grid(data_file, ax,
                                              5 + (p % 2),
                                              p/2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='o',size=40,alpha=0.7)

        if contour:
            if plotdata:
                plot_contour_on_grid(ax,data1,data2,color='b')
            if plotfits:
                plot_contour_on_grid(ax,fits1,fits2,color='r')

        ax.set_xlabel(ab[p%2])
        ax.set_ylabel(o[p/2])
        # ax.set_xlim(0.87,1.05)
        # ax.set_ylim(0.728,1.197)
        # ax.set_xlim(-0.060,0.021)
        # ax.set_ylim(-0.138,0.078)
        if p < 4:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        if p % 2 == 1:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    fig.subplots_adjust(hspace=0.0001,wspace=0.0001)
    fig.suptitle('{}\nGenerated on {}'.format(output,time.asctime()))
    
    if not contour:
        cb = fig.colorbar(scat,ax=axes)
        cb.set_label('|z| [kpc]')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
    return 

def plot_cuts(output, x='Mgb', y='Fe', basedir='.', exclude=excl, zcuts=[0.4], rcuts=[3,8], 
              spy=False, err=True, grid=False, line=False, plotbreak=True, plotlabels=True,
              isochrones=False):

    band_d = {'Hb': {'label': r'$H\beta$', 'num': 0, 'lim': [-10,5.4]},
              'HdA': {'label': r'$H\delta_A$', 'num': 1, 'lim': [-4.3,8.4], 'spynum': 2}, #break = 2
              'HgA': {'label': r'$H\gamma_A$', 'num': 2, 'lim': [-8,8.4]},
              'HdF': {'label': r'$H\delta_F$', 'num': 3, 'lim': [-2,7.4]},
              'HgF': {'label': r'$H\gamma_F$', 'num': 4, 'lim': [-5,5.4]},
              'Fe': {'label': r'<Fe>', 'num': 5, 'lim': [0,3.4], 'break': 1.6, 'spynum':6},
              'MgFe': {'label': r'[MgFe]', 'num': 6, 'lim': [-0.5,4.8], 'ticks': [0,1,2,3,4],'spynum':8}, #break = 2
              'Mgb': {'label': r'Mg$b$', 'num': 7, 'lim': [0,5.4], 'break': 2.4, 'spynum': 4}}

    fig = plt.figure()
    lax = fig.add_subplot(111, label='biglabel') #label is necessary if len(zcuts) == len(rcuts) == 0
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_xlabel(band_d[x]['label'])
    lax.set_ylabel(band_d[y]['label'])
    lax.tick_params(axis='both',pad=20,length=0)

    bigz = [0] + zcuts + [2.6]
    bigr = [0] + rcuts + [11]
    
    i = 1
    for z in range(len(zcuts) + 1):
        zc = [bigz[-z-2], bigz[-z-1]]
        for r in range(len(rcuts) + 1):
            rc = [bigr[r], bigr[r+1]]
            print zc, rc
            ax = fig.add_subplot(len(zcuts)+1,len(rcuts)+1,i)
            for p in range(6):
                if spy:
                    data_file = '{}/NGC_891_P{}_bin30.msoz.spy.dat'.format(basedir,p+1)
                    band1 = band_d[x]['spynum']
                    band2 = band_d[y]['spynum']
                else:
                    data_file = '{}/NGC_891_P{}_bin30.msoz.bands.dat'.format(basedir,p+1)
                    band1 = band_d[x]['num']
                    band2 = band_d[y]['num']

                print data_file
                scat = plot_quick_on_grid(data_file, ax, band1, band2, spy=spy, 
                                          exclude=exclude[p], nocolor=True, err=err,
                                          marker='o', size=40, plot_r=False, 
                                          zcut=zc, rcut=rc, basedir=basedir)
            ax.set_ylim(*band_d[y]['lim'])
            ax.set_xlim(*band_d[x]['lim'])
            
            if plotbreak:
                try:
                    ax.axvline(band_d[x]['break'],alpha=0.6,ls='--',color='k')
                except KeyError:
                    pass
                try:
                    ax.axhline(band_d[y]['break'],alpha=0.6,ls='--',color='k')
                except KeyError:
                    pass
            if line:
                ax.plot([-0.2,8.8],[0,6],'--',color='k', alpha=0.6)

            try:
                ax.set_xticks(band_d[x]['ticks'])
            except KeyError:
                pass

            if i % (len(rcuts) + 1) == 0:
                tax = ax.twinx()
                tax.set_ylabel('${}\leq |z| <{}$ kpc'.format(*zc))
                #rotation='horizontal',labelpad=20)
                tax.set_ylim(*ax.get_ylim())
                tax.set_yticklabels([])

            if i <= len(rcuts) + 1:
                tax = ax.twiny()
                tax.set_xlabel('${}\leq |r| <{}$ kpc'.format(*rc))
                #rotation='horizontal',labelpad=20)
                tax.set_xlim(*ax.get_xlim())
                tax.set_xticklabels([]) 
            
            if i <= (len(zcuts)) * (len(rcuts) + 1):
                ax.set_xticklabels([])
            if len(rcuts) > 0 and i % (len(rcuts)+1) != 1:
                ax.set_yticklabels([])
            if grid:
                if spy:
                    model_file = '{}/BC03_spy.fits'.format(basedir)
                    band1 = band_d[x]['spynum']/2
                    band2 = band_d[y]['spynum']/2
                else:
                    model_file = '{}/BC03_bands.fits'.format(basedir)
                    band1 = band_d[x]['num']
                    band2 = band_d[y]['num']
                plot_model_grid(model_file,ax,band1,band2,alpha=0.5,
                                plotlabels=plotlabels, isochrones=isochrones)
            if spy and not err and i == (len(zcuts)+1) * (len(rcuts)+1) - len(rcuts):
                res = np.loadtxt(data_file)
                repxerr = np.nanmedian(res[:,band_d[x]['spynum']+1])
                repyerr = np.nanmedian(res[:,band_d[y]['spynum']+1])
                xpos = band_d[x]['lim'][1] - repxerr*2
                ypos = band_d[y]['lim'][1] - repyerr*2
                ax.errorbar(xpos,ypos,xerr=repxerr,yerr=repyerr,fmt='none',capsize=0,ecolor='k')

            i += 1

    fig.subplots_adjust(hspace=0.00001,wspace=0.0001)
    
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return fig

def plot_z_index(output, basedir='.', exclude=excl):
    
    fig = plt.figure()
    Feax = fig.add_subplot(311)
    Feax.set_ylabel(r'<Fe>')
    Feax.set_xticklabels([])
    MgFeax = fig.add_subplot(312)
    MgFeax.set_ylabel('[MgFe]')
    MgFeax.set_xticklabels([])
    Mgbax = fig.add_subplot(313)
    Mgbax.set_xlabel('|$z$| [kpc]')
    Mgbax.set_ylabel('Mg$b$')

    for p in range(6):
        datafile = glob('{}/NGC_891_P{}_bin30.*bands.dat'.format(basedir,p+1))[0]
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
        res = quick_eat(datafile)
        r, z = np.loadtxt(loc,usecols=(4,5),unpack=True)
        r = np.abs(r)
        z = np.abs(z)
        
        exar = np.array(exclude[p]) - 1
        res = np.delete(res,exar,axis=0)
        z = np.delete(z,exar)
        r = np.delete(r,exar)
        
        Feax.plot(z, res[:,5], '.', color='k')
        MgFeax.plot(z, res[:,6], '.', color='k')
        Mgbax.plot(z, res[:,7], '.', color='k')

    Feax.set_ylim(-0.2,3)
    MgFeax.set_ylim(-0.2,3.7)
    Mgbax.set_ylim(0.5,5.2)
    Feax.set_xlim(*Mgbax.get_xlim())
    MgFeax.set_xlim(*Mgbax.get_xlim())

    fig.subplots_adjust(hspace=0.00001)
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return

def plot_all_pointing_smallgrid(output, plotdata=True, plotfits=False, 
                                ma11=False, exclude=excl, contour=False,
                                zcut=[-99,99]):

    if ma11:
        model_file = 'MA11_data.fits'
        outpre = 'ma11'
    else:
        model_file = 'BC03_data.fits'
        outpre = 'bc03'

    fig = plt.figure(figsize=(9,4))
    
    bandlist = ['Hb','HdA','HgA','HdF','HgF','Fe','MgFe','Mgb']

    ab = ['<MgFe>','Mgb']
    o = [r'$H_{\delta,A}$']

    # ab = ['<MgFe>','Mgb']
    # o = [r'$H_{\delta,A}$','<Fe>']

    axes = []
    
    for p in range(2):
        ax = fig.add_subplot(1,2,p+1)
        axes.append(ax)

        bid1 = 6 + 1*(p % 2)
        bid2 = 1
        
        plot_model_grid(model_file,ax,
                        bid1,
                        bid2, ma11=ma11)

        print bandlist[bid1], bandlist[bid2]
        
        if contour:
            data1 = np.array([])
            data2 = np.array([])
            fits1 = np.array([])
            fits2 = np.array([])

        for pointing in range(6):

            if plotfits:
                bestfits = 'NGC_891_P{}_bin30_allz2.fitz.bands.dat'.format(pointing+1)
                if contour:
                    f1, f2 = prep_contour_data(bestfits,
                                               bid1,
                                               bid2, zcut=zcut,
                                               exclude=exclude[pointing])
                    fits1 = np.r_[fits1,f1]
                    fits2 = np.r_[fits2,f2]

                else:
                    scat = plot_quick_on_grid(bestfits, ax,
                                              bid1,
                                              bid2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='s',size=20,alpha=0.7)
            
            if plotdata:
                data_file = 'NGC_891_P{}_bin30.msoz.bands.dat'.format(pointing+1)
                if contour:
                    d1, d2 = prep_contour_data(data_file,
                                               bid1,
                                               bid2, zcut=zcut,
                                               exclude=exclude[pointing])
                    data1 = np.r_[data1, d1]
                    data2 = np.r_[data2, d2]
                    
                    scat = plot_quick_on_grid(data_file, ax,
                                              bid1,
                                              bid2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='o',size=40,alpha=0.2)


                else:
                    scat = plot_quick_on_grid(data_file, ax,
                                              bid1,
                                              bid2, zcut=zcut,
                                              exclude=exclude[pointing],
                                              marker='o',size=40,alpha=0.7)

        if contour:
            if plotdata:
                plot_contour_on_grid(ax,data1,data2,color='b')
            if plotfits:
                plot_contour_on_grid(ax,fits1,fits2,color='r')

        ax.set_xlabel(ab[p])
        ax.set_ylabel(o[0])
        # ax.set_xlim(0.87,1.05)
        # ax.set_ylim(0.728,1.197)
        # ax.set_xlim(-0.060,0.021)
        # ax.set_ylim(-0.138,0.078)
        # if p == 1:
        #     ax.set_xlim(1,4)
        #     ax.set_ylim(0.5,3)
        # if p < 4:
        #     ax.set_xticklabels([])
        #     ax.set_xlabel('')
        if p % 2 == 1:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    fig.subplots_adjust(hspace=0.0001,wspace=0.0001)
    fig.suptitle('{}\nGenerated on {}'.format(output,time.asctime()))
    
    if not contour:
        cb = fig.colorbar(scat,ax=axes)
        cb.set_label('|z| [kpc]')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
    return 

def plot_all_grids(usefits=False, ma11=False):

    if ma11:
        model_file = 'MA11_data.fits'
        outpre = 'ma11'
    else:
        model_file = 'BC03_data.fits'
        outpre = 'bc03'

    for i in range(6):
        print i+1
        if usefits:
            bf = 'NGC_891_P{}_bin30_allz2.fitz.bands.dat'.format(i+1)
        else:
            bf = None
            
        plot_index_grid(model_file,
                        'NGC_891_P{}_bin30.msoz.bands.dat'.format(i+1),
                        '{}_P{}_balmsum_tau_grid.pdf'.format(outpre,i+1),
                        exclude=excl[i],
                        ma11=ma11,
                        bestfits=bf)

    return

def do_all420(outname, ma11=False):

    if ma11:
        model_file = 'MA11_data.fits'
        outpre = 'ma11'
    else:
        model_file = 'BC03_data.fits'
        outpre = 'bc03'

    make_galaxies(ma11=ma11)
    prep_all_data()
    prep_all_fits()
    run_sbands('*z.fits')
    run_sbands('*tau.fits')
    combine_sbands(model_file,ma11=ma11)

    plot_all_pointing_grid('{}_{}_taugrid_data.pdf'.format(outpre,outname),
                           ma11=ma11,plotdata=True,plotfits=False,
                           contour=False,exclude=excl)
    plot_all_pointing_grid('{}_{}_taugrid_fits.pdf'.format(outpre,outname),
                           ma11=ma11,plotdata=False,plotfits=True,
                           contour=False,exclude=excl)
    plot_all_pointing_grid('{}_{}_taugrid_test.pdf'.format(outpre,outname),
                           ma11=ma11,plotdata=True,plotfits=True,
                           contour=True,exclude=excl)
    plt.close('all')
    return

def model_compare(dataloc, modellist, output, prep=True, sband=True):

    pp = PDF(output)
    
    bandlist = ['Hb','HdA','HgA','HdF','HgF','Fe','MgFe']
    plotbands = [0,  1,     2,                5,   6]
    
    for p in range(6):
        
        fig = plt.figure()
        lax = fig.add_subplot(111)
        lax.set_axis_off()
        lax.set_ylabel('Data')
        lax.set_title('P{}\n{}'.format(p+1,time.asctime()),fontsize=8)
        
        #First, get the actual data
        datafile = glob('{}/NGC_891_P{}*.ms.fits'.format(dataloc,p+1))[0]
        velfile = glob('{}/NGC_891_P{}*velocities.dat'.format(dataloc,p+1))[0]
        preppedfile = '{}/NGC_891_P{}_bin30.ms.bp.fits'.format(dataloc,p+1)
    
        print 'Grabbing', datafile

        if prep:
            data_prep(datafile, velfile, preppedfile)
        
        if sband:
            try:
                run_sbands(preppedfile, 'LICK.bands'.format(dataloc))
            except Exception as e:
                print e
                run_sbands(preppedfile, 'LICK.bands'.format(dataloc))

        bandfile = '{}/NGC_891_P{}_bin30.ms.bp.bands.dat'.format(dataloc,p+1)
        data = quick_eat(bandfile)

        #Now the different models
        linelist = []
        leglist = []
        axlist = []
        for axn, band in enumerate(plotbands):
            ax = fig.add_subplot(3,2,axn+2)
            ax.text(0.1,0.85,bandlist[band],transform=ax.transAxes)
            ax.tick_params(axis='both',labelsize=8)
            axlist.append(ax)

        for mnum, mloc in enumerate(modellist):
            lax.set_xlabel(mloc)
            
            print '{}/NGC_891_P{}*.fit.fits'.format(mloc,p+1)
            modelfile = glob('{}/NGC_891_P{}*.fit.fits'.format(mloc,p+1))[0]
            Mvelfile = glob('{}/NGC_891_P{}*velocities.dat'.format(mloc,p+1))[0]
            Mpreppedfile = '{}/NGC_891_P{}_bin30_allz2.fit.bp.fits'.format(mloc,p+1)

            print 'Grabbing', modelfile
            
            if prep:
                data_prep(modelfile, Mvelfile, Mpreppedfile, isdata=False)

            if sband:
                try:
                    run_sbands(Mpreppedfile, 'LICK.bands'.format(dataloc))
                except Exception as e:
                    #Fuck you, IRAF
                    print e
                    run_sbands(Mpreppedfile, 'LICK.bands'.format(dataloc))
            
            Mbandfile = '{}/NGC_891_P{}_bin30_allz2.fit.bp.bands.dat'.format(mloc,p+1)
            model = quick_eat(Mbandfile)
            
            print data.shape
            print model.shape

            for axn, band in enumerate(plotbands):
                l = axlist[axn].plot(data[:,band],model[:,band],'.')[0]
                tau, pval = ss.kendalltau(data[:,band], model[:,band])
                CI = (1-ss.norm.cdf(ss.norm.isf(tau)))*100
                axlist[axn].text(0.9,0.9-mnum*0.05,'{:4.2f}'.format(tau),color=l.get_color(),
                                 transform=axlist[axn].transAxes,fontsize=5)
                # if axn != 0 and axn %2 == 0:
                #     ax.set_yticklabels([])
                # if axn < 3:
                #     ax.set_xticklabels([])

            linelist.append(l)
            leglist.append(os.path.basename(mloc))

        for ax in axlist:
            ax.autoscale(False)
            ax.plot([-100,100],[-100,100],':k',alpha=0.3)
                        

        kax = fig.add_axes([0.2,0.7,0.2,0.2])
        kax.set_xlabel('Data')
        kax.set_ylabel('Best Fit Model')
        kax.set_xticklabels([])
        kax.set_yticklabels([])
        kax.legend(linelist, leglist, loc='center',fontsize=8, numpoints=1,frameon=False)
        # fig.subplots_adjust(hspace=0.0001, wspace=0.0001)
        fig.tight_layout(h_pad=0.2,w_pad=1.2)
        pp.savefig(fig)
        plt.close(fig)

    pp.close()
