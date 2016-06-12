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

    # if not emcorr:
    #     bs = '_balmsub'
    # else:
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
               bands='/d/monk/eigenbrot/WIYN/14B-0456/anal/D4000.bands', 
               clobber=True):
    
    inputlist = glob(findstr)
    
    for data in inputlist:
        output = data.split('.fits')[0]+'.Dn4000.dat'
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

    numbands = 2
    if ma11:
        fraclist = ma11_fraclist
    else:
        fraclist = bc03_fraclist
    fraclist = np.sort(fraclist)

    results = np.zeros((numaps, fraclist.size, numbands))
    for f, frac in enumerate(fraclist):
        
        if ma11:
            fracfile = 'MA11_Z{:04n}_tau.Dn4000.dat'.format(frac*1000)
        else:
            fracfile = 'BC03_Z{:04n}_tau.Dn4000.dat'.format(frac*1000)
        data = np.loadtxt(fracfile,usecols=(0,1,5,6),
                          dtype={'names':('aps','bands','index','eqwidth'),
                                 'formats':('S50','S11','f4','f4')})

        aps = np.unique(data['aps'])
        sar = [int(s.split(',')[-1].split(']')[0]) for s in aps]
        sidx = np.argsort(sar)
        print aps[sidx]
        for a, ap in enumerate(aps[sidx]):
            idx = np.where(data['aps'] == ap)
            results[a,f,:] = eat_index(data[idx])

    
    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','aperture')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','HdA')
    outhdu.header.update('i1','Dn4000')
    outhdu.writeto(output,clobber=True)

    return

def quick_eat(datafile, numbands=2):

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
        results[a,:] = eat_index(data[idx])
    
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
    #1 Dn4000
    #
    # Output is:
    #0 HdA
    #1 Dn4000

    HdA = index['eqwidth'][0]
    Dn4 = index['index'][1]
    
    return np.array([HdA,Dn4])

def plot_model_grid(model_data_file, ax, band1, band2, ma11 = False):

    if ma11:
        fraclist = ma11_fraclist
    else:
        fraclist = bc03_fraclist
    fraclist = np.sort(fraclist)
    
    mlwa_list = np.array(deftmlwa)

    modeldata = pyfits.open(model_data_file)[0].data
    numtau, numZ, numindex = modeldata.shape


    for t in range(numtau):
        ax.plot(modeldata[t,:,band1],
                modeldata[t,:,band2],
                '-k')
        ax.text(modeldata[t,-1,band1],
                modeldata[t,-1,band2],
                '{:4.1f} Gyr'.format(mlwa_list[t]),fontsize=14,ha='left')

    for z in range(numZ):
        ax.plot(modeldata[:,z,band1],
                modeldata[:,z,band2],
                ':k')
        ax.text(modeldata[-1,z,band1],
                modeldata[-1,z,band2],
                '{:4.2f} Z/Z$_{{\odot}}$'.format(fraclist[z]),fontsize=6,ha='center',va='top')

    return

def plot_quick_on_grid(datafile, ax, band1, band2, exclude=[], basedir='.', plot_r = False, 
                       nocolor=False, zcut=[-99,99], rcut=[-99,99], size=40, marker='o', alpha=0.7):
    
    res = quick_eat(datafile)
    pointing = int(re.search('_P([1-6])_',datafile).groups()[0])
    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc,usecols=(4,5),unpack=True)
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
    if nocolor:
        d = 'k'

    scat = ax.scatter(res[:,band1], res[:,band2], s=size, linewidths=0,
                      marker=marker, vmin=-0.1, vmax=vmx,
                      c=d, alpha=alpha, cmap=plt.cm.gnuplot2)

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
    
    fig = plt.figure()
    
    bandlist = ['HdA','Dn4000']

    ab = ['<Fe>','<MgFe>']
    o = [r'$H\beta$',r'$H_{\delta,A}$','$H_{\gamma,A}$']

    ax = fig.add_subplot(111)
    plot_model_grid(model_data_file,ax,0,1)

    scat = plot_quick_on_grid(data_file, ax,
                              0,1,
                              exclude=exclude,
                              marker='o',size=40)
    ax.set_xlabel(r'$H_{\delta,A}$')
    ax.set_ylabel('Dn4000')
    
    cb = fig.colorbar(scat,ax=axes)
    cb.set_label('|z| [kpc]')
    fig.suptitle('{}\nGenerated on {}'.format(data_file,time.asctime()))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
    return 

def plot_D4000(data_file, output, exclude=[], zcut=[-99,99]):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    scat = plot_quick_on_grid(data_file, ax, 1, 0, exclude=exclude,
                              marker='o', size=40, zcut=zcut)

    ax.set_ylabel(r'$H_{\delta,A}$')
    ax.set_xlabel('Dn4000')
    cb = fig.colorbar(scat)
    cb.set_label('|z| [kpc]')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def plot_all_pointing_D4000(output, exclude=excl, r=False, zcut=[-99,99], rcut=[-99,99]):

    fig = plt.figure()
    ax = fig.add_subplot(111)
#    plot_model_grid('BC03_Dn4000.fits',ax,1,0)
    for p in range(6):
    
        data_file = 'NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(p+1)
        print data_file
        scat = plot_quick_on_grid(data_file, ax, 1, 0, exclude=exclude[p],
                                  marker='o', size=40, plot_r=r, zcut=zcut, rcut=rcut)

    ax.set_ylabel(r'$H_{\delta,A}$')
    ax.set_xlabel('Dn4000')
    ax.set_xlim(0.8,2.6)
    ax.set_ylim(-6,10)
    cb = fig.colorbar(scat)
    if r:
        cb.set_label('|r| [kpc]')
    else:
        cb.set_label('|z| [kpc]')
    fig.suptitle('{}\nGenerated on {}'.format(output,time.asctime()))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def plot_cuts_D4000(output, basedir='.', exclude=excl, zcuts=[0.4], rcuts=[3,8]):

    fig = plt.figure()
    lax = fig.add_subplot(111)
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_xlabel('Dn4000')
    lax.set_ylabel(r'$H_{\delta,A}$')
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
                data_file = '{}/NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(basedir,p+1)
                print data_file
                scat = plot_quick_on_grid(data_file, ax, 1, 0, exclude=exclude[p], nocolor=True,
                                          marker='o', size=40, plot_r=False, zcut=zc, rcut=rc, basedir=basedir)
            ax.text(2.5,8,'${}\leq |z| <{}$ kpc\n${}\leq |r| <{}$ kpc'.format(*(zc+rc)),ha='right',va='center')
            ax.set_ylim(-4,9.7)
            ax.set_xlim(0.82,2.66)
            if i < 4:
                ax.set_xticklabels([])
            if len(rcuts) > 0 and i % (len(rcuts)+1) != 1:
                ax.set_yticklabels([])
            i += 1

    fig.subplots_adjust(hspace=0.00001,wspace=0.0001)
    
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return fig
