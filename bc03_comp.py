from glob import glob
import os
import re
import numpy as np
import pyfits
from yanny import yanny
from pyraf import iraf
import time
import prep_balmer as pb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()
iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

def bc03_prep():

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]

    for i, model in enumerate(modellist):
        h = pyfits.open(model)[1]
        linwave = np.arange(h.data['WAVE'].min(),h.data['WAVE'].max(),2.1)
        numage = h.data['FLUX'].shape[1]
        out = np.zeros((numage,linwave.size))
        for j in range(numage):
            out[j,:] = np.interp(linwave,h.data['WAVE'][0],
                                 h.data['FLUX'][0,j])
        
        outhdu = pyfits.PrimaryHDU(out)
        outhdu.header.update('Z',fraclist[i])
        outhdu.header.update('CTYPE1','LINEAR')
        outhdu.header.update('CRPIX1',1)
        outhdu.header.update('CRVAL1',linwave.min())
        outhdu.header.update('CDELT1',np.mean(np.diff(linwave)))
        outhdu.header.update('CTYPE2','LINEAR')
        outhdu.header.update('CRPIX2',1)
        outhdu.header.update('CRVAL2',numage)
        outhdu.header.update('CDELT2',1)
        outhdu.writeto('BC03_Z{:04n}.fits'.format(fraclist[i]*1000),clobber=True)

    return

def data_prep(datafile, velocity, output, isdata=True):
    
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
        pointing = int(re.search('_P([1-6])_',datafile).groups()[0])
        fitfile = '{}_allz2.fit.fits'.format(datafile.split('.')[0])
        contsub = '{}_contsub.ms.fits'.format(datafile.split('.')[0])
        pb.prep_spectra(datafile, fitfile, contsub, velocity)
        pb.do_fitprof(contsub, pointing)
        emline = pyfits.open('P{}_HB_fits.fits'.format(pointing))[0].data
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

    if isdata:
#        emline = emline[:,idx]
        print shift.shape
        print emline.shape
        shift -= emline/1e17

    header.update('CRVAL1', 3800.)
    pyfits.PrimaryHDU(shift,header).writeto(output,clobber=True)

    return

def prep_all_data():
    
    for i in range(6):
        
        output = 'NGC_891_P{}_bin30.msz.fits'.format(i+1)
        data_prep('NGC_891_P{}_bin30.ms.fits'.format(i+1),
                  'NGC_891_P{}_bin30_velocities.dat'.format(i+1),output)

    return

def run_sbands(findstr, bands):
    
    inputlist = glob(findstr)
    
    for data in inputlist:
        output = data.split('.fits')[0]+'.bands.dat'
        print '{} -> {}'.format(data,output)
        iraf.sbands(data,output,bands,
                    normali=True,
                    mag=False,
                    verbose=True)
        
    return

def combine_sbands(output,numaps=10):

    numbands = 7
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)

    results = np.zeros((numaps, fraclist.size, numbands))
    for f, frac in enumerate(fraclist):
        
        fracfile = 'BC03_Z{:04n}.bands.dat'.format(frac*1000)
        data = np.loadtxt(fracfile,usecols=(0,1,5,6),
                          dtype={'names':('aps','bands','index','eqwidth'),
                                 'formats':('S26','S11','f4','f4')})

        aps = np.unique(data['aps'])
        sar = [int(s.split(',')[-1].split(']')[0]) for s in aps]
        sidx = np.argsort(sar)
        print aps[sidx]
        for a, ap in enumerate(aps[sidx]):
            idx = np.where(data['aps'] == ap)
            results[a,f,:] = eat_index(data['index'][idx])

    
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
        bands = eat_index(galdata['index'][idx])
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

    FeAvg = 0.5*(index[5] + index[6])
    MgFe = np.sqrt(index[4]*(0.72*index[5] + 0.28*index[6]))
    
    return np.array([index[7], index[0], index[1], 
                     index[2], index[3], 
                     FeAvg, MgFe])

def plot_bc03_grid(bc03_data_file, ax, band1, band2):

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    
    agelist = np.array([  5.01186000e+06,   2.51188000e+07,   1.01518000e+08,
                          2.86119008e+08,   6.40542976e+08,   9.04792000e+08,
                          1.43400000e+09,   2.50000000e+09,   5.00000000e+09,
                          1.00000000e+10])/1e9

    bc03data = pyfits.open(bc03_data_file)[0].data
    numages, numZ, numindex = bc03data.shape

    idx = np.where(agelist >= 0.3)[0]

    for a in range(numages - 3):
        # ax.plot(np.log10(bc03data[a+3,:,band1]),
        #         np.log10(bc03data[a+3,:,band2]),
        #         '-k')
        # ax.text(np.log10(bc03data[a+3,-1,band1]),
        #         np.log10(bc03data[a+3,-1,band2]),
        #         '{:4.1f} Gyr'.format(agelist[a+3]),fontsize=6,ha='right')
        ax.plot(bc03data[a+3,:,band1],
                bc03data[a+3,:,band2],
                '-k')
        ax.text(bc03data[a+3,-1,band1],
                bc03data[a+3,-1,band2],
                '{:4.1f} Gyr'.format(agelist[a+3]),fontsize=8,ha='right')
    
    for z in range(numZ):
        # ax.plot(np.log10(bc03data[idx,z,band1]),
        #         np.log10(bc03data[idx,z,band2]),
        #         ':k')
        # ax.text(np.log10(bc03data[idx[-1],z,band1]),
        #         np.log10(bc03data[idx[-1],z,band2]),
        #         '{:4.1f} Z/Z$_{{\odot}}$'.format(fraclist[z]),fontsize=6)
        ax.plot(bc03data[idx,z,band1],
                bc03data[idx,z,band2],
                ':k')
        ax.text(bc03data[idx[-1],z,band1],
                bc03data[idx[-1],z,band2],
                '{:4.2f} Z/Z$_{{\odot}}$'.format(fraclist[z]),fontsize=8)

    return

def plot_yanny_on_grid(parfile, ax, band1, band2):

    par = yanny(parfile,np=True)
    
    scat = ax.scatter(par['APINFO'][band1], par['APINFO'][band2],
                      c=np.abs(par['APINFO']['z']), s=40, linewidths=0,
                      alpha=0.7, cmap=plt.cm.gnuplot2)
    
    return scat

def plot_index_grid(bc03_data_file,data_file,output):
    
    fig = plt.figure(figsize=(12,11))
    
    bandlist = ['Hb','HdA','HgA','HdF','HgF','Fe','MgFe']

    ab = ['<Fe>','<MgFe>']
    o = [r'$H\beta$',r'$H_{\delta,A}$','$H_{\gamma,A}$']

    axes = []
    for p in range(6):
        ax = fig.add_subplot(3,2,p+1)
        axes.append(ax)
        plot_bc03_grid(bc03_data_file,ax,
                       5 + (p % 2),
                       p/2)
        scat = plot_yanny_on_grid(data_file,ax,
                                  bandlist[5 + (p % 2)],
                                  bandlist[p/2])
        ax.set_xlabel(ab[p%2])
        ax.set_ylabel(o[p/2])
        ax.set_xlim(0.87,1.05)
        ax.set_ylim(0.728,1.197)
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
    cb.set_label('z [kpc]')
    fig.suptitle('{}\nGenerated on {}'.format(data_file,time.asctime()))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
    return 
