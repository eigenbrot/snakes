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

def make_galaxies(tausf_list = [10,8,5,1]):

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    
    for z in range(len(modellist)):
        #get the length of the wavelength vector
        tmp = tm.make_galaxy('tmp',SSPs=modellist[z],makeplot=False,writeoutput=False)
        nwave = tmp['wave'].size
        output = np.zeros((len(tausf_list), nwave))
        outhdu = pyfits.PrimaryHDU()
        for i, t in enumerate(tausf_list):
            gal = tm.make_galaxy('tmp',tau_sf = t,
                                 SSPs=modellist[z],makeplot=False,
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
        outhdu.writeto('BC03_Z{:04n}_tau.fits'.format(fraclist[z]*1000),clobber=True)
        
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
        
        output = 'NGC_891_P{}_bin30.msoz.fits'.format(i+1)
        data_prep('NGC_891_P{}_bin30.mso.fits'.format(i+1),
                  'NGC_891_P{}_bin30_velocities.dat'.format(i+1),output)

    return

def run_sbands(findstr, bands, clobber=True):
    
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

def combine_sbands(output,numaps=10):

    numbands = 7
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)

    results = np.zeros((numaps, fraclist.size, numbands))
    for f, frac in enumerate(fraclist):
        
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

def quick_eat(datafile, numbands=7):

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

    FeAvg = 0.5*(index[5] + index[6])
    MgFe = np.sqrt(index[4]*(0.72*index[5] + 0.28*index[6]))
    
    return np.array([index[7], index[0], index[1], 
                     index[2], index[3], 
                     FeAvg, MgFe])

def plot_bc03_grid(bc03_data_file, ax, band1, band2):

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    
    tausf_list = np.array([10, 8, 5, 1])

    bc03data = pyfits.open(bc03_data_file)[0].data
    numtau, numZ, numindex = bc03data.shape


    for t in range(numtau):
        ax.plot(bc03data[t,:,band1],
                bc03data[t,:,band2],
                '-k')
        ax.text(bc03data[t,-1,band1],
                bc03data[t,-1,band2],
                '{:4.1f} Gyr'.format(tausf_list[t]),fontsize=6,ha='left')
    
    for z in range(numZ):
        ax.plot(bc03data[:,z,band1],
                bc03data[:,z,band2],
                ':k')
        ax.text(bc03data[-1,z,band1],
                bc03data[-1,z,band2],
                '{:4.2f} Z/Z$_{{\odot}}$'.format(fraclist[z]),fontsize=6,ha='center',va='top')

    return

def plot_yanny_on_grid(parfile, ax, band1, band2):

    par = yanny(parfile,np=True)
    
    scat = ax.scatter(par['APINFO'][band1], par['APINFO'][band2],
                      c=np.abs(par['APINFO']['z']), s=40, linewidths=0,
                      alpha=0.7, cmap=plt.cm.gnuplot2)
    
    return scat

def plot_quick_on_grid(datafile, ax, band1, band2):
    
    res = quick_eat(datafile)

    scat = ax.scatter(res[:,band1], res[:,band2], s=40, linewidths=0,
                      alpha=0.7)

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
        scat = plot_quick_on_grid(data_file, ax,
                                  5 + (p % 2),
                                  p/2)
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
    
    # cb = fig.colorbar(scat,ax=axes)
    # cb.set_label('z [kpc]')
    fig.suptitle('{}\nGenerated on {}'.format(data_file,time.asctime()))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)
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
