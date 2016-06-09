import numpy as np
import re
import pyfits
import plot_simple as ps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

exclude = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]

def data_prep(datafile, velocity, output):
    
    wavemin=3800.
    wavemax=6000.

    vel = np.loadtxt(velocity,usecols=(1,),unpack=True)

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    header = hdu.header    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval
    idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
    
    wave = wave[idx]
    data = data[:,idx]

    shift = np.vstack([np.interp(wave,wave*(1 - vel[i]/3e5),data[i,:]) for i in range(data.shape[0])])

    header.update('CRVAL1', 3800.)
    pyfits.PrimaryHDU(shift,header).writeto(output,clobber=True)

    return

def prep_all_data():
    
    for i in range(6):
        
        output = 'NGC_891_P{}_bin30.msoz.fits'.format(i+1)
        data_prep('NGC_891_P{}_bin30.mso.fits'.format(i+1),
                  'NGC_891_P{}_bin30_velocities.dat'.format(i+1),output)

    return

def cut_pointing(pointing, zcut=0.5):
    """Simply split and combine all spectra into two bins: above and below 
    zcut
    """
    
    loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)
    log = 'NGC_891_P{}_bin30.log'.format(pointing)
    dfile = 'NGC_891_P{}_bin30.msoz.fits'.format(pointing)
    print loc, log, dfile
    z = np.loadtxt(loc, usecols=(5,), unpack=True)
    z = np.abs(z)

    ss = re.compile(r'(\d+): \[(.*)\], SNR: (\d+\.\d+)')    
    SN = np.zeros(z.size)
    with open(log,'r') as lff:
        lf = lff.readlines()
    for i, l in enumerate(lf):
        grp = ss.search(l).groups()
        SN[i] = float(grp[2])

    hdu = pyfits.open(dfile)[0]
    data = hdu.data
    header = hdu.header    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval

    exar = np.array(exclude[pointing - 1]) - 1
    z = np.delete(z,exar)
    SN = np.delete(SN,exar)
    data = np.delete(data,exar,axis=0)

    zhidx = np.where(z >= zcut)[0]
    zlidx = np.where(z < zcut)[0]

    low_d = np.sum(data[zlidx,:]*SN[zlidx,None], axis=0)/np.sum(SN[zlidx])
    high_d = np.sum(data[zhidx,:]*SN[zhidx,None], axis=0)/np.sum(SN[zhidx])

    return wave, low_d, high_d

def get_all_data(exclude=exclude):
    
    biglist = []
    bigr = []
    bigz = []
    bigSN = []
    ss = re.compile(r'(\d+): \[(.*)\], SNR: (\d+\.\d+)')    
    
    for p in range(6):
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(p+1)
        log = 'NGC_891_P{}_bin30.log'.format(p+1)
        dfile = 'NGC_891_P{}_bin30.msoz.fits'.format(p+1)
        print loc, log, dfile
        r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
        z = np.abs(z)
        r = np.abs(r)

        SN = np.zeros(z.size)
        with open(log,'r') as lff:
            lf = lff.readlines()
        for i, l in enumerate(lf):
            grp = ss.search(l).groups()
            SN[i] = float(grp[2])
        
        hdu = pyfits.open(dfile)[0]
        data = hdu.data

        exar = np.array(exclude[p]) - 1
        r = np.delete(r, exar)
        z = np.delete(z, exar)
        SN = np.delete(SN, exar)
        data = np.delete(data, exar, axis=0)
        
        biglist.append(data)
        bigr.append(r)
        bigz.append(z)
        bigSN.append(SN)

    header = hdu.header    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval

    return wave, np.hstack(bigr), np.hstack(bigz), np.hstack(bigSN), np.vstack(biglist)

def plot_cut_all(wave, r, z, sn, data, ax = None, zcut=[0,2.7], rcut=[0,11]):
    
    if ax is None:
        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Wavelength [$\AA$]')
        ax.set_ylabel('Normalized Flux')

    idx = np.where((z >= zcut[0]) & (z < zcut[1])
                   & (r >= rcut[0]) & (r < rcut[1]))[0]

    avg = np.nansum(data[idx,:]*sn[idx,None],axis=0)/np.nansum(sn[idx])
    
    avg /= avg[810] #Approx 5500

    ax.plot(wave, avg, 'k')

    return ax, avg

def plot_cut_grid(output, zcuts=[0.4], rcuts=[3,8], exclude=exclude):

    fig = plt.figure()
    lax = fig.add_subplot(111)
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_xlabel('Wavelength [$\AA$]')
    lax.set_ylabel('Normalized Flux')
    lax.tick_params(axis='both',pad=20,length=0)

    bigz = [0] + zcuts + [2.6]
    bigr = [0] + rcuts + [11]

    w, r, z, sn, data = get_all_data(exclude=exclude)

    i = 1
    for zz in range(len(zcuts) + 1):
        zc = [bigz[-zz-2], bigz[-zz-1]]
        for rr in range(len(rcuts) + 1):
            rc = [bigr[rr], bigr[rr+1]]
            print zc, rc
            ax = fig.add_subplot(len(zcuts)+1,len(rcuts)+1,i)
            plot_cut_all(w,r,z,sn,data,zcut=zc,rcut=rc,ax=ax)

            ax.text(4000,1.5,'${}\leq |z| <{}$ kpc\n${}\leq |r| <{}$ kpc'.format(*(zc+rc)),ha='left',va='center')
            ax.set_ylim(0,1.9)
            ax.set_xlim(3600,5900)
            if i < 4:
                ax.set_xticklabels([])
            if i % 3 != 1:
                ax.set_yticklabels([])
            i += 1

    fig.subplots_adjust(hspace=0.00001,wspace=0.0001)
    
    pp = PDF(output)
    pp.savefig(fig)
    pp.close() 
    
    return

def plot_pointing(pointing, zcut=0.5, ax = None):

    if ax is None:
        ax = plt.figure().add_subplot(111)
    ax.set_xlabel('Wavelength [$\AA$]')
    ax.set_ylabel('Normalized Flux + offset')
    ax.set_xlim(3500, 6500)
    ax.set_ylim(0,3)
    ax.set_title('P{}'.format(pointing))
    
    wave, low, high = cut_pointing(pointing, zcut)

    #810 ~= 5500 AA
    low = low/low[810]
    high = high/high[810]

    high += 1.0 #np.max(low) - np.min(high) + 0.1

    ax.plot(wave, low, 'k')
    ax.plot(wave, high, 'k')
    ax.text(6000, high[-1], 'z $\geq$ {} kpc'.format(zcut))
    ax.text(6000, low[-1], 'z $<$ {} kpc'.format(zcut))

    return ax

def plot_all_pointings_grid(output, zcut = 0.5):

    al = []
    plist = [6,3,4,2,1,5]
    fig = plt.figure(figsize=(20,6))
    for p in range(6):
        ax = fig.add_subplot(1,6,p+1)
        plot_pointing(plist[p],zcut=zcut,ax=ax)
        if p > 0:
            ax.set_yticklabels([])
            ax.set_ylabel('')

    fig.subplots_adjust(wspace=0.00001)
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
#    plt.close('all')

    return fig

def plot_all_pointings(outpre, zcut=0.5):

    if not isinstance(zcut, list):
        zcut = [zcut] * 6

    for p in range(6):
        output = '{}_P{}.pdf'.format(outpre,p+1)
        pp = PDF(output)
        pp.savefig(plot_pointing(p+1,zcut[p]).figure)
        pp.close()

    plt.close('all')

    return
