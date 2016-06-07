import numpy as np
import pyfits
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
    dfile = 'NGC_891_P{}_bin30.msoz.fits'.format(pointing)
    print loc, dfile
    z = np.loadtxt(loc, usecols=(5,), unpack=True)
    z = np.abs(z)
    hdu = pyfits.open(dfile)[0]
    data = hdu.data
    header = hdu.header    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval

    exar = np.array(exclude[pointing - 1]) - 1
    z = np.delete(z,exar)
    data = np.delete(data,exar,axis=0)

    zhidx = np.where(z >= zcut)[0]
    zlidx = np.where(z < zcut)[0]

    low_d = np.mean(data[zlidx,:], axis=0)
    high_d = np.mean(data[zhidx,:], axis=0)

    return wave, low_d, high_d

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
