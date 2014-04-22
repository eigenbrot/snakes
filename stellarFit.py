import numpy as np
import pyfits
import glob
import ADEUtils as ADE
import matplotlib.pyplot as plt

glob = glob.glob

def make_SSP_fits(input_file,output_fits):

    ages, metallicities, wavelengths, fluxes = np.loadtxt(input_file,
                                                          unpack=True)

    hdulist = []

    for age in np.unique(ages):
        aidx = np.where(ages == age)
        for Z in np.unique(metallicities[aidx]):
            print age, Z
            Zidx = np.where(metallicities[aidx] == Z)
            
            wave = wavelengths[aidx][Zidx]
            flux = fluxes[aidx][Zidx]
            hdu = pyfits.ImageHDU(flux)
            hdu.header.update('CDELT1',np.mean(np.diff(wave)))
            hdu.header.update('CRVAL1',np.min(wave))
            hdu.header.update('CRPIX1',1)
            hdu.header.update('Z',Z,comment='Metallicity')
            hdu.header.update('AGE',age,comment='Age [Gyr]')
            hdulist.append(hdu)

    
    primary = pyfits.PrimaryHDU()
    primary.header.update('INPUT',input_file,
                          comment='File used to construct FITS')

    pyfits.HDUList([primary] + hdulist).writeto(output_fits,clobber=True)

    return

def test(spectrum,template,line=10, order=5):

    data_hdu = pyfits.open(spectrum)
    wave = data_hdu[0].header['CRVAL1'] + \
        np.arange(data_hdu[0].data[line].size)*data_hdu[0].header['CDELT1']
    galaxy = data_hdu[0].data[line]
    
    print wave.size, galaxy.size

    temp_hdu = pyfits.open(template)
    twave = temp_hdu[1].header['CRVAL1'] + \
        np.arange(temp_hdu[1].data.size)*temp_hdu[1].header['CDELT1']
    template = temp_hdu[1].data

    mask = np.where((twave > wave.min()) & (twave < wave.max()))
    twave = twave[mask]
    template = template[mask]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(wave,galaxy)
    datafit = ADE.polyclip(wave,galaxy,order)
    ax.plot(wave,datafit(wave))
    galaxy -= datafit(wave)
    ax2 = fig.add_subplot(212)
    ax2.plot(wave,galaxy)

    fig2 = plt.figure()
    ax3 = fig2.add_subplot(211)
    ax3.plot(twave,template)
    tmpfit = ADE.polyclip(twave,template,order)
    ax3.plot(twave,tmpfit(twave))
    template -= tmpfit(twave)
    ax4 = fig2.add_subplot(212)
    ax4.plot(twave,template)
    
    fig.show()
    fig2.show()

    c = 299792.458
    velscale = np.log(wave[1]/wave[0])*c
