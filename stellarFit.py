import numpy as np
import pyfits
import glob
import ADEUtils as ADE
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from ppxf import ppxf
import ppxf_util as pputil
from scipy import ndimage

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

def test(spectrum,template,line=10, order=5, cut=0.75):

    data_hdu = pyfits.open(spectrum)
    wave = data_hdu[0].header['CRVAL1'] + \
        np.arange(data_hdu[0].data[line].size)*data_hdu[0].header['CDELT1']
    galaxy = data_hdu[0].data[line][106:]
    wave = wave[106:]

    temp_hdu = pyfits.open(template)
    twave = temp_hdu[1].header['CRVAL1'] + \
        np.arange(temp_hdu[1].data.size)*temp_hdu[1].header['CDELT1']
    template = temp_hdu[1].data

    mask = np.where((twave > wave.min()) & (twave < wave.max()))
    # twave = twave[mask]
    # template = template[mask]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(wave,galaxy)
    datafit = ADE.polyclip(wave,galaxy,order)
    ax.plot(wave,datafit(wave))
    galaxy /= datafit(wave)
    ax2 = fig.add_subplot(212)
#    ax2.plot(wave,galaxy)

    fig2 = plt.figure()
    ax3 = fig2.add_subplot(211)
    ax3.plot(twave,template)
    tmpfit = ADE.polyclip(twave,template,order)
    ax3.plot(twave,tmpfit(twave))
    template /= tmpfit(twave)
    ax4 = fig2.add_subplot(212)
#    ax4.plot(twave,template)
    
    print wave.shape, twave.shape

    c = 299792.458
    FWHM_gal = 0.95 #AA
    FWHM_tmp = 0.55 #AA
    FWHM_diff = np.sqrt(FWHM_gal**2 - FWHM_tmp**2)
    sigma = FWHM_diff/2.355/data_hdu[0].header['CDELT1']

    lamRange1 = np.array([wave.min(),wave.max()])
    loggalaxy, logLam1, velscale = pputil.log_rebin(lamRange1,galaxy)

    lamRange2 = np.array([twave.min(),twave.max()])
    template = ndimage.gaussian_filter1d(template,sigma)
    ax4.plot(twave,template)
    logtmp, logLam2, velscale_tmp = pputil.log_rebin(lamRange2,
                                                     template,
                                                     velscale=velscale)
    templates = np.empty((logtmp.size,len(temp_hdu)-1))

    for i in range(len(temp_hdu)-1):
        hdu = temp_hdu[i+1]
        twave = temp_hdu[1].header['CRVAL1'] + \
        np.arange(temp_hdu[1].data.size)*temp_hdu[1].header['CDELT1']
        td = hdu.data
        tfit = ADE.polyclip(twave,td,order)
        td /= tfit(twave)
        td = ndimage.gaussian_filter1d(td,sigma)
        logtd, logLam2, velscale = pputil.log_rebin(lamRange2,td,
                                                    velscale=velscale)
        templates[:,i] = logtd/np.median(logtd)
    
    print templates.shape

    dv = (logLam2[0] - logLam1[0])*c
    vel = c*1.008246
#    goodpixels = pputil.determine_goodpixels(np.log(twave),lamRange,vel)
#    goodpixels = np.where(loggalaxy > 0.8)[0]
    start = [0,180.]
#    template /= np.median(template)
    loggalaxy /= np.median(loggalaxy)
    goodpixels = np.where(np.abs(loggalaxy - np.mean(loggalaxy)) < \
                              cut*np.std(loggalaxy))[0]
    collection = collections.BrokenBarHCollection.span_where(
        np.arange(loggalaxy.size),
        ymin=0,ymax=4,
        where=np.abs(loggalaxy - np.mean(loggalaxy)) >= cut*np.std(loggalaxy), 
        facecolor='red',alpha=0.5)
    ax2.plot(np.arange(loggalaxy.size),loggalaxy)
    ax2.add_collection(collection)
        

    fig.show()
    fig2.show()

    pp = ppxf(templates,loggalaxy, np.ones(loggalaxy.shape), 
              velscale,start, 
              vsyst=dv,plot=True, goodpixels=goodpixels)

    return pp
    
