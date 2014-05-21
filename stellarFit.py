import numpy as np
import pyfits
import glob
import ADEUtils as ADE
import matplotlib.pyplot as plt
import matplotlib.collections as collections
from ppxf import ppxf
import ppxf_util as pputil
from scipy import ndimage
from pyraf import iraf
from matplotlib.backends.backend_pdf import PdfPages as PDF
iraf.noao.onedspec()

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

def logify(infits, outfits, w1='INDEF',w2='INDEF',
           dw='INDEF',nw='INDEF',report=False):
    
    iraf.noao.onedspec.dispcor(infits,outfits,
                               linearize=True,log=False,flux=True,
                               dw=dw, w1=w1, w2=w2, nw=nw)
    if report:
        hdu = pyfits.open(outfits)[0]
        dw = hdu.header['CDELT1']
        nw = hdu.header['NAXIS1']
        w1 = hdu.header['CRVAL1']
        w2 = w1 + nw*dw

    return w1, w2, dw, nw

def prep_templates(template_fits, outputfits, dw):

    thdu = pyfits.open(template_fits)[0]
    tw1 = thdu.header['CRVAL1']
    tn2 = thud.header['NAXIS1']
    

def fitms(spectrum,template_list, out_prefix, order=5, cut=0.75, pre=225):

    pd = PDF(out_prefix+'.pdf')

    #make the galaxy wavelength vector
    data_hdu = pyfits.open(spectrum)
    ddw = data_hdu[0].header['CDELT1']
    wave = data_hdu[0].header['CRVAL1'] + \
        np.arange(data_hdu[0].shape[1])*ddw
#    wave = wave[pre:]
    lamRange1 = np.array([wave.min(),wave.max()])
    galaxy = data_hdu[0].data[0]#[pre:]
    loggalaxy, logLam1, velscale = pputil.log_rebin(lamRange1,galaxy)

    print data_hdu[0].shape

    c = 299792.458
    FWHM_gal = 0.95 #AA
    FWHM_tmp = 0.55 #AA
    FWHM_diff = np.sqrt(FWHM_gal**2 - FWHM_tmp**2)
    sigma = FWHM_diff/2.355/ddw

    temp_hdu = pyfits.open(template_list[0])
    twave = temp_hdu[1].header['CRVAL1'] + \
        np.arange(temp_hdu[1].data.size)*temp_hdu[1].header['CDELT1']
    mask = np.where((twave > 4650.) & (twave < 5500.))
    template = temp_hdu[1].data[mask]
    twave = twave[mask]
    lamRange2 = np.array([twave.min(),twave.max()])
    template = ndimage.gaussian_filter1d(template,sigma)
    logtmp, logLam2, velscale_tmp = pputil.log_rebin(lamRange2,
                                                     template,
                                                     velscale=velscale)

    bestfits = np.empty((1,wave.size))
    templates = np.empty((logtmp.size))
    print templates.shape

    for template_file in template_list:
        temp_hdu = pyfits.open(template_file)
        for i in range(len(temp_hdu)-1):
            hdu = temp_hdu[i+1]
            twave = temp_hdu[i+1].header['CRVAL1'] + \
                np.arange(temp_hdu[i+1].data.size)*temp_hdu[i+1].header['CDELT1']
            mask = np.where((twave > 4650.) & (twave < 5500.))
            twave = twave[mask]
            td = hdu.data[mask]
            tfit = ADE.polyclip(twave,td,order)
            td /= tfit(twave)
            td = ndimage.gaussian_filter1d(td,sigma)
            logtd, logLam2, velscale = pputil.log_rebin(lamRange2,td,
                                                        velscale=velscale)
            templates = np.vstack((templates,logtd/np.median(logtd)))
    
    templates = templates[1:].T
    print templates.shape

    for i in range(data_hdu[0].data.shape[0]):

        fig = plt.figure()

        galaxy = data_hdu[0].data[i]#[pre:]
#        datafit = ADE.polyclip(wave,galaxy,order)

        ax = fig.add_subplot(211)
        ax.plot(wave,galaxy)
#        ax.plot(wave,datafit(wave))
#        galaxy /= datafit(wave)

        loggalaxy, logLam1, velscale = pputil.log_rebin(lamRange1,galaxy)
#        pyfits.PrimaryHDU(templates).writeto('templates.fits')

        dv = (logLam2[0] - logLam1[0])*c
        vel = c*1.008246
        start = [0,180.]
        loggalaxy /= np.median(loggalaxy)
        goodpixels = np.where(np.abs(loggalaxy - np.mean(loggalaxy)) < \
                                  cut*np.std(loggalaxy))[0]
        goodpixels = goodpixels[goodpixels > pre]
 
        collection = collections.BrokenBarHCollection.span_where(
            np.arange(loggalaxy.size),
            ymin=0,ymax=4,
            where=np.abs(loggalaxy - np.mean(loggalaxy)) \
                >= cut*np.std(loggalaxy), 
            facecolor='red',alpha=0.5)
        collection2 = collections.BrokenBarHCollection.span_where(
            np.arange(loggalaxy.size),
            ymin=0,ymax=4,
            where=np.arange(loggalaxy.size) <= pre, 
            facecolor='red',alpha=0.5)

        ax2 = fig.add_subplot(212)
        ax2.plot(np.arange(loggalaxy.size),loggalaxy)
        ax2.add_collection(collection)
        ax2.add_collection(collection2)
        
        pp = ppxf(templates,loggalaxy, np.ones(loggalaxy.shape)*1., 
                  velscale,start,bias=0,degree=25,mdegree=0,
                  vsyst=dv,plot=True, goodpixels=goodpixels)

        print bestfits.shape, pp.bestfit.size
        bestfits = np.vstack((bestfits,pp.bestfit))

        ax2.plot(np.arange(loggalaxy.size),pp.bestfit)
        ax2.set_ylim(0.8,1.2)
        pd.savefig(fig)

    bestfits = bestfits[1:]
    pd.close()
    bfh = pyfits.PrimaryHDU(np.array(bestfits,dtype=np.float32))
    bfh.header.update('CDELT1',ddw)
    bfh.header.update('CRVAL1',wave.min())
    bfh.header.update('CRPIX1',1)
    bfh.header.update('CTYPE1','LINEAR')
    bfh.writeto(out_prefix+'.fits',clobber=True)
    iraf.noao.onedspec.dispcor(out_prefix+'.fits',out_prefix+'_lin.fits',
                               linearize=True,log=False,flux=False,
                               dw=ddw,w1=wave.min(),w2=wave.max(),nw='INDEF',
                               samedis=True)


    return
    
