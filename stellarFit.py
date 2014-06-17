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
import os

glob = glob.glob

def make_SSP_fits(input_file,output_fits):

    ages, metallicities, wavelengths, fluxes = np.loadtxt(input_file,
                                                          unpack=True)

    header = pyfits.header.Header()
    stack = []
    i = 0
    pd = 0
    pm = 0
    for age in np.unique(ages):
        aidx = np.where(ages == age)
        for Z in np.unique(metallicities[aidx]):
            print age, Z
            Zidx = np.where(metallicities[aidx] == Z)   
            wave = wavelengths[aidx][Zidx]
            flux = fluxes[aidx][Zidx]
            stack.append(flux/np.mean(flux))
            header.update('Z{}'.format(i),Z,comment='Metallicity')
            header.update('AGE{}'.format(i),age,comment='Age [Gyr]')
            d = np.mean(np.diff(wave))
            m = np.min(wave)
            if pd != d or pm != m:
                print 'WARNING, {} != {} or {} != {}'.format(pd,d,pm,m)
            pd = d
            pm = m
            i += 1

    header.update('CDELT1',np.mean(np.diff(wave)))
    header.update('CRVAL1',np.min(wave))
    header.update('CRPIX1',1)
    
    data = np.vstack(stack)
    header.update('INPUT',input_file,
                  comment='File used to construct FITS')
    primary = pyfits.PrimaryHDU(data,header)
    primary.writeto(output_fits,clobber=True)

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
    

def fitms(spectrum,error,template_list, out_prefix, cut=0.75, pre=0, mdegree=25, degree=4):

    pd = PDF(out_prefix+'.pdf')

    #make the galaxy wavelength vector
    data_hdu = pyfits.open(spectrum)
    error_hdu = pyfits.open(error)
    ddw = data_hdu[0].header['CDELT1']
    lambda0 = data_hdu[0].header['CRVAL1']
    wave =  lambda0+ \
        np.arange(data_hdu[0].shape[1])*ddw
#    wave = wave[pre:]
    lamRange1 = np.array([wave.min(),wave.max()])#/1.008246
    galaxy = data_hdu[0].data[0]#[pre:]
    loggalaxy, logLam1, velscale = pputil.log_rebin(lamRange1,galaxy)

    masklow = wave.min() - 500.
    maskhigh = wave.max() + 500.

    print data_hdu[0].shape
    print masklow,maskhigh

    c = 299792.458
    FWHM_gal = 0.95/1.008246 #AA
    FWHM_tmp = 0.55 #AA
    FWHM_diff = np.sqrt(FWHM_gal**2 - FWHM_tmp**2)
    sigma = FWHM_diff/2.355/ddw

    temp_hdu = pyfits.open(template_list[0])
    twave = temp_hdu[0].header['CRVAL1'] + \
        np.arange(temp_hdu[0].data.shape[1])*temp_hdu[0].header['CDELT1']
    mask = np.where((twave > masklow) & (twave < maskhigh))
    template = temp_hdu[0].data[0][mask]
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
        temp_hdu = pyfits.open(template_file)[0]
        data = temp_hdu.data
        header = temp_hdu.header
        twave = header['CRVAL1'] + np.arange(data.shape[1])*header['CDELT1']
        mask = (twave > masklow) & (twave < maskhigh)
        twave = twave[mask]
        for i in range(data.shape[0]):
            td = data[i][mask]
            # tfit = ADE.polyclip(twave,td,order)
            # td /= tfit(twave)
            td = ndimage.gaussian_filter1d(td,sigma)
            logtd, logLam2, velscale = pputil.log_rebin(lamRange2,td,
                                                        velscale=velscale)
            templates = np.vstack((templates,logtd/np.median(logtd)))
    
    templates = templates[1:].T
    print templates.shape, twave.shape

    for i in range(data_hdu[0].data.shape[0]):

        fig = plt.figure()

        galaxy = data_hdu[0].data[i]#[pre:]
        noise = error_hdu[0].data[i]
#        datafit = ADE.polyclip(wave,galaxy,order)

        # ax = fig.add_subplot(212)
        # ax.plot(wave,galaxy)
#        ax.plot(np.arange(templates[:,0].size),templates[:,0])
#        ax.set_ylim(0,2)
#        ax.plot(wave,datafit(wave))
#        galaxy /= datafit(wave)

        loggalaxy, logLam1, velscale = pputil.log_rebin(lamRange1,galaxy)
        logerror, _, _ = pputil.log_rebin(lamRange1,noise)
#        pyfits.PrimaryHDU(templates).writeto('templates.fits')
        print velscale

        dv = (logLam2[0] - logLam1[0])*c
        print dv
        vel = c*1.008246
        start = [0.,2.]
        loggalaxy /= np.median(loggalaxy)
        logerror /= np.median(loggalaxy)
        
        goodfit = ADE.polyclip(np.arange(loggalaxy.size),loggalaxy,2)
        tmpgood = np.where(np.abs(loggalaxy - goodfit(np.arange(loggalaxy.size))) < \
                                  cut*np.std(loggalaxy))[0]
        tmpgood = tmpgood[tmpgood > pre]

        goodsmooth = np.zeros(loggalaxy.size)
        goodsmooth[tmpgood] = 1.
        goodsmooth = ndimage.gaussian_filter1d(goodsmooth,11)
        goodpixels = np.where(goodsmooth > 0.9)[0]

        pp = ppxf(templates,loggalaxy, logerror, 
                  velscale,start,bias=None,degree=degree,mdegree=mdegree,
                  vsyst=dv,plot=True, goodpixels=goodpixels)
        print bestfits.shape
        bestfits = np.vstack((bestfits,pp.bestfit))
        
        ###PLOT###
        plot_px = np.exp(logLam1)
        plot_gal = ndimage.gaussian_filter1d(loggalaxy,3)
        collection = collections.BrokenBarHCollection.span_where(
            plot_px,
            ymin=0,ymax=4,
            where=goodsmooth < 0.9,
            facecolor='red',alpha=0.5)
        collection2 = collections.BrokenBarHCollection.span_where(
            plot_px,
            ymin=0,ymax=4,
            where=np.arange(loggalaxy.size) <= pre, 
            facecolor='red',alpha=0.5)


        ax2 = fig.add_subplot(411)
        ax2.plot(np.exp(logLam1),plot_gal)
        ax2.add_collection(collection)
        ax2.add_collection(collection2)

        ax2.plot(plot_px,pp.bestfit)
        ax2.set_xlim(plot_px[0],plot_px[-1])
        ax5 = ax2.twiny()
        ax5.set_xlim(0,loggalaxy.size)
        ax2.set_ylim(0.6,1.4)
        bbox2 = ax2.get_position().get_points().flatten()

        plt.setp(ax2.get_xticklabels(),visible=False)
        ax3 = fig.add_subplot(412,sharex=ax2)
        bbox3 = ax3.get_position().get_points().flatten()
        newpos = [bbox3[0],
                  bbox2[1] - (bbox3[3] - bbox3[1]),
                  bbox3[2] - bbox3[0],
                  bbox3[3] - bbox3[1]]
        ax3.set_position(newpos)
        ax3.plot(plot_px,loggalaxy - pp.bestfit,'r',lw=0.7)     
        ax3.set_ylim(-0.2,0.2)
        ax3.set_xlim(plot_px[0],plot_px[-1])

        ax = fig.add_subplot(212)
        pidx = np.arange(1400,1600,1)
        ax.plot(plot_px[pidx],plot_gal[pidx])
        ax.plot(plot_px[pidx],pp.bestfit[pidx])
        ax4 = ax.twiny()
        velx = np.arange(pidx.size)*velscale
        velx -= velx[-1]/2.
        ax4.set_xlim(velx[0],velx[-1])

        pd.savefig(fig)

    bestfits = bestfits[1:]
    pd.close()
    bfh = pyfits.PrimaryHDU(np.array(bestfits,dtype=np.float32))
    bfh.header.update('CDELT1',ddw)
    bfh.header.update('CRVAL1',wave.min())
    bfh.header.update('CRPIX1',1)
    bfh.header.update('CTYPE1','LINEAR')
    bfh.writeto(out_prefix+'.ms.fits',clobber=True)
    iraf.noao.onedspec.dispcor(out_prefix+'.ms.fits',out_prefix+'_lin.ms.fits',
                               linearize=True,log=False,flux=False,
                               dw=ddw,w1=wave.min(),w2=wave.max(),nw='INDEF',
                               samedis=True)


    return
    
def flatten_spectra(input_file,output_file, error=None):

    iraf.fit1d(input_file,'tmp.fits',type='fit',interact=False,
               functio='legendre',order=2)
    iraf.imarith(input_file,'/','tmp.fits','tmp2.fits')

    hdu = pyfits.open('tmp2.fits')[0]
    header = hdu.header
    data = hdu.data
    t = np.r_[data.T,data.T,data.T]
    pyfits.PrimaryHDU(t.T,header).writeto('tmp3.fits')
    
    iraf.fit1d('tmp3.fits','tmp4.fits',type='fit',interact=False,
               functio='spline3',order=100,low_rej=2.,high_rej=2.)
    hdu2 = pyfits.open('tmp4.fits')[0]
    data2 = hdu2.data
    pyfits.PrimaryHDU(data2[:,data2.shape[1]/3.:data2.shape[1]*2./3.],
                      header).writeto('tmp5.fits')
    iraf.imarith('tmp5.fits','*','tmp.fits','tmp6.fits')
    d6 = pyfits.open('tmp6.fits')[0].data
    # d7 = d6/np.mean(d6,axis=1)[:,np.newaxis]
    # pyfits.PrimaryHDU(d7,header).writeto('tmp7.fits')
    iraf.imarith(input_file,'/','tmp6.fits',output_file)
#    os.system('rm tmp*.fits')
    
    return np.mean(d6,axis=1)[:,np.newaxis]

def runtest(input_spectra,output_prefix):

    templates = ['ELODIE_krz004.norm.ms.fits',
                 'ELODIE_krz002.norm.ms.fits',
                 'ELODIE_krz001.norm.ms.fits']

    # First flatten the data
    # print 'Normalizing input spectra...'
    input_prefix = os.path.basename(input_spectra).split('.ms.fits')[0]
    # flat_name = '{}.norm.ms.fits'.format(input_prefix)
    # means = flatten_spectra(input_spectra,flat_name)
    
    # Now fit stellar templates
    print 'Fitting templates...'
    error = '{}_error.{}'.format(input_spectra.split('.')[0],
                                '.'.join(input_spectra.split('.')[1:]))
    fitms(input_spectra,error,templates,output_prefix,cut=0.4,mdegree=0,degree=120)

    #grab the best fit, scale it by the slit function, and subtract it from
    #the OG spectra
    print 'Subtracting template...'
    bestfit_file = '{}_lin.ms.fits'.format(output_prefix)
    bfdata = pyfits.open(bestfit_file)[0].data
    
    spech = pyfits.open(input_spectra)[0]
    header = spech.header
    spectra = spech.data

    sub_spectra = (spectra - bfdata)

    final_name = '{}.norm.sub.ms.fits'.format(input_prefix)
    pyfits.PrimaryHDU(sub_spectra,header).writeto(final_name,clobber=True)

    return
