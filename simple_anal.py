import sys
import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import scipy.ndimage as spnd
plt.ioff()

def fiber_rms(coefs, wave, models, flux, ax = None):

    frac = coefs['LIGHT_FRAC']
    fracerr = coefs['LIGHT_FRAC_ERR']
    
    print frac
    
    plist = []
    plist.append(np.linspace(coefs['TAUV'] - coefs['TAUV_ERR'],
                             coefs['TAUV'] + coefs['TAUV_ERR'],
                             10))
#    plist.append(coefs['TAUV'])
    for i in range(frac.size):
        if frac[i] == 0:
            plist.append(np.array([0.]))
        else:
            plist.append(np.linspace(frac[i] - fracerr[i],
                                     frac[i] + fracerr[i],
                                     10))

    exrange = zip(*[i.flatten() for i in np.meshgrid(*plist)])

    M = np.zeros(wave.size)
    S = np.zeros(wave.size)
    Mean = np.zeros(wave.size)
    RMS = np.zeros(wave.size)

    rows, cols = os.popen('stty size', 'r').read().split()
    cols = int(cols) - 2

    for i, f in enumerate(exrange):
        sys.stdout.write('\r')
        sys.stdout.write('[{{:<{}}}]'.format(cols).format('='*((i+1)*cols/len(exrange))))
        sys.stdout.flush()

        tmodel = mcombine(wave,models,np.array(f))
        if ax is not None:
            ax.plot(wave,tmodel,'r-',lw=0.2,alpha=0.02)

        Mean += tmodel
        RMS += tmodel**2
        newM = M + (tmodel - M)/(i+1)
        newS = S + (tmodel - M)*(tmodel - newM)

        M = newM
        S = newS

        # try:
        #     bigmodel = np.vstack((bigmodel,tmodel))
        # except UnboundLocalError:
        #     bigmodel = tmodel
            
    Std = np.sqrt(S/i)
    Mean /= i+1

#    RMS = np.sqrt(np.sum((bigmodel - np.mean(bigmodel,axis=0))**2,axis=0)/(i+1))

    rms = np.sqrt(Mean**2 + Std**2)
    rms2 = np.sqrt(RMS/(i+1))
    if ax is not None:
        ax.plot(wave,flux,'k',lw=0.5,alpha=0.6)
    return Mean, Std

def mcombine(wave, models, coefs):
                 
    tauv = coefs[0]

    y = np.sum(models * coefs[1:,None], axis=0)

    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tauv*klam)
    y *= e_tau_lam

    return y

def do_fibers(galaxy, datafile, model):

    dhdu = pyfits.open(galaxy)[0]
    header = dhdu.header
    data = dhdu.data*1e19
    numfibers, wavesize = data.shape
    wave = (np.arange(wavesize) - header['CRPIX1']) * header['CDELT1'] \
           +  header['CRVAL1']
    idx = np.where((wave >= 3750) & (wave <= 6800))
    wave = wave[idx]

    m = pyfits.open(model)[1].data
    nmodels = m['AGE'][0].size
    npix = wave.size

    vdisp = 691.0
    bc03_pix = 70.0
    bc03_vdisp = 75.0
    vdisp_add = np.sqrt(vdisp**2 - bc03_vdisp**2)
    sigma_pix = vdisp_add / bc03_pix

    custom_lib = np.zeros((nmodels, npix))

    for i in range(nmodels):
        cflux = spnd.filters.gaussian_filter1d(m['FLUX'][0][i], sigma_pix)
        custom_lib[i] = np.interp(wave,m['WAVE'][0],cflux)

    fitcoefs = pyfits.open(datafile)[1].data

    for f in [74]:

        flux = data[f,idx][0]

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Wavelength [Angstroms]')
        ax.set_ylabel('Flux')
        mean, std = fiber_rms(fitcoefs[f], wave, custom_lib, flux, ax=ax)
        
        sax = plt.figure().add_subplot(111)
        sax.plot(wave,std/mean)
        sax.set_xlabel('Wavelength [Angstroms]')
        sax.set_ylabel('Fit $\sigma$')
        sax.figure.show()
        ax.figure.show()

        return mean, std
