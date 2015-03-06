import numpy as np
import pyfits
import matplotlib.pyplot as plt
import scipy.ndimage as spnd
plt.ioff()

def fiber_rms(coefs, wave, models, flux):

    frac = coefs['LIGHT_FRAC']
    fracerr = coefs['LIGHT_FRAC_ERR']
    
    print frac
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Wavelenght [Angstroms]')
    ax.set_ylabel('Flux')
    ax.set_title('Fiber 75')
    for i in range(frac.size):
        print i
        inf = frac.copy()
        if frac[i] == 0:
            print 'skipped '+str(i)
            continue

        exrange = np.linspace(frac[i] - fracerr[i],
                              frac[i] + fracerr[i],
                              10)
        for f in exrange:
            
            inf[i] = f
            tmodel = mcombine(wave,models,coefs['TAUV'],inf)
            ax.plot(wave,tmodel,'r-',lw=0.3,alpha=0.4)
            try:
                bigmodel = np.vstack((bigmodel,tmodel))
            except UnboundLocalError:
                print 'Creating'
                bigmodel = tmodel


    print wave.shape, flux.shape
    ax.plot(wave,flux,'k',lw=0.5,alpha=0.6)
    fig.show()
    return bigmodel

def mcombine(wave, models, tauv, coefs):

    y = np.sum(models * coefs[:,None], axis=0)

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

        return fiber_rms(fitcoefs[f], wave, custom_lib, flux)
        

    
