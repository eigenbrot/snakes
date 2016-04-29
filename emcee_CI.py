import time
import numpy as np
import pyfits
import emcee
import triangle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def main(datafile, errorfile, location, coeffile, output,
         model='/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_allZ_vardisp.fits',
         plot=True, wavemin=3800., wavemax=6800., fitregion=[3850.,6650.],
         lightmin=5450., lightmax=5550., emmaskw=500.,
         nsample=1000, burn=100, nwalkers=256, threads=1):

    #read the models
    m = pyfits.open(model)[1].data[0] #[0] because IDL creates an extra dim

    #read the data
    dhdu= pyfits.open(datafile)[0]
    data = dhdu.data
    header = dhdu.header
    error = pyfits.open(errorfile)[0].data
    numfibers, wavesize = data.shape

    #read LM best fits
    coef_arr = pyfits.open(coeffile)[1].data

    CDELT = header['CDELT1']
    CRVAL = header['CRVAL1']
    CRPIX = header['CRPIX1']
    print 'CDELT = ', CDELT
    print 'CRVAL1 = ', CRVAL
    print 'CRPIX1 = ', CRPIX
    print data.shape
    wave = (np.arange(wavesize) - CRPIX - 1) * CDELT + CRVAL

    fiber_radii = np.loadtxt(location, usecols=(1,), unpack=True)
    sizeidx = [0.937,1.406,1.875,2.344,2.812]

    idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
    lightidx = np.where((wave >= lightmin) & (wave <= lightmax))[0]

    wave = wave[idx]

    agearr = m['AGE']/1e9
    numages = agearr.shape[0]

    # f = open(output,'w')
    # f.write('{:11}'.format('# Fiber Num'))
    # f.write((numages*'{:9.3f} Gyr').format(*agearr))
    # f.write('{:>10}{:>10}'.format('Tau_V, Vsys'))
    # f.write('\n#\n')

    if plot:
        plotname = output.split('.')[0] + '_tri.pdf'
        pp = PDF(plotname)

    flux_factor = 1e17
    
    outputarr = np.zeros(numfibers,dtype={'names':['VSYS',
                                                   'FIXEDVBOOL',
                                                   'TAUV',
                                                   'TAUV_ERR',
                                                   'emmaskw',
                                                   'LIGHT_FRAC',
                                                   'LIGHT_FRAC_ERR',
                                                   'SNR',
                                                   'MLWA','MLWA_L','MLWA_H',
                                                   'MMWA',
                                                   'MLWZ',
                                                   'MMWZ',
                                                   'chisq','redchi','bluechi','hkchi'],
                                          'formats':['f4',
                                                     'i8',
                                                     'f4',
                                                     'f4',
                                                     'f4',
                                                     '{}f4'.format(numages),
                                                     '{}f4'.format(numages),
                                                     'f4',
                                                     'f4','f4','f4',
                                                     'f4',
                                                     'f4',
                                                     'f4',
                                                     'f4','f4','f4','f4']})

    yfitfile = output.split('.')[0] + '.emceefit.fits'
    yfitarr = np.zeros((numfibers, wave.size))

    for i in [14]:#range(numfiber):

        print "Doing fiber {}".format(i+1)
        flux = data[i,idx]*flux_factor
        err = error[i,idx]*flux_factor

        vdidx = np.where(sizeidx == fiber_radii[i])[0][0]
        LMcoefs = coef_arr[i]
        
        MCcoefs, yfit, S = EMfit(m, wave, flux, err, vdidx, LMcoefs,
                                 fitregion=fitregion,emmaskw=emmaskw,
                                 lightidx=lightidx,nsample=nsample,
                                 burn=burn,nwalkers=nwalkers,threads=threads)

        outputarr[i] = MCcoefs
        yfitarr[i,:] = yfit/flux_factor

        if plot:
            fig = triangle.corner(S.flatchain*np.r_[1,np.ones(numages)*100.],
                                  labels=[r'$\tau_V$'] + ['$w_{}$'.format(i) for i in range(numages)],
                                  truths=np.r_[MCcoefs['TAUV'],MCcoefs['LIGHT_FRAC']])
            pp.savefig(fig)
            plt.close(fig)

    if plot:
        pp.close()

    pyfits.BinTableHDU(outputarr).writeto(output,clobber=True)
    fithdu = pyfits.PrimaryHDU(yfitarr)
    fithdu.header.update('CDELT1',CDELT)
    fithdu.header.update('CRPIX1',1)
    fithdu.header.update('CRVAL1',wavemin)
    fithdu.writeto(yfitfile,clobber=True)

    return MCcoefs, S

def EMfit(model, wave, flux, err, vdidx, LMcoefs, fitregion=[3850.,6650.],
          emmaskw=500., lightidx=None, threads=1,
          nsample = 1000, burn=100, nwalkers=256):

    nmodels = model['AGE'].shape[0]
    npix = wave.size
    quality = np.ones(npix)
    light_factor = 100.

    bad = np.where((np.isnan(flux) == True) | (np.isnan(err) == True) |\
                   (err == 0))[0]
    quality[bad] = 0

    em= np.array([4959.91, # OIII
                  5006.84, # OIII
                  6562.82, # Ha
                  6583.41, # NII
                  6548.04, # NII
                  6716.44, # SII
                  6730.81]) # SII
    # bad sky lines
    sk = [6300., 5890., 5577.]
    
    em *= (LMcoefs['VSYS']/3e5 + 1)
    dz = emmaskw/3e5
    dzsk = 1500./3e5

    for e in em:
        maskout = np.where((wave > e*(1-dz)) & (wave < e*(1+dz)))[0]
        if maskout.size > 0:
            quality[maskout] = 0
    for s in sk:
        maskout = np.where((wave > s*(1-dzsk)) & (wave < s*(1+dzsk)))[0]
        if maskout.size > 0:
            quality[maskout] = 0

    okidx = np.where(quality == 1)[0]
    
    custom_lib = np.zeros((nmodels,npix))
    for i in range(nmodels):
        custom_lib[i,:] = np.interp(wave,model['WAVE'],model['FLUX'][i,:,vdidx])
    
    fitidx = np.where((wave[okidx] >= fitregion[0]) &\
                      (wave[okidx] <= fitregion[1]))[0]
    fitflux = flux[okidx[fitidx]]
    fiterr = err[okidx[fitidx]]
    fitwave = wave[okidx[fitidx]]
    fitlib = custom_lib[:,okidx[fitidx]]

    #Set up the emcee stuff
    ndim = nmodels + 1
    LMtheta = np.r_[LMcoefs['TAUV'],LMcoefs['LIGHT_FRAC']/light_factor]
    print LMtheta
    p0 = [LMtheta + 1e-1*LMtheta*np.random.randn(ndim) for i in range(nwalkers)]
    S = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=threads, args=(fitwave,
                                                                             fitflux,
                                                                             fiterr,
                                                                             fitlib,
                                                                             LMcoefs['VSYS']))
    #run
    print 'Burning...'
    t1 = time.time()
    pos, prob, state = S.run_mcmc(p0,burn)
    print 'Burned {} steps in {} s'.format(burn,time.time() - t1)
    
    S.reset()
    print 'Running...'
    t1 = time.time()
    pos, prob, state = S.run_mcmc(pos,nsample)
    print 'Ran {} steps in {} s'.format(nsample,time.time() - t1)

    maxidx = np.argmax(S.flatlnprobability)
    EMtheta = S.flatchain[maxidx]
    print '###', EMtheta
    stds = np.mean(S.flatchain,axis=0)

    cents, lowCI, highCI = compute_theta_CI(S)
    print '$$$', cents
    print lowCI
    print highCI

    MLWA, MLWA_L, MLWA_H = compute_MLWA_CI(S, wave, custom_lib, 
                                           model['AGE'][:,0]/1e9, lightidx)
    print MLWA, MLWA_L, MLWA_H

    yfit = mcombine(EMtheta,wave,custom_lib,LMcoefs['VSYS'])
    chisq = lnprob(EMtheta, wave, flux, err, custom_lib, LMcoefs['VSYS'])

    fitcoefs = np.zeros(1,dtype={'names':['VSYS',
                                          'FIXEDVBOOL',
                                          'TAUV',
                                          'TAUV_ERR',
                                          'emmaskw',
                                          'LIGHT_FRAC',
                                          'LIGHT_FRAC_ERR',
                                          'SNR',
                                          'MLWA','MLWA_L','MLWA_H',
                                          'MMWA',
                                          'MLWZ',
                                          'MMWZ',
                                          'chisq','redchi','bluechi','hkchi'],
                                'formats':['f4',
                                           'i8',
                                           'f4',
                                           'f4',
                                           'f4',
                                           '{}f4'.format(nmodels),
                                           '{}f4'.format(nmodels),
                                           'f4',
                                           'f4','f4','f4',
                                           'f4',
                                           'f4',
                                           'f4',
                                           'f4','f4','f4','f4']})

    fitcoefs['VSYS'] = LMcoefs['VSYS']
    fitcoefs['emmaskw'] = emmaskw
    fitcoefs['TAUV'] = EMtheta[0]
    fitcoefs['TAUV_ERR'] = stds[0]
    fitcoefs['LIGHT_FRAC'] = EMtheta[1:]*light_factor
    fitcoefs['LIGHT_FRAC_ERR'] = stds[1:]*light_factor

    fitcoefs['FIXEDVBOOL'] = 1
    fitcoefs['SNR'] = LMcoefs['SNR']
    fitcoefs['MLWA'] = MLWA
    fitcoefs['MLWA_L'] = MLWA_L
    fitcoefs['MLWA_H'] = MLWA_H
    fitcoefs['MMWA'] = LMcoefs['MMWA']
    fitcoefs['MLWZ'] = LMcoefs['MLWZ']
    fitcoefs['MMWZ'] = LMcoefs['MMWZ']
    fitcoefs['chisq'] = chisq
    fitcoefs['redchi'] = LMcoefs['redchi']
    fitcoefs['bluechi'] = LMcoefs['bluechi']
    fitcoefs['hkchi'] = LMcoefs['hkchi']

    return fitcoefs[0], yfit, S

def compute_theta_CI(S, bins=20):

    numpar = S.flatchain.shape[1]
    lowCI = np.zeros(numpar)
    highCI = np.zeros(numpar)
    cent = np.zeros(numpar)

    for p in range(numpar):
        hist, b = np.histogram(S.flatchain[:,p],bins=bins)
        cdf = np.cumsum(1.0*hist/np.sum(hist))
        bcent = 0.5*(b[1:] + b[:-1])
        cent[p] = np.interp(0.5,cdf,bcent)
        lowCI[p] = np.interp(0.32,cdf,bcent)
        highCI[p] = np.interp(0.68,cdf,bcent)

    return cent, cent-lowCI, highCI-cent

def compute_MLWA_CI(S, wave, mlib, agearr, lightidx, nsample=None, bins=20):

    if nsample is None:
        nsample = int(S.flatchain.shape[0]/10.)
    
    MLWA_samples = np.zeros(nsample)
    ids = np.random.randint(0,S.flatchain.shape[0],nsample)

    for i in range(nsample):
        MLWA_samples[i] = compute_MLWA(S.flatchain[ids[i],:],
                                       wave, mlib, agearr, lightidx)

    hist, b = np.histogram(MLWA_samples,bins=bins)
    cdf = np.cumsum(1.0*hist/np.sum(hist))
    bcent = 0.5*(b[1:] + b[:-1])
    
    mid = np.interp(0.5,cdf,bcent)
    low = np.interp(0.32,cdf,bcent)
    high = np.interp(0.68,cdf,bcent)

    return mid, mid - low, high - mid

def compute_MLWA(theta, wave, mlib, agearr, lightidx):
    
    tauV = theta[0]
    weights = theta[1:] * 100.

    redd = np.exp(-1*tauV*(wave[lightidx]/5500.)**(-0.7))
    light_weight = np.mean(mlib[:,lightidx] * redd, axis=1) * weights
    MLWA = np.sum(light_weight * agearr)/np.sum(light_weight)

    return MLWA

def mcombine(theta, wave, mlib, vel):
    
    light_factor = 100.
    tauV = theta[0]
    weights = np.array(theta[1:])

    y = np.sum(mlib * weights[:,None] * light_factor, axis=0)
    
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tauV*klam)
    y *= e_tau_lam

    redwave = wave * (vel/3e5 + 1)
    y = np.interp(wave,redwave,y)

    return y

def lnprob(theta, wave, data, err, mlib, vel):
    
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    
    model = mcombine(theta, wave, mlib, vel)
    return lp - np.sum(((data - model)/err)**2)

def lnprior(theta):
    tauV = theta[0]
    weights = np.array(theta[1:])
    
    if tauV < 0 or tauV > 20:
        return -np.inf
    for w in weights:
        if w < 0 or w > 1e3:
            return -np.inf

    return 0.0

def plot_traces(S,output):
    
    labels = [r'$\tau_V$'] + ['$w_{}$'.format(i) for i in range(S.chain.shape[2]-1)]

    pp = PDF(output)
    subplot = 1
    
    for p in range(S.chain.shape[2]):
        if subplot == 1:
            tax = plt.figure().add_subplot(211)
        else:
            tax = tax.figure.add_subplot(212)
        tax.plot(S.chain[:,:,p].T,color='k',alpha=0.2)
        if labels is not None:
            tax.set_ylabel(labels[p])
        tax.set_xlabel('steps')
        if subplot == 2:
            pp.savefig(tax.figure)
            plt.close(tax.figure)
            subplot = 1
        else:
            subplot = 2

    if subplot == 2:
        pp.savefig(tax.figure)
        plt.close(tax.figure)

    pp.close()

    return
