import time
from glob import glob
import numpy as np
import pyfits
import emcee
import triangle
import matplotlib.pyplot as plt
import h5py
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def main(datafile, errorfile, location, coeffile, outputpre,
         model='/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_allZ_vardisp.fits',
         plot=True, wavemin=3800., wavemax=6800., fitregion=[3850.,6650.],
         lightmin=5450., lightmax=5550., emmaskw=500., fitaps=None, justMLWA=False,
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

    wave = wave[idx]
    lightidx = np.where((wave >= lightmin) & (wave <= lightmax))[0]

    agearr = m['AGE']/1e9
    numages = agearr.shape[0]

    # f = open(outputpre,'w')
    # f.write('{:11}'.format('# Fiber Num'))
    # f.write((numages*'{:9.3f} Gyr').format(*agearr))
    # f.write('{:>10}{:>10}'.format('Tau_V, Vsys'))
    # f.write('\n#\n')

    if plot:
        plotname = outputpre + '_tri.pdf'
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

    yfitfile = outputpre + '.emceefit.fits'
    yfitarr = np.zeros((numfibers, wave.size))

    h5file = h5py.File(outputpre + '_emcee.h5','a')

    if fitaps is None:
        fitaps = range(numfibers)
    else:
        #because ap numbers start at 1
        fitaps = [i-1 for i in fitaps]

    for i in fitaps:

        print "Doing fiber {}".format(i+1)
        flux = data[i,idx]*flux_factor
        err = error[i,idx]*flux_factor

        vdidx = np.where(sizeidx == fiber_radii[i])[0][0]
        LMcoefs = coef_arr[i]
        
        if justMLWA:
            print 'computing MLWA from prevchain'
            grp = h5file['Ap{}'.format(i+1)]
            chain = grp['chain']
            flatchain = np.reshape(chain,(chain.shape[0]*chain.shape[1],chain.shape[2]))
            MLWA_S = EMfit(m, wave, flux, err, vdidx, LMcoefs,
                           fitregion=fitregion,emmaskw=emmaskw,
                           lightidx=lightidx,nsample=nsample,
                           prevchain= flatchain,
                           burn=burn,nwalkers=nwalkers,threads=threads)
            grp.create_dataset('MLWA_S2',data=MLWA_S,compression='gzip',compression_opts=0)
            h5file.flush()
            continue

        else:
            MCcoefs, yfit, S, MLWA_S = EMfit(m, wave, flux, err, vdidx, LMcoefs,
                                             fitregion=fitregion,emmaskw=emmaskw,
                                             lightidx=lightidx,nsample=nsample,
                                             burn=burn,nwalkers=nwalkers,threads=threads)

        outputarr[i] = MCcoefs
        yfitarr[i,:] = yfit/flux_factor

        grp = h5file.create_group('Ap{}'.format(i+1))
        grp.create_dataset('chain',data=S.chain,compression='gzip',compression_opts=9)
        grp.create_dataset('lnprob',data=S.lnprobability,compression='gzip',compression_opts=9)
        grp.create_dataset('MLWA_S',data=MLWA_S,compression='gzip',compression_opts=9)
        h5file.flush()

        if plot:
            try:
                fig = triangle.corner(S.flatchain*np.r_[1,np.ones(numages)*100.],
                                      labels=[r'$\tau_V$'] + ['$w_{{{}}}$'.format(i) for i in range(numages)],
                                      truths=np.r_[MCcoefs['TAUV'],MCcoefs['LIGHT_FRAC']])
                                      #extents=[(0,3)] + numages*[(0,200)])
                pp.savefig(fig)
                plt.close(fig)
            except:
                pass
    if plot:
        pp.close()

    h5file.close()

    if justMLWA:
        return

    pyfits.BinTableHDU(outputarr).writeto(outputpre + '.emceecoef.fits',clobber=True)
    fithdu = pyfits.PrimaryHDU(yfitarr)
    fithdu.header.update('CDELT1',CDELT)
    fithdu.header.update('CRPIX1',1)
    fithdu.header.update('CRVAL1',wavemin)
    fithdu.writeto(yfitfile,clobber=True)

    return MCcoefs, S, MLWA_S

def EMfit(model, wave, flux, err, vdidx, LMcoefs, fitregion=[3850.,6650.],
          emmaskw=500., lightidx=None, threads=1, prevchain=None,
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

    if prevchain is not None:
        MLWA_samples = np.zeros(prevchain.shape[0])
        for s in range(prevchain.shape[0]):
            if s % 1000 == 0:
                print s*1.0/prevchain.shape[0]
            MLWA_samples[s] = compute_MLWA(prevchain[s,:],
                                           fitwave, fitlib,
                                           model['AGE'][:,0]/1e9, lightidx)
        return MLWA_samples

    #Set up the emcee stuff
    ndim = nmodels + 1
    LMtheta = np.r_[LMcoefs['TAUV'],LMcoefs['LIGHT_FRAC']/light_factor]
    print LMtheta
    p0 = [LMtheta + 1e1*LMtheta*np.random.ranf(ndim) for i in range(nwalkers)]
    S = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=threads, args=(fitwave,
                                                                             fitflux,
                                                                             fiterr,
                                                                             fitlib,
                                                                             LMcoefs['VSYS'],
                                                                             model['AGE'][:,0]/1e9, 
                                                                             lightidx))

    #run
    print 'Burning...'
    t1 = time.time()
    pos, prob, state, _ = S.run_mcmc(p0,burn)
    print 'Burned {} steps in {} s'.format(burn,time.time() - t1)
    
    S.reset()
    print 'Running...'
    t1 = time.time()
    pos, prob, state, _ = S.run_mcmc(pos,nsample)
    print 'Ran {} steps in {} s'.format(nsample,time.time() - t1)

    maxidx = np.argmax(S.flatlnprobability)
    EMtheta = S.flatchain[maxidx]
    print '###', EMtheta
    stds = np.mean(S.flatchain,axis=0)

    cents, lowCI, highCI = compute_theta_CI(S.flatchain)
    print '$$$', cents
    print lowCI
    print highCI

    MLWA_samples = np.vstack(S.blobs).T
    print MLWA_samples.shape, S.chain.shape
    flat_MLWA = np.reshape(MLWA_samples,(nwalkers*nsample))
    print flat_MLWA.shape

    MLWA, MLWA_L, MLWA_H  = compute_MLWA_CI(flat_MLWA)
    print MLWA, MLWA_L, MLWA_H

    yfit = mcombine(EMtheta,wave,custom_lib,LMcoefs['VSYS'])
    chisq, best_MLWA = lnprob(EMtheta, wave, flux, err, custom_lib, LMcoefs['VSYS'],
                              model['AGE'][:,0]/1e9,lightidx)
    print best_MLWA

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

    return fitcoefs[0], yfit, S, MLWA_samples

def compute_theta_CI(flatchain, bins=20):

    numpar = flatchain.shape[1]
    lowCI = np.zeros(numpar)
    highCI = np.zeros(numpar)
    cent = np.zeros(numpar)

    for p in range(numpar):
        hist, b = np.histogram(flatchain[:,p],bins=bins)
        cdf = np.cumsum(1.0*hist/np.sum(hist))
        bcent = 0.5*(b[1:] + b[:-1])
        cent[p] = np.interp(0.5,cdf,bcent)
        lowCI[p] = np.interp(0.32,cdf,bcent)
        highCI[p] = np.interp(0.68,cdf,bcent)

    return cent, cent-lowCI, highCI-cent

def combine_S_CI(Slist):
    
    #find the best
    bestln = -np.inf
    bestS = None
    for S in Slist:
        if np.max(S.flatlnprobability) > bestln:
            bestln = np.max(S.flatlnprobability)
            bestS = S
    
    #Compute the ln prob CI
    hist, b = np.histogram(bestS.flatlnprobability,bins=100)
    cdf = np.cumsum(1.0*hist/np.sum(hist))
    bcent = 0.5*(b[1:] + b[:-1])
    limit = np.interp(0.68,cdf,bcent)

    newlnprob = np.array([])
    newflatchain = np.zeros(bestS.flatchain.shape[1])

    for i, S in enumerate(Slist):
        idx = np.where(S.flatlnprobability > limit)[0]
        if idx.size > 0:
            print idx.shape
            print i
            newlnprob = np.r_[newlnprob, S.flatlnprobability[idx]]
            newflatchain = np.vstack((newflatchain,S.flatchain[idx,:]))

    return newlnprob, newflatchain[1:,:]

def compute_MLWA_CI(MLWA_samples, bins=200):

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

def lnprob(theta, wave, data, err, mlib, vel, agearr, lightidx):
    
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, 0
    
    model = mcombine(theta, wave, mlib, vel)
    MLWA = compute_MLWA(theta, wave, mlib, agearr, lightidx)

    return lp - np.sum(((data - model)/err)**2), MLWA

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

def plot_distributions(file_list, apnum, distname, bins=100, lims=[0,10]):

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(distname)
    ax.set_ylabel('N')

    for f in file_list:
        
        ff= h5py.File(f,'r')
        d = ff['Ap{}'.format(apnum)][distname]
        if distname == 'lnprob':
            d = np.reshape(d,(d.shape[0]*d.shape[1]))
            d /= -1*(1334 - 1 - ff['Ap{}'.format(apnum)]['chain'].shape[2])
            ax.set_xlabel(r'$\chi_{\nu}^2$')
            d = d[np.where((d > lims[0]) & (d < lims[1]))]

        ax.hist(d,bins=bins,histtype='stepfilled',alpha=0.5,label=f)
        ff.close()

    ax.legend(loc=0)

    return ax

def reject_metal(folder_list, pointing, apnum, lims=[-1,10]):

    ax = plt.figure().add_subplot(111)
    ax.set_title('P{}.{}\n{}'.format(pointing,apnum,time.asctime()))
    ax.set_xlabel(r'$\Delta\chi_{\nu}^2$')
    ax.set_ylabel('N')
    ax.set_xscale('log')

    file_list = [glob('{}/*P{}*.h5'.format(f,pointing))[0] for f in folder_list]

    dlist = []

    #first, find the minimum
    minchi = np.inf
    for f in file_list:
        ff = h5py.File(f,'r')
        d = ff['Ap{}'.format(apnum)]['lnprob']
        d = np.reshape(d,(d.shape[0]*d.shape[1]))
        d /= -1*(1334 - 1 - ff['Ap{}'.format(apnum)]['chain'].shape[2])
        dlist.append(d)
        if np.min(d) < minchi:
            minchi = np.min(d)

    for i, d in enumerate(dlist):
        d -= minchi
        d = d[np.where((d > lims[0]) & (d < lims[1]))]
        ax.hist(d,bins=200,normed=True,histtype='stepfilled',
                alpha=0.5,label=folder_list[i])

    ax.axvline(1,ls=':',alpha=0.6,color='k')
    ax.legend(loc=0)

    return ax
    
