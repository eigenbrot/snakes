import matplotlib
matplotlib.use('Agg')
import pyfits
from matplotlib.backends.backend_pdf import PdfPages as PDF
import numpy as np
import scipy.ndimage as spnd
import scipy.optimize as spo
import os
import sys
#import lmfit
#import pymc
import emcee
import triangle
#from sklearn.mixture import GMM
import time
import matplotlib.pyplot as plt
plt.ioff()

def do_simple(datafile, errorfile, output, 
              model = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
              plot = True, wavemin = 3750., wavemax = 6800., 
              lightmin = 5450., lightmax = 5550.,
              nsample = 10000, burn = 1000, nwalkers = 256):

    m = pyfits.open(model)[1].data
    
    dhdu = pyfits.open(datafile)[0]
    header = dhdu.header
    data = dhdu.data
    error = pyfits.open(errorfile)[0].data

    try:
        numfibers, wavesize = data.shape
    except ValueError:
        numfibers = 1
        wavesize = data.size
        print 'Found one fiber with length {}'.format(wavesize)

    wave = (np.arange(wavesize) - header['CRPIX1']) * header['CDELT1'] \
           +  header['CRVAL1']
    idx = np.where((wave >= wavemin) & (wave <= wavemax))
    wave = wave[idx]

    vdisp = np.array([493., 589., 691., 796., 966.])/2.355
    size_borders = [19, 43, 62, 87, 109]
    size_switch = 0

    lightidx = np.where((wave >= lightmin) & (wave <= lightmax))[0]
    
    agearr = m['AGE'][0]/1e9
    numages = agearr.size

    f = open(output, 'w')
    f.write('{:11}'.format('# Fiber Num'))
    f.write((numages*'{:9.3f} Gyr').format(*agearr))
    f.write((5*'{:>13}').\
            format('MMWA [Gyr]','MLWA [Gyr]','Tau_V','S/N','Chisq'))
    f.write('\n#\n')

    if plot:
        plotname = output.split('.')[0] + '.pdf'
        pp = PDF(plotname)

    L_sun = 3.826e33 #erg / s
    dist_mpx = 10.062
    flux_factor = 1e17
    tau = 2*np.pi

    for i in [48]:#range(numfibers):

        print 'Doing fiber {}'.format(i+1)
        if numfibers == 1:
            flux = data[idx]*flux_factor
            err = error[idx]*flux_factor
        else:
            flux = data[i,idx][0]*flux_factor
            err = error[i,idx][0]*flux_factor

        if i == size_borders[0]:
            size_switch += 1
            size_borders = size_borders[1:]

        coef, fig, S = superfit(m, wave, flux, err, vdisp[size_switch],
                                plotlabel='Fiber {}'.format(i+1),
                                nsample = nsample, burn = burn)
    
        pp.savefig(fig)
        plt.close(fig)

        SNR = np.sqrt(
            np.sum((flux[lightidx]/err[lightidx])**2)/lightidx.size)

        MMWA = np.sum(agearr*coef['light_frac']*\
                      m['M_REMAINING'][0]/\
                      m['NORM'][0])/\
            np.sum(coef['light_frac']*m['M_REMAINING'][0]/m['NORM'][0])

        redd = np.exp(-1*coef['tauv']*(wave[lightidx]/5500)**(-0.7))
        light_weight = np.mean(m['FLUX'][0][:,lightidx] * redd, axis=1)*\
                       coef['light_frac']

        MLWA = np.sum(light_weight * agearr)/np.sum(light_weight)

        f.write('{:11}'.format(i+1))
        f.write((numages*'{:13.3e}').format(*coef['light_frac']))#/m['NORM'][0]))
        f.write('{:13.7f}{:13.7f}{:13.3f}{:13.3f}{:13.3e}\n'.\
                format(MMWA,MLWA,coef['tauv'],SNR,coef['chisq']))

    pp.close()
    f.close()

    return S

def superfit(model, restwl, flux, err, vdisp, 
             emmaskw = 400.0, plotlabel = '', 
             nsample = 10000, burn = 1000, nwalkers = 256):

    nmodels = model['AGE'][0].size
    npix = restwl.size
    quality = np.ones(npix)

    if np.min(restwl) < np.min(model['WAVE'][0]):
        print "Warning, wavelengths are too small"
    if np.max(restwl) > np.max(model['WAVE'][0]):
        print "Warning, wavelengths are too big"

    bad = np.where((np.isnan(flux) == True) | (np.isnan(err) == True) |\
                   (err == 0))[0]
    quality[bad] = 0

    em= [3726.03, # OII
         3728.82, # OII
         3889.05, # H8
         3869.06, # NeIII
         4101.73, # Hg
         4340.46, # Hd
         4861.33, # Hb
         4959.91, # OIII
         5006.84, # OIII
         5875.67, # He I
         6300.30, # OI
         6548.04, # NII
         6562.82, # Ha
         6583.41, # NII
         6716.44, # SII
         6730.81] # SII
    # bad sky lines
    sk = [5569., 5882.6]

    dz = emmaskw / 3e5
    dzsk = 1000. / 3e5

    for line in em:
        maskout = np.where((restwl > line*(1-dz)) & (restwl < line*(1+dz)))
        quality[maskout] = 0

    for sky in sk:
        maskout = np.where((restwl > sky*(1-dzsk)) & (restwl < sky*(1+dzsk)))
        quality[maskout] = 0

    ok = np.where(quality == 1)[0]
    

    bc03_pix = 70.0
    bc03_vdisp = 75.0

    if vdisp < bc03_vdisp: 
        vdisp_add = 0
    else:
        vdisp_add = np.sqrt(vdisp**2 - bc03_vdisp**2)
    sigma_pix = vdisp_add / bc03_pix

    custom_lib = np.zeros((nmodels, npix))

    for i in range(nmodels):
        cflux = spnd.filters.gaussian_filter1d(model['FLUX'][0][i], sigma_pix)
        custom_lib[i] = np.interp(restwl,model['WAVE'][0],cflux)
    
    #Setup MCMC stuff
    # Find best fit
    nll = lambda *args: -lnprob(*args)
    x0 = np.zeros(nmodels+1)
    x0[1:] += 1e5
    result = spo.minimize(nll, x0, method='Powell', args=(restwl[ok], custom_lib[:,ok], flux[ok], err[ok]))
    #result = spo.fmin(nll, x0, xtol=1e-5, ftol=1e-5, args=(restwl[ok], custom_lib[:,ok], flux[ok], err[ok]))
    #xf = result
    xf = result['x']

    lmy = mcombine(xf, restwl, custom_lib)
    ax = plt.figure().add_subplot(111)
    ax.plot(restwl,flux,'k')
    ax.plot(restwl,lmy,'r')
    ax.figure.savefig('minimized_fit.png')
#    print xf
#    raw_input('minimized fit')

    ndim = nmodels + 1
    p0 = [xf + 1e-1*xf*np.random.randn(ndim) for i in range(nwalkers)]

    S = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(restwl[ok], custom_lib[:,ok], flux[ok], err[ok]))
    #run MCMC
    print 'Burn'
    rows, cols = os.popen('stty size', 'r').read().split()
    cols = int(cols)*2/3
    t1 = time.time()
    tsum = 0
    pos, prob, state = S.run_mcmc(p0,1)
    for i in range(burn-1):
        tt1 = time.time()        
        pos, prob, state = S.run_mcmc(pos,1,rstate0=state,lnprob0=prob)
        tt2 = time.time()
        tsum += tt2 - tt1
        delta = tsum/(i+1)
        sys.stdout.write('\r')
        sys.stdout.write('[{{:<{}}}] in {{:5.2f}}s ({{:5.0f}}s remaining)'.format(cols).\
                         format('='*(int((i+1)*cols/(burn-1))),time.time() - t1, delta*(burn-i)))
        sys.stdout.flush()
    
    print '\nBurned {} steps in {} s'.format(burn,time.time() - t1)
    print 'reset'
    S.reset()
    print 'Run'
    t1 = time.time()
    tsum = 0
    for i in range(nsample):
        tt1 = time.time()
        pos, prob, state = S.run_mcmc(pos,1,rstate0=state,lnprob0=prob)
        tt2 = time.time()
        tsum += tt2 - tt1
        delta = tsum/(i+1)
        sys.stdout.write('\r')
        sys.stdout.write('[{{:<{}}}] in {{:5.2f}}s ({{:5.0f}}s remaining)'.format(cols).\
                         format('='*(int((i+1)*cols/nsample)),time.time() - t1, delta*(nsample-i)))
        sys.stdout.flush()

    print 
    #Collect results
#    fitcoefs = np.mean(S.flatchain,axis=0)
    fitcoefs = np.zeros(nmodels+1)
    for t in range(nmodels+1):
        hist, bins = np.histogram(S.flatchain[:,t],bins=100)
        bins = 0.5*(bins[1:] + bins[:-1])
        fitcoefs[t] = bins[np.argmax(hist)]

    fiterrs = np.std(S.flatchain,axis=0)
#     labels=[r'$\tau_V$'] + ['$w_{{{}}}$'.format(ii+1) for ii in range(nmodels)]
#     fig = triangle.corner(S.flatchain,
#                           labels=labels,
#                           truths=fitcoefs)
#     fig.suptitle(plotlabel)
#     print 'emc.pdf'
#     pS = PDF('emc.pdf')
#     pS.savefig(fig)
#     pS.close()
#     pT = PDF('trace.pdf')
#     print 'trace'
#     for p in range(nmodels+1):
#         axt = plt.figure().add_subplot(111)
#         axt.plot(S.chain[:,:,p].T,color='k',alpha=0.3)
#         axt.set_xlabel('step')
#         axt.set_ylabel(labels[p])
# #        axt.set_ylim(fitcoefs[p] - fiterrs[p]*2, fitcoefs[p] + fiterrs[p]*2)
#         pT.savefig(axt.figure)
#     pT.close()

    t2 = time.time()
    yfit = mcombine(fitcoefs, restwl, custom_lib)
    chisq = np.sum((yfit - flux)**2/err**2)

    print '\tRan {:n} steps in {:4.2f} seconds'.format(nsample,t2-t1)
    print '\tComputed chisq = {}'.format(chisq)

    coefs = {'tauv': fitcoefs[0], 'light_frac': fitcoefs[1:],
             'model_age': model['AGE'][0], 'chisq': chisq}
    
    fig = plt.figure(figsize=(10,8))
    pax = fig.add_axes([0.1,0.3,0.85,0.64])
    rax = fig.add_axes([0.1,0.15,0.85,0.15])
    
    ymax = np.max(yfit) * 1.1
    xmin = np.min(restwl) * 0.98
    xmax = np.max(restwl) * 1.02
    pax.set_xticklabels([])
    pax.set_ylabel('Flux')
    pax.set_title(plotlabel)
    pax.set_xlim(xmin,xmax)
    pax.set_ylim(-1,ymax)

    galfit = np.ones(flux.size) + np.nan
    galfit[ok] = flux[ok]
    masked = np.copy(flux)
    masked[ok] = np.nan

    pax.plot(restwl, masked, 'g')
    pax.plot(restwl, galfit, 'k')
    pax.plot(restwl, yfit, 'r')
    pax.fill_between(restwl, flux-err, flux+err, color='k', alpha=0.1)
    pax.text(0.1,0.9,'Tau V = {:4.2f}'.format(fitcoefs[0]),
             transform=pax.transAxes)
    
    lidx = np.where((restwl >= 5450.) & (restwl <= 5550.))
    for i in range(nmodels):
        yi = fitcoefs[i+1] * custom_lib[i,:] *\
             np.exp(-1*fitcoefs[0]*(restwl/5500.)**(-0.7))
        pax.plot(restwl,yi,'b--')
        pax.text(0.1, 0.87 - 0.022*i, 'f_{:02} = {:10.3e}'.\
                 format(i,np.mean(yi[lidx])),
                 transform=pax.transAxes,fontsize=9)
        pax.text(0.8, 0.45 - 0.022*i,'{:6.3f}: {:10.3e}'.\
                 format(model['AGE'][0][i]/1e9, 
                        coefs['light_frac'][i]/model['NORM'][0][i]),
                 transform=pax.transAxes,fontsize=9)

    rax.set_ylabel('Residuals/error')
    rax.set_xlabel('Wavelength [Angstroms]')
    rax.set_xlim(xmin,xmax)
    rax.set_ylim(-5,5)
    rax.plot(restwl, (galfit - yfit)/err,'k')

    return coefs, fig, S

def plot_emcee(sampler,outputprefix,labels=None,truths=None):

    if labels is None:
        labels = [r'$\tau_V$'] + ['$w_{{{}}}$'.format(ii+1) for ii in range(sampler.flatchain.shape[1])]
    if truths is None:
#        truths = np.mean(sampler.flatchain,axis=0)
        fitcoefs = np.zeros(sampler.flatchain.shape[1])
        for t in range(sampler.flatchain.shape[1]):
            hist, bins = np.histogram(sampler.flatchain[:,t],bins=50)
            bins = 0.5*(bins[1:] + bins[:-1])
            fitcoefs[t] = bins[np.argmax(hist)]
        truths = fitcoefs

    tracename = '{}_trace.pdf'.format(outputprefix)
    tracepp = PDF(tracename)
    subplot = 1
    for p in range(sampler.chain.shape[2]):
        if subplot == 1:
            tax = plt.figure().add_subplot(211)
        else:
            tax = tax.figure.add_subplot(212)
        tax.plot(sampler.chain[:,:,p].T,color='k',alpha=0.2)
        if labels is not None:
            tax.set_ylabel(labels[p])
        tax.set_xlabel('steps')
        if subplot == 2:
            tracepp.savefig(tax.figure)
            plt.close(tax.figure)
            subplot = 1
        else:
            subplot = 2

    if subplot == 2:
        tracepp.savefig(tax.figure)
        plt.close(tax.figure)

    tracepp.close()

    triname = '{}_tri.pdf'.format(outputprefix)
    tripp = PDF(triname)
    trifig = triangle.corner(sampler.flatchain,labels=labels,truths=truths)
    tripp.savefig(trifig)
    plt.close(trifig)
    tripp.close()

    return

def plot_age(sampler, datafile, sift=1,
             model = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
             lightmin = 5450., lightmax = 5550., wavemin = 3750., wavemax = 6800.):

    m = pyfits.open(model)[1].data
    agearr = m['AGE'][0]/1e9
    
    dhdu = pyfits.open(datafile)[0]
    header = dhdu.header
    data = dhdu.data

    try:
        numfibers, wavesize = data.shape
    except ValueError:
        numfibers = 1
        wavesize = data.size
        print 'Found one fiber with length {}'.format(wavesize)

    wave = (np.arange(wavesize) - header['CRPIX1']) * header['CDELT1'] \
           +  header['CRVAL1']
    idx = np.where((wave >= wavemin) & (wave <= wavemax))
    wave = wave[idx]
    
    lightidx = np.where((wave >= lightmin) & (wave <= lightmax))[0]

    redd = np.exp(-1*sampler.flatchain[::sift,0][:,None]*(wave[lightidx]/5500)**(-0.7))
    # print redd.shape, redd[:,None].shape
    # print m['FLUX'][0][:,lightidx][None,:].shape
    # print (m['FLUX'][0][:,lightidx][None,:] * redd[:,None]).shape
    # print np.mean(m['FLUX'][0][:,lightidx][None,:] * redd[:,None], axis=2).shape
    # print sampler.flatchain[:,1:].shape
    light_weight = np.mean(m['FLUX'][0][:,lightidx][None,:] * redd[:,None], axis=2)*\
                   sampler.flatchain[::sift,1:]

    # print light_weight.shape
    # print (light_weight * agearr).shape
    MLWA = np.sum(light_weight * agearr,axis=1)/np.sum(light_weight,axis=1)

    return MLWA

def mcombine(theta, wave, mlib):

    tauv = theta[0]
    weights = np.array(theta[1:])
    
    y = np.sum(mlib * weights[:,None],axis=0)
    
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tauv*klam)
    y *= e_tau_lam

    return y

def lnprob(theta, wave, mlib, data, err):

    lp = lnprior(theta)

    model = mcombine(theta, wave, mlib)
    # val = lp - np.sum(((data - model)/err)**2)
    
    # if np.isnan(val):
    #     return -np.inf

    return lp - np.sum(((data - model)/err)**2)

def lnprior(theta):

    supersmall = -np.inf
    tauv = theta[0]
    weights = np.array(theta[1:])
    
    if tauv < -5.0 or tauv > 5.0:
        return supersmall

    for w in weights:
        if w < 0.0 or w > 1e11:
            return supersmall

    return 0.0
