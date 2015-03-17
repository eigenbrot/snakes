import pyfits
from matplotlib.backends.backend_pdf import PdfPages as PDF
import numpy as np
import scipy.ndimage as spnd
#import scipy.optimize as spo
#import lmfit
import pymc
#from sklearn.mixture import GMM
import matplotlib.pyplot as plt
import time
plt.ioff()

def do_simple(datafile, errorfile, output, 
              model = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
              plot = True, wavemin = 3750., wavemax = 6800., 
              lightmin = 5450., lightmax = 5550.,
              nsample = 10000, burn = 1000):

    m = pyfits.open(model)[1].data
    
    dhdu = pyfits.open(datafile)[0]
    header = dhdu.header
    data = dhdu.data
    error = pyfits.open(errorfile)[0].data

    numfibers, wavesize = data.shape

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

    for i in [77]:#range(numfibers):

        print 'Doing fiber {}'.format(i+1)
        flux = data[i,idx][0]*flux_factor
        err = error[i,idx][0]*flux_factor

        if i == size_borders[0]:
            size_switch += 1
            size_borders = size_borders[1:]

        coef, fig = superfit(m, wave, flux, err, vdisp[size_switch],
                             plotlabel='Fiber {}'.format(i+1),
                             nsample = nsample, burn = burn)
    
        pp.savefig(fig)
#        plt.close(fig)

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
        f.write('{:13.7f}{:13.7f}{:13.3f}{:13.3f}{:13.3f}\n'.\
                format(MMWA,MLWA,coef['tauv'],SNR,coef['chisq']))

    pp.close()
    f.close()

    return 

def superfit(model, restwl, flux, err, vdisp, 
             emmaskw = 400.0, plotlabel = '', 
             nsample = 10000, burn = 1000):

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
    
    #define the priors
    tauvprior = pymc.Uniform('tauv',0.0,5.0)
    params = {'wave': restwl[ok],
              'mlib': custom_lib[:,ok],
              'tauv': tauvprior}
    MCmodel = {'tauv': tauvprior}
    for p in range(nmodels):
        parname = 'w_{}'.format(p)
        prior = pymc.Uniform(parname,0,1e5)
        params[parname] = prior
        MCmodel[parname] = prior

    galaxy_model = pymc.Deterministic(eval = mcombine,
                                      name = 'Galaxy model',
                                      doc = 'Galaxy model from SSPs',
                                      parents = params,
                                      trace = True,
                                      verbose = 0,
                                      dtype = float,
                                      plot = True)

    galaxy = pymc.Normal('galaxy',
                         mu = galaxy_model,
                         tau = 1/(err[ok]**2),
                         value = flux[ok],
                         observed = True)
    MCmodel['galaxy'] = galaxy

    #run MCMC
    t1 = time.time()
    S = pymc.MCMC(MCmodel)
    S.sample(nsample,burn)

    #collect results
    fitcoefs = np.zeros(nmodels+1)
    fiterrs = np.zeros(nmodels+1)
    
    for i, k in enumerate(['tauv'] + \
                          ['w_{}'.format(p) for p in range(nmodels)]):
        trace = S.trace(k)[:]
        # print 'fitting ', k
        # amp, mu, std = fit_pdf(trace)
        # px = np.linspace(-1*trace.max(),trace.max(),1000)
        # pp = amp*np.exp(-1*((px - mu)/(2*std))**2)
        mu = np.mean(trace)
        # hist, bins = np.histogram(trace,bins=nsample/100)
        # bins = 0.5*(bins[1:] + bins[:-1])
        # hist = spnd.filters.gaussian_filter1d(hist,3)
        # mu = bins[np.argmax(hist)]
        std = np.std(trace)
        fitcoefs[i] = mu
        fiterrs[i] = std
        params[k] = mu
#         tax = plt.figure().add_subplot(111)
#         tax.set_xlabel(k)
#         tax.set_ylabel('N')
#         #tax.plot(bins,hist,color='k')
#         tax.hist(trace,bins=nsample/100,histtype='step',color='k')
#         tax.axvline(x=mu)
# #        tax.plot(px,pp)
# #        tax.set_xlim(-100,trace.max()*1.1)
#         tax.figure.show()
        
    t2 = time.time()
    pymc.Matplot.plot(S)
    params['wave'] = restwl
    params['mlib'] = custom_lib
    yfit = mcombine(**params)
    chisq = np.sum((yfit - flux)**2/err**2)

    print '\tRan {:n} steps in {:4.2f} seconds'.format(nsample,t2-t1)
    print '\tComputed chisq = {}'.format(chisq)

    fitcoefs = np.zeros(nmodels+1)
    fitcoefs[0] = params['tauv']
    for p in range(10):
        fitcoefs[p+1] = params['w_{}'.format(p)]

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
    pax.set_ylim(1,ymax)

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
#    plt.show()

#    return S
    return coefs, fig

def mcombine(**pardict):

    weights = np.zeros(len(pardict) - 3)
    for k in pardict.keys():
        if 'w_' in k:
            # if pardict[k] < 0:
            #     return np.ones(pardict['wave'].size)*9e77
            n = int(k.split('_')[1])
            weights[n] = pardict[k]
    
    y = np.sum(pardict['mlib'] * weights[:,None],axis=0)
    
    klam = (pardict['wave'] / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*pardict['tauv']*klam)
    y *= e_tau_lam

    return y

def fit_pdf(pdf):

    prob, bins = np.histogram(pdf,bins=1000)
    x = 0.5*(bins[1:] + bins[:-1])
    mx = np.linspace(-1*x.max(),x.max(),1000)

    params = lmfit.Parameters()
    params.add_many(('mu', np.mean(pdf), True, 0, 1e5, None),
                ('sigma', np.std(pdf), True, 0, 200, None),
                ('amp', 1, True, 0, None, None))

    status = lmfit.minimize(fit_pdf_func, params,
                            method='leastsq',
                            ftol=1e-9,#maxfev=1000000000,
                            args=(x,prob,mx))
    
    # fit = params['amp'].value*np.exp(
    #     -1*((x - params['mu'].value)/(2*params['sigma'].value))**2)

    return params['amp'].value, params['mu'].value, params['sigma'].value

def fit_pdf_func(params,x,y,mx):

    by = params['amp'].value*np.exp(
        -1*((mx - params['mu'].value)/(2*params['sigma'].value))**2)

    by[mx < 0] = 0.0
    my = np.interp(x,mx,by)

    return y - my
