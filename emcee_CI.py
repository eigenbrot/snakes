import time
from glob import glob
import numpy as np
import pyfits
import emcee
import triangle
import matplotlib.pyplot as plt
import h5py
import gc
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
                                                   'MLWZ','MLWZ_L','MLWZ_H',
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
                                                     'f4','f4','f4',
                                                     'f4',
                                                     'f4','f4','f4','f4']})

    yfitfile = outputpre + '.emceefit.fits'
    yfitarr = np.zeros((numfibers, wave.size))
    h5name = outputpre + '_emcee.h5'

    if fitaps is None:
        fitaps = range(numfibers)
    else:
        #because ap numbers start at 1
        fitaps = [i-1 for i in fitaps]

    for i in fitaps:

        #First, check if ap has already been run
        h5file = h5py.File(h5name,'a')
        if 'Ap{}'.format(i+1) in h5file.keys():
            print "Fiber {} has already been run, moving on...".format(i+1)
            h5file.close()
            continue
            
        #Close this file often to keep memory free
        h5file.close()

        print "Starting fiber {} on {}".format(i+1,time.asctime())
        flux = data[i,idx]*flux_factor
        err = error[i,idx]*flux_factor

        vdidx = np.where(sizeidx == fiber_radii[i])[0][0]
        LMcoefs = coef_arr[i]
        
        if justMLWA:
            print 'computing MLWA from prevchain'
            h5file = h5py.File(h5name,'a')
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
            h5file.close()
            continue

        else:
            MCcoefs, yfit, S, MLWA_S, MLWZ_S = EMfit(m, wave, flux, err, vdidx, LMcoefs,
                                                     fitregion=fitregion,emmaskw=emmaskw,
                                                     lightidx=lightidx,nsample=nsample,
                                                     burn=burn,nwalkers=nwalkers,threads=threads)

        outputarr[i] = MCcoefs
        yfitarr[i,:] = yfit/flux_factor

        print "writing hdf5..."
        h5file = h5py.File(h5name,'a')
        grp = h5file.create_group('Ap{}'.format(i+1))
        grp.create_dataset('chain',data=S.chain,compression='gzip',compression_opts=9)
        grp.create_dataset('lnprob',data=S.lnprobability,compression='gzip',compression_opts=9)
        grp.create_dataset('MLWA_S',data=MLWA_S,compression='gzip',compression_opts=9)
        grp.create_dataset('MLWZ_S',data=MLWZ_S,compression='gzip',compression_opts=9)
        h5file.flush()
        h5file.close()
        print "clearing memory"
        del grp
        del MLWA_S
        del MLWZ_S
        del S
        print gc.collect()

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

    try:
        h5file.close()
        print "H5 file was not closed for some reason"
    except Exception as e:
        #File is already closed
        print e
        print "Probably because you deleted it"
        pass

    if justMLWA:
        return

    pyfits.BinTableHDU(outputarr).writeto(outputpre + '.emceecoef.fits',clobber=True)
    fithdu = pyfits.PrimaryHDU(yfitarr)
    fithdu.header.update('CDELT1',CDELT)
    fithdu.header.update('CRPIX1',1)
    fithdu.header.update('CRVAL1',wavemin)
    fithdu.writeto(yfitfile,clobber=True)

    return

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
    
    if threads > 1:
        print "Creating pool with {} processes".format(threads)
        pool = emcee.interruptible_pool.InterruptiblePool(processes=threads)
    else:
        pool = None

    S = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=(fitwave,
                                                                       fitflux,
                                                                       fiterr,
                                                                       fitlib,
                                                                       LMcoefs['VSYS'],
                                                                       model['AGE'][:,0]/1e9,
                                                                       model['Z'][:,0],
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
    # print '$$$', cents
    # print lowCI
    # print highCI

    blobstack = np.dstack(S.blobs)
    MLWA_samples = blobstack[:,0,:]
    MLWZ_samples = blobstack[:,1,:]
    flat_MLWA = np.reshape(MLWA_samples,(nwalkers*nsample))
    flat_MLWZ = np.reshape(MLWZ_samples,(nwalkers*nsample))

    MLWA, MLWA_L, MLWA_H  = compute_MLW_CI(flat_MLWA)

    MLWZ, MLWZ_L, MLWZ_H  = compute_MLW_CI(flat_MLWZ)

    yfit = mcombine(EMtheta,wave,custom_lib,LMcoefs['VSYS'])
    chisq, best_blob  = lnprob(EMtheta, wave, flux, err, custom_lib, LMcoefs['VSYS'],
                               model['AGE'][:,0]/1e9,model['Z'][:,0],lightidx)
    best_MLWA = best_blob[0]
    best_MLWZ = best_blob[1]

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
                                          'MLWZ','MLWZ_L','MLWZ_H',
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
                                           'f4','f4','f4',
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
    fitcoefs['MLWZ_L'] = MLWZ_L
    fitcoefs['MLWZ_H'] = MLWZ_H    
    fitcoefs['MMWZ'] = LMcoefs['MMWZ']
    fitcoefs['chisq'] = chisq
    fitcoefs['redchi'] = LMcoefs['redchi']
    fitcoefs['bluechi'] = LMcoefs['bluechi']
    fitcoefs['hkchi'] = LMcoefs['hkchi']

    if threads > 1:
        print "cleaning up pool..."
        pool.close()
        pool.join()
        print "done"

    return fitcoefs[0], yfit, S, MLWA_samples, MLWZ_samples

def compute_theta_CI(flatchain, bins=20):

    numpar = flatchain.shape[1]
    lowCI = np.zeros(numpar)
    highCI = np.zeros(numpar)
    cent = np.zeros(numpar)

    for p in range(numpar):
        sortchain = np.sort(flatchain[:,p])
        cdf = np.arange(sortchain.size*1.0)/(sortchain.size - 1)
        cent[p] = np.interp(0.5,cdf,sortchain)
        lowCI[p] = np.interp(0.32,cdf,sortchain)
        highCI[p] = np.interp(0.68,cdf,sortchain)

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

def compute_MLW_CI(MLW_samples):

    fidx = np.isfinite(MLW_samples)
    MLW_samples = MLW_samples[fidx]
    
    sS = np.sort(MLW_samples)
    cdf = np.arange(sS.size*1.0)/(sS.size - 1)

    mid = np.interp(0.5,cdf,sS)
    low = np.interp(0.32,cdf,sS)
    high = np.interp(0.68,cdf,sS)

    return mid, mid - low, high - mid

def compute_MLWA(theta, wave, mlib, agearr, lightidx):
    
    tauV = theta[0]
    weights = theta[1:] * 100.

    redd = np.exp(-1*tauV*(wave[lightidx]/5500.)**(-0.7))
    light_weight = np.mean(mlib[:,lightidx] * redd, axis=1) * weights
    MLWA = np.sum(light_weight * agearr)/np.sum(light_weight)

    return MLWA

def compute_MLWZ(theta, wave, mlib, Zarr, lightidx):
    
    tauV = theta[0]
    weights = theta[1:] * 100.

    redd = np.exp(-1*tauV*(wave[lightidx]/5500.)**(-0.7))
    light_weight = np.mean(mlib[:,lightidx] * redd, axis=1) * weights
    MLWZ = np.sum(light_weight * Zarr)/np.sum(light_weight)

    return MLWZ

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

def lnprob(theta, wave, data, err, mlib, vel, agearr, Zarr, lightidx):
    
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf, [0, 0]
    
    model = mcombine(theta, wave, mlib, vel)
    MLWA = compute_MLWA(theta, wave, mlib, agearr, lightidx)
    MLWZ = compute_MLWZ(theta, wave, mlib, Zarr, lightidx)

    return lp - np.sum(((data - model)/err)**2), [MLWA, MLWZ]

def lnprior(theta):
    tauV = theta[0]
    weights = np.array(theta[1:])
    
    if tauV < 0 or tauV > 20:
        return -np.inf
    for w in weights:
        if w < 0 or w > 1e3:
            return -np.inf

    return 0.0

def make_coeffile(pointing,suffix='allZ'):
    #I'll assume all the names are defaults here, because I'm lazy
    
    prevcoef = 'NGC_891_P{}_bin30_allz2.coef.fits'.format(pointing)
    emceefile = 'em_P{}_{}_emcee.h5'.format(pointing,suffix)
    newcoef = 'em_P{}_{}.emceecoef.fits'.format(pointing,suffix)

    print prevcoef
    print emceefile
    print newcoef

    f = h5py.File(emceefile,'r')
    numap = len(f.keys())
    numages = f['Ap1']['chain'].shape[2] - 1
    
    LMdata = pyfits.open(prevcoef)[1].data

    light_factor = 100.

    dtype = {'names':['VSYS',
                      'FIXEDVBOOL',
                      'TAUV',
                      'TAUV_ERR',
                      'emmaskw',
                      'LIGHT_FRAC',
                      'LIGHT_FRAC_ERR',
                      'SNR',
                      'MLWA','MLWA_L','MLWA_H',
                      'MMWA',
                      'MLWZ','MLWZ_L','MLWZ_H',
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
                        'f4','f4','f4',
                        'f4',
                        'f4','f4','f4','f4']}

    outputarr = np.zeros(numap,dtype=dtype)

    for i in range(numap):
        print i+1
        fitcoefs = np.zeros(1,dtype=dtype)
        LMcoefs = LMdata[i]

        grp = f['Ap{}'.format(i+1)]
        print "\treshaping..."
        flatlnprob = np.reshape(grp['lnprob'],grp['lnprob'].size)
        flatchain = np.reshape(grp['chain'],
                               (grp['chain'].shape[0]*grp['chain'].shape[1],
                                grp['chain'].shape[2]))
        
        flat_MLWA = np.reshape(grp['MLWA_S'],grp['MLWA_S'].size)
        flat_MLWZ = np.reshape(grp['MLWZ_S'],grp['MLWZ_S'].size)

        maxidx = np.argmax(flatlnprob)
        EMtheta = flatchain[maxidx]
        stds = np.mean(flatchain,axis=0)
        
        print "\tcomputing CI's..."
        cents, lowCI, highCI = compute_theta_CI(flatchain)
        MLWA, MLWA_L, MLWA_H  = compute_MLW_CI(flat_MLWA)
        MLWZ, MLWZ_L, MLWZ_H  = compute_MLW_CI(flat_MLWZ)

        fitcoefs['VSYS'] = LMcoefs['VSYS']
        fitcoefs['emmaskw'] = LMcoefs['emmaskw']
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
        fitcoefs['MLWZ_L'] = MLWZ_L
        fitcoefs['MLWZ_H'] = MLWZ_H    
        fitcoefs['MMWZ'] = LMcoefs['MMWZ']
        fitcoefs['chisq'] = -1*np.max(flatlnprob)
        fitcoefs['redchi'] = LMcoefs['redchi']
        fitcoefs['bluechi'] = LMcoefs['bluechi']
        fitcoefs['hkchi'] = LMcoefs['hkchi']

        outputarr[i] = fitcoefs[0]

    pyfits.BinTableHDU(outputarr).writeto(newcoef,clobber=True)

    return

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
            #d /= -1*(1428 - 1 - ff['Ap{}'.format(apnum)]['chain'].shape[2])
            ax.set_xlabel(r'$\chi_{\nu}^2$')
            d = d[np.where((d > lims[0]) & (d < lims[1]))]

        ax.hist(d,bins=bins,histtype='stepfilled',alpha=0.5,label=f)
        ff.close()

    ax.legend(loc=0)

    return ax

def reject_metal(folder_list, pointing, apnum, lims=[-1,10], thin=1000):

    clist = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d']
    ax = plt.figure().add_subplot(111)
    ax.set_title('P{}.{}\n{}'.format(pointing,apnum,time.asctime()))
    ax.set_ylabel('lnprob')
    ax.set_xlabel('MLWA')

    file_list = [glob('{}/*P{}*.h5'.format(f,pointing))[0] for f in folder_list]

    dlist = []
    MLWAlist = []

    #first, find the minimum
    maxprob = -np.inf
    maxi = 0
    for i, f in enumerate(file_list):
        print f
        ff = h5py.File(f,'r')
        d = ff['Ap{}'.format(apnum)]['lnprob']
        d = np.reshape(d,(d.shape[0]*d.shape[1]))
        d /= (1428 - 1 - ff['Ap{}'.format(apnum)]['chain'].shape[2])
        dlist.append(d)
        MLWA = ff['Ap{}'.format(apnum)]['MLWA_S']
        MLWA = np.reshape(MLWA,(MLWA.size))
        MLWAlist.append(MLWA)
        if np.max(d) > maxprob:
            maxi = i
            maxprob = np.max(d)
            bestMLWA = MLWA[np.argmax(d)]
    
    #Compute CI
    print 'Found most likely model in', file_list[maxi]
    best_dist, bbcent, bcdf = clean_cdf(dlist[maxi])
    CI_cut = np.interp(0.68,bcdf,bbcent)
    ax.axhline(CI_cut,ls=':',alpha=0.6,color='k')

    goodMLWA = np.array([])

    for i in range(len(dlist)):
        dist = dlist[i]
    
        gidx = np.where(dist > CI_cut)[0]
        if gidx.size > 0:
            goodMLWA = np.r_[goodMLWA, MLWAlist[i][gidx]]

        dist, bcent, cdf = clean_cdf(dist)

        limit68 = np.interp(0.68,cdf,bcent)
        limit10 = np.interp(0.1,cdf,bcent)
        limit1 = np.interp(0.01,cdf,bcent)

        d = dlist[i][::thin]
        M = MLWAlist[i][::thin]
        idx = np.where(d > limit1)[0]
        if idx.size > 0:
            ax.scatter(M[idx],d[idx],alpha=0.5,color=clist[i],edgecolors='none', label=folder_list[i])

        # ax.hist(d,bins=200,normed=True,histtype='stepfilled',
        #         alpha=0.5,label=folder_list[i])


    ax.legend(loc=0)
    ax.axvline(np.max(goodMLWA),ls='--',color='k',alpha=0.4)
    ax.axvline(np.min(goodMLWA),ls='--',color='k',alpha=0.4)
    ax.axvline(bestMLWA,ls='-.',color='k',alpha=0.4)

    print 'MLWA = {:4.3} + {:4.3} - {:4.3}'.format(bestMLWA,
                                                   np.max(goodMLWA)-bestMLWA,
                                                   bestMLWA-np.min(goodMLWA))

    return ax
    
def clean_cdf(dist, bins=5000):

    hist, b = np.histogram(dist, bins=bins)
    cdf = np.cumsum(1.0*hist/np.sum(hist))
    bcent = 0.5*(b[1:] + b[:-1])

    low_limit = np.interp(0.1,cdf,bcent)
    bidx = np.where(bcent > low_limit)[0]
    print bcent.size, bidx.size

    # cax = plt.figure().add_subplot(111)
    # cax.plot(bcent,cdf)
    # cax.axvline(low_limit,ls=':')
    # cax.figure.show()

    while bidx.size < bins/10.:
        idx = np.where(dist > low_limit)[0]
        dist = dist[idx]
        hist, b = np.histogram(dist, bins=bins)
        cdf = np.cumsum(1.0*hist/np.sum(hist))
        bcent = 0.5*(b[1:] + b[:-1])
        
        low_limit = np.interp(0.05,cdf,bcent)
        bidx = np.where(bcent > low_limit)[0]
        print bidx.size

        # cax = plt.figure().add_subplot(111)
        # cax.plot(bcent,cdf)
        # cax.axvline(low_limit,ls=':')
        # cax.figure.show()

    return dist, bcent, cdf

def A_and_Z(inputfile,ap,numpoints=100):

    import scipy.interpolate as spi

    f = h5py.File(inputfile)
    grp = f['Ap{}'.format(ap)]
    
    lnprob = grp['lnprob']
    MLWA = grp['MLWA_S']
    MLWZ = grp['MLWZ_S']
    
    numsamp = MLWA.size
    
    lnprob = np.reshape(lnprob,numsamp)
    MLWA = np.reshape(MLWA,numsamp)
    MLWZ = np.reshape(MLWZ,numsamp)

    # ax = plt.figure().add_subplot(111)
    # ax.scatter(MLWZ,lnprob)
    # ax.figure.show()

    fidx = np.isfinite(lnprob)
    lnprob = lnprob[fidx]
    MLWA = MLWA[fidx]
    MLWZ = MLWZ[fidx]

    numsamp = MLWA.size

    slnprob = np.sort(lnprob)
    cdf = 1.0 * np.arange(numsamp)/(numsamp -1)

    lowlim = np.interp(0.05,cdf,slnprob)
    idx = np.where(lnprob > lowlim)
    MLWA = MLWA[idx]
    MLWZ = MLWZ[idx]
    lnprob = lnprob[idx]

    Ai = np.linspace(MLWA.min(),MLWA.max(),numpoints)
    Zi = np.linspace(MLWZ.min(),MLWZ.max(),numpoints)

    lni = spi.griddata((MLWA,MLWZ),
                       lnprob,
                       (Ai[None,:],Zi[:,None]),
                       method='nearest',
                       fill_value=np.nan)

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.6,0.6])
    # ax = fig.add_subplot(111)
    ax.set_xlabel('MLWA')
    ax.set_ylabel('MLWZ')
    
    im = ax.imshow(-lni,origin='lower',extent=(Ai.min(),Ai.max(),
                                              Zi.min(),Zi.max()),
                   cmap='gray',aspect='auto')

    levels = [np.interp(i,cdf,slnprob) for i in [0.32,0.68,0.9]]
    c = ax.contour(Ai, Zi, lni, levels, colors='k', linestyles='solid')

    for dist, margdim, ls, rect in zip([Ai,Zi],
                                       [0,1],
                                       [':','--'],
                                       [[0.1,0.71,0.6,0.22],[0.71,0.1,0.22,0.6]]):

        hist = np.nansum(lni,axis=margdim)
        shist = np.sort(hist)
        macdf = 1.0 * np.arange(shist.size)/(shist.size - 1)
        CIlim = np.interp(0.32,macdf,shist)
        CIidx = np.where(hist > CIlim)
        L = np.min(dist[CIidx])
        H = np.max(dist[CIidx])
        M = dist[np.argmax(hist)]
        
        maax = fig.add_axes(rect)
        maax.set_xticklabels([])
        maax.set_yticklabels([])
        
        if margdim == 0:
            maax.plot(dist,hist,drawstyle='steps-mid',color='k')
            maax.axvline(L,ls=ls,color='k')
            maax.axvline(H,ls=ls,color='k')
            maax.axvline(M,ls='-',color='k')
            maax.set_xlim(ax.get_xlim()) 
            maax.set_ylabel('lnprob')
            
            MLWA_L = L
            MLWA_H = H
            MLWA_M = M
        else:
            maax.plot(hist,dist,drawstyle='steps-mid',color='k')
            maax.axhline(L,ls=ls,color='k')
            maax.axhline(H,ls=ls,color='k')
            maax.axhline(M,ls='-',color='k')
            maax.set_ylim(ax.get_ylim()) 
            maax.set_xlabel('lnprob')
            
            MLWZ_L = L
            MLWZ_H = H
            MLWZ_M = M
            
    min_age, max_age = ax.get_xlim()
    min_Z, max_Z = ax.get_ylim()

    ax.scatter([MLWA_M],[MLWZ_M],s=40,c='none',edgecolors='b',linewidths=2)
    ax.set_xlim(min_age,max_age)
    ax.set_ylim(min_Z,max_Z)

    MLWA_string = 'MLWA = {:4.3f} + {:4.3f} - {:4.3f} Gyr'.format(MLWA_M,
                                                          MLWA_M - MLWA_L,
                                                          MLWA_H - MLWA_M)

    MLWZ_string = 'MLWZ = {:4.3f} + {:4.3f} - {:4.3f} Zsol'.format(MLWZ_M,
                                                           MLWZ_M - MLWZ_L,
                                                           MLWZ_H - MLWZ_M)

    fig.text(0.71,0.93,
             '{} Ap# {}\n{}'.format(inputfile,ap,time.asctime()),
             fontsize=8,va='top')
    fig.text(0.71,0.8, MLWA_string, fontsize=8)
    fig.text(0.71,0.75, MLWZ_string, fontsize=8)

    print MLWA_string
    print MLWZ_string

    f.close()

    return ax, c

def AZ_pointing(inputfile, output):
    
    pp = PDF(output)

    i = 1
    while True:
        print i
        try:
            ax, c = A_and_Z(inputfile,i)
            pp.savefig(ax.figure)
            plt.close(ax.figure)
            i += 1
        except KeyError:
            break

    pp.close()

    return

def plot_ages(output,
              plot_LM = False,
              exclude = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]):

    axlist = []
    bigax = plt.figure().add_subplot(111)
    bigax.set_xlabel('|Height [kpc]|')
    bigax.set_ylabel('MLWA')

    plist = [6,3,4,2,1,5]
    color_list = ['blue','seagreen','sienna','orange','yellowgreen','darkturquoise']
    style_list = ['-','-','-','--','--','--']

    for i in range(6):
        pointing = plist[i]
        color = color_list[i]
        style = style_list[i]

        print 'em_P{}_allZ*.emceecoef.fits'.format(pointing)
        coef = glob('em_P{}_allZ*.emceecoef.fits'.format(pointing))[0]
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)

        print coef
        print loc

        c_data = pyfits.open(coef)[1].data
        MLWA = c_data['MLWA']
        MLWA_L = c_data['MLWA_L']
        MLWA_H = c_data['MLWA_H']
        err = np.vstack((MLWA_L,MLWA_H))

        if plot_LM:
            LM_dat = glob('NGC_891_P{}_*CI*dat'.format(pointing))[0]
            print LM_dat
            LM, LM_L, LM_H = np.loadtxt(LM_dat, usecols=(1,2,3), unpack=True)
            LM_err = np.vstack((LM_L,LM_H))
            # LM = np.loadtxt(LM_dat, usecols=(26,), unpack=True)
            # LM_err = np.zeros((2,LM.size))

        r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
        avgr = np.mean(r)

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('|Height [kpc]|')
        ax.set_ylabel('MLWA')
        ax.set_title('{:}\nP{:} ({:3.1f} kpc)'.format(time.asctime(),pointing,avgr))
        
        exarr = np.array(exclude[pointing-1]) - 1
        MLWA = np.delete(MLWA,exarr)
        err = np.delete(err,exarr,axis=1)
        r = np.delete(r,exarr)
        z = np.delete(z,exarr)
        if plot_LM:
            LM = np.delete(LM, exarr)
            LM_err = np.delete(LM_err, exarr, axis=1)

        bigax.errorbar(np.abs(z), MLWA, yerr=err, fmt='.', label='P{}'.format(pointing),
                       color=color,capsize=0)
        ax.errorbar(np.abs(z), MLWA, yerr=err, fmt='.', color=color, capsize=0)
        if plot_LM:
            blme = bigax.errorbar(np.abs(z), LM, yerr=LM_err, fmt='s', mfc='none', mec=color,
                                  lw=0.1,ms=3,color=color, capsize=0)
            lme = ax.errorbar(np.abs(z), LM, yerr=LM_err, fmt='s', mfc='none', mec=color, 
                              lw=0.2,ms=5, color=color, capsize=0)
            blme[-1][0].set_linestyle(':')
            lme[-1][0].set_linestyle(':')

        ax.set_xlim(-0.1,2.6)
        ax.set_ylim(0,11)
        axlist.append(ax)

    bigax.set_title(time.asctime())
    bigax.legend(loc=0, numpoints=1, scatterpoints=1)
    bigax.set_xlim(-0.1,2.6)
    bigax.set_ylim(0,11)

    axlist = [bigax] + axlist

    pp = PDF(output)
    [pp.savefig(a.figure) for a in axlist]
    pp.close()
    plt.close('all')

    return
