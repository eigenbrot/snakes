import pyfits
from matplotlib.backends.backend_pdf import PdfPages as PDF
import numpy as np
import scipy.ndimage as spnd
import scipy.optimze as spo

def do_simple(datafile, errorfile, output, 
              model = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
              plot = True, wavemin = 3750. wavemax = 6800., 
              lightmin = 5450., ligthmax = 5550.):

    m = pyfits.open(model)[1]
    
    dhdu = pyfits.open(datafile)[0]
    header = hdud.header
    data = dhdu.data
    error = pyfits.open(errorfile)[0].data

    numfibers, wavesize = data.shape

    wave = (np.arange(wavesize) - header['CRPIX1']) * header['CDELT1'] \
           +  header['CRVAL1']
    idx = np.where((wave >= wavemin) & (wave <= wavemax))
    wave = wave[idx]

    vdisp = [493., 589., 691., 796., 966.]/2.355
    size_borders = [19, 43, 62, 87, 109]
    size_switch = 0

    lightidx = np.where((wave >= lightmin) & (wave <= lightmax))
    
    agearr = m['AGE']/1e9
    numages = agearr.size

    f = open(output, 'w')
    f.write('{:11}'.format('# Fiber Num'))
    f.write((numages*'{:9} Gyr').format(*agearr))
    f.write((5*'{:13}').\
            format('MMWA [Gyr]','MLWA [Gyr]','Tau_V','S/N/','Chisq'))
    f.write('\n#\n')

    if plot:
        plotname = output.split('.')[0] + '.pdf'
        pp = PDF(plotname)

    L_sun = 3.826e33 #erg / s
    dist_mpx = 10.062
    flux_factor = 1e19
    tau = 2*np.pi

    for i in range(numfibers):

        print 'Doing fiber {}'.format(i+1)
        flux = data[i,idx]*flux_factor
        err = error[i,idx]*flux_factor

        if i == size_boarders[0]:
            size_switch += 1
            size_borders = size_borders[1:]

        #coef = superfit()
    
        SNR = np.sqrt(
            np.sum((flux[lightidx]/err[lightidx])**2)/lightidx.size)

        MMWA = np.sum(agearr*coef['light_frac']*m['M_REMAINING']/m['NORM'])/\
               np.sum(coef['light_frac']*m['M_REMAINING']/m['NORM'])

        redd = np.exp(-1*coef['tauv']*(wave[lightidx]/5500)**(-0.7))
        light_weight = np.mean(m['FLUX'][0][:,lightidx] * redd, axis=1)*\
                       coef['light_frac']

        MLWA = np.sum(light_weight * agearr)/np.sum(light_weight)

        f.write('{:11}'.format(i+1))
        f.write((numages*'{:13.3e}').format(*coef['light_frac']/m['NORM']))
        f.write('{:13.7f}{:13.7f}{:13.3f}{:13.3f}{:13.3f}\n'.\
                format(MMWA,MLWA,coef['tauv'],SNR,coef['chisq']))

    pp.close()
    f.close()

    return 

def superfit(model, restwl, flux, error, vdisp, 
             emmaskw = 400.0, plotlable = ''):

    nmodels = model['AGE'].size
    npix = restwl.size
    quality = np.ones(npix)

    if np.min(resetwl) < np.min(model['WAVE']):
        print "Warning, wavelengths are too small"
    if np.max(resetwl) > np.max(model['WAVE']):
        print "Warning, wavelengths are too big"

    bad = np.where((np.isnan(flux) == True) | (np.isnan(err) == True) |\
                   (err == 0))
    
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

    ok = np.where(quality == 1)
    

    bc03_pix = 70.0
    bc03_vdisp = 75.0

    if vdisp < bc03_vdisp: 
        vdisp_add = 0
    else:
        vdisp_add = np.sqrt(vdisp**2 - bc03_vdisp**2)
    sigma_pix = vdisp_add / bc03_pix

    custom_lib = np.zeros(nmodels, npix)

    for i in range(nmodels):
        cflux = spnd.filters.gaussian_filter1d(m['FLUX'][0][i], sigma_pix)
        custom_lib[i] = np.interp(restwl,model['WAVE'],cflux)

    
    x0 = np.zeros(nmodels+1)
    fitcoefs = spo.leastsq(mcombine, x0, 
                           args=(restwl[ok],flux[ok],err[ok],custom_lib,False))
    
    yfit = mcombine(fitcoefs,restwl,flux,err,custom_lib,True)
    chisq = np.sum((yfit - flux)**2/err**2)

    coefs = {'tauv': fitcoefs[0], 'light_frac': fitcoefs[1:],
             'model_age': m['AGE'], 'chisq': chisq}
    
    fig = plt.figure(figsize=(10,8))
    pax = fig.add_axes([0.15,0.3,0.80,0.69])
    
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
    pax.text(0.2,0.9,'Tau V = {:4.2f}'.format(fitcoefs[0]),
             transform=pax.TransAxes)
    
    lidx = np.where((restwl >= 5450.) & (restwl <= 5550.))
    for i in range(nmodels):
        yi = fitcoefs[i+1] * custom_lib[i,:] *\
             np.exp(-1*fitcoefs[0]*(restwl/5500.)**(-0.7))
        pax.plot(restwl,yi,'b:')
        pax.text(0.2, 0.88 - 0.02*i, 'f_{:02} = {:10.3e}'.\
                 format(i,np.mean(yi[lidx])))
        pax.text(0.8, 0.55 - 0.02*i,'{:6.3f}: {:10.3f}'.\
                 format(model['AGE'][i]/1e9, 
                        coefs['light_frac'][i]/m['NORM'][i]))

    rax = fig.add_axes([0.15,0.15,0.80,0.15])
    rax.set_ylabel('Residuals/error')
    rax.set_xlabel('Wavelength [Angstroms]')
    rax.set_xlim(xmin,xmax)
    rax.set_ylim(-5,5)
    rax.plot(restwl, (galfit - yfit)/err)


    return coefs, fig

def mcombine(X, wave, flux, err, mlib, final):

    if np.any(X[1:] < 0) or np.any(X[1:] > 1):
        return np.ones(flux.size)*1e13

    y = np.sum(mlib * X[1:][:,None],axis=0)
    
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*X[0]*klam)
    y *= e_tau_lam

    if final:
        return y
    else:
        return (flux - y)/err
