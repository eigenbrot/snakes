import sys
import time
import numpy as np
import pyfits
import scipy.ndimage as spnd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def plot_bc(coeffile, fitfile, datafile, errorfile, model,
            output=None, location=None, wavemin=3800., wavemax=6800.,
            plotblue = False):

    flux_factor = 1e17
    smoothkern = 2
    
    dhdu = pyfits.open(datafile)[0]
    data = dhdu.data
    head = dhdu.header
    error = pyfits.open(errorfile)[0].data
    
    numfibers, wavesize = data.shape

    cdelt = head['CDELT1']
    crval = head['CRVAL1']
    crpix = head['CRPIX1']

    print 'CDELT1 = ', cdelt
    print 'CRVAL1 = ', crval
    print 'CRPIX1 = ', crpix

    wave = (np.arange(wavesize) - crpix) * cdelt + crval

    idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
    restwl = wave[idx]

    if location is not None:
        fiber_radii, rkpc, zkpc = np.loadtxt(location, 
                                             usecols=(1,4,5),
                                             unpack=True)
        sizeidx = [0.937,1.406,1.875,2.344,2.812]

    m = pyfits.open(model)[1].data[0]
    nmodels = m['FLUX'].shape[0]

    if output is None:
        if plotblue:
            output = fitfile.split('.')[0]+'.fit.blue.pdf'
        else:
            output = fitfile.split('.')[0]+'.fit.pdf'

    pp = PDF(output)

    coef_arr = pyfits.open(coeffile)[1].data
    yfits = pyfits.open(fitfile)[0].data

    size_borders = [19, 43, 62, 87, 109] # The last one is needed to prevent indexing errors
    size_switch = 0
    
    redidx = np.where(restwl >= 5400)
    blueidx = np.where(restwl < 5400)
    hklow = 3920
    hkhigh = 4000
    hkidx = np.where((restwl > hklow) & (restwl < hkhigh))
    npix = restwl.size

    for i in range(numfibers):
        
        flux = data[i,idx]*flux_factor
        err = error[i,idx]*flux_factor

        if location is not None:
            vdidx = np.where(sizeidx == fiber_radii[i])[0][0]
            plotlabel = 'Aperture {:n}, r={:6.2f}, z={:5.2f}'.\
                        format(i+1,rkpc[i],zkpc[i])
        else:
            if i == size_borders[0]:
                size_switch += 1
                size_borders = size_borders[1:]
            vdidx = size_switch
            plotlabel = 'Fiber {:n}'.format(i+1)
    
        print plotlabel

        quality = np.ones(npix)
        bad = np.where((err == 0) | (np.isnan(flux)) | (np.isnan(err)))
        quality[bad] = 0
    
        outside_model = np.where((restwl < m['WAVE'].min()) |\
                                 (restwl > m['WAVE'].max()))
        quality[outside_model] = 0

        sk =    np.array([6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0])
        sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
        
        sk2 = np.array([6300., 5890., 5577.])
        em2 = [6563.8,  6716.0, 6583.41, 6548.04]
        
        em = np.array([6563.8,  6716.0, 6583.41, 6548.04, 4959., 5006.8])
        emnam = [r'H$\alpha$', 'S2', 'NII', 'NII', '[OIII]', '[OIII]']
        
        ab = np.array([3820.4, 3835.4,      3889.0,     3933.7, 3968.5, 3970.18,         4304.4,   4341.,       5175.3, 5894.0, 4861.,  4102., 3820.4])
        absnam = ['L',   r'H$\eta$', r'H$\zeta$', 'K',   'H'   , r'H$\epsilon$',    'G',     r'H$\gamma$',  'Mg',   'Na',   r'H$\beta$',   r'H$\delta$',  'L']
        
        dz = 1500. / 3e5
        dzsk = 1600. / 3e5
        
        for ee in em2:
            maskout = np.where((restwl > ee*(1-dz)) & (restwl < ee*(1+dz)))
            quality[maskout] = 0
    
        for ss in sk2:
            maskout = np.where((restwl > ss*(1-dzsk)) & (restwl < ss*(1+dzsk)))
            quality[maskout] = 0

        ok = quality == 1

        # Interpolate models to data wavelength grid
    
        custom_lib = np.zeros((nmodels, npix))
        for ii in range(nmodels):
            custom_lib[ii,:] = np.interp(restwl, 
                                         m['WAVE'], m['FLUX'][ii,:,vdidx])
        custom_lib[:,outside_model] = 0

        yfit = yfits[i,:] * flux_factor
        coefs = coef_arr[i]

    
        ############################
        ############################
        
        if plotblue:
            pidx = np.where(restwl < 4500.)
        else:
            pidx = np.where(restwl == restwl)
            
        pwave = restwl[pidx]

        xmin = pwave.min() * 0.98
        xmax = pwave.max() * 1.02
        
        fig = plt.figure(figsize=(11,8))
        fax = fig.add_axes([0.1,0.25,0.85,0.69])
        fax.set_xticklabels([])
        fax.set_ylabel('Log Flux + 17 [erg/s/cm$^2$/$\AA$]')
        fax.set_xlim(xmin,xmax)
        fax.set_ylim(-0.49, 2.6)
    
        fax.axvline(x=5400, color='k', ls=':', alpha=0.3)
        fax.axvline(x=hklow, color='k', ls=':', alpha=0.3)
        fax.axvline(x=hkhigh, color='k', ls=':', alpha=0.3)
        
        for j in range(coefs['LIGHT_FRAC'].size):
            xred = restwl * (coefs['VSYS']/3e5 + 1)
            yi = coefs['LIGHT_FRAC'][j] * custom_lib[j,:] *\
                 np.exp(-1 * coefs['TAUV']*(restwl/5500.)**(-0.7))
            yi = np.interp(restwl, xred, yi)
            fax.plot(pwave,
                     np.log10(spnd.filters.gaussian_filter1d(yi,smoothkern))[pidx],
                     color='b',alpha=0.4)

        galfit = np.zeros(flux.size) + yfit
        galfit[ok] = flux[ok]
        plotgal = spnd.filters.gaussian_filter1d(galfit,smoothkern)
        galfit[~ok] = np.NAN
        plotgal[~ok] = np.NAN
        masked = spnd.filters.gaussian_filter1d(flux,smoothkern)
        masked[ok] = np.NAN
        plotfit = spnd.filters.gaussian_filter1d(yfit,smoothkern)

        fax.plot(pwave,np.log10(plotgal)[pidx],color='k')
        fax.plot(pwave,np.log10(masked)[pidx],color='c',lw=4)
        fax.plot(pwave,np.log10(plotfit)[pidx],color='r')

        fax.fill_between(pwave,
                         np.log10(plotgal - 
                                  spnd.filters.gaussian_filter1d(
                                      err,smoothkern*2))[pidx],
                         np.log10(plotgal +
                                  spnd.filters.gaussian_filter1d(
                                      err,smoothkern*2))[pidx],
                         color='k', alpha=0.2)

        ###################################
        ###################################
        Z = coefs['VSYS']/3e5 + 1

        ypos = 1
        for s, sn in zip(sk, sknam):
            if s > 5500. and plotblue: continue
            tidx = np.where((restwl >= s - 10) & (restwl <= s + 10.))
            try:
                ypos = np.log10(np.nanmax(np.r_[flux[tidx],yfit[tidx]])) + 0.1
            except ValueError:
                pass
            fax.text(s, ypos, sn, fontsize=8, 
                     ha='center', va='center')
            if plotblue:
                fax.axvline(s, color='k', ls=':', alpha=0.7)
            else:
                fax.plot((s,s), (ypos-0.04,ypos-0.1), color='k', alpha=0.8)

            
        prevy = 99
        for a, an in zip(ab*Z, absnam):
            if a > 5500. and plotblue: continue
            tidx = np.where((restwl >= a - 10) & (restwl <= a + 10.))
            try:
                ypos = np.log10(np.nanmin(np.r_[flux[tidx],yfit[tidx]])) - 0.1
            except ValueError:
                pass
            if (an == r'H$\gamma$' or 
                an == r'H$\eta$' or
                an == r'H$\epsilon$') and np.abs(ypos - prevy) <= 0.04:
                ypos += 0.06
            prevy = ypos
            if np.isnan(ypos) or ypos < -0.4:
                ypos = -0.4
            fax.text(a, ypos, an, color='r', fontsize=8, 
                     ha='center', va='center')
            if plotblue:
                fax.axvline(a, color='r', ls=':', alpha=0.7)
            else:
                fax.plot((a,a), (ypos+0.04,ypos+0.1), color='r', alpha=0.8)

        for e, en in zip(em*Z, emnam):
            if e > 5500. and plotblue: continue
            tidx = np.where((restwl >= e - 10) & (restwl <= e + 10.))
            try:
                ypos = np.log10(np.nanmax(np.r_[flux[tidx],yfit[tidx]])) + 0.1
            except ValueError:
                pass
            fax.text(e, ypos, en, color='b', fontsize=8, 
                     ha='center', va='center')
            if plotblue:
                fax.axvline(e, color='b', ls=':', alpha=0.7)
            else:
                fax.plot((e,e), (ypos-0.04,ypos-0.1), color='b', alpha=0.8)


        ############################################
        ############################################

        chivec = (galfit - yfit)/err
        chivec[~ok] = np.nanmean(chivec)
        plotchi = spnd.filters.gaussian_filter1d(chivec,smoothkern)
        plotchi[~ok] = np.NAN
        eax = fig.add_axes([0.1,0.1,0.85,0.15])
        eax.set_xlabel('Wavelength [$\AA$]')
        eax.set_ylabel('Residuals/error')
        eax.set_xlim(fax.get_xlim())
        eax.set_ylim(-6,6)
        eax.set_yticks([-5,0,5])
        eax.plot(pwave, plotchi[pidx], color='k')

        ###########################################
        ###########################################
    
        fs = 10

        fig.text(0.15, 0.92, plotlabel, fontsize=fs)

        fig.text(0.15, 0.89, 'SNR = {:8.2f}'.format(coefs['SNR']), fontsize=fs)
        fig.text(0.15, 0.87, 'V = {:8.2f} km/s'.format(coefs['VSYS']), 
                 fontsize=fs)
        fig.text(0.15, 0.85, 'V_disp = {}'' fiber'.format(vdidx+2),
                 fontsize=fs)
        fig.text(0.15, 0.83, r'$\tau_V$ = {:8.2f}'.format(coefs['tauv']), 
                 fontsize=fs)
        
        fig.text(0.35, 0.89, 'MLWA = {:8.2f}'.format(coefs['MLWA']),
                 fontsize=fs)
        fig.text(0.35, 0.87, 'MMWA = {:8.2f}'.format(coefs['MMWA']),
                 fontsize=fs)
        fig.text(0.35, 0.85, 'MLWZ = {:8.2f}'.format(coefs['MLWZ']),
                 fontsize=fs)
        fig.text(0.35, 0.83, 'MMWZ = {:8.2f}'.format(coefs['MMWZ']),
                 fontsize=fs)
        
        fig.text(0.55, 0.89, r'$\chi^2$ = {:8.2f}'.format(coefs['chisq']), 
                 fontsize=fs)
        fig.text(0.55, 0.87, r'$\chi^2_\mathrm{{red}}$ = {:8.2f}'.\
                 format(coefs['redchi']), fontsize=fs)
        fig.text(0.55, 0.85, r'$\chi^2_\mathrm{{blue}}$ = {:8.2f}'.\
                 format(coefs['bluechi']), fontsize=fs)
        fig.text(0.55, 0.83, r'$\chi^2_\mathrm{{HK}}$ = {:8.2f}'.\
                 format(coefs['hkchi']), fontsize=fs)
        
        fig.suptitle(time.asctime())

        pp.savefig(fig)
        plt.close(fig)

    pp.close()
    return 

def parse_input(inputlist):
    
    coeffile = inputlist[0]
    fitfile = inputlist[1]
    datafile = inputlist[2]
    errorfile = inputlist[3]
    modelfile = inputlist[4]

    kwar = {}

    i = 5
    while i < len(inputlist):
        
        if inputlist[i] == '-o':
            kwar['output'] = inputlist[i+1]
            i += 1
            
        if inputlist[i] == '-l':
            kwar['location'] = inputlist[i+1]
            i += 1
            
        if inputlist[i] == '-w':
            kwar['wavemin'] = inputlist[i+1]
            kwar['wavemax'] = inputlist[i+2]
            i += 2

        if inputlist[i] == '-b':
            kwar['plotblue'] = True

        if inputlist[i] == '-n':
            import nice_plots
            nice_plots.format_plots(False)

        i += 1

    return [coeffile, fitfile, datafile, errorfile, modelfile], kwar

if __name__ == '__main__':
    
    args, kwargs = parse_input(sys.argv[1:])
    print args
    print kwargs
    plot_bc(*args, **kwargs)
