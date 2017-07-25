import os
import time
import numpy as np
import pyfits
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.ndimage as spnd
from matplotlib.backends.backend_pdf import PdfPages as PDF

def make_galaxy(output,
                SSPs = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
                sigma = '/d/monk/eigenbrot/WIYN/14B-0456/anal/pre_ind/SN_realization/Noise.fits',
                vdisp = 209.34,
                moddisp = 75.0,
                tau_sf = 8.0,
                t_form = 12.0,
                tau_V = 1.5,
                Mtot = 12e9,
                SN = np.inf,
                SNmin = 5450.,
                SNmax = 5550.,
                lightmin = 5450.,
                lightmax = 5550.,
                makeplot=True,
                writeoutput=True):

    '''Makes a galaxy by adding up SSPs with different ages based on a tau star
    formation model:

    \psi(t) = \psi_0 exp(-t/tau_sf),

    where \psi(t) is the star formation rate.
    '''

    model = pyfits.open(SSPs)[1].data
    flux = model['FLUX'][0]
    wave = model['WAVE'][0]
    ssp_age = model['AGE'][0]/1e9 #in Gyr
    norm = model['NORM'][0]
    flux *= norm[:,None]
    
    logt = np.log10(np.r_[1e-99,ssp_age,t_form])
    tdiff = np.diff(logt)
    borders = 10**(logt[1:] - tdiff/2.)
    borders[0] = 1e-99
    borders[-1] = t_form

    while np.diff(borders)[-1] < 0:
        useidx = np.where(np.diff(borders) >= 0)
        useidx = np.r_[useidx[0],useidx[0][-1]+1]
        borders = borders[useidx]
        borders[-1] = t_form
        flux = flux[(useidx[:-1],)]
    
    mass = tau_sf*(np.exp(borders[1:]/tau_sf) 
                   - np.exp(borders[:-1]/tau_sf))

    weighted_age = tau_sf*np.log(0.5*(np.exp(borders[1:]/tau_sf)
                                      + np.exp(borders[:-1]/tau_sf)))
    
    plot_age = np.linspace(0,t_form,1000)
    plot_psi = np.exp(plot_age/tau_sf)
    psi = np.exp(weighted_age/tau_sf)

    psi0 = Mtot/np.sum(mass)
    mass *= psi0
    psi *= psi0
    plot_psi *= psi0
    MMWA = np.sum(weighted_age*mass)/np.sum(mass)

    galaxy = np.sum(flux*mass[:,None],axis=0)
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tau_V*klam)
    galaxy *= e_tau_lam
    idx = np.where((wave >= lightmin) & (wave <= lightmax))[0]

    light_weight = np.mean(flux[:,idx] * e_tau_lam[idx],axis=1)*mass
    MLWA = np.sum(light_weight * weighted_age)/np.sum(light_weight)

#    linwave = np.linspace(wave.min(),wave.max(),wave.size)
    linwave = np.arange(wave.min(),wave.max(),2.1)
    # lingal = np.interp(linwave, wave, galaxy)

    npix = np.round((np.log10(wave.max()) - np.log10(wave.min())) / 1e-4)
    logwave = 10**(np.arange(npix) * 1e-4 + np.log10(wave.min()))
    loggal = np.interp(logwave, wave, galaxy)

    mpix = np.mean(np.diff(logwave)/logwave[1:]*3e5)

    if vdisp > moddisp:
        vdisp_add = np.sqrt(vdisp**2 - moddisp**2)
        sigma_pix = vdisp_add / mpix
        loggal = spnd.filters.gaussian_filter1d(loggal,sigma_pix)

    lingal = np.interp(linwave, logwave, loggal)
        
    if np.isfinite(SN):
        linwave, lingal, error = add_noise(linwave, lingal, sigma, SN,
                                           SNmin = SNmin, 
                                           SNmax = SNmax)
    else:
        error = np.ones(linwave.size)

    lingal *= 1e-17
    error *= 1e-17

    if makeplot:
        fig = plt.figure(figsize=(10,8))
        fig.suptitle('MLWA = {:5.3} Gyr'.format(MLWA))
        axg = fig.add_subplot(212)
        axs = fig.add_subplot(211)

        axg.plot(linwave,lingal,'k')
        if np.isfinite(SN):
            axg.fill_between(linwave, lingal-error, lingal + error, 
                             color='k', alpha=0.6)
        axg.set_xlabel('Wavelength [Angstroms]')
        axg.set_ylabel('Flux [Arbitrary]')
        axg.text(0.1,0.9,'$S/N =  {}$'.format(SN),transform=axg.transAxes)

        axs.plot(weighted_age,np.log10(psi),'.k')
        axs.plot(weighted_age,np.log10(mass),'.g', label='Log(Mass [$M_{\odot}$])')
        for i in range(weighted_age.size):
            axs.fill_between(plot_age,np.log10(plot_psi),alpha=0.3,
                             where=(plot_age > borders[:-1][i]) & 
                             (plot_age < borders[1:][i]))

        axs.set_xlabel('Lookback time [Gyr]')
        axs.set_ylabel('Log( $\psi(t)$ [$M_{\odot}$/Gyr] )')
        axs.set_xlim(t_form + 1,-1)
        axs.set_ylim(5,11)
        axs.text(0.1,0.7,r'$\tau_{{SF}} = {:3.1f}$'.format(tau_sf),
                 transform=axs.transAxes)
        axs.text(0.1,0.65,r'$\tau_V = {:3.1f}$'.format(tau_V),
                 transform=axs.transAxes)
        axs.text(0.1,0.6,r'$M_{{tot}} = {:4.1e} M_{{\odot}}$'.format(Mtot),
                 transform=axs.transAxes)
        axs.text(0.1,0.55,r'$\Rightarrow\psi_0 = {:4.1e} M_{{\odot}}/Gyr$'.\
                 format(psi0),transform=axs.transAxes)
        axs.legend(loc=0,frameon=True,numpoints=1,fontsize=8)

        pp = PDF('{}_galaxy.pdf'.format(output))
        pp.savefig(fig)
        pp.close()

    if writeoutput:
        outname = '{}.ms_lin.fits'.format(output)
        hdu = pyfits.PrimaryHDU(lingal)
        hdu.header.update('CTYPE1','LINEAR')
        hdu.header.update('CRPIX1',1)
        hdu.header.update('CRVAL1',linwave[0])
        hdu.header.update('CDELT1',np.mean(np.diff(linwave)))
        hdu.writeto(outname,clobber=True)

        errname = '{}.me_lin.fits'.format(output)
        ehdu = pyfits.PrimaryHDU(error)
        ehdu.header.update('CTYPE1','LINEAR')
        ehdu.header.update('CRPIX1',1)
        ehdu.header.update('CRVAL1',linwave[0])
        ehdu.header.update('CDELT1',np.mean(np.diff(linwave)))
        ehdu.writeto(errname,clobber=True)

        f = open('{}_model.dat'.format(output),'w')
        f.write(str('# Generated on {}\n' + 
                    '# tau_sf = {:}\n' +
                    '# tau_V = {:}\n' +
                    '# vdisp = {:}\n' +
                    '# Mtot = {:3.1e}\n').format(time.asctime(),tau_sf,tau_V,vdisp,Mtot))
        f.write(str('# {:>11}'+weighted_age.size*'{:>9.3f} Gyr'+'{:>13}{:>13}\n#\n').\
                format('Num',*weighted_age.tolist()+['MMWA [Gyr]','MLWA [Gyr]']))
        f.write(str('{:13}'+'{:13.3e}'*(weighted_age.size+2)).format(1,*mass.tolist()+[MMWA,MLWA]))
        f.write('\n')
        f.close()

    return {'wave':linwave, 'flux':lingal, 'err':error,
            'MMWA': MMWA, 'MLWA': MLWA, 'age': weighted_age, 'mass': mass}
    
def add_noise(wave, spectrum, sigmafile, 
              desSN, SNmin = 5450., SNmax = 5550.):

    hdu = pyfits.open(sigmafile)[0]
    sigma = hdu.data
    swave = (np.arange(sigma.size) + hdu.header['CRPIX1'])*hdu.header['CDELT1'] + hdu.header['CRVAL1']
    bidx = np.where((wave >= swave.min()) & (wave <= swave.max()))
    wave = wave[bidx]
    spectrum = spectrum[bidx]
    sigma = np.interp(wave, swave, sigma)
    
    idx = np.where((wave >= SNmin) & (wave <= SNmax))[0]

    s = np.mean(spectrum[idx]/sigma[idx])/desSN
    noise = s * sigma

    r = np.random.randn(wave.size)
    fp = spectrum + r * s * sigma
    
    print np.mean(spectrum[idx]/(noise[idx]))
    print np.mean(fp[idx]/(noise[idx]))

    return wave, fp, noise

def tau_age(output, taulist = None, t_formlist = None,
            mintau = 0.1, maxtau = 12, numtau = 10):

    
    if taulist is None:
        taulist = np.linspace(mintau, maxtau, numtau)

    ax = plt.figure().add_subplot(111)
    ax1 = plt.figure().add_subplot(111)
    for t_form in t_formlist:
        MMWA_list = np.zeros(taulist.size)
        MLWA_list = np.zeros(taulist.size)
        for i, tau in enumerate(taulist):
            name = 'tau_{:3.1f}_{:3.1f}'.format(tau,t_form)
            print name
            d = make_galaxy(name,tau_sf=tau,t_form=t_form)
            MMWA_list[i] = d['MMWA']
            MLWA_list[i] = d['MLWA']

        #    ax.plot(taulist,MMWA_list,'.',label='MMWA')
        ax.plot(taulist,MLWA_list,label=t_form)
        ax1.plot(taulist,MMWA_list - MLWA_list,'.')

    ax.set_xlabel(r'$\tau_{sf}$ [Gyr]')
    ax.set_ylabel('MLWA [Gyr]')
    ax.legend(loc=0, frameon=False, numpoints=1, scatterpoints=1, title='t$_{\mathrm{form}}$')

    ax1.set_xlabel(r'$\tau_{sf}$ [Gyr]')
    ax1.set_ylabel('MMWA - MLWA')
    
    pp = PDF('{}_age.pdf'.format(output))
    pp.savefig(ax.figure)
#    pp.savefig(ax1.figure)
    pp.close()
    
    return taulist, MMWA_list

# def make_monte(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60],
#                N = 100):

#     ssp = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_008Z_ChabIMF.fits'
#     for SN in SNlist:
#         direc = 'SN{}'.format(SN)
#         if not os.path.exists(direc):
#             os.makedirs(direc)
#         f = open('{}/run.pro'.format(direc),'w')
#         for tau in taulist:
#             for i in range(N):
#                 name = '{}/SN{:02}_t{:03}_N{:03}'.format(direc,SN,tau,i+1)
#                 print name
#                 make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau)
#                 f.write("do_simple, '{0:}.ms_lin.fits', '{0:}.me_lin.fits', '{0:}_fit.dat', wavemin=3750., wavemax=6800., model='{1:}', /plot\n".format(os.path.basename(name),ssp))
#         f.close()
#     return

# def compare_SN(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60], N = 10):

#     ratios = np.zeros((len(SNlist),len(taulist)))
#     errs = np.zeros(ratios.shape)

#     for s, SN in enumerate(SNlist):
#         direc = 'SN{}'.format(SN)
#         for t, tau in enumerate(taulist):
#             tmp = np.zeros(N)
#             for i in range(N):
#                 model_name = '{}/SN{:02}_t{:03}_N{:03}_model.dat'.\
#                              format(direc,SN,tau,i+1)
#                 fit_name = '{}/SN{:02}_t{:03}_N{:03}_fit.dat'.\
#                            format(direc,SN,tau,i+1)
#                 print model_name
#                 mMMWA, mMLWA = np.loadtxt(model_name,
#                                           usecols=(11,12),unpack=True)
#                 fMMWA, fMLWA = np.loadtxt(fit_name,usecols=(11,12),unpack=True)
#                 tmp[i] = 1 - fMMWA/mMMWA
#                 ratios[s,t] = np.mean(tmp)
#                 errs[s,t] = np.std(tmp)

#     ax = plt.figure().add_subplot(111)
#     for i in range(len(taulist)):
#         ax.errorbar(SNlist, ratios[:,i], yerr=errs[:,i], label=taulist[i])
# #        ax.plot(SNlist,errs[:,i], label=taulist[i])

#     ax.set_xlim(0,70)
#     ax.set_xlabel('SNR')
# #    ax.set_ylabel(r'$\sigma$(1 - MLWA$_{fit}$/MLWA$_{true}$)')
#     ax.set_ylabel('1 - MMWA$_{fit}$/MMWA$_{true}$')
#     ax.legend(loc=0,numpoints=1,scatterpoints=1,frameon=False,title=r'$\tau_{sf}$')

#     return ratios, errs, SNlist, taulist, ax
