import numpy as np
import pyfits
import matplotlib.pyplot as plt
import scipy.ndimage as spnd
from matplotlib.backends.backend_pdf import PdfPages as PDF

def make_galaxy(output,
                SSPs = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
                vdisp = 209.34,
                tau_sf = 8.0,
                tau_V = 1.5,
                Mtot = 12e9,
                SN = np.inf):

    '''Makes a galaxy by adding up SSPs with different ages based on a tau star
    formation model:

    \psi(t) = \psi_0 exp(-t/tau_sf),

    where \psi(t) is the star formation rate.
    '''
    form_age = 12.0 #in Gyr

    model = pyfits.open(SSPs)[1].data
    flux = model['FLUX'][0]
    wave = model['WAVE'][0]
    ssp_age = model['AGE'][0]/1e9 #in Gyr
    norm = model['NORM'][0]
    flux *= norm[:,None]

    logt = np.log10(np.r_[1e-99,ssp_age,form_age])
    tdiff = np.diff(logt)
    borders = 10**(logt[1:] - tdiff/2.)
    borders[0] = 1e-99
    borders[-1] = form_age
    
    plot_age = np.linspace(0,form_age,1000)
    plot_psi = np.exp(plot_age/tau_sf)
    psi = np.exp(ssp_age/tau_sf)

    mass = tau_sf*(np.exp(borders[1:]/tau_sf) 
                   - np.exp(borders[:-1]/tau_sf))
    psi0 = Mtot/np.sum(mass)
    mass *= psi0
    psi *= psi0
    plot_psi *= psi0

    galaxy = np.sum(flux*mass[:,None],axis=0)
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tau_V*klam)
    galaxy *= e_tau_lam

    linwave = np.linspace(wave.min(),wave.max(),wave.size)
    lingal = np.interp(linwave, wave, galaxy)
    
    bc03_pix = 70.0
    bc03_vdisp = 75.0

    vdisp_add = np.sqrt(vdisp**2 - bc03_vdisp**2)
    sigma_pix = vdisp_add / bc03_pix
    lingal = spnd.filters.gaussian_filter1d(lingal,sigma_pix)
    
    if np.isfinite(SN):
        error, noise = add_noise(linwave, lingal, SN)
        lingal += noise
    else:
        error = np.ones(linwave.size)


    lingal *= 1e-17
    error *= 1e-17

    fig = plt.figure(figsize=(10,8))
    axg = fig.add_subplot(212)
    axs = fig.add_subplot(211)

    axg.plot(linwave,lingal,'k')
    if np.isfinite(SN):
        axg.fill_between(linwave, lingal-error, lingal + error, 
                         color='k', alpha=0.6)
    axg.set_xlabel('Wavelength [Angstroms]')
    axg.set_ylabel('Flux [Arbitrary]')
    axg.text(0.1,0.9,'$S/N =  {}$'.format(SN),transform=axg.transAxes)

    axs.plot(ssp_age,np.log10(psi),'.k')
    axs.plot(ssp_age,np.log10(mass),'.g', label='Log(Mass [$M_{\odot}$])')
    for i in range(ssp_age.size):
        axs.fill_between(plot_age,np.log10(plot_psi),alpha=0.3,
                         where=(plot_age > borders[:-1][i]) & 
                         (plot_age < borders[1:][i]))

    axs.set_xlabel('Lookback time [Gyr]')
    axs.set_ylabel('Log( $\psi(t)$ )')
    axs.set_xlim(13,-1)
    axs.text(0.1,0.7,r'$\tau_{{SF}} = {:3.1f}$'.format(tau_sf),
             transform=axs.transAxes)
    axs.text(0.1,0.65,r'$\tau_V = {:3.1f}$'.format(tau_V),
             transform=axs.transAxes)
    axs.text(0.1,0.6,r'$M_{{tot}} = {:4.1e} M_{{\odot}}$'.format(Mtot),
             transform=axs.transAxes)
    axs.text(0.1,0.55,r'$\Rightarrow\psi_0 = {:4.1e} M_{{\odot}}/Gyr$'.\
             format(psi0),transform=axs.transAxes)
    axs.legend(loc=0,frameon=False,numpoints=1,fontsize=8)


    pp = PDF('{}_galaxy.pdf'.format(output))
    pp.savefig(fig)
    pp.close()

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

    return linwave, lingal, error
    
def add_noise(wave, spectrum, desSN, lightmin = 5450., lightmax = 5550.):

    idx = np.where((wave >= lightmin) & (wave <= lightmax))[0]
    error = np.random.rand(wave.size)*0.01
    error = np.random.rand(wave.size)
    mult = 0
    SN = np.inf
    while SN >= desSN:
        mult += 0.001
        tmp = error + mult
        SN = np.sqrt(np.sum((spectrum[idx]/tmp[idx])**2)/idx.size)

    specnoise = np.mean(error + mult)*np.random.randn(wave.size)
        
    print SN

    return error + mult, specnoise
