import numpy as np
import pyfits
import matplotlib.pyplot as plt

def make_galaxy(output,
                SSPs = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_solarZ_ChabIMF.fits',
                tau_sf = 8.0,
                tau_V = 1.5,
                SN = np.inf):

    '''Makes a galaxy by adding up SSPs with different ages based on a tau star
    formation model:

    \psi(t) = \psi_0 exp(-t/tau_sf),

    where \psi(t) is the star formation rate.
    '''

    model = pyfits.open(SSPs)[1].data
    flux = model['FLUX'][0]
    wave = model['WAVE'][0]
    ages = model['AGE'][0]/1e9 #in Gyr

    psi = 10*np.exp(-ages/tau_sf)
    
    galaxy = np.sum(flux*psi[:,None],axis=0)
    klam = (wave / 5500.)**(-0.7)
    e_tau_lam = np.exp(-1*tau_V*klam)
    galaxy *= e_tau_lam

    fig = plt.figure()
    axg = fig.add_subplot(211)
    axs = fig.add_subplot(212)

    axg.plot(wave,galaxy,'k')
    axg.set_xlabel('Wavelength [Angstroms]')
    axg.set_ylabel('Flux [Arbitrary]')
    axs.plot(ages,np.log(psi),'k')
    axs.bar(ages,np.log(psi),align='center',alpha=0.5)
    axs.set_xlabel('$t$ [Gyr]')
    axs.set_ylabel('$\psi(t)$')
    axs.text(0.8,0.8,r'$\tau_{{SF}} = {:3.1f}$'.format(tau_sf),
             transform=axs.transAxes)
    axs.text(0.8,0.7,r'$\tau_V = {:3.1f}$'.format(tau_V),
             transform=axs.transAxes)

    fig.show()

    return galaxy
    
