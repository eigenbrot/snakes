import numpy as np
import fsps

def make_single_age(age, Z=1.0, wavemin=3000, wavemax=9000):

    z_arr = np.array([0.0002, 0.0003, 0.0004, 
                      0.0005, 0.0006, 0.0008, 
                      0.0010, 0.0012, 0.0016, 
                      0.0020, 0.0025, 0.0031, 
                      0.0039, 0.0049, 0.0061, 
                      0.0077, 0.0096, 0.0120, 
                      0.0150, 0.0190, 0.0240, 0.0300])/0.0190

    
    zidx = np.argmin(np.abs(Z - z_arr))
    print '\t',zidx
    if zidx.size == 0:
        print 'Could not find Z = {}'.format(Z)
        return
        
    sp = fsps.StellarPopulation(add_stellar_remnants=False, 
                                add_dust_emission=False, 
                                add_agb_dust_model=False, 
                                dust_type=0, 
                                dust_index=0, 
                                imf_type=1, 
                                sf_start=age, 
                                zmet=zidx, 
                                sfh=0)

    wave, spec = sp.get_spectrum(tage=13.8, peraa=True)

    idx = np.where((wave > wavemin) & (wave < wavemax))[0]

    return wave[idx], spec[idx]

def make_many(ages, output, Z=1.0):

    with open(output, 'w') as f:
        for a in ages:
            print a
            wave, spec = make_single_age(a, Z=Z)
            for i in range(wave.size):
                f.write('{:13.6f}{:10.3f}{:10.1f}{:13.4e}\n'.format(a,
                                                                Z*0.0190,
                                                                wave[i],
                                                                spec[i]))

    return

def match_bc03_ages(output, Z=1.0):

    dfk_file = '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/eps_0.07k_kmeangrp_avgnoweightsplit_newdt_26Apr14_age_sspweight.dat'

    ages = np.loadtxt(dfk_file, usecols=(0,), unpack=True)
    
    make_many(ages/1e9, output, Z=Z)

    return
