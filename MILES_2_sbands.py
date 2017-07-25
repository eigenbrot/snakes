import os
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
import scipy.ndimage as spnd

full_age = [00.03, 00.04, 00.05, 00.06, 00.07, 00.08, 00.09, 00.10,
00.15, 00.20, 00.25, 00.30, 00.35, 00.40, 00.45, 00.50, 00.60, 00.70,
00.80, 00.90, 01.00, 01.25, 01.50, 01.75, 02.00, 02.25, 02.50, 02.75,
03.00, 03.25, 03.50, 03.75, 04.00, 04.50, 05.00, 05.50, 06.00, 06.50,
07.00, 07.50, 08.00, 08.50, 09.00, 09.50, 10.00, 10.50, 11.00, 11.50,
12.00, 12.50, 13.00, 13.50, 14.00]

def read_fits(fitsfile):

    hdu = pyfits.open(fitsfile)[0]
    header = hdu.header

    flux = hdu.data
    wave = (np.arange(flux.shape[0]) - (header['CRPIX1'] - 1))*header['CDELT1'] + header['CRVAL1']

    return wave, flux

def convert_MILES(mwave, mflux, vdisp):

    linwave = np.arange(mwave.min(), mwave.max(), 2.1)

    npix = np.round((np.log10(mwave.max()) - np.log10(mwave.min())) / 1e-4)
    logwave = 10**(np.arange(npix) * 1e-4 + np.log10(mwave.min()))
    logflux = np.interp(logwave, mwave, mflux)

    mpix = np.mean(np.diff(logwave)/logwave[1:]*3e5)
    mdisp = 58.4

    if vdisp > mdisp:
        disp_add = np.sqrt(vdisp**2 - mdisp**2)
        disp_pix = disp_add / mpix
        logflux = spnd.filters.gaussian_filter1d(logflux, disp_pix)

    newflux = np.interp(linwave, logwave, logflux)

    return linwave, newflux

def combine_Z(MH, vdisp, agelst = [0.35, 0.6, 1, 3, 4, 8, 12], alpha=0.0):

    basedir = '/Users/Arthur/Documents/School/891_research/MILES'

    if MH < 0:
        pm = 'm'
        strmetal = MH * -1
    else:
        pm = 'p'
        strmetal = MH
    
    final_list = []
    outdirname = ''.join('{:6.2f}'.format(vdisp).split('.'))
    outnamebase = '{}/MILES_MH{}{}_Ep{}'.format(outdirname,pm,strmetal,alpha)

    if not os.path.exists(outdirname):
        os.system('mkdir {}'.format(outdirname))

    pp = PDF(outnamebase+'.pdf')
    
    for age in agelst:
        model_file = '{0:}/MILES_BASTI_CH_Ep{1:4.2f}/Mch1.30Z{2:}{3:4.2f}T{4:07.4f}_iTp{1:4.2f}_Ep{1:4.2f}.fits'.\
           format(basedir, alpha, pm, strmetal, age)
        print model_file

        wave, flux = read_fits(model_file)
        conv_wave, conv_flux = convert_MILES(wave, flux, vdisp)

        final_list.append(conv_flux)

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        ax.plot(conv_wave, conv_flux)
        pp.savefig(ax.figure)
        plt.close(ax.figure)
        
    output = np.vstack(final_list)
    print output.shape

    pp.close()
    
    HDU = pyfits.PrimaryHDU(output)
    HDU.header.update(CRPIX1=1, CRVAL1=conv_wave.min(), CDELT1=np.mean(np.diff(conv_wave)))
    
    HDU.writeto(outnamebase+'.fits',clobber=True)

    return

def convert_to_IDL_fmt(MH, agelst=full_age, alpha=0.0):

    basedir = '/Users/Arthur/Documents/School/891_research/MILES'

    if MH < 0:
        pm = 'm'
        strmetal = MH * -1
    else:
        pm = 'p'
        strmetal = MH

    Zdict = {-2.27: 0.0001, -1.79: 0.0003, -0.66:0.004, -0.35:0.008, 0.06:0.0198, 0.4:0.04}
    Zname = Zdict[MH]/0.02
        
    flux_list = []
    norm_list = []
    outnamebase = 'MILES_IDL_{}Z_E{}'.format(''.join(str(Zname).split('.')),alpha)

    pp = PDF(outnamebase+'.pdf')
    
    for age in agelst:
        model_file = '{0:}/MILES_BASTI_CH_Ep{1:4.2f}/Mch1.30Z{2:}{3:4.2f}T{4:07.4f}_iTp{1:4.2f}_Ep{1:4.2f}.fits'.\
           format(basedir, alpha, pm, strmetal, age)
        print model_file

        wave, flux = read_fits(model_file)
        nwl = np.where((wave > 5450) & (wave < 5550))[0]

        norm = np.median(flux[nwl])
        norm_list.append(norm)
        flux_list.append(flux / norm)

        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Normalized flux')
        ax.plot(wave, flux/norm)
        pp.savefig(ax.figure)
        plt.close(ax.figure)
        
    pp.close()
    flux = np.vstack(flux_list)
    norm = np.hstack(norm_list)
    age = np.array(agelst)
    npix = wave.size
    nage = age.size

    flux = flux.reshape(1,nage,npix)
    wave = wave.reshape(1,npix)
    norm = norm.reshape(1,nage)
    age = age.reshape(1,nage)*1e9
    
    print npix, nage, flux.shape, norm.shape, age.shape, wave.shape
    
    w = pyfits.Column(name='WAVE', format='{:n}E'.format(npix), array=wave)
    f = pyfits.Column(name='FLUX', format='{:n}E'.format(npix*nage), array=flux, dim='({}, {})'.format(npix,nage))
    n = pyfits.Column(name='NORM', format='{:n}D'.format(nage), array=norm)
    a = pyfits.Column(name='AGE', format='{:n}D'.format(nage), array=age)

    BT = pyfits.BinTableHDU.from_columns([w,f,a,n])

    BT.writeto(outnamebase+'.fits', clobber=True)

    return

def all_IDL(MHlst = [-2.27, -1.79, -0.66, -0.35, 0.06, 0.4]):

    for MH in MHlst:
        convert_to_IDL_fmt(MH)

    return

def do_all(vdisplst = [293.77,241.92,199.58,187.38,179.47,
                       391.80,343.56,257.51,249.95,244.71,
                       470.14,428.06,313.04,302.67,297.51],
           MHlst = [-0.35, 0.4], alphalst = [0.0, 0.4],
           agelst = [0.35, 0.6, 1, 3, 4, 8, 12]):

    for vd in vdisplst:
        for mh in MHlst:
            for a in alphalst:
                combine_Z(mh, vd, alpha=a, agelst = agelst)


    return

def spy_to_mab():

    MHlst = ['m0.35','p0.4']
    #ages = np.array([0.35, 0.6, 1, 3, 4, 8, 12])
    ages = np.array(full_age)
    
    Z04_lst = []
    Z2_lst = []
    
    for g in range(3):
        solar = pyfits.open('MILES_E0.0_group{}_spy.fits'.format(g+1))[0].data
        enhanced = pyfits.open('MILES_E0.4_group{}_spy.fits'.format(g+1))[0].data

        Z04_lst.append(solar[:,0,2])
        Z04_lst.append(solar[:,0,3])
        Z04_lst.append(enhanced[:,0,2])
        Z04_lst.append(enhanced[:,0,3])

        Z2_lst.append(solar[:,1,2])
        Z2_lst.append(solar[:,1,3])
        Z2_lst.append(enhanced[:,1,2])
        Z2_lst.append(enhanced[:,1,3])

    Z04_data = np.vstack(Z04_lst).T
    Z2_data = np.vstack(Z2_lst).T

    print Z04_data.shape
    print Z2_data.shape

    print [ages[0]]+Z04_data[0].tolist()
    with open('mgfe_MILES.dat_Z04','w') as f:
        for i in range(ages.size):
            f.write('0.4 '+str('{:8.3f}'*13+'\n').format(*([ages[i]]+Z04_data[i].tolist())))

    with open('mgfe_MILES.dat_Z2','w') as f:
        for i in range(ages.size):
            f.write('2.0 '+str('{:8.3f}'*13+'\n').format(*([ages[i]]+Z2_data[i].tolist())))

    return solar            
