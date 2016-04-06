import pyfits
import re
import os
import numpy as np
import time
import ir2py_lines as i2p
import matplotlib.pyplot as plt
from pyraf import iraf
from matplotlib.backends.backend_pdf import PdfPages as PDF
import plot_simple as ps
import plot_allZ2 as pa2

plt.ioff()
iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

llist = ['HB','Ha']
#rlist = [np.array([4793., 4923.]), np.array([6485., 6685.])]
rlist = [np.array([4847., 4872.]), np.array([6510., 6625.])]
centlist = [[4861.], np.array([6563., 6549., 6585.])]

def prep_spectra(datafile, fitfile, output, velocity, smooth=3.):

    wavemin=3800.
    wavemax=6800.

    vel = np.loadtxt(velocity,usecols=(1,),unpack=True)

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    header = hdu.header
    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']
    
    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval
    idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
    
    wave = wave[idx]
    data = data[:,idx]

    fits = pyfits.open(fitfile)[0].data

    print fits.shape[1], wave.size

    diff = data - fits + 1e-16

    shift = np.vstack([np.interp(wave,wave*(1 - vel[i]/3e5),diff[i,:]) for i in range(diff.shape[0])])
    
    if smooth > 0:
        smoothed = np.vstack([np.convolve(shift[i,:],np.ones(smooth)/smooth,'same') for i in range(diff.shape[0])])
    else:
        smoothed = shift

    header.update('CRVAL1', 3800.)
    #os.system('rm tmp'+output)
    os.system('rm '+output)
    pyfits.PrimaryHDU(smoothed*1e17,header).writeto(output,clobber=True)

    # iraf.continuum('tmp'+output,output,
    #                type='ratio',
    #                interac='no',
    #                functio='spline3',
    #                order=33,
    #                low_rej=2,
    #                high_rej=2)

    return

def do_fitprof(datafile, pointing):
    """Run IRAF fitprof routine to measure line centers and widths.

    For each line specified in the module header generate the necessary inputs
    to fitprof and execture the routine. The fitting regions are also
    specified in the module header.

    Parameters
    ----------
    datafile : str
        Name of a multispec fits file to pass as input to fitprof. Must have a WCS solution in the header.

    Returns
    -------
    None
        Nothing is returned. Instead, IRAF writes the results to a file.

    """
    
    dirn = os.path.dirname(datafile)
    if len(dirn) > 0:
        dirn += '/'

    for l, cl, r in zip(llist, centlist, rlist):
        with open('{}P{}_{}.lines'.format(dirn,pointing,l),'w') as f:
            for c in cl:
                f.write('{:4.2f} INDEF g 15\n'.format(c))

        proffile = '{}P{}_{}.fitp'.format(dirn,pointing,l)
        if os.path.exists(proffile):
            print 'Removing ' + proffile
            os.system('rm '+proffile)
            os.system('rm {}P{}_{}_fits.fits'.format(dirn,pointing,l))
        try:
            iraf.fitprofs(datafile,
                          region='{:4.0f} {:4.0f}'.format(*r),
                          positio='{}P{}_{}.lines'.format(dirn,pointing,l),
                          logfile='{}P{}_{}.fitp'.format(dirn,pointing,l),
                          fitpositions='single',
                          fitgfwhm='single',
#                          plotfile='test_{}.pt'.format(l),
                          output='{}P{}_{}_fits.fits'.format(dirn,pointing,l),
                          option='fit',
                          nerrsam=100,
                          sigma0=0.3,
                          invgain=0.)
        except Exception as e:
            print 'IRAF died:'
            print e
            pass

    return

def get_results(pointing, output):
    """Parse fitprof output and display results

    The line centers are taken from the output of fitprof. For each line specified in the module header the average offset and stddev across all fibers in the IFU is computed. Output is a textfile and a plot of accuracy and stochasticity as a function of wavelength.

    Parameters
    ----------
    output : str
        Name of the output text file. This file will contain the mean, offset, stddev, and number of rejected apertures.
    threshold : float, optional
        Threshold value for iterative sigma clipping in mean across IFU. The total number of rejected fibers will be recorded in the output file.

    Returns
    -------
    None :
       The result is a text file containing the results and a plot containing the accuracy as a function of wavelength.
    
    Notes
    -----
    Right now the plot is hardcoded to be writting to WLC.png

    """
    print 'Consolidating measurements'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel('Av')
    ax.set_xlabel('Aperture number')
    ax.set_title(time.asctime())

    tHad = i2p.parse_fitprofs('P{}_Ha.fitp'.format(pointing),3)
    Haf = tHad[4][:,0]
    Hafe = tHad[5][:,0]
    tHBd = i2p.parse_fitprofs('P{}_HB.fitp'.format(pointing),1)
    HBf = tHBd[4][:,0]
    HBfe = tHBd[5][:,0]

    #Correct for underfit [NII] line
    Hacorr = tHad[4][:,1] - tHad[4][:,2]/3.
    Haf += Hacorr

    ratio = Haf/HBf
    ratio_e = np.sqrt((Hafe/HBf)**2 + (HBfe*Haf/HBf**2)**2)

    EBV = 1.97*np.log10(ratio/2.86)
    Av = 4.05*EBV
    tauV = Av/1.086
    Av_e = 1.97*4.05*ratio_e/(1.086*ratio*np.log(10))

    print 'Plotting'
    ax.errorbar(np.arange(ratio.size)+1,Av,yerr=Av_e,fmt='.')
    fig.savefig(output+'.png')

    with open(output+'.txt','w') as f:
        f.write('# Written on {}\n'.format(time.asctime()))
        f.write('#\n#{:>3}{:>10}{:>10}\n\n'.format('Ap','Tau_V','TauV_err'))
        for i in range(Av.size):
            f.write('{:4}{:10.2f}{:10.2f}\n'.format(i+1,Av[i],Av_e[i]))

    return
    
def make_balmer_model(balmerD, HBfitp, location, velocities, output,
                      dispdata = '/d/monk/eigenbrot/WIYN/14B-0456/anal/disp/GP_disp_batch_avg_int.fits'):
    
    tauV = np.loadtxt(balmerD, usecols=(1,), unpack=True)
    sizes = np.loadtxt(location, usecols=(1,), unpack=True)
    vels = np.loadtxt(velocities, usecols=(1,), unpack=True)
    numap = tauV.size
    HB_data = i2p.parse_fitprofs(HBfitp,1)
    HB_peak = HB_data[7][:,0]/1e17 #- 1e-16

    #It's been long enough I think I can just use these magic numbers
    wave = np.arange(2011)*2.1 + 3340.
    #wave = np.arange(1428)*2.1 + 3800.
    m_wave = np.logspace(np.log10(3340), np.log10(7600), 16000)
    mpix = np.mean(np.diff(m_wave)/m_wave[1:]*3e5) #size of 1 model pixel in km/s
    balmer_flux = np.zeros((numap,m_wave.size))

    #Read in disp data
    disph = pyfits.open(dispdata)[0]
    disp_head = disph.header
    disp_data = disph.data
    disp_numwave = disp_data.shape[1]
    disp_wave = np.arange(disp_numwave)*disp_head['CDELT1'] + disp_head['CRVAL1']
    disp_arr = np.zeros((5,m_wave.size))
    for d in range(5):
        disp_arr[d,:] = np.interp(m_wave,disp_wave,disp_data[d,:])/2.355

    sized = {0.937: 0,
             1.406: 1,
             1.875: 2,
             2.344: 3,
             2.812: 4}

    cents = [6563, 4861, 4341, 4102, 3970]
    ratios = [2.86, 1, 0.47, 0.25, 0.16]
    
    for c, r in zip(cents, ratios):
        idx = np.argmin(np.abs(m_wave - c))
        balmer_flux[:,idx] = 1
        
    final_model = np.zeros((numap,wave.size))
    for i in range(numap):
        print i
        wave_red = wave * (1 + vels[i]/3e5)
        if np.isnan(tauV[i]):
            red = np.ones(m_wave.size)
        else:
            red = np.exp(-1 * tauV[i]*(m_wave/5500.)**(-0.7))
        balmer_flux[i,:] *= red
        did = sized[sizes[i]]
        sigma_pix = disp_arr[did,:]/mpix
        balmer_flux[i,:] = mconv(balmer_flux[i,:],sigma_pix)
        HBid = np.argmin(np.abs(m_wave - cents[1]))
        # HBid = np.where((m_wave > cents[1] - 30) & (m_wave < cents[1] + 30))[0]
        # mHB_flux = np.sum(balmer_flux[i,HBid])
        mHB_peak = balmer_flux[i,HBid]
        print HB_peak[i], mHB_peak
        balmer_flux[i,:] *= HB_peak[i]/mHB_peak
        final_model[i,:] = np.interp(wave,m_wave,balmer_flux[i,:])
        final_model[i,:] = np.interp(wave,wave_red,final_model[i,:])

    final_hdu = pyfits.PrimaryHDU(final_model)
    final_hdu.header.update('CRPIX1',1)
    final_hdu.header.update('CRVAL1',3340.)
    final_hdu.header.update('CDELT1',2.1)
    
    final_hdu.writeto(output,clobber=True)

def mconv(y, sig):

    d1 = y.size
    x = np.arange(d1)
    yp = np.zeros(d1)
    norm = np.sqrt(2*np.pi)*sig

    for i in range(d1):
        div = ((x-i)/sig)**2
        close = np.where(div < 50)
        kern = np.exp(-0.5*div[close])/norm[close]
        yp[i] = np.sum(y[close]*kern)/np.sum(kern)
        
    return yp

def do_all(datafile, fitfile, velocity, location, smooth=3., balmer=False):

    pointing = int(re.search('_P([1-6])_',datafile).groups()[0])
    prepped = 'P{}_contsub.ms.fits'.format(pointing)
    output = 'P{}_balmerD'.format(pointing)

    prep_spectra(datafile, fitfile, prepped, velocity, smooth=smooth)
    do_fitprof(prepped,pointing)
    get_results(pointing,output)

    if balmer:
        bal_out = 'NGC_891_P{}_bin30_balmer_model.ms.fits'.format(pointing)
        sub_out = 'NGC_891_P{}_bin30_balmsub.mso.fits'.format(pointing)
        make_balmer_model(output+'.txt','P{}_HB.fitp'.format(pointing),location, velocity, bal_out)
        data = pyfits.open(datafile)[0]
        balm = pyfits.open(bal_out)[0]
        data.data -= balm.data
        data.writeto(sub_out,clobber=True)

    return

def Pbatch(balmer=False):

    for i in range(6):
        
        data = 'NGC_891_P{}_bin30.mso.fits'.format(i+1)
        fit = 'NGC_891_P{}_bin30_allz2.fit.fits'.format(i+1)
        vel = 'NGC_891_P{}_bin30_velocities.dat'.format(i+1)
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(i+1)

        print data
        print fit
        print vel
        
        do_all(data,fit,vel,loc,balmer=balmer)

    return

def compare_tau(offset=0):
    #offset = 57 for ma11

    for i in range(6):
        
        balmerdat = 'P{}_balmerD.txt'.format(i+1)
        SSPdat = 'NGC_891_P{}_bin30_allz2.dat'.format(i+1)
        print balmerdat
        print SSPdat

        balmT = np.loadtxt(balmerdat, usecols=(1,), unpack=True)
        SSPT = np.loadtxt(SSPdat, usecols=(66-offset,), unpack=True)

        diff = (balmT - SSPT)/SSPT

        with open('P{}_Tdiff.txt'.format(i+1),'w') as f:
            f.write('# Written on {}\n'.format(time.asctime()))
            f.write('#\n# Tdiff = Balmer Tau - SSP Tau\n')
            f.write('#\n#{:>3}{:>8}\n\n'.format('Ap','Tdiff'))
            for i in range(diff.size):
                f.write('{:4}{:8.2f}\n'.format(i+1,diff[i]))\

    return

def make_plots():

    exclude = [[5,34],
               [1,2,35],
               [59],
               [2,8],
               [1,2,3,27,28,29],
               [35,36,38]]

    ps.all_maps('TauV_diff_map.pdf',col=1,inputprefix='',
                inputsuffix='Tdiff.txt',minval=-5,maxval=5,
                binned=True,exclude=exclude,
                label=r'$\frac{\tau_{V,Balm} - \tau_{V,SSP}}{\tau_{V,SSP}}$')
    
    ps.all_maps('TauV_balmer_map.pdf',col=1,inputprefix='',
                inputsuffix='balmerD.txt',
                minval=-1,maxval=6,binned=True,exclude=exclude,
                label=r'$\tau_{V,Balm}$')

    al = pa2.simple_plot(inputsuffix='balmerD.txt',
                         col=1,label=r'$\tau_{V,Balm}$',
                         ylims=[-1,6],exclude=exclude)
    al2 = pa2.simple_plot(inputsuffix='Tdiff.txt',col=1,
                          label=r'$\frac{\tau_{V,Balm} - \tau_{V,SSP}}{\tau_{V,SSP}}$',
                          ylims=[-5,5],exclude=exclude)

    pp = PDF('TauV_balm_heights.pdf')
    [pp.savefig(a.figure) for a in al]
    pp.close()
    
    pp2 = PDF('TauV_diff_heights.pdf')
    [pp2.savefig(a.figure) for a in al2]
    pp2.close()
    
    plt.close('all')
    return

def plot_fits(datafile, fitfile, fitp, output, numlines):

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    header = hdu.header
    
    cdelt = header['CDELT1']
    crpix = header['CRPIX1']
    crval = header['CRVAL1']
    
    print cdelt

    wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval
    
    fits = pyfits.open(fitfile)[0].data
    
    d = i2p.parse_fitprofs(fitp,numlines)

    if numlines == 3:
        x = np.arange(rlist[1][0], rlist[1][1])
    else:
        x = np.arange(rlist[0][0], rlist[0][1])

    idx = np.where((wave >= x.min()) & (wave <= x.max()))[0]
    
    wave = wave[idx]
    data = data[:,idx]
    fits = fits[:,idx]

    pp = PDF(output)
    for ap in range(d[1].shape[0]):
        print ap+1
        y = np.zeros(wave.size)
        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Flux')
        ax.set_title('Ap {}'.format(ap+1))
        back = []
        for i in range(numlines):
            back.append(d[6][ap,i])
            tmp = np.exp(-0.5*(wave - d[1][ap,i])**2/(d[2][ap,i]/2.355)**2)
            tmp *= d[7][ap,i]/np.max(tmp)
#            tmp *= d[4][ap,i]/np.sum(tmp)
            ax.plot(wave,tmp+back[-1],'r:')
            y += tmp
    
        if numlines == 3:
            ratio = d[7][ap,2]/d[7][ap,1]
            ax.text(0.8,0.8,'ratio = {:4.3f}'.format(ratio),transform=ax.transAxes)

        ax.plot(wave,y+np.mean(back),'b',alpha=0.7)
        ax.plot(wave,data[ap,:],'k')
        ax.plot(wave,fits[ap,:]+np.mean(back),'r')
        pp.savefig(ax.figure)
        
    pp.close()
    plt.close('all')

    return

def fit_batch():
    
    for i in range(6):

        plot_fits('P{}_contsub.ms.fits'.format(i+1),
                  'P{}_Ha_fits.fits'.format(i+1),
                  'P{}_Ha.fitp'.format(i+1),
                  'P{}_Ha_profs.pdf'.format(i+1),3)

        plot_fits('P{}_contsub.ms.fits'.format(i+1),
                  'P{}_HB_fits.fits'.format(i+1),
                  'P{}_HB.fitp'.format(i+1),
                  'P{}_HB_profs.pdf'.format(i+1),1)

    return
