#! /usr/bin/env python

import pyfits
import numpy as np
import numpy.ma as ma
import sys
import os
from scipy.interpolate import interp1d
import pyspeckit as psk
import matplotlib.pyplot as plt
import scipy.optimize as spo
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages as PDF
import time
import ADEUtils as ADE
import bottleneck as bn

centlambda = [4901.416,5048.126]
tau = np.pi*2.

def apextract(filename, errorimage, apcenters, nrows):
    '''extracts apertures from a rectified SALT image. This is designed to
    perform the same function as IRAF's apextract routine but without all
    the fluff of tracing etc. because a rectified SALT image should not
    need to be traced (that's what recitification is for).

    Inputs:
        filename - the name of the rectified SALT image (or any image) you want
        to extract apertures from

        errorimage - an image the same size as the data image that contains
        error estimates for each pixel. See mkerr for for info

        apcenters - a python list containing the rows corresponding to the
        center of each aperture. The number of apertures extracted will be
        equal to len(apcenters).

        nrows - int. The number of rows to sum over when extracting apertures

    Output:
        A FITS file that _should_ be identical to IRAF's own .ms.fits files.
    '''
    
    hdu = pyfits.open(filename)[0]
    head = hdu.header
    data = hdu.data
    error = pyfits.open(errorimage)[0].data

    apertures = []
    erraps = []


    apnum = 1
    for r in apcenters:
        r1 = r - int(nrows/2.0)
        r2 = r1 + nrows

        print "Extracing from rows {} to {}".format(r1,r2)
        
        head.update('APNUM'+str(apnum),'{} {} {} {}'.format(apnum, apnum, r1, r2))
        apnum += 1

        print np.mean(data[r1:r2+1,:],axis=0)
        apertures.append(np.mean(data[r1:r2+1,:],axis=0))
        erraps.append(np.sqrt(np.sum(np.abs(error[r1:r2+1,:]),axis=0))/nrows)

    data_output_list = np.vstack(apertures)
    error_output_list = np.vstack(erraps)

    if data_output_list.shape[0] == 1:
        data_output_list = np.squeeze(data_output_list)
        error_output_list = np.squeeze(error_output_list)

    outname = filename.split('.fits')[0]+'.ms.fits'
    erroutput = filename.split('.fits')[0]+'_error.ms.fits'
    pyfits.PrimaryHDU(error_output_list,head).writeto(erroutput,clobber=True)

    head.update('SEPERR',True,comment='Error vectors are in a separate file')
    pyfits.PrimaryHDU(data_output_list,head).writeto(outname,clobber=True)

def meatspin(specfile,inguess,tied=None,interact=False,fig_path='./specfigs'):
    '''Takes a multi-spectrum fits file (.ms.fits, probably produced by
    apextract) and fits emission lines to each spectrum in the file. This is
    really an interactive wrapper to PySpecKit's specfit functions. 

    Inputs:
        specfile - the multi-spectrum fits file you want to fit
        
        inguess - a python list that has initial guesses for the parameters
        you want to fit. Each line has three parameters corresponding to
        increasing moments:
            0 - amplitude
            1 - centroid
            2 - width
        For example, if you want to fit two lines then your inguess should look
        something like:
            [amp1,cent1,sig1,amp2,cent2,sig2]

        tied - None or a python list. If you want any parameters to depend on
        each other you need to put a STRING CONTAINING THE DEPENDENCY as an
        element in this list. Use p[x] to describe the x'th parameter.
        i.e. if amp2 should be 1/3 amp1 (like [OIII]) then tied will be 
        (assuming you're fitting only two lines):
            ['','','','p[0]/3.','','']

        interact - boolean. If True you will have an opportunity to fine tune
        the fit for each aperture.

        figpath - path to directory where the fit for each individual aperture
        is saved. This directory MUST already exist.
    '''
    
    guesses = inguess[:]

    if tied==None: tied = ['']*len(guesses)

    fitmin = min(guesses[1::3]) - 30.
    fitmax = max(guesses[1::3]) + 30.

    hdus = pyfits.open(specfile)[0]

    header = hdus.header
    try:
        seperr = header['SEPERR']
        errorfile = specfile.split('.ms.fits')[0] + '_error.ms.fits'
        print "taking errors from {}".format(errorfile)
        ehdus = pyfits.open(errorfile)[0]
    except KeyError:
        seperr = False

    print seperr
    if seperr:
        numspec = hdus.data.shape[0]
    else:
        numspec = hdus.data.shape[0]/2

    fit_pars = np.zeros((numspec,len(inguess)))
    fit_errs = np.zeros(fit_pars.shape)
    moments = np.zeros((numspec,3))

    for i in np.arange(numspec):
        
        if seperr: 
            # Dear authors of pyspeckit,
            # I SHOULDN'T HAVE TO DO THIS!
            #
            # p.s. write some legible code
            dv = header['CDELT1']
            v0 = header['CRVAL1']
            p3 = header['CRPIX1']
            xarr = (np.arange(hdus.data.shape[1]) - p3 + 1)*dv + v0
            spec = psk.Spectrum(xarr = xarr,
                                data = hdus.data[i],
                                error = ehdus.data[i],
                                xarrkwargs = {'unit':'angstroms'},
                                header = header)
        else:
            spec = psk.Spectrum(specfile,specnum=i*2,errspecnum=i*2+1)
        
        infostr = header['APNUM'+str(i+1)].split(' ')
        print "\n  Fitting aperture {0} (lines {1} to {2})".format(infostr[0],
                                                                   infostr[2],
                                                                   infostr[3])
        fig_name = fig_path+'/{:02n}'.format(int(infostr[0]))+'_'+\
            str(guesses[1])[0:7]+'.pdf'
        pp = PDF(fig_name)

        #changed this as of 12.7 to baselineorder = 2
        spec.plotter(figure=0,xmin=fitmin,xmax=fitmax,errstyle='fill',linestyle='-')
        spec.plotter.figure.show()
        spec.baseline(order=3,fit_plotted_area=True)
        spec.specfit(guesses=guesses,tied=tied,negamp=False,fit_plotted_area=True)
        spec.plotter(xmin=fitmin,xmax=fitmax,errstyle='fill',linestyle='')
#        spec.specfit(guesses=guesses,tied=tied,negamp=False)
        spec.specfit.plot_fit(linestyle='-')

        
        if np.abs(np.average(np.array(guesses[1::3]) - np.array(spec.specfit.modelpars[1::3]))) > 2.0:
            scratch = 'd'
            print 'BAD!'
        else:
            scratch = ''

        if interact:
            scratch = raw_input("'q' moves to next line\n'g' redefines guesses\n")
            
            while scratch != 'q':
                if scratch == 'g': 
                    while True:
                        print "Guesses are:\n"+''.join('{:^9n}'.format(j) for j in range(len(guesses)))
                        print '-'*9*len(guesses)
                        print ''.join("{:^9n}".format(k) for k in guesses)+'\n'
                        cidx = raw_input("Change guess # ('r' to refit):  ")
                        if cidx == 'r': break
                        cval = float(raw_input("To:  "))
                        guesses[int(cidx)] = cval
                            
                spec.specfit(guesses=guesses,tied=tied,negamp=False,fit_plotted_area=True)
                            
                if scratch == 'd':
                    break
                    
                if scratch == 'Q':
                    interact=False
                    break

                scratch = raw_input("'q' moves to next line\n'g' redefines guesses\n")

        if scratch != 'd': guesses = spec.specfit.modelpars

        fit_pars[i] = spec.specfit.modelpars
        fit_errs[i] = spec.specfit.modelerrs
        
        center = spec.specfit.modelpars[1]
        std = spec.specfit.modelpars[2]
        
        spec_moments = meat_moment(spec)
        moments[i] = spec_moments
       
        ax = spec.plotter.figure.gca()
        ax.set_title('Aperture {0:n} in {1} on\n'.format(int(infostr[0]),specfile)+datetime.now().isoformat(' '))
        ax.set_xlim(center-10.*std,center+10.*std)
        ax.text(0.05,0.95,
                '$\mu$= {1:4.4f}\n$\mu_2$= {2:4.4f} $\Rightarrow\sigma$= {0:4.4f}\n$\mu_3$= {3:4.4f}'\
                    .format(spec_moments[1]**0.5,*spec_moments),transform=ax.transAxes,ha='left',va='top')

        pp.savefig(spec.plotter.figure)
        pp.close()

    return (fit_pars, fit_errs, moments)

def meat_moment(spec):
    
    center = spec.specfit.modelpars[1]
    std = spec.specfit.modelpars[2]

    cdf_minidx = spec.xarr.x_to_pix(center - 10.*std)
    cdf_maxidx = spec.xarr.x_to_pix(center + 10.*std) + 1 #b/c we want to include this point
    
    cropped_spec = spec.specfit.spectofit[cdf_minidx:cdf_maxidx]
    speccdf = np.cumsum(cropped_spec)
    speccdf /= speccdf.max()
    
    try:
        moment_minidx = np.where(speccdf <= 0.05)[0][-1]
        moment_maxidx = np.where(speccdf >= 0.95)[0][0]
    # moment_minidx = int(np.interp(0.05,speccdf,np.arange(speccdf.size)))
    # moment_maxidx = int(np.interp(0.95,speccdf,np.arange(speccdf.size))) + 1
    except IndexError:
        return np.array([0,0,0])

    ax = spec.plotter.figure.gca()
    ax.axvline(x=spec.xarr[cdf_minidx:cdf_maxidx][moment_minidx],linestyle=':')
    ax.axvline(x=spec.xarr[cdf_minidx:cdf_maxidx][moment_maxidx],linestyle=':')

    moment_spec = cropped_spec[moment_minidx:moment_maxidx]
    moment_lambda = np.array(spec.xarr[cdf_minidx:cdf_maxidx][moment_minidx:moment_maxidx])
    return ADE.ADE_moments(moment_lambda,moment_spec,threshold=np.inf)
    

def make_curve(specimage, radii,guesses,outputfile,tied=[],\
                   interact=False):
    '''Takes a rectified SALT image and extracts some apertures and fits
    some lines. It will try to fit all lines you give it simultaneously and
    so should be used cautiously.
    It has been largely replaced by slayer (see below).
    '''
    

    apextract(specimage,radii,radii[1] - radii[0])

    specfile = specimage.split('.fits')[0]+'.ms.fits'

    fit_pars, fit_errs = meatspin(specfile,guesses,tied=tied,interact=interact)
    
    phead = pyfits.PrimaryHDU(None)
    datahead = pyfits.ImageHDU(fit_pars)
    errorhead = pyfits.ImageHDU(fit_errs)
    datahead.header.update('EXTNAME','FIT')
    errorhead.header.update('EXTNAME','ERROR')
    
    pyfits.HDUList([phead,datahead,errorhead]).writeto(outputfile,clobber=True)

    return

def slayer(specimage,errimage,radii,guesses,outputfile,nrows=False,interact=False):
    '''Takes a rectified SALT image and extracts some apertures and fits some
    lines. Unlike make_curve, each line is fit seperately which is nice when
    some of you lines suck. This is currently the perfered method.
    '''

    specfile = specimage.split('.fits')[0]+'.ms.fits'

    if not nrows: 
        nrows = radii[1] - radii[0]

    apextract(specimage,errimage,radii,nrows)

    total_results = []
    total_errs = []
    total_moments = []

    numlines = len(guesses)/3

    for l in range(numlines):

        lineguess = guesses[l*3:(l+1)*3]
        print lineguess

        linepars, lineerrs, linemoments = meatspin(specfile,lineguess,interact=interact,tied=['','',''])

        total_results.append(linepars)
        total_errs.append(lineerrs)
        total_moments.append(linemoments)

    phead = pyfits.PrimaryHDU(None)
    datahead = pyfits.ImageHDU(np.hstack(total_results))
    errorhead = pyfits.ImageHDU(np.hstack(total_errs))
    momenthead = pyfits.ImageHDU(np.hstack(total_moments))
    radiihead = pyfits.ImageHDU(np.array(radii))
    phead.header.update('SLAYDATE',datetime.now().isoformat(' '))
    phead.header.update('SPECIM',specimage,comment='Input 2D spectrum')
    phead.header.update('ERRIM',errimage,comment='Input error spectrum')
    phead.header.update('NROWS',nrows,comment='Number of rows summed per aperture')
    for l in range(numlines):
        phead.header.update('L{}M0'.format(l),guesses[l*3],comment=\
                                'Line {} moment 0 guess'.format(l))
        phead.header.update('L{}M1'.format(l),guesses[l*3+1],comment=\
                                'Line {} moment 1 guess'.format(l))
        phead.header.update('L{}M2'.format(l),guesses[l*3+2],comment=\
                                'Line {} moment 2 guess'.format(l))
        
    datahead.header.update('EXTNAME','FIT')
    errorhead.header.update('EXTNAME','ERROR')
    momenthead.header.update('EXTNAME','MOMENTS')
    radiihead.header.update('EXTNAME','RADII')
    pyfits.HDUList([phead,datahead,errorhead,radiihead,momenthead]).writeto(outputfile,clobber=True)

    return

def gravity_gun(specimage,errimage,template_image,outputfile,combinedfile,interact=False,addonly=False):
    '''for combining multiple frames of the same data. I.e. the same scale
    height taken on different nights.  As of now it does not support rising
    and setting tracks being combined. Sorry!
    '''

    thdus = pyfits.open(template_image)
    tradii = thdus['RADII'].data.tolist()
    tfit = thdus['FIT'].data
    tguess = thdus['FIT'].data[int(len(tradii)/2)].tolist()
    terr = thdus['ERROR'].data

    if not addonly: slayer(specimage,errimage,tradii,tguess,outputfile,interact=interact)
    
    nhdus = pyfits.open(outputfile)
    nfit = nhdus['FIT'].data
    nerr = nhdus['ERROR'].data

    phead = pyfits.PrimaryHDU(None)
    fithead = pyfits.ImageHDU(np.mean(np.dstack((nfit,tfit)),axis=2))
    comberr = ((nerr/2.)**2 + (terr/2.)**2)**0.5
    errhead = pyfits.ImageHDU(comberr)
    radiihead = pyfits.ImageHDU(np.array(tradii))
    phead.header.update('SLAYDATE',datetime.now().isoformat(' '))
    phead.header.update('URIMAGE',template_image.split('/')[-1])
    phead.header.update('NEWIMAGE',specimage)
    fithead.header.update('EXTNAME','FIT')
    errhead.header.update('EXTNAME','ERROR')
    radiihead.header.update('EXTNAME','RADII')
    pyfits.HDUList([phead,fithead,errhead,radiihead]).writeto(combinedfile,clobber=True)

    return

def plot_curve(datafile,central_lambda=[4901.416,5048.126],flip=False,ax=False,label=None,hr=1):
    '''Takes a slayer output file and plots the rotation curve associated with
    the lines that were fit. Has lots of neat options for plotting.
    '''

    kpcradii, avg_centers, std_centers = openslay(datafile,central_lambda=central_lambda,flip=flip)    

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else: fig = False

    # ax2 = ax.twiny()
    # ax2.set_xlim(arcsecradii[0],arcsecradii[-1])
    # ax2.set_xlabel('Arcsec from center of galaxy')
    if hr == 1: ax.set_xlabel('Radius [kpc]')
    else: ax.set_ylabel('Radius [$r/h_r$]')
    ax.set_ylabel('LOS velocity [km/s]')

    ax.errorbar(kpcradii/hr,avg_centers,yerr=std_centers,fmt='.',label=label)

    ax.axvline(x=0,ls='--',color='k',alpha=0.3)
    ax.axhline(y=0,ls='--',color='k',alpha=0.3)

    ax.set_xlim(-50,50)
    ax.set_ylim(-500,500)

    ax.set_title(datafile+'\n'+datetime.now().isoformat(' '))
    
    if fig: fig.show()
    

def mkerr(image,stdimg,outimage):
    
    data = ma.array(pyfits.open(stdimg)[0].data)
    signal = pyfits.open(image)[0].data

    for i in range(10):
        std = np.array([np.std(data,axis=1)]).T
        badidx = np.where(np.abs(data) > 3*std)
        
        data[badidx] = ma.masked

    std = np.array([np.std(data,axis=1)]).T

    error = np.sqrt(signal + std**2)
    error[np.isnan(error)] = 999.

    pyfits.PrimaryHDU(error).writeto(outimage,clobber=True)

def dirtydown(img,outimg,HDU=0,axis=0):
    
    hdus = pyfits.open(img)
    data = hdus[HDU].data
    means = np.mean(data,axis=axis)
    data /= means
    pyfits.PrimaryHDU(data,hdus[HDU].header).writeto(outimg)


def openslay(datafile,central_lambda=[4901.416,5048.126],flip=False,moments=False):

    hdus = pyfits.open(datafile)

    pars = hdus[1].data
    errs = hdus[2].data
    pxradii = hdus[3].data

    centers = pars[:,1::3]
    amps = pars[:,::3]
    centerr = errs[:,1::3]

    velocenters = (centers-central_lambda)/central_lambda*3e5
    veloerrors = centerr/central_lambda*3e5

    avg_centers = np.sum(amps*velocenters,axis=1)/np.sum(amps,axis=1)
    std_centers = np.std(velocenters,axis=1)

    offset = helpoff(pxradii,avg_centers)
    #offset = 271.750855446
    print "Offset is "+str(offset)
    
    kpcradii = pxradii - offset
    kpcradii *= 0.118*8. # 0.118 "/px (from KN 11.29.12) times 8x binning
    kpcradii *= 34.1e3/206265. # distance (34.1Mpc) / 206265"
    
    if flip: kpcradii *= -1.
    badidx = np.where(std_centers > 500.)
    kpcradii = np.delete(kpcradii,badidx)
    avg_centers = np.delete(avg_centers,badidx)
    std_centers = np.delete(std_centers,badidx)
    
    if moments:
        all_moments = hdus[4].data
        m1 = (all_moments[:,::3] - central_lambda)/central_lambda*3e5
        m2 = all_moments[:,1::3]
        m3 = all_moments[:,2::3]
#        m1 = np.delete(m1,badidx)
#        m2 = np.delete(m2,badidx)
#        m3 = np.delete(m3,badidx)
        return (kpcradii,avg_centers,std_centers,m1,m2,m3)

    else: return (kpcradii,avg_centers,std_centers)

def helpoff(radii,centers):

    x0 = np.array([np.median(radii)])
    
    xf = spo.fmin(offunc,x0,args=(radii,centers),disp=False)
    
    radii = radii.copy() - xf[0]
    pidx = np.where(radii >= 0.0)
    nidx = np.where(radii < 0.0)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(radii[pidx],np.abs(centers[pidx]))
    # ax.plot(np.abs(radii[nidx]),np.abs(centers[nidx]))
    # fig.show()

    return xf[0]

def offunc(x,radii,centers):
    
    radii = radii.copy() - x[0]
    pidx = np.where(radii >= 0.0)
    nidx = np.where(radii < 0.0)

    if pidx[0].size <= 1 or nidx[0].size <= 1: return 999.
    
    pcent = centers[pidx]
    ncent = centers[nidx]

    pcoef = np.polyfit(radii[pidx],np.abs(centers[pidx]),4)
    ncoef = np.polyfit(np.abs(radii[nidx][::-1]),np.abs(centers[nidx][::-1]),4)

    pf = np.poly1d(pcoef)
    nf = np.poly1d(ncoef)

    minr = max(radii[pidx].min(),np.abs(radii[nidx]).min())
    maxr = min(radii[pidx].max(),np.abs(radii[nidx]).max())

    r = np.linspace(minr,maxr,100)

    chisq = np.sum((pf(r) - nf(r))**2)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(r,pf(r))
    # ax.plot(r,nf(r))
    # fig.show()
    # raw_input('ssds')

#    print chisq
    
    return chisq

def plot_line(datafile,radius,wavelength=5048.126,ax=False,
              central_lambda=[4901.416,5048.126],flip=False,
              plot=True,window=20):
    """ Plots a single line from a .ms file. It also needs a corresponding
    .slay file to get the pixel -> kpc radius conversion.

    datafile - str. can be the name either of a .slay or .ms file. This
    function requires both files to be present in the current dirctory. If you
    extracted lines with slayer() then everything should be good automatically

    Radius - in kpc. The distance from the center of the galaxy where you want
    to plot a spectrum. Can be negative or positve (for different sides of the
    galaxy)

    wavelength - in Angstroms, the central wavelength of the plotted region

    window - in Angstroms, the range of the plotted region

    ax - a matplotlilb.axes.AcesSubplot object that will be plotted on if
    provided. If ax is provided then no additional labels or titles will be
    added to it
    
    central_lambda - python list. The rest-frame wavelengths of the line
    centers that are contained in the .slay file. Used by openslay()

    plot - boolean. Set to True to actually plot something. This is here so
    that this function can be called as a helper to just return the output
    arrays

    flip - boolean. Passed to openslay(). Do you want to flip the galaxy
    around r=0? Changing this will change the positive/negative convention of
    the radius input

    """

    if '.slay.' in datafile:
        datafile = ''.join(datafile.split('.')[:-2] + ['.ms.fits'])

    slayfile = ''.join(datafile.split('.')[:-2] + ['.slay.fits'])

    kpcradii, _, _ = openslay(slayfile,central_lambda=central_lambda,
                              flip=flip,moments=False)
    pxradii = pyfits.open(slayfile)[3].data

    row = np.where(np.abs(kpcradii-radius) == np.min(np.abs(kpcradii-radius)))[0][0]
    print "using pixel value {} where radius is {} kpc".format(pxradii[row],kpcradii[row])

    hdu = pyfits.open(datafile)[0]
    CRVAL = hdu.header['CRVAL1']
    Cdelt = hdu.header['CDELT1']
    
    # We use '=f8' to force the endianess to be the same as the local
    # machine. This is so the precompiled bottleneck (bn) functions don't
    # complain
    spectrum = np.array(hdu.data[row*2],dtype='=f8')
    error = hdu.data[row*2 + 1]
    wave = np.arange(spectrum.size)*Cdelt + CRVAL
    
    idx = np.where((wave >= wavelength - window/2.) & (wave <= wavelength + window/2.))

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wavelength [Angstroms]')
        ax.set_ylabel('ADU/s')
        ax.set_title(datetime.now().isoformat(' '))
    
    ax.errorbar(wave[idx],spectrum[idx],yerr=error[idx])
    fig = ax.figure
    
    if plot:
        fig.show()

    return wave[idx], spectrum[idx], error[idx]

def plot_row(msfile,rownum,smooth=False,ax=False):
    """
    A very simple function to plot a single row (aperture) from a .ms.fits
    file. It also allows you to do a moving-median smooth on the data to make
    plots look a little nicer. msfile is the file in question and rownum is
    whatever row you want. rownum + 1 is assumed to be the error.

    """

    hdu = pyfits.open(msfile)[0]
    exptime = hdu.header['EXPTIME']
    CRVAL = hdu.header['CRVAL1']
    Cdelt = hdu.header['CDELT1']
    
    # We use '=f8' to force the endianess to be the same as the local
    # machine. This is so the precompiled bottleneck (bn) functions don't
    # complain
    spectrum = np.array(hdu.data[rownum],dtype='=f8')
    error = hdu.data[rownum + 1]
    wave = np.arange(spectrum.size)*Cdelt + CRVAL

    if smooth:
        spectrum = bn.move_median(spectrum,smooth)

    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Wavelength [Angstroms]')
        ax.set_ylabel('ADU/s')
        ax.set_title(datetime.now().isoformat(' '))

    ax.errorbar(wave,spectrum,yerr=error)
    fig = ax.figure
    
    fig.show()

    return ax
