import pyfits
import numpy as np
import sys
import os
from scipy.interpolate import interp1d
import scipy.optimize as spo
import matplotlib.pyplot as plt
import scipy.optimize as spo
import NaCl as na
import gc
import ADEUtils as ADE
from PyGalPot import PyGalPot as PGP
from datetime import datetime
#from matplotlib.backends.backend_pdf import PdfPages as PDF
import time
from ADESALT import openslay, plot_curve

centlambda = [4901.416,5048.126]
tau = np.pi*2.

def gen_models(zlist,fraclist,SALTdata,hrot=5.45):
    '''Produces a model XV diagram for each combination of height and disk
    fraction specified by the user. Model tangential velocity curves are 
    constructed with GalPot and then "observed" to produce the XV diagrams.
    Information on V_c and h_rot are taken from SALT data'''

    SALTr, SALTv, SALTerr = openslay(SALTdata)

    V_c = find_Vc(SALTr, SALTv, SALTerr)

    print "V_c is {:4.3f} km/s".format(V_c)

    modelr = np.linspace(0,np.abs(SALTr).max(),100)
    modelv = V_c*np.tanh(modelr/hrot)

    if type(fraclist[0]) != float:
        frac_results = fraclist[:]
    
    else:
        frac_results = []
        print "Generating galaxies..."
        for frac in fraclist:
        
            galpars = na.make_galaxy(frac,rot_curve=(modelr,modelv))
        #        galpars = [3.363e8, 7.171e6,22.759]
            frac_results.append(PGP(galpars))

            print "done."
    plotlist = []
    print 'Generating plots...'
    for height in zlist:
        
        plotlist.append(plt.figure())
        ax = plotlist[-1].add_subplot(111)
        ax.set_title('Height = {:4.3f} kpc'.format(height))
        ax.set_xlabel('$r$ [kpc]')
        ax.set_ylabel('$V(r)$ [km/s]')
        for galaxy, diskfraction in zip(frac_results,[0.4,0.6,0.8]):
            zr, zv = galaxy.get_TVC(height)
            xvr, xvv, _ = simcurve(1001,height,zv[-1],hrot,0.0,5,0.5,np.pi,
                                   kappa_0=1.652,z_d=0.245,
                                   scale=(SALTr.max()-SALTr.min())/1001.,
                                   rot_curve=(zr,zv))
            ax.plot(xvr,xvv,label=str(diskfraction))
            ax.plot(zr,zv,'b:')
            ax.plot(modelr,modelv,'k:')
        ax.set_xlim(0,35)
        ax.set_ylim(0,270)
        ax.legend(loc=0,title='Disk maximality')
#        plotlist[-1].show()
    print 'done'

    return frac_results, plotlist
        
def find_Vc(r, v, err, back=True):

    frontidx = [0,1] ##
    frontstd = np.std(v[frontidx]) ##
    fronterr = 999. ##

    '''commented out to deal with one-sided line_profile data.
    Uncomment to use on actual data'''
    backidx = [-1,-2]
    backstd = np.std(v[backidx])
    backerr = 999.

    while frontstd < 2*fronterr: ###
        frontidx.append(frontidx[-1]+1) ##
        frontstd = np.std(v[frontidx]) ##
        fronterr = ((np.sum(err[frontidx]**2))**0.5)/len(frontidx) ##

    '''this should be 2*backerr for real data'''
    while backstd < 2*backerr:        
        backidx.append(backidx[-1]-1)
        backstd = np.std(v[backidx])
        backerr = ((np.sum(err[backidx]**2))**0.5)/len(backidx)
    
    if back:
        finalv = (fronterr*np.abs(np.mean(v[frontidx])) + backerr*np.abs(np.mean(v[backidx])))/(fronterr + backerr)
    # finalv = np.mean(v[backidx])
    else:
        finalv = np.mean(v[frontidx])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(r,v,yerr=err,fmt='.')
    ax.axvline(x=r[frontidx[-1]]) ##
    ax.axvline(x=r[backidx[-1]])
    ax.axhline(y=finalv)
    ax.axhline(y=-1*finalv) ##
    fig.show()
    
    return finalv


def simcurve(size,Z,v_r,h_rot,
             ax=False,scale=1.,
             kappa_0=1.62,z_d=0.245,label='',rot_label=False,
             p=False,rot_curve=False,output='test.fits',
             spiralpars=None,flarepars=None,ringpars=None,warppars=None):

    #Z is in kpc, not scale heights


    #first generate the radial distance array. Taken from ADEUtils
    vecvec = np.indices((size,size),dtype=np.float32)
    distances = scale*((size/2. - vecvec[0,:])**2 + (size/2. - vecvec[1,:])**2)**0.5

    #The angle array will be used to figure out the los velocity for each bin.
    # taken from ADEUtils
    angles = np.abs(np.arctan2(size/2. - vecvec[0],size/2. - vecvec[1]) + np.pi)

    #some values from Xilouris '99 for ESO 435-G25 (IC 2531)
#    kappa_0 = 0.652 #=\tau/2z_d units are kpc^-1
    h_d = 8.43 #kpc
    h_z = 0.43 #kpc scale height
#    z_d = 0.23 #kpc dust scale height

    #N needs to be an integer
    # if N > 7.0: N = 8
    # elif N < 5.0: N = 4
    # else: N = 6

    #This array holds the opacity of each bin
    kaparray = np.exp(-1*(distances/h_d))
    kaparray *= kappa_0 * np.exp(-1*(np.abs(Z)/z_d)) / kaparray[size/2,size/2]
    # And this one has the surface brightness of each bin, assuming a
    # doubly-exp disc 
    # total normalization is irrelevant
    Iarray = np.exp(-1*(distances/h_d))
    Iarray *= np.exp(-1*(np.abs(Z)/h_z))

    # Now add whatever morphological extras the user desires
    if spiralpars:
        spiral = LSP(distances,angles,**spiralpars)
        kaparray *= spiral
        Iarray *= spiral

    if flarepars:
        kaparray /= np.exp(-1*(np.abs(Z)/h_z))
        Iarray /= np.exp(-1*(np.abs(Z)/h_z))
        if flarepars['ftype'] == 'exp':
            flare = disco(distances,Z,h_z,**flarepars)
        elif flarepars['ftype'] == 'linear':
            flare = quickmatch(distances,Z,h_z,**flarepars)
        else:
            print "Flare type not recognized. Accepted types are 'exp' or 'linear'."
            flare = 1.
        kaparray *= flare
        Iarray *= flare
    if ringpars:
        ring = bigben(distances,**ringpars)
        totalI = np.sum(Iarray)
        totalkap = np.sum(kaparray)
        Iarray *= (1 - ringpars['r_w'])
        Iarray += ring*totalI
        kaparray *= (1 - ringpars['r_w'])
        kaparray += ring*totalkap

    if warppars:
        kaparray /= np.exp(-1*(np.abs(Z)/h_z))
        Iarray /= np.exp(-1*(np.abs(Z)/h_z))
        warp = star_trek(distances,angles,Z,h_z,**warppars)
        Iarray *= warp
        kaparray *= warp

    #for each sight line, we'll now figure out the relative light
    # contributions from each bin along the sight line
    kapcum = np.cumsum(kaparray,axis=0)
    tauarray = scale*kapcum
       
    #This array will hold the relative contribution from each bin
    fracarray = Iarray*np.exp(-1*tauarray)

    #Compute the rotation curve either with a tanh model, or from a provided
    # 1D curve
    if not rot_curve: 
        TVC = v_r*np.tanh(distances/h_rot)
    else:
        def comp_TVC(rr): return np.interp(rr,rot_curve[0],rot_curve[1])
        vcomp = np.vectorize(comp_TVC)
        TVC = vcomp(distances)
    
    #Now let's compute the projected los velocity of each bin
    v_sarray = TVC*np.cos(angles)

    #finally, the light-weighted contribution to total LOS velocity
    LOSfracarray = (v_sarray*fracarray)/np.sum(fracarray,axis=0)

    if ax:
        radii_vec = scale*(np.arange(size)-size/2)
        rot_curve_r = radii_vec[np.where(radii_vec > 0)]
        rot_curve_vec = v_r*np.tanh(rot_curve_r/h_rot) ## this isn't right
        ax.plot(radii_vec,np.sum(LOSfracarray,axis=0),label=label)
        ax.plot(rot_curve_r,rot_curve_vec,':',label=rot_label)
        ax.set_xlim(-50,50)
        ax.set_ylim(-500,500)
        if p:
            s = '$Z/h_z$: {0:3.3f}\n$V_r$: {1:5.3f}\n$h_{{rot}}$: {2:4.3f}\n$w$: {3:4.3f}\n$N$: {4:4.3f}\n$p$: {5:4.3f}\n$\\theta_{{view}}$: {6:4.3f}\n$\kappa_0$: {7:4.3f}\n$z_d$: {8:4.3f}'\
                .format(Z/h_z, v_r, h_rot, w, N, pitch,view_ang,kappa_0,z_d)
            ax.text(15, -20, s,horizontalalignment='left',va='top',fontsize=12)
        

    frachdu = pyfits.PrimaryHDU(fracarray)
    tauhdu = pyfits.ImageHDU(tauarray)
    kaphdu = pyfits.ImageHDU(kaparray)
    Ihdu = pyfits.ImageHDU(Iarray)
    vshdu = pyfits.ImageHDU(v_sarray)
    LOShdu = pyfits.ImageHDU(LOSfracarray)
    dhdu = pyfits.ImageHDU(distances)
    ahdu = pyfits.ImageHDU(angles)
    rothdu = pyfits.ImageHDU(TVC)
    
    frachdu.header.update('Z',Z,comment='kpc')
    frachdu.header.update('h_rot',h_rot,comment='kpc')
    frachdu.header.update('v_r',v_r,comment='v_r')
    frachdu.header.update('z_d',z_d,comment='Dust scale height in kpc')
    frachdu.header.update('h_z',h_z,comment='Gas scale height in kpc')
    frachdu.header.update('h_d',h_d,comment='Gas and dust scale length in kpc')
    frachdu.header.update('scale',scale,comment='pixel scale in kpc/px')

    hdulist = [frachdu, tauhdu, kaphdu, Ihdu, 
               vshdu, LOShdu, dhdu, ahdu, rothdu]

    if spiralpars:
        spiralhdu = pyfits.ImageHDU(spiral)
        spiralhdu.header.update('EXTNAME','SPIRAL')
        hdulist.append(spiralhdu)
        frachdu.header.update('w',spiralpars['w'],comment='Spiral weight')
        frachdu.header.update('N',spiralpars['N'],comment='Number of spiral arms')
        frachdu.header.update('pitch',spiralpars['pitch'],comment='Spiral wind degree')
        frachdu.header.update('VIEWANG',spiralpars['view_ang'],comment='Viewing angle')

    if flarepars:
        flarehdu = pyfits.ImageHDU(flare)
        flarehdu.header.update('EXTNAME','FLARE')
        hdulist.append(flarehdu)
        frachdu.header.update('FTYPE',flarepars['ftype'],comment='Type of flare')
        frachdu.header.update('h_zR',flarepars['h_zR'],comment='Scale height scale length [kpc]')

    if ringpars:
        ringhdu = pyfits.ImageHDU(ring)
        ringhdu.header.update('EXTNAME','RING')
        hdulist.append(ringhdu)
        frachdu.header.update('r_R',ringpars['r_R'],comment='Ring radius [kpc]')
        frachdu.header.update('r_sig',ringpars['r_sig'],comment='Ring width [kpc]')
        frachdu.header.update('r_w',ringpars['r_w'],comment='Ring strength')

    if warppars:
        warphdu = pyfits.ImageHDU(warp)
        warphdu.header.update('EXTNAME','WARP')
        hdulist.append(warphdu)
        frachdu.header.update('warpfac',warppars['warp_factor'],comment='Warp factor')
        frachdu.header.update('warpang',warppars['warp_ang'],comment='Warp angle')

    # Add WCS coordinates (kpc) to headers
    for HDU in hdulist:

        HDU.header.update('CRPIX1',size/2,comment='WCS: X reference pixel')
        HDU.header.update('CRPIX2',size/2,comment='WCS: Y reference pixel')
        HDU.header.update('CRVAL1',0.0,
                              comment='WCS: X reference coordinate value')
        HDU.header.update('CRVAL2',0.0,
                              comment='WCS: Y reference coordinate value')
        HDU.header.update('CDELT1',scale,comment='WCS: X pixel size')
        HDU.header.update('CDELT2',scale,comment='WCS: Y pixel size')
        HDU.header.update('CTYPE1','LINEAR',comment='X type')
        HDU.header.update('CTYPE2','LINEAR',comment='Y type')

    if rot_curve: frachdu.header.update('rotcurve','yes')
    else: frachdu.header.update('rotcurve','False')

    frachdu.header.update('EXTNAME','FRAC')
    tauhdu.header.update('EXTNAME','TAU')
    kaphdu.header.update('EXTNAME','KAP')
    Ihdu.header.update('EXTNAME','SB')
    vshdu.header.update('EXTNAME','V_S')
    LOShdu.header.update('EXTNAME','LOSFRAC')
    dhdu.header.update('EXTNAME','DIST')
    ahdu.header.update('EXTNAME','ANG')
    rothdu.header.update('EXTNAME','ROT')

    pyhdus = pyfits.HDUList(hdulist)
    pyhdus.writeto(output,clobber=True)
    pyhdus.close()

    return scale*(np.arange(size)-size/2), np.sum(LOSfracarray,axis=0), TVC

def LSP(distances, angles, w, N, pitch, view_ang):
    '''Generates an array that redistributes light into spiral arms via the
    prescription of ASR 2012
    '''

    prodlist = []
    
    for n in np.arange(2,N+1,2):
        sinarray = (n*w)/(n-1) * np.sin( np.log(distances)/np.tan(pitch) - angles +\
                                             view_ang)**N
        prodlist.append(sinarray)


    return 1 - w + np.array(prodlist).prod(axis=0)
    
def disco(distances, Z, h_z, **flarepars):
    '''Generates an array that creates a galaxy where the scale height depends
    on r via  h_z(R) = exp(R/h_zR), i.e. a flare
    '''
    
    h_zR = flarepars['h_zR']
    h_zprime = np.exp(distances/h_zR)
    ideal_flare = np.exp(-1*Z / h_zprime) # The ideal flare formulation

    return ideal_flare * h_z / h_zprime
    

def quickmatch(distances, Z, h_z, **flarepars):
    '''Generates an array that is used to redistribute the light in a galaxy
    into a linear flare. A linear flare has a linearly increasing scale
    height, as opposed to the flares produced by disco, which exponentially
    increase.

    In this case h_z(R) = h_z + h_zR*R
    '''
    
    h_zR = flarepars['h_zR']
    h_zprime = h_z + h_zR*distances
    ideal_flare = np.exp(-1*Z / h_zprime)

    return ideal_flare * h_z / h_zprime

def bigben(distances, **ringpars):
    '''Generates an array that can be used to redistribute light into a
    ring. The ring is parametrized as a gaussian in r such that:

    I(r) = exp(-(r - r_R)**2/(2*r_sig**2))
    '''

    r_R = ringpars['r_R']
    r_sig = ringpars['r_sig']
    r_w = ringpars['r_w']

    ring = np.exp(-1*(distances - r_R)**2/(2*r_sig**2))
    total = np.sum(ring)

    return ring*r_w/total

def star_trek(distances, angles, Z, h_z, **warppars):
    '''Generates an array that can be used to redistribute light into a
    warp. The warp is parametrized as a change in Z.

    '''
    
    warp_factor = warppars['warp_factor']
    warp_ang = warppars['warp_ang']
    
    Zprime = Z - (warp_factor * np.cos(angles - warp_ang) * distances**3)
    
    return np.exp(-1*(Zprime)/h_z)

def fit_curve(datafile,central_lambda=[4901.416,5048.126],flip=False,ax=False,label='',\
                  rot_label='rotation_curve',pars=np.array([0,230,5.5,0.8,6.,0.36,np.pi/2.,0.652]),fixed=[],p=False):

    kpcradii, avg_centers, std_centers = openslay(datafile,central_lambda=central_lambda,flip=flip)

    x0 = [pars[i] for i in range(len(pars)) if i not in fixed]

    xf = spo.fmin(func,x0,args=(kpcradii,avg_centers,std_centers,fixed,pars),disp=False)
    
    pid = [i for i in range(len(pars)) if i not in fixed]
    k = 0
    for j in pid:
        pars[j] = xf[k]
        k += 1
    
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else: fig = False
    
    radii_width = kpcradii.max() - kpcradii.min()

#    plot_curve(datafile,ax=ax,central_lambda=central_lambda,flip=flip)
    model_r, model_v, _ = simcurve(1001,pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],ax=ax,
                                   kappa_0=pars[7],scale=radii_width/1001.,label=label,rot_label=rot_label,p=p)
    
    if p:
        ax.text(12,-30,('{}\n'*8).format(*['$*$' if i in fixed else '' for i in range(len(pars))]),color='r',
                horizontalalignment='left',va='top',fontsize=12)
        ax.text(15,-400,'$*$ = fixed',color='r')
    

    if fig: fig.show()

    return (pars, xf)

def profile_curve(fitsfile,in_radii,Iwidth=17,fig=False,sub=False,title=''):
    '''makes tangential velocity and PV curves based of a fits file that
    contains the output of simcurve. The user specifies what columns in in the
    fits file to use as radii coordinates. The TVC and "observed" velocity are
    computed via line profile fitting provided by line_profile'''

    hdus = pyfits.open(fitsfile)
    dists = hdus['DIST'].data
    V_S = hdus['V_S'].data

    radii = np.array([])
    fitvelos = np.array([])
    fitsig = np.array([])
    velos = np.array([])

    if not fig: fig = plt.figure()
    if not sub: ax = fig.add_subplot(111)
    else: ax = fig.add_subplot(sub)
    ax.set_title(title)
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Signal [arbitrary]')

    for r in in_radii:
        
        radii = np.append(radii,r)
        column = radius_to_column(dists,r)
        v_c = V_S[:,column].max()
        velos = np.append(velos,v_c)
        
        v,l, pars = line_profile(fitsfile,r,Iwidth=Iwidth,plot=False)
        gauss = gaussfunc(v,*pars)
        fitvelos = np.append(fitvelos,pars[1])
        fitsig = np.append(fitsig,pars[2])

        l = ax.plot(v,l,label='{:4.3f}'.format(r))[0]
        color = l.get_color()
        ax.plot(v,gauss,':',color=color)
        ax.axvline(x=v_c,color=color,linestyle='--')
        print r,pars

    ax.legend(loc=0,title='r [kpc]')

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(radii,velos,label='True $V_c$')
    ax1.errorbar(radii,fitvelos,yerr=fitsig,label='Fit velocity')
    ax1.set_xlabel('Radius [kpc]')
    ax1.set_ylabel('Velocity [km/s]')
    ax1.legend(loc=0)

#    fig.show()
#    fig1.show()

    return radii, velos, fitvelos, fitsig

def line_profile(fitsfile,radius,Iwidth=17.,
                 width=1.,plot=True,fit=True,
                 observe=True):
    """ Radius is in kpc"""


    hdus = pyfits.open(fitsfile)
    vs = hdus['V_S'].data
    frac = hdus[0].data
    dist = hdus['DIST'].data

    col1 = radius_to_column(dist,radius - width/2)
    col2 = radius_to_column(dist,radius + width/2)

    '''in this case, the desired width was less than the pixel resolution of
    the simulation'''
    if col1 == col2:
        col2 += 1
        print 'Desired width is less than simulation resolution. Using a width of {} kpc instead'.format(hdus[0].header['SCALE'])

    numsamp = 1000

    print 'building histogram'
    vhist, bins = np.histogram(
        np.mean(vs[:,col1:col2],axis=1),
        bins=20,
        weights=np.mean(frac[:,col1:col2],axis=1),
        density=True)
    bincent = 0.5*(bins[1:]+bins[:-1])
    v = np.linspace(bincent.min()-200.,bincent.max()+200,numsamp)
    ihist = np.interp(v,bincent,vhist,right=0.0,left=0.0)
    gauss_arr = np.empty(numsamp)

    if plot:
        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        ax0.set_xlabel('Velocity [km/s]')
        ax0.set_ylabel('Flux')
        ax0.plot(v,ihist,'--')
        fig0.show()
        raw_input('asdas')

    scale = np.mean(np.diff(v))
    print scale
    Iwidthpx = Iwidth/scale
    _, kernel = ADE.ADE_gauss(numsamp,numsamp/2-1,0,FWHM=Iwidthpx,NORM=True)
    print kernel.sum()
#    ADE.eplot(np.arange(kernel.size)*scale,kernel)
    lineshape = np.convolve(kernel,ihist,'same')

#     print 'building {} gaussians'.format(ihist.size)
#     for i in range(ihist.size):
# #        print ihist[i]
#         print i
#         r, gauss = ADE.ADE_gauss(numsamp,i,0.,PEAK_VAL=ihist[i]+0.0000001,FWHM=Iwidth/scale,NORM=False)
#         gauss_arr = np.vstack((gauss_arr,gauss))
#         if i % 25 == 0 and plot: 
#             ax0.plot(v,gauss)

    if plot: 
        fig0.show()
    # gauss_arr = gauss_arr[1:]

    # gauss_arr /= np.sum(gauss_arr)

    # lineshape2 = np.sum(gauss_arr,axis=0)

    if observe:
        print 'simulating effects of RSS'
        v, lineshape = observify(v,lineshape)
#        ov, lineshape2 = observify(v,lineshape2)
    # else:
    #     ov = v
        
    if fit:
        print 'fitting'
        fitpars = spo.curve_fit(gaussfunc,v,lineshape,
                                p0=(lineshape.max(),
                                    v[np.where(lineshape == \
                                                   lineshape.max())[0][0]],
                                    30.))[0]
        
        mgauss = gaussfunc(v, *fitpars)
    
    else:
        fitpars = False

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    #    ax.plot(bincent,vhist,'.')
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Flux')
        ax.plot(v,lineshape/lineshape.sum())
#        ax.plot(ov,lineshape2/lineshape2.sum(),':')
        if fit:
            ax.plot(v,mgauss)

        fig.show()

    hdus.close()
    del gauss_arr
    return v, lineshape, fitpars

def gaussfunc(x,peak,center,width): return peak*np.exp(-1*(x - center)**2/(2*width**2))

def radius_to_column(dist,radius):
    '''a small helper function to convert a radius (in kpc) into a column
    index that can be used to access various arrays produced by simcurve'''
    
    cidx = np.where(dist[int(dist.shape[0]/2),:] >= np.abs(radius))[0]
    if radius > 0: 
        pcidx = np.where(cidx >= dist.shape[1]/2)[0]
    else: 
        pcidx = np.where(cidx < dist.shape[1]/2)[0]
    scidx = np.argsort(dist[int(dist.shape[0]/2),cidx[pcidx]])

    try: column = cidx[pcidx[scidx]][0]
    except IndexError:
        column = dist.shape[1] - 1
        print "Invalid radius, using max radius of {} instead".format(dist[int(dist.shape[0]/2),column])
#    print "Using column {} where distance is {}".format(column,dist[int(dist.shape[0]/2),column])

    return column

def observify(velocity, flux, resolution=54., binsize=12.06):
    '''Designed to be a helper function to line_profile. It takes a line
    profile and simulates the effects of RSS on the data. It does this by
    first broadening the profile by the instrumental resolution (resolution)
    and then resampling it with pixels the same size as RSS pixles
    (binsize). All resolutions and binsizes are in km/s
    '''

    '''Assume resolution is the FWHM'''
    sigma = resolution/2.35482
    length = 2*4*sigma # 4 sigma should be enough to get the wings
    x = np.arange(length) - length/2.
    _, kernel = ADE.ADE_gauss(x, 0, sigma, NORM=True)
    smeared_flux = np.convolve(flux,kernel,'same')

    resampled_velo = np.arange(velocity.min(), velocity.max(), binsize)
    resampled_flux = np.interp(resampled_velo, velocity, smeared_flux)

    return resampled_velo, resampled_flux

def fit_rot(datafile, pars=np.array([230,2.5])):
    
    r, v = np.loadtxt(datafile,unpack=True)
    r *= 34.1e3/206265
    v -= 2.455e3

    xf = spo.fmin(rotfunc,pars,args=(r,v))

    # radii_width = r.max() - r.min()
    # rm, _, vm = simcurve(1000,0,xf[0],xf[1],scale = radii_width/1000.)

    rm = np.linspace(r.min(),r.max(),1000)
    vm = xf[0]*np.tanh(rm/xf[1])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(r,v,label='HI data')
    ax.plot(rm,vm,label='Model')
    ax.set_xlabel('distance from center [kpc]')
    ax.set_ylabel('Velocity [km/s]')
    fig.show()

    return xf

def rotfunc(x,r,v):

    rm = np.linspace(r.min(),r.max(),1000)
    vm = x[0]*np.tanh(rm/x[1])

    interp_v = np.interp(r, rm, vm)

    chisq = np.sum((v - interp_v)**2)/(r.size - x.size - 1)

    print chisq
    return chisq
    
def bootstrap(datafile, numtry):
    
    kpcradii, avg_centers, std_centers = openslay(datafile)

    fig = plt.figure()

    xlist = []

    x0 = np.array([3.86,230,0.652])

    for i in range(numtry):
#        fig.clf()
        ax = fig.add_subplot(111)
        sampleidx = np.random.randint(kpcradii.size, size = kpcradii.size)
        
        print sampleidx

        tmpradii = kpcradii[sampleidx]
        tmpcents = avg_centers[sampleidx]
        tmperr = std_centers[sampleidx]

        ax.set_xlabel('Radius [kpc]')
        ax.set_ylabel('LOS velocity [km/s]')
        ax.errorbar(tmpradii,tmpcents,yerr=tmperr,fmt='.') 
        fig.show()

        xf = spo.fmin(func,x0,args=(kpcradii,avg_centers,std_centers))
        
        radii_width = tmpradii.max() - tmpradii.min()
        model_r, model_v, _ = simcurve(1000,xf[0],xf[1],xf[2],ax=ax,kappa_0=xf[3],scale=radii_width/1000.)
        fig.show()
        
        xlist.append(xf)

        x0 = xf + xf*np.random.randn(xf.size)*0.3
    
    return np.vstack(xlist)

def func(x,kpcradii,avg_centers,std_centers,fixed,par0):
    
#    if len(x) == 0: return 1.0000
    xl = list(x)
    pars = [par0[j] if j in fixed else xl.pop(0) for j in range(len(par0))]
    
    if pars[3] < 0: return 9999999.
    
    '''N needs to be an integer'''
    # if pars[4] > 6.0: pars[4] = 7
    # if pars[4] < 6.0: pars[4] = 5
    
    radii_width = kpcradii.max() - kpcradii.min()

    model_r, model_v, _ = simcurve(501,pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],kappa_0=pars[7],scale=radii_width/501.)
    
    interp_v = np.interp(kpcradii, model_r, model_v)

    chisq = np.sum((avg_centers - interp_v)**2/std_centers**2)/(kpcradii.size - x.size - 1)

    # if np.isnan(chisq):
    #     print interp_v, pars
    #     raw_input('sda')

    print chisq
    return chisq

def bin_data(distvec,valuevecs,bin_amount):

    if np.log10(bin_amount)/np.log10(2) % 1 != 0:
        print "Can't bin by {}. ".format(bin_amount),
        bin_amount = int(2.**round(np.log10(bin_amount)/np.log10(2)))
        print "Binning by {} instead.".format(bin_amount)

    if type(valuevecs) != list or type(valuevecs) != tuple: valuevecs = [valuevecs]

    binned_dist = (distvec[bin_amount-1::bin_amount] + distvec[::bin_amount])/2.
    
    binned_vals = []
    for valvec in valuevecs:
        binned_value = np.array(valvec)
        for i in range(int(np.log10(bin_amount)/np.log10(2))):
            binned_value = binned_value[:-1:2] + binned_value[1::2]

        binned_vals.append(binned_value)

    return binned_dist, binned_vals

def bin_sim(sim_file,bin_amount,output):
    '''Takes the multi-HDU FITS produced by sim_curve and bins all the data
    by whatever amount the user specifies'''

    hdus = pyfits.open(sim_file)
    
    if np.log10(bin_amount)/np.log10(2) % 1 != 0:
        print "Can't bin by {}. ".format(bin_amount),
        bin_amount = int(2.**round(np.log10(bin_amount)/np.log10(2)))
        print "Binning by {} instead.".format(bin_amount)
    
    # '''first let's take care of the distance HDU. We need to treat this
    # differently b/c we average the distances, instead of summing them'''
    # dist_data = hdus['DIST'].data
    # dist_header = hdus['DIST'].header
    # if dist_data.shape[1] % 2 == 1: dist_data = dist_data[:,:-1]
    # binned_dist = (dist_data[:,bin_amount-1::bin_amount] + dist_data[:,::bin_amount])/2.

    # hdulist = [pyfits.ImageHDU(binned_dist)]
    # hdulist[0].header = dist_header

    hdulist = []

    for hdu in hdus:
        # if hdu.header['EXTNAME'] == 'DIST': continue
        
        print hdu.header['EXTNAME']

        binned_data = hdu.data
        if binned_data.shape[1] % 2 == 1: binned_data = binned_data[:,:-1]

        for i in range(int(np.log10(bin_amount)/np.log10(2))):
            '''so now (12.17.12) we're averaging all of the HDUs. We can't do
            them in the same way we did DIST b/c that assumed that the values
            were evenly spaced and increasing monatomically, which the other
            HDUs might not be'''
            binned_data = (binned_data[:,:-1:2] + binned_data[:,1::2])/2.

        if hdu.header['EXTNAME'] != 'FRAC':
            hdulist.append(pyfits.ImageHDU(binned_data))
            hdulist[-1].header = hdu.header
        else:
            hdulist.insert(0,pyfits.PrimaryHDU(binned_data))
            hdulist[0].header = hdu.header


    pyfits.HDUList(hdulist).writeto(output,clobber=True)

    return


def profile_Vz():
    '''just a little script to do line profile stuff at a few different
    heights and compare to Vc vs. z data'''

    cols = range(0,35,5)
    r0,v0,pv0,e0 = profile_curve('simZ0.fits',cols)
    r96,v96,pv96,e96 = profile_curve('simZ96.fits',cols)
    r19,v19,pv19,e19 = profile_curve('simZ19.fits',cols)
    r38,v38,pv38,e38 = profile_curve('simZ38.fits',cols)

    vc0 = find_Vc(r0,pv0,e0)
    vc96 = find_Vc(r96,pv96,e96)
    vc19 = find_Vc(r19,pv19,e19)
    vc38 = find_Vc(r38,pv38,e38)

    z = np.array([0,0.965,1.93,3.86])*0.43 #kpc
    Vc = np.array([vc0,vc96,vc19,vc38])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Height [kpc]')
    ax.set_ylabel('Tangential Velocity [km/s]')
    ax.set_xlim(-0.5,2)
    ax.plot(z,Vc,marker='s',linestyle='')

    fig.show()

    return z, Vc

def dust_or_nodust():
    
    radii = range(0,50,2)
    
    rnd0,vnd0,pvnd0,end0 = profile_curve('simZ0_changerc_nodust.fits',radii)
    rnd96,vnd96,pvnd96,end96 = profile_curve('simZ96_changerc_nodust.fits',radii)
    rnd19,vnd19,pvnd19,end19 = profile_curve('simZ19_changerc_nodust.fits',radii)
    rnd38,vnd38,pvnd38,end38 = profile_curve('simZ38_changerc_nodust.fits',radii)

    rd0,vd0,pvd0,ed0 = profile_curve('simZ0_changerc.fits',radii)
    rd96,vd96,pvd96,ed96 = profile_curve('simZ96_changerc.fits',radii)
    rd19,vd19,pvd19,ed19 = profile_curve('simZ19_changerc.fits',radii)
    rd38,vd38,pvd38,ed38 = profile_curve('simZ38_changerc.fits',radii)

    figlist = []
    
    for (rnd, pvnd, rd, pvd, vnd, z) in zip((rnd0,rnd96,rnd19,rnd38),
                                            (pvnd0,pvnd96,pvnd19,pvnd38),
                                            (rd0,rd96,rd19,rd38),
                                            (pvd0,pvd96,pvd19,pvd38),
                                            (vnd0,vnd96,vnd19,vnd38),
                                            (0,0.96,1.93,3.86)):
        fig = plt.figure()
        figlist.append(fig)
        ax = fig.add_subplot(111)
        ax.plot(rnd,vnd,label='TVC')
        ax.plot(rd,pvd,label='With dust')
        ax.plot(rnd,pvnd,label='No dust')
        ax.set_xlabel('Radius [kpc]')
        ax.set_ylabel('Radius [km/s]')
        ax.legend(loc=0)
        ax.set_title('$z/h_z = ${:4.2f}'.format(z))
        ax.set_ylim(0,250)
        ax.set_xlim(0,50)

    for fig in figlist: fig.show()

    return

def make_sims0():

    SALTr, SALTv, SALTerr = openslay('tiESO_z0_MgI.slay.fits') # should be .slay.gg.fits ??

    V_c = find_Vc(SALTr, SALTv, SALTerr)

    modelr = np.linspace(0,np.abs(SALTr).max(),100)
    modelv = V_c*np.tanh(modelr/5.45)

    for ff in [0,0.96,1.93,3.86]:
        if ff == 0.96:
            out = 'simZ96_nodust.fits'
            out2 = 'simZ96.fits'
        elif ff == 3.86:
            out = 'simZ38_nodust.fits'
            out2 = 'simZ38.fits'
        else:
            out = 'simZ{:1.0f}_nodust.fits'.format(ff*10)
            out2 = 'simZ{:1.0f}.fits'.format(ff*10)
        simcurve(1001,ff*0.43,0.0,5.45,rot_curve=(modelr,modelv),
                 scale=0.0999,kappa_0=0.0,output=out)
        simcurve(1001,ff*0.43,0.0,5.45,rot_curve=(modelr,modelv),
                 scale=0.0999,output=out2)

    return


def make_sims(galaxy):

    for ff in [0,0.96,1.93,3.86]:
        if ff == 0.96:
            out = 'simZ96_changerc_nodust.fits'
            out2 = 'simZ96_changerc.fits'
        elif ff == 3.86:
            out = 'simZ38_changerc_nodust.fits'
            out2 = 'simZ38_changerc.fits'
        else:
            out = 'simZ{:1.0f}_changerc_nodust.fits'.format(ff*10)
            out2 = 'simZ{:1.0f}_changerc.fits'.format(ff*10)
        simcurve(1001,ff*0.43,0.0,5.45,rot_curve=galaxy.get_TVC(ff*0.43),
                 scale=0.0999,kappa_0=0.0,output=out)
        simcurve(1001,ff*0.43,0.0,5.45,rot_curve=galaxy.get_TVC(ff*0.43),
                 scale=0.0999,output=out2)

    return
