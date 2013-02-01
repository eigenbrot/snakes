import pyfits
import numpy as np
#import numpy.ma as ma
import sys
import os
from scipy.interpolate import interp1d
#import pyspeckit as psk
import matplotlib.pyplot as plt
import scipy.optimize as spo
from datetime import datetime
from matplotlib.backends.backend_pdf import PdfPages as PDF
import time
from ADESALT import openslay, plot_curve

centlambda = [4901.416,5048.126]
tau = np.pi*2.

def simcurve(size,Z,v_r,h_rot,w,N,pitch,view_ang,ax=False,scale=1.,\
                 kappa_0=0.652,z_d=0.23,label='',rot_label=False,p=False):
    
    #first generate the radial distance array. Taken from ADEUtils
    vecvec = np.indices((size,size),dtype=np.float64)
    distances = scale*((size/2. - vecvec[0,:])**2 + (size/2. - vecvec[1,:])**2)**0.5

    #The angle array will be used to figure out the los velocity for each bin.
    # taken from ADEUtils
    angles = np.abs(np.arctan2(size/2. - vecvec[0],size/2. - vecvec[1]) + np.pi)

    #some values from Xilouris '99 for ESO 435-G25 (IC 2531)
#    kappa_0 = 0.652 #=\tau/2z_d units are kpc^-1
    h_d = 8.43 #kpc
#    z_d = 0.23 #kpc

    #N needs to be an integer
    if N > 7.0: N = 8
    elif N < 5.0: N = 4
    else: N = 6

    #This array holds the opacity of each bin
    kaparray = np.exp(-1*(distances/h_d + np.abs(Z)/z_d))*LSP(distances, angles, w*0.8,N,pitch,view_ang + np.pi/6)
    kaparray *= kappa_0/kaparray.max()

    #And this one has the surface brightness of each bin, assuming a doubly-exp disc
    # total normalization is irrelevant
    Iarray = np.exp(-1*(distances/h_d + Z/z_d))
    lumpy = LSP(distances, angles, w,N,pitch,view_ang)
    Iarray *= lumpy/Iarray.max()


    #for each sight line, we'll now figure out the relative light contributions from each
    # bin along the sight line   
    kapcum = np.cumsum(kaparray,axis=0)
    tauarray = scale*kapcum
       
    #This array will hold the relative contribution from each bin
    fracarray = Iarray*np.exp(-1*tauarray)

    #We'll start by computing the rotation curve suggested by MAB
    rot_curve = v_r*np.tanh(distances/h_rot)
    
    #Now let's compute the projected los velocity of each bin, assuming the entire disc
    # rotates with v_r
    v_sarray = rot_curve*np.cos(angles)

    #finally, the light-weighted contribution to total LOS velocity
    LOSfracarray = (v_sarray*fracarray)/np.sum(fracarray,axis=0)

    if ax:
        radii_vec = scale*(np.arange(size)-size/2)
        rot_curve_r = radii_vec[np.where(radii_vec > 0)]
        rot_curve_vec = v_r*np.tanh(rot_curve_r/h_rot)
        ax.plot(radii_vec,np.sum(LOSfracarray,axis=0),label=label)
        ax.plot(rot_curve_r,rot_curve_vec,':',label=rot_label)
        ax.set_xlim(-50,50)
        ax.set_ylim(-500,500)
        if p:
            s = '$Z/h_z$: {0:3.3f}\n$V_r$: {1:5.3f}\n$h_{{rot}}$: {2:4.3f}\n$w$: {3:4.3f}\n$N$: {4:4.3f}\n$p$: {5:4.3f}\n$\\theta_{{view}}$: {6:4.3f}\n$\kappa_0$: {7:4.3f}\n$z_d$: {8:4.3f}'\
                .format(Z/0.43, v_r, h_rot, w, N, pitch,view_ang,kappa_0,z_d)
            ax.text(15, -20, s,horizontalalignment='left',va='top',fontsize=12)
        

    frachdu = pyfits.PrimaryHDU(fracarray)
    tauhdu = pyfits.ImageHDU(tauarray)
    kaphdu = pyfits.ImageHDU(kaparray)
    Ihdu = pyfits.ImageHDU(Iarray)
    vshdu = pyfits.ImageHDU(v_sarray)
    LOShdu = pyfits.ImageHDU(LOSfracarray)
    dhdu = pyfits.ImageHDU(distances)
    ahdu = pyfits.ImageHDU(angles)
    rothdu = pyfits.ImageHDU(rot_curve)
    

    frachdu.header.update('EXTNAME','FRAC')
    tauhdu.header.update('EXTNAME','TAU')
    kaphdu.header.update('EXTNAME','KAP')
    Ihdu.header.update('EXTNAME','SB')
    vshdu.header.update('EXTNAME','V_S')
    LOShdu.header.update('EXTNAME','LOSFRAC')
    dhdu.header.update('EXTNAME','DIST')
    ahdu.header.update('EXTNAME','ANG')
    rothdu.header.update('EXTNAME','ROT')

#    pyfits.HDUList([frachdu,tauhdu,kaphdu,Ihdu,vshdu,LOShdu,dhdu,ahdu,rothdu]).writeto('test.fits',clobber=True)

    return scale*(np.arange(size)-size/2), np.sum(LOSfracarray,axis=0), rot_curve

def LSP(distances, angles, w, bigN, p, vtheta):
    
    prodlist = []
    
    for n in np.arange(2,bigN+1,2):
        sinarray = (n*w)/(n-1) * np.sin( np.log(distances)/np.tan(p) - angles + vtheta)**bigN
        prodlist.append(sinarray)


    return 1 - w + np.array(prodlist).prod(axis=0)
    
# def plot_curve(datafile,central_lambda=[4901.416,5048.126],flip=False,ax=False,label=None):
#     '''Takes a slayer output file and plots the rotation curve associated with
#     the lines that were fit. Has lots of neat options for plotting.
#     '''

#     kpcradii, avg_centers, std_centers = openslay(datafile,central_lambda=central_lambda,flip=flip)    

#     if not ax:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#     else: fig = False

#     # ax2 = ax.twiny()
#     # ax2.set_xlim(arcsecradii[0],arcsecradii[-1])
#     # ax2.set_xlabel('Arcsec from center of galaxy')
#     ax.set_xlabel('Radius [kpc]')
#     ax.set_ylabel('LOS velocity [km/s]')
#     ax.errorbar(kpcradii,avg_centers,yerr=std_centers,fmt='.',label=label)
    
#     ax.axvline(x=0,ls='--',color='k',alpha=0.3)
#     ax.axhline(y=0,ls='--',color='k',alpha=0.3)

#     ax.set_xlim(-50,50)
#     ax.set_ylim(-500,500)

#     ax.set_title(datafile+'\n'+datetime.now().isoformat(' '))
    
#     if fig: fig.show()

def simhelp(zs):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for z in zs:
        simcurve(400,z,200,ax=ax,scale=0.2)

    ax.legend(loc=0)
    fig.show()

def simgrid():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_curve('tiESO_z05_MgI.slay.fits',flip=True,ax=ax)
    for v in [0,tau/8.,tau/4.,3*tau/8.,tau/2.]:
        _,_,_ = simcurve(1001,0.965*0.43,222.57,5.45,0.8,6,0.378,v,ax=ax,kappa_0=0.8,scale=0.0799,label=str(v/tau),rot_label=False,p=True)
    
    ax.legend(loc=0,title='Viewing angle/$\\tau$')

    fig.show()

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

# def openslay(datafile,central_lambda=[4901.416,5048.126],flip=False):

#     hdus = pyfits.open(datafile)

#     pars = hdus[1].data
#     errs = hdus[2].data
#     pxradii = hdus[3].data

#     centers = pars[:,1::3]
#     amps = pars[:,::3]
#     centerr = errs[:,1::3]

#     velocenters = (centers-central_lambda)/central_lambda*3e5
#     veloerrors = centerr/central_lambda*3e5

#     avg_centers = np.sum(amps*velocenters,axis=1)/np.sum(amps,axis=1)
#     std_centers = np.std(velocenters,axis=1)

#     offset = helpoff(pxradii,avg_centers)
#     print "Offset is "+str(offset)
    
#     kpcradii = pxradii - offset
#     kpcradii *= 0.1267*8.
#     kpcradii *= 0.178
    
#     if flip: kpcradii *= -1.

#     badidx = np.where(std_centers > 500.)
#     kpcradii = np.delete(kpcradii,badidx)
#     avg_centers = np.delete(avg_centers,badidx)
#     std_centers = np.delete(std_centers,badidx)

#     return (kpcradii,avg_centers,std_centers)

# def helpoff(radii,centers):

#     x0 = np.array([np.median(radii)])
    
#     xf = spo.fmin(offunc,x0,args=(radii,centers),disp=False)
    
#     radii = radii.copy() - xf[0]
#     pidx = np.where(radii >= 0.0)
#     nidx = np.where(radii < 0.0)
    
#     # fig = plt.figure()
#     # ax = fig.add_subplot(111)
#     # ax.plot(radii[pidx],np.abs(centers[pidx]))
#     # ax.plot(np.abs(radii[nidx]),np.abs(centers[nidx]))
#     # fig.show()

#     return xf[0]

# def offunc(x,radii,centers):
    
#     radii = radii.copy() - x[0]
#     pidx = np.where(radii >= 0.0)
#     nidx = np.where(radii < 0.0)

#     if pidx[0].size <= 1 or nidx[0].size <= 1: return 999.
    
#     pcent = centers[pidx]
#     ncent = centers[nidx]

#     pcoef = np.polyfit(radii[pidx],np.abs(centers[pidx]),4)
#     ncoef = np.polyfit(np.abs(radii[nidx][::-1]),np.abs(centers[nidx][::-1]),4)

#     pf = np.poly1d(pcoef)
#     nf = np.poly1d(ncoef)

#     minr = max(radii[pidx].min(),np.abs(radii[nidx]).min())
#     maxr = min(radii[pidx].max(),np.abs(radii[nidx]).max())

#     r = np.linspace(minr,maxr,100)

#     chisq = np.sum((pf(r) - nf(r))**2)

#     # fig = plt.figure()
#     # ax = fig.add_subplot(111)
#     # ax.plot(r,pf(r))
#     # ax.plot(r,nf(r))
#     # fig.show()
#     # raw_input('ssds')

# #    print chisq
    
#     return chisq


def dosims(filename,fitfix):
    '''to be run in /d/monk/eigenbrot/SALT/2011-3-UW_002/rotation'''

    pp = PDF(filename)

    print "Getting UrPars..."
    fig09 = plt.figure()
    ax09 = fig09.add_subplot(111)
    plot_curve('tiESO_z05_MgI.slay.fits',flip=True,ax=ax09)
    urpars = fit_curve('tiESO_z05_MgI.slay.fits',flip=True,ax=ax09,pars=[0.965*0.43,242.992,5.45,1.0,5,0.5,tau/4.,1.358],fixed=fitfix,p=True,rot_label='Ur-rotation curve',label='Urfit')[0]
    ax09.legend(loc=0)
#    fig09.show()
    
    fig19 = plt.figure()
    ax19 = fig19.add_subplot(111)
    plot_curve('tiESO_z1_MgI.slay.fits',ax=ax19,flip=True)
    print "Running z=1.93 step 1...",
    fit_curve('tiESO_z1_MgI.slay.fits',ax=ax19,flip=True,pars=[1.93*0.43] + urpars[1:],fixed=[0,1,2,3,4,5,6,7],label='Urfit extrapolation',rot_label='Ur-rotation curve')
    print "2..."
    pars19 = fit_curve('tiESO_z1_MgI.slay.fits',ax=ax19,flip=True,pars=[1.93*0.43,238.512] + urpars[2:7] + [1.374],fixed=fitfix,label='Rotation curve fit to data',p=True,rot_label='Fit rotation curve')[0]
    ax19.legend(loc=0)
#    fig19.show()

    fig38 = plt.figure()
    ax38 = fig38.add_subplot(111)
    plot_curve('tiESO_z2_MgI.slay.fits',ax=ax38,flip=True)
    print "Running z=3.98 step 1...",
    fit_curve('tiESO_z2_MgI.slay.fits',ax=ax38,flip=True,pars=[3.86*0.43] + urpars[1:],fixed=[0,1,2,3,4,5,6,7],label='Urfit extrapolation',rot_label='Ur-rotation curve')
    print "2..."
    pars38 = fit_curve('tiESO_z2_MgI.slay.fits',ax=ax38,flip=True,pars=[3.86*0.43,286.819] + urpars[2:7] + [11.240],fixed=fitfix,label='Rotation curve fit to data',p=True,rot_label='Fit rotation curve')[0]
    ax38.legend(loc=0)
#    fig38.show()

    fig0 = plt.figure()
    ax0 = fig0.add_subplot(111)
    plot_curve('tiESO_z0_MgI.slay.fits',ax=ax0)
    print "Running z=3.86 step 1...",
    fit_curve('tiESO_z0_MgI.slay.fits',ax=ax0,pars=[0] + urpars[1:],fixed=[0,1,2,3,4,5,6,7],label='Urfit extrapolation',rot_label='Ur-rotation curve')
    print "2..."
    pars0 = fit_curve('tiESO_z0_MgI.slay.fits',ax=ax0,pars=[0,256.641] + urpars[2:7] + [1.731],fixed=fitfix,label='Rotation curve fit to data',p=True,rot_label='Fit rotation curve')[0]
    ax0.legend(loc=0)
#    fig0.show()
    
    pp.savefig(fig0)
    pp.savefig(fig09)
    pp.savefig(fig19)
    pp.savefig(fig38)
    pp.close()

    return (pars0,urpars,pars19,pars38)

def uberboot(numiter,pars=[0.000,256.292,5.45,
                   0.965*0.43,244.636,5.45,
                   1.930*0.43,241.189,5.45,
                   3.860*0.43,237.407,5.45,
                   0.0,5,0.5,tau/4,1.62,0.254],
             fixed=[0,2,3,5,6,8,9,11,12,13,14,15,16,17]):

    r0, v0, e0 = openslay('tiESO_z0_MgI.slay.gg.fits')
    r09, v09, e09 = openslay('tiESO_z05_MgI.slay.fits',flip=True)
    r19, v19, e19 = openslay('tiESO_z1_MgI.slay.fits',flip=True)
    r38, v38, e38 = openslay('tiESO_z2_MgI.slay.fits',flip=True)

    bigarray = np.empty(len(pars) - len(fixed))
    
    for n in range(numiter):
        
        idx0 = np.random.randint(r0.size, size = r0.size)
        idx09 = np.random.randint(r09.size, size = r09.size)
        idx19 = np.random.randint(r19.size, size = r19.size)
        idx38 = np.random.randint(r38.size, size = r38.size)

        rr0, rv0, re0 = r0[idx0], v0[idx0], e0[idx0]
        rr09, rv09, re09 = r0[idx09], v0[idx09], e0[idx09]
        rr19, rv19, re19 = r0[idx19], v0[idx19], e0[idx19]
        rr38, rv38, re38 = r0[idx38], v0[idx38], e0[idx38]
        
        data = dict((('r0',rr0),('v0',rv0),('e0',re0),
                 ('r09',rr09),('v09',rv09),('e09',re09),
                 ('r19',rr19),('v19',rv19),('e19',re19),
                 ('r38',rr38),('v38',rv38),('e38',re38)))

        x0 = [pars[i] for i in range(len(pars)) if i not in fixed]

        print x0

        print "fitting # {:}".format(n)
        t1 = time.time()
        xf = spo.fmin(uberfunc,x0,args=(data,fixed,pars),xtol=0.1)#,epsfcn=1)[0]
        t2 = time.time()
        print "   fitting took {} seconds".format(t2-t1)
        
        bigarray = np.vstack((bigarray,xf))

    bigarray = bigarray[1:]

    return bigarray


def uberfit(filename,pars=[0.000,256.292,5.45,
                           0.965*0.43,244.636,5.45,
                           1.930*0.43,241.189,5.45,
                           3.860*0.43,237.407,5.45,
                           0.0,5,0.5,tau/4,1.62,0.1225],fixed=[]):

    r0, v0, e0 = openslay('tiESO_z0_MgI.slay.gg.fits')
    r09, v09, e09 = openslay('tiESO_z05_MgI.slay.fits',flip=True)
    r19, v19, e19 = openslay('tiESO_z1_MgI.slay.fits',flip=True)
    r38, v38, e38 = openslay('tiESO_z2_MgI.slay.fits',flip=True)

    data = dict((('r0',r0),('v0',v0),('e0',e0),
                 ('r09',r09),('v09',v09),('e09',e09),
                 ('r19',r19),('v19',v19),('e19',e19),
                 ('r38',r38),('v38',v38),('e38',e38)))

    x0 = [pars[i] for i in range(len(pars)) if i not in fixed]

    print x0

    print "fitting"
    t1 = time.time()
    xf = spo.fmin(uberfunc,x0,args=(data,fixed,pars))#,epsfcn=1)[0]
    t2 = time.time()
    print "fitting took {} seconds".format(t2-t1)
    
    print "pre-plot"
    pid = [i for i in range(len(pars)) if i not in fixed]
    k = 0
    for j in pid:
        pars[j] = xf[k]
        k += 1

    pp=PDF(filename)

    pf0 = [0,1,2,12,13,14,15,16,17]
    pf09 = [3,4,5,12,13,14,15,16,17]
    pf19 = [6,7,8,12,13,14,15,16,17]
    pf38 = [9,10,11,12,13,14,15,16,17]

    print "plotting"
    ###########
    fig0 = plt.figure()
    ax0 = fig0.add_subplot(111)
    rw0 = r0.max() - r0.min()
    f0 = [i for i in fixed if i in pf0]
    plot_curve('tiESO_z0_MgI.slay.gg.fits',ax=ax0)
    mr0, mv0, _ = simcurve(1001,pars[0],pars[1],pars[2],pars[12],pars[13],pars[14],pars[15],kappa_0=pars[16],\
                           z_d=pars[17],scale=rw0/1001.,ax=ax0,p=True,label='Global fit',rot_label='Rotation curve')
    ax0.text(12,-30,('{}\n'*9).format(*['$*$' if i in f0 else '' for i in pf0]),color='r',
            horizontalalignment='left',va='top',fontsize=12)
    ax0.text(15,-450,'$*$ = fixed',color='r')
    ax0.legend(loc=0)
    iv0 = np.interp(r0, mr0, mv0)

    ###########
    # fig0_1 = plt.figure()
    # ax0_1 = fig0_1.add_subplot(111)
    # rw0_1 = r0_1.max() - r0_1.min()
    # f0_1 = [i for i in fixed if i in pf0_1]
    # plot_curve('tiESO_z0_MgI.slay.2_16.fits',ax=ax0_1)
    # _, _, _ = simcurve(1001,pars[3],pars[4],pars[5],pars[12],pars[13],pars[14],pars[15],kappa_0=pars[16],\
    #                        z_d=pars[20],scale=rw0_1/1001.,ax=ax0_1,p=True,label='Global fit',rot_label='Rotation curve')
    # ax0_1.text(12,-30,('{}\n'*9).format(*['$*$' if i in f0_1 else '' for i in pf0_1]),color='r',
    #         horizontalalignment='left',va='top',fontsize=12)
    # ax0_1.text(15,-450,'$*$ = fixed',color='r')
    # ax0_1.legend(loc=0)

    ###########
    fig09 = plt.figure()
    rw09 = r09.max() - r09.min()
    f09 = [i for i in fixed if i in pf09]
    ax09 = fig09.add_subplot(111)
    plot_curve('tiESO_z05_MgI.slay.fits',ax=ax09,flip=True)
    mr09, mv09, _ = simcurve(1001,pars[3],pars[4],pars[5],pars[12],pars[13],pars[14],pars[15],kappa_0=pars[16],\
                       scale=rw09/1001.,z_d=pars[17],ax=ax09,p=True,label='Global fit',rot_label='Rotation curve')
    ax09.text(12,-30,('{}\n'*9).format(*['$*$' if i in f09 else '' for i in pf09]),color='r', 
              horizontalalignment='left',va='top',fontsize=12)
    ax09.text(15,-450,'$*$ = fixed',color='r')
    ax09.legend(loc=0)
    iv09 = np.interp(r09, mr09, mv09)

    ###########
    fig19 = plt.figure()
    rw19 = r19.max() - r19.min()
    f19 = [i for i in fixed if i in pf19]
    ax19 = fig19.add_subplot(111)
    plot_curve('tiESO_z1_MgI.slay.fits',flip=True,ax=ax19)
    mr19, mv19, _ = simcurve(1001,pars[6],pars[7],pars[8],pars[12],pars[13],pars[14],pars[15],kappa_0=pars[16],\
                           scale=rw19/1001.,z_d=pars[17],ax=ax19,p=True,label='Global fit',rot_label='Rotation curve')
    ax19.text(12,-30,('{}\n'*9).format(*['$*$' if i in f19 else '' for i in pf19]),color='r', 
              horizontalalignment='left',va='top',fontsize=12)
    ax19.text(15,-450,'$*$ = fixed',color='r')
    ax19.legend(loc=0)
    iv19 = np.interp(r19, mr19, mv19)

    ############
    fig38 = plt.figure()
    rw38 = r38.max() - r38.min()
    f38 = [i for i in fixed if i in pf38]
    ax38 = fig38.add_subplot(111)
    plot_curve('tiESO_z2_MgI.slay.fits',flip=True,ax=ax38)
    mr38, mv38, _ = simcurve(1001,pars[9],pars[10],pars[11],pars[12],pars[13],pars[14],pars[15],kappa_0=pars[16],\
                           scale=rw38/1001.,z_d=pars[17],ax=ax38,p=True,label='Global fit',rot_label='Rotation curve')
    ax38.text(12,-30,('{}\n'*9).format(*['$*$' if i in f38 else '' for i in pf38]),color='r', 
              horizontalalignment='left',va='top',fontsize=12)
    ax38.text(15,-450,'$*$ = fixed',color='r')
    ax38.legend(loc=0)
    iv38 = np.interp(r38, mr38, mv38)
    
    #############
    
    bigr = np.concatenate((r0,r09,r19,r38))
    bigiv = np.concatenate((iv0,iv09,iv19,iv38))
    bigv = np.concatenate((v0,v09,v19,v38))
    bige = np.concatenate((e0,e09,e19,e38))
    
    chisq = np.sum((bigv - bigiv)**2/bige**2)/(bigr.size - xf.size - 1)

    ax0.text(-40,200,'Red. $\chi^2$: {:4.3f}'.format(chisq))
    ax09.text(-40,200,'Red. $\chi^2$: {:4.3f}'.format(chisq))
    ax19.text(-40,200,'Red. $\chi^2$: {:4.3f}'.format(chisq))
    ax38.text(-40,200,'Red. $\chi^2$: {:4.3f}'.format(chisq))

    pp.savefig(fig0)
    pp.savefig(fig09)
    pp.savefig(fig19)
    pp.savefig(fig38)
    pp.close()

    return pars

def uberfunc(x,data,fixed,par0):
    
    print x
    xl = list(x)
    pars = [par0[j] if j in fixed else xl.pop(0) for j in range(len(par0))]
    
#    if pars[2] != 5.45: print x,pars
#    if pars[17] < 0 or pars[14] > 4: return np.zeros(10) + 9e9

    ##z=0
    r0 = data['r0']
    v0 = data['v0']
    e0 = data['e0']
    rw0 = r0.max() - r0.min()
    mr0, mv0, _ = simcurve(501,pars[0],pars[1],pars[2],pars[12],pars[13],pars[14],pars[15],\
                               kappa_0=pars[16],z_d=pars[17],scale=rw0/501.)
    iv0 = np.interp(r0, mr0, mv0)
    
    # ##z=0, part 2
    # r0_1 = data['r0_1']
    # v0_1 = data['v0_1']
    # e0_1 = data['e0_1']
    # rw0_1 = r0_1.max() - r0_1.min()
    # mr0_1, mv0_1, _ = simcurve(501,pars[3],pars[4],pars[5],pars[12],pars[13],pars[14],pars[15],\
    #                            kappa_0=pars[19],z_d=pars[20],scale=rw0/501.)
    # iv0_1 = np.interp(r0_1, mr0_1, mv0_1)

    ##z=0.96
    r09 = data['r09']
    v09 = data['v09']
    e09 = data['e09']
    rw09 = r09.max() - r09.min()
    mr09, mv09, _ = simcurve(501,pars[3],pars[4],pars[5],pars[12],pars[13],pars[14],pars[15],\
                                 kappa_0=pars[16],z_d=pars[17],scale=rw09/501.)
    iv09 = np.interp(r09, mr09, mv09)

    ##z=1.98
    r19 = data['r19']
    v19 = data['v19']
    e19 = data['e19']
    rw19 = r19.max() - r19.min()
    mr19, mv19, _ = simcurve(501,pars[6],pars[7],pars[8],pars[12],pars[13],pars[14],pars[15],\
                                 kappa_0=pars[16],z_d=pars[17],scale=rw19/501.)
    iv19 = np.interp(r19, mr19, mv19)

    ##z=3.86
    r38 = data['r38']
    v38 = data['v38']
    e38 = data['e38']
    rw38 = r38.max() - r38.min()
    mr38, mv38, _ = simcurve(501,pars[9],pars[10],pars[11],pars[12],pars[13],pars[14],pars[15],\
                                 kappa_0=pars[16],z_d=pars[17],scale=rw38/501.)
    iv38 = np.interp(r38, mr38, mv38)

    bigr = np.concatenate((r0,r09,r19,r38))
    bigiv = np.concatenate((iv0,iv09,iv19,iv38))
    bigv = np.concatenate((v0,v09,v19,v38))
    bige = np.concatenate((e0,e09,e19,e38))

    chisq = np.sum((bigv - bigiv)**2/bige**2)/(bigr.size - x.size - 1)
    
    print chisq
#    return (bigv - bigiv)/bige
    return chisq


