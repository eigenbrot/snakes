import time
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from axisartist import angle_helper
from axisartist.grid_finder import MaxNLocator
from axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
plt.ioff()

def compute_rphi(location, velocity, Vsys=528., Vc=225., rflat=3.3,
                 dvdz = 15.789, output=False, dV=23., chidV_file=None):

    rho, z = np.loadtxt(location, usecols=(4,5), unpack=True) #kpc
    V = np.loadtxt(velocity, usecols=(1,), unpack=True) #km/s
    if chidV_file is not None:
        chidV = pyfits.open(chidV_file)[1].data['VSYS_ERROR']
        print 'loaded chi vel err from ', chidV_file
    else:
        chidV = 0.

    Vcvec = np.interp(np.abs(rho),[0,rflat,15],[0,Vc,Vc])
    Vcvec -= dvdz*np.abs(z)

    V -= Vsys

    idx = np.where(np.abs(V/Vcvec) > 1)
    V[idx] = Vcvec[idx]*np.sign(V[idx])

    r = np.abs(Vcvec/V) * np.abs(rho)
    phi = np.arccos(-1*V/Vcvec) #-1 b/c we define phi=0 to be the approaching
                                #side

    dr = np.abs(Vcvec/V**2) * np.abs(rho) * np.sqrt(dV**2 + chidV**2)

    if output:
        with open(output,'w') as f:
            f.write("""# Generate on {}
# from:
#     {}
#     {}
#
""".format(time.asctime(),location,velocity))
            f.write('#{:>6}{:>10}{:>10}{:>10}\n#\n'.format('Apnum','r (kpc)', 'phi (deg)', 'dr (kpc)'))
            for i in range(r.size):
                f.write('{:7n}{:10.3f}{:10.3f}{:10.3f}\n'.format(i,r[i],phi[i]*180./np.pi,dr[i]))


    return r, phi*180./np.pi, dr

def compute_rphi_tau(location, coeffile, scale=30./1001.):

    rho, z = np.loadtxt(location, usecols=(4,5), unpack=True) #kpc
    coefs = pyfits.open(coeffile)[1].data
    tau = coefs['TAUV']
    
    R = 15.

    size = int(R*2/scale) #size of array, in pixels
    if size % 2 == 0:
        size += 1

    idvec = np.indices((size,size), dtype=np.float32)
    radii = scale*((size/2. - idvec[0,:])**2 + (size/2. - idvec[1,:])**2)**0.5

    hz = 0.29 #For dust,
    hr = 7.68 # from Xilouris '99
    kappa0 = 0.85/2./hz #from Xilouris '99

    kaparray = np.exp(-1*(radii/hr))
    kaparray *= kappa0 / np.max(kaparray)
    
    print kappa0, kaparray[size/2,size/2]
    
    r = np.zeros(rho.size)
    phi = np.zeros(rho.size)

    for i in range(rho.size):
        tkap = kaparray*np.exp(-1*np.abs(z[i])/hz)
        kapcum = np.cumsum(tkap,axis=0) #still in pixels
        tauarray = scale*kapcum #times the ds of each pixel = integral!
        col = int(size/2 + rho[i]/scale)
        tau_loc = np.interp(tau[i],
                            tauarray[:,col],
                            np.arange(size)*scale)
        # print tauarray[:,col], tau[i], tau_loc
        if tau_loc > R:
            phi[i] = 270 - np.arctan(rho[i]/(tau_loc - R))*180./np.pi #+ 180 because 0 is the approaching side
        else:
            phi[i] = 90 + np.arctan(rho[i]/(R - tau_loc))*180./np.pi #+ 90 because 0 is the approaching side
        r[i] = np.sqrt(rho[i]**2 + (R - tau_loc)**2)

    return r, phi

def plot_rphi(r,phi,data,rlim=(0,12),thlim=(0,180)):

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='polar')
    # ax.set_theta_offset(np.pi)
    # ax.set_ylim(*rlim)
    # scat = ax.scatter(phi*np.pi/180., r, c=data, s=40, cmap=plt.cm.gnuplot2, 
    #                   edgecolors='none', alpha=0.6)

    rmax = rlim[1]
    ax = fractional_polar_axes(fig,thlim=thlim,rlim=rlim,step=(45,round(rmax/5)),
                               ticklabels=False, thlabel='', rlabel='')
    for i in range(4):
        rloc = round(rmax/4)*(i+1)
        ax.text(90,rloc,'{:3.0f}'.format(rloc), ha='left', fontsize=10)
    
    scat = ax.scatter(phi, r, c=data, s=40, cmap=plt.cm.gnuplot,
                      edgecolors='none', alpha=0.6, vmin=0, vmax=2.6)

    return ax, scat

def plot_galaxy(rlim=(0,12),thlim=(0,180),tau=False,basedir='.',
                chidV_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/var_disp/final_disp/chisq_vel/4166',                
                componentfile=None, exclude=[[],[],[],[],[],[]]):

    rr = np.array([])
    pphi = np.array([])
    zz = np.array([])
    plims = []

    for i in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir, i+1)
        vel = '{}/NGC_891_P{}_bin30_velocities.dat'.format(basedir, i+1)
        chidV_file = '{}/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(chidV_dir, i+1)
        rho, z = np.loadtxt(loc,usecols=(4,5),unpack=True)
        z = np.abs(z)
        if tau:
            coef = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir, i+1)
            r, phi = compute_rphi_tau(loc,coef)
        else:
            r, phi, _ = compute_rphi(loc,vel,chidV_file=chidV_file)

        exarr = np.array(exclude[i]) - 1
        rr = np.r_[rr,np.delete(r,exarr)]
        pphi = np.r_[pphi, np.delete(phi,exarr)]
        zz = np.r_[zz,np.delete(z,exarr)]

        plims.append([-1*rho.min(),-1*np.mean(rho),-1*rho.max()])

    ax, scat = plot_rphi(rr,pphi,zz,rlim=rlim,thlim=thlim)
    
    ax.set_ylim(0,20)
    rmax = rlim[1]
    if componentfile:
        add_IRAF_component(componentfile, ax)

    for i in range(6):
        print plims[i]
        l = Line2D([90. - 90.*np.sign(plims[i][1]),
                    np.arccos(plims[i][1]/rmax)*180/np.pi],
                   [np.abs(plims[i][1]),rmax],alpha=0.3)
        lmin = Line2D([90. - 90.*np.sign(plims[i][0]),
                       np.arccos(plims[i][0]/rmax)*180/np.pi],
                      [np.abs(plims[i][0]),rmax],linestyle=':',alpha=0.3)
        lmax = Line2D([90. - 90.*np.sign(plims[i][2]),
                       np.arccos(plims[i][2]/rmax)*180/np.pi],
                      [np.abs(plims[i][2]),rmax],linestyle=':',alpha=0.3)
#        ax.add_line(l)
        ax.add_line(lmin)
        ax.add_line(lmax)        

        ax.text(90. - 90.*np.sign(plims[i][1]), np.abs(plims[i][1]),
                r'$\mathrm{{P{}}}$'.format(i+1), ha='center',va='bottom')

    cbax = ax.figure.add_axes([0.16,0.865,0.7,0.05])
    cb = ax.figure.colorbar(scat, cax=cbax, orientation='horizontal')
    cb.set_alpha(1)
    cb.draw_all()
    cb.set_ticks(np.arange(0,3,0.5))
    cbax.tick_params(axis='x', labelsize=12)
    cbax.text(0.5,1.3,r'$|z|\mathrm{\ [kpc]}$',va='bottom', ha='center',
              transform=cbax.transAxes, fontsize=20)

    ax.text(0,rmax,r'$\phi = 0^{\circ}$',va='center', ha='right',fontsize=20)
    ax.text(180,rmax*1.02,r'$180^{\circ}$',va='center', ha='left',fontsize=20)
    ax.text(90,rmax*1.09,r'$r\mathrm{\ [kpc]}$',ha='center',fontsize=20)
        
    # ax.set_xlabel('True radius [kpc]')
    # ax.set_ylabel(r'$\phi$')
    # ax.set_rgrids([5,10,15,20],angle=20,fontsize=9)
    # ax.text(np.pi*20./180., rmax*1.05, 'radius [kpc]', ha='right', fontsize=10)
    # ax.text(0,rmax*1.2, r'$\phi$ = ', va='center', ha='right')

    # ax.figure.show()
    
    return ax

def add_IRAF_component(componentfile, ax):
    
    r, phi = np.loadtxt(componentfile, usecols=(5,6), unpack=True)

    ax.scatter(phi, r, color='k', marker='s', s=40, facecolor='none', alpha=0.7)

    return

def create_rphi_files():
    
    for p in range(6):
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(p+1)
        vel = 'NGC_891_P{}_bin30_velocities.dat'.format(p+1)
        chidV_file='/d/monk/eigenbrot/WIYN/14B-0456/anal/var_disp/final_disp/chisq_vel/4166/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(p+1)
        out = 'NGC_891_P{}_bin30_rphi.dat'.format(p+1)
        compute_rphi(loc,vel,output=out,chidV_file=chidV_file)

    return

def fractional_polar_axes(f, thlim=(0, 180), rlim=(0, 1), step=(30, 0.2),
                          thlabel='theta', rlabel='r', ticklabels=True):
    ## From https://github.com/neuropy/neuropy/blob/master/neuropy/scripts/polar_demo.py
    """Return polar axes that adhere to desired theta (in deg) and r limits. steps for theta
    and r are really just hints for the locators. Using negative values for rlim causes
    problems for GridHelperCurveLinear for some reason"""
    th0, th1 = thlim # deg
    r0, r1 = rlim
    thstep, rstep = step

    # scale degrees to radians:
    tr_scale = Affine2D().scale(np.pi/180., 1.)
    tr_rot = Affine2D().rotate(np.pi)
    tr = tr_scale + PolarAxes.PolarTransform() + tr_rot
    theta_grid_locator = angle_helper.LocatorDMS((th1-th0) // thstep)
    r_grid_locator = MaxNLocator((r1-r0) // rstep)
    theta_tick_formatter = angle_helper.FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extremes=(th0, th1, r0, r1),
                                        grid_locator1=theta_grid_locator,
                                        grid_locator2=r_grid_locator,
                                        tick_formatter1=theta_tick_formatter,
                                        tick_formatter2=None)

    a = FloatingSubplot(f, 111, grid_helper=grid_helper)
    f.add_subplot(a)

    # adjust x axis (theta):
    a.axis["bottom"].set_visible(False)
    a.axis["top"].set_axis_direction("bottom") # tick direction
    a.axis["top"].toggle(ticklabels=ticklabels, label=bool(thlabel))
    a.axis["top"].major_ticklabels.set_axis_direction("top")
    a.axis["top"].label.set_axis_direction("top")

    # adjust y axis (r):
    a.axis["left"].set_axis_direction("bottom") # tick direction
    a.axis["right"].set_axis_direction("top") # tick direction
    a.axis["left"].toggle(ticklabels=ticklabels, label=bool(rlabel))

    # add labels:
    a.axis["top"].label.set_text(thlabel)
    a.axis["left"].label.set_text(rlabel)

    # create a parasite axes whose transData is theta, r:
    auxa = a.get_aux_axes(tr)
    # make aux_ax to have a clip path as in a?:
    auxa.patch = a.patch 
    # this has a side effect that the patch is drawn twice, and possibly over some other
    # artists. So, we decrease the zorder a bit to prevent this:
    a.patch.zorder = -2

    # add sector lines for both dimensions:
    thticks = grid_helper.grid_info['lon_info'][0]
    rticks = grid_helper.grid_info['lat_info'][0]
    for th in thticks[1:-1]: # all but the first and last
        auxa.plot([th, th], [r0, r1], '--', c='grey', zorder=-1)
    for ri, r in enumerate(rticks):
        # plot first r line as axes border in solid black only if it isn't at r=0
        if ri == 0 and r != 0:
            ls, lw, color = 'solid', 2, 'black'
        else:
            ls, lw, color = 'dashed', 1, 'grey'
        # From http://stackoverflow.com/a/19828753/2020363
        auxa.add_artist(plt.Circle([0, 0], radius=r, ls=ls, lw=lw, color=color, fill=False,
                        transform=auxa.transData._b, zorder=-1))
    return auxa
