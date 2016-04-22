import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.floating_axes import GridHelperCurveLinear, FloatingSubplot
plt.ioff()

def compute_rphi(location, velocity, Vsys=528., Vc=225., rflat=3.3):

    rho, z = np.loadtxt(location, usecols=(4,5), unpack=True) #kpc
    V = np.loadtxt(velocity, usecols=(1,), unpack=True) #km/s

    Vcvec = np.interp(np.abs(rho),[0,rflat,15],[0,Vc,Vc])

    V -= Vsys

    idx = np.where(np.abs(V/Vcvec) > 1)
    V[idx] = Vc*np.sign(V[idx])

    r = Vcvec/V * rho
    phi = np.arccos(-1*V/Vcvec) #-1 b/c we define phi=0 to be the approaching
                                #side

    return r, phi*180./np.pi

def plot_rphi(r,phi,data,rlim=(0,10)):

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='polar')
    # ax.set_theta_offset(np.pi)
    # scat = ax.scatter(phi, r, c=data, s=40, cmap=plt.cm.gnuplot2, 
    #                   edgecolors='none', alpha=0.6)

    rmax = rlim[1]
    ax = fractional_polar_axes(fig,rlim=rlim,step=(45,round(rmax/4)),
                               ticklabels=False, thlabel='', rlabel='')
    ax.text(0,rmax,r'$\phi = 0^{\circ}$',va='center', ha='right')
    ax.text(180,rmax,r'$180^{\circ}$',va='center', ha='left')
    ax.text(90,rmax*1.05,'r [kpc]',ha='center')
    for i in range(4):
        rloc = round(rmax/4)*(i+1)
        ax.text(90,rloc,'{:3.0f}'.format(rloc), ha='left', fontsize=10)
    
    scat = ax.scatter(phi, r, c=data, s=40, cmap=plt.cm.gnuplot2,
                      edgecolors='none', alpha=0.6)

    return ax, scat

def plot_galaxy(rlim=(0,12)):

    rr = np.array([])
    pphi = np.array([])
    zz = np.array([])
    plims = []

    for i in range(6):
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(i+1)
        vel = 'NGC_891_P{}_bin30_velocities.dat'.format(i+1)
        rho, z = np.loadtxt(loc,usecols=(4,5),unpack=True)
        r, phi = compute_rphi(loc,vel)

        rr = np.r_[rr,r]
        pphi = np.r_[pphi, phi]
        zz = np.r_[zz,z]

        plims.append([-1*rho.min(),-1*np.mean(rho),-1*rho.max()])

    ax, scat = plot_rphi(rr,pphi,zz,rlim=rlim)
    
    ax.set_ylim(0,20)
    rmax = rlim[1]

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
                'P{}'.format(i+1), ha='center',va='bottom')

    cbax = ax.figure.add_axes([0.15,0.87,0.7,0.05])
    cb = ax.figure.colorbar(scat, cax=cbax, orientation='horizontal')
    cb.set_label('z [kpc]')
        
    # ax.set_xlabel('True radius [kpc]')
    # ax.set_ylabel(r'$\phi$')
    # ax.set_rgrids([5,10,15,20],angle=20,fontsize=9)
    # ax.text(np.pi*20./180., rmax*1.05, 'radius [kpc]', ha='right', fontsize=10)
    # ax.text(0,rmax*1.2, r'$\phi$ = ', va='center', ha='right')

    ax.figure.show()
    
    return ax

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
