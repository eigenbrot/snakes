import numpy as np
import scipy.interpolate as spi
from scipy.misc import imread
import matplotlib.pyplot as plt

def main(order):

    r, v = np.loadtxt('Yim_Vc.dat', usecols=(0,1), unpack=True)

    # fit = np.poly1d(np.polyfit(r,v,order))
    fit = spi.UnivariateSpline(r,v,k=order)
    
    bigr = np.linspace(0,20,100)
    
    ax = plt.figure().add_subplot(111)
    ax.plot(r,v,'.')
    ax.plot(bigr,fit(bigr),'-r')
    ax.set_ylim(0,500)
    
    ax.figure.show()

    return fit

def fit_YimOos(order, savepick=False):

    Yr, Yv = np.loadtxt('Yim_Vc.dat', usecols=(0,1), unpack=True)
    Or, Ov = np.loadtxt('Oos_Vc.dat', usecols=(0,2), unpack=True)

    Yidx = np.where(Yr < 5.)
    Oidx = np.where(Or >= 5.)

    r = np.r_[Yr[Yidx],Or[Oidx]]
    v = np.r_[Yv[Yidx],Ov[Oidx]]

    fit = spi.UnivariateSpline(r,v,k=order)

    bigr = np.linspace(0,20,100)
    
    ax = plt.figure().add_subplot(111)
    ax.plot(r,v,'.')
    ax.plot(bigr,fit(bigr),'-r')
    ax.set_ylim(0,500)
    
    ax.figure.show()

    if savepick:
        import pickle
        with open(savepick, 'wb') as f:
            pickle.dump(fit,f)
    
    return fit
    

def plot_Vc(fit, output, oosfig=False, swatfig=False, yimfig=False):

    r = np.linspace(0,20,50)
    Vc0 = fit(r)
    corr = np.max(np.vstack((np.zeros(r.size) + 15,
                             45. - 35.*r/15.)),axis=0)
    print corr    
    zlist = np.linspace(0,3.0,5)

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('r [kpc]')
    ax.set_ylabel('Vc [km/s]')

    for z in zlist:
        Vc = Vc0 - corr*z
        ax.plot(r,Vc,label='{:3.1f}'.format(z))

    add_oos(ax)
    add_yim(ax)

    if oosfig:
        add_oosfig(ax)

    if swatfig:
        add_swatfig(ax)

    if yimfig:
        add_yimfig(ax)
    
    ax.legend(loc=0, numpoints=1, scatterpoints=1, title='z [kpc]')

    ax.figure.savefig(output)
    
    return

def add_oos(ax):

    r, Vc0, Vc15 = np.loadtxt('Oos_Vc.dat', usecols=(0,2,5), unpack=True)

    idx = np.where(r < 13.)
    
    ax.plot(r,Vc0,'xb')
    ax.plot(r[idx],Vc15[idx],'xr')

    return

def add_yim(ax):

    r, Vc = np.loadtxt('Yim_Vc.dat', usecols=(0,1), unpack=True)

    ax.plot(r,Vc,'.k')

    return

def add_oosfig(ax):

    im = imread('Oos15.png')
    ax.set_ylim(0, 250)
    ax.set_xlim(0, 21.05)
    ax.imshow(im, zorder=0, aspect='auto', extent=[0,21.05,0,250])

    return

def add_swatfig(ax):

    im = imread('Swat5.png')
    ax.set_ylim(0,290)
    ax.set_xlim(0,34.74)
    ax.imshow(im, zorder=0, aspect='auto', extent=[0,34.74,0,290])

    return

def add_yimfig(ax):

    im = imread('Yim5.png')
    ax.set_ylim(0,300)
    ax.set_xlim(-0.485,14.544)
    ax.imshow(im, zorder=0, aspect='auto', extent=[-0.485,14.544,0,300])

    return
