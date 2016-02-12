#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import sys
from scipy.misc import imread

plt.ioff()

havel = np.array([651.52, 556.33, 350.69, 428.43, 674.46, 331.71])

def plot_PV(prefix='',
            EA=True,
            output='Swat_comp.png', col=1, chivel=False,
            colorcol=None, colorlabel='Height [kpc]',
            offset=73.7):

    swatfig = '/d/monk/eigenbrot/WIYN/14B-0456/anal/Swat.png'

    rr = np.array([])
    cc = np.array([])
    vv = np.array([])
    avV = np.array([])
    Har = np.array([])

    for i in range(6):

        dat = glob('{}*P{}*xcor.dat'.format(prefix,i+1))[0]
        print dat
        loc = glob('{}*P{}*locations.dat'.format(prefix,i+1))[0]
        print loc
        
        v = np.loadtxt(dat,usecols=(col,),unpack=True)

        if chivel:
            vel = glob('{}*P{}*vel.dat'.format(prefix,i+1))[0]
            print vel
            cv = np.loadtxt(vel,usecols=(1,),unpack=True)
            v += cv

        v += offset

        try:
            loc = loc.split('/')[1]
        except IndexError:
            pass
        r, c = np.loadtxt(loc,usecols=(2,5),unpack=True)
        
        if colorcol is not None:
            c = np.loadtxt(dat,usecols=(colorcol,),unpack=True)

        rr = np.r_[rr,r]
        cc = np.r_[cc,c]
        vv = np.r_[vv,v]
        avV = np.r_[avV,np.mean(v)]
        Har = np.r_[Har,np.mean(r)]

    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(rr,zz,vv)
    ax = fig.add_subplot(111)
    ax.set_xlabel('Distance along axis [arcmin]')
    ax.set_ylabel('Velocity [km/s]')
    ax.set_xticks([-8, -4, 0, 4, 8, 12])
    ax.set_yticks([400,600,800])
    im = imread(swatfig)
    ax.axvline(0,ls=':',alpha=0.7)
    ax.axhline(528,ls=':',alpha=0.7)
    s = ax.scatter(rr/60.,vv,c=cc,cmap='CMRmap',edgecolors='none',s=15,zorder=1)
    ax.scatter(Har/60.,havel, c='g', s=25)
    ax.scatter(Har/60.,avV, c='r', s=25)

    idx = np.where(Har < 0)
    ax.scatter(Har[idx]/60*-1,np.abs(havel[idx] - 528) + 528, c='g', alpha=0.7, s=20, marker='s')
    ax.scatter(Har[idx]/60*-1,np.abs(avV[idx] - 528) + 528, c='r', alpha=0.7, marker='s',s=20)

    print havel - avV
    print np.mean(havel - avV), np.std(havel - avV)
    
    ymin = 215
    ymax = 860
    xmin = -8.33
    xmax = 12.47
    ax.imshow(im,zorder=0, extent=[xmin,xmax,ymin,ymax], aspect='auto')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    cb = fig.colorbar(s)
    cb.set_label(colorlabel)

    if EA:
        Ez, Er, Ev = np.loadtxt('/d/monk/eigenbrot/WIYN/14B-0456/anal/PV_Ha_ngc891.txt',unpack=True)
        Ez *= 9.8/10
        Er *= 206265/9.8e3/60.
        ax.scatter(Er,Ev,c=Ez,cmap='CMRmap',vmin=cc.min(),vmax=cc.max(),edgecolors='none',marker='*',s=20,alpha=0.6)

    kax = ax.twiny()
    kax.set_xlim(ax.get_xlim()*np.array([2.91,2.91]))
    kax.set_xlabel('Radius [kpc]')


    fig.savefig(output,dpi=200)
    plt.close('all')

    return

if __name__ == '__main__':

    plot_PV(output=sys.argv[1])
