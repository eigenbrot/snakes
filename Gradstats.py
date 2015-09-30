import numpy as np
import pyfits
import matplotlib.pyplot as plt
plt.ioff()
from glob import glob
import re
import os
import sys
import time

wavemin = 4500.
wavemax = 5500.

def go(pointing='*',prefix='NGC_891'):

    skylist = glob('{}_P{}*otdcs*ms_lin.fits'.format(prefix,pointing))
    skylist.sort()
    sslist = [s.replace('ms_lin','ms_s_lin') for s in skylist]
    print skylist
    print sslist

    p = re.compile('\d_(\d{3})_')
    fnum = np.array([int(p.search(i).group(1)) for i in skylist])
    
    mjd = np.array([float(pyfits.open(s)[0].header['JD']) for s in skylist])
    mjd -= np.min(mjd)
    mjd *= 24.

    sky = get_sky(skylist)
    source_sky = get_source(sslist)
    rms = get_rms(skylist)

    fig = plt.figure()
    skax = fig.add_subplot(311)
    skax.set_ylabel('$\\frac{B_S}{<B_S>_{\mathrm{night}}}$',fontsize=14)
    
    soax = fig.add_subplot(312)
    soax.set_ylabel('$\\frac{F - S}{<F - S>_{\mathrm{night}}}$',fontsize=14)

    rmax = fig.add_subplot(313)
    rmax.set_xlabel('Hours since first exposure')
    rmax.set_ylabel('rms($S - <S>_{\mathrm{size}}$)',fontsize=10)

    skax.errorbar(mjd, sky[0], yerr=sky[1],ls='',marker='o',color='k')
    soax.errorbar(mjd, source_sky[0], yerr=source_sky[1],ls='',marker='o',color='k')
    rmax.errorbar(mjd, rms[0], yerr=rms[1],ls='',marker='o',color='k')

    for i in range(mjd.size):
        skax.text(mjd[i]+0.08, sky[0,i],'{:4.2f}'.format(sky[0,i]),
                  va='center',fontsize=7)
        soax.text(mjd[i]+0.08, source_sky[0,i],'{:4.2f}'.format(source_sky[0,i]),
                  va='center',fontsize=7)
        rmax.text(mjd[i]+0.08, rms[0,i],'{:4.2f}'.format(rms[0,i]),
                  va='center',fontsize=7)

    skax.set_xticklabels([])
    soax.set_xticklabels([])
    
    utcax = skax.twiny()
    utcax.set_xticks(mjd)
    utcax.set_xticklabels(fnum)
    utcax.set_xlabel('Exposure number')

    skax.set_xlim(-0.5,np.max(mjd)+0.5)
    soax.set_xlim(skax.get_xlim())
    rmax.set_xlim(skax.get_xlim())
    utcax.set_xlim(skax.get_xlim())

    skax.set_ylim(skax.get_ylim()*np.array([0.999,1.]))
    soax.set_ylim(soax.get_ylim()*np.array([0.999,1.1]))
    rmax.set_ylim(rmax.get_ylim()*np.array([1.,1.1]))

    fig.subplots_adjust(hspace=0.0001)

    skax.text(0.1,0.9,'Error bars are $\sigma$\nacross all fibers',
              transform=skax.transAxes, va='top',
              fontsize=9)
    soax.text(0.1,0.9,'Error bars are $\sigma$\nacross all fibers',
              transform=soax.transAxes, va='top',
              fontsize=9)
    rmax.text(0.1,0.9,'Error bars are $\sigma$\nof RMS across each\nfiber size',
              transform=rmax.transAxes, va='top',
              fontsize=9)

    skax.text(-0.1,1.3,
              '$B_S$ = Sky surface brightness',
              transform=skax.transAxes, va='top',
              fontsize=7)
    skax.text(-0.1,1.23,
              '$F$ = Galaxy flux',
              transform=skax.transAxes, va='top',
              fontsize=7)
    skax.text(-0.1,1.16,
              '$S$ = Sky flux',
              transform=skax.transAxes, va='top',
              fontsize=7)

    skax.text(0.85,1.3,'{}_P{}\n{}'.format(prefix,pointing,time.asctime()),
              transform=skax.transAxes, va='top',
              fontsize=7)

    return fig

def openfits(image):

    h = pyfits.open(image)[0]
    data = h.data
    wave = np.arange(data.shape[1])*h.header['CDELT1'] + h.header['CRVAL1']

    return wave, data

def get_sky(skylist):

    skidx = np.array([1,2,18,19,20,31,32,43,44,53,54,62,63,71,79,87,88,95,102,109]) - 1
    print 'Sky:'
    tmp = []
    for frame in skylist:
        print '\t', frame
        wave, data = openfits(frame)
        idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        d = data[skidx]
        d = d[:,idx]
        d[0:4] /= 2.
        d[4:8] /= 3.
        d[8:12] /= 4.
        d[12:16] /= 5.
        d[16:20] /= 6.
        med = np.median(d,axis=1)
        tmp.append(med)
        
    d = np.vstack(tmp)
    skym = np.mean(d,axis=0)
    sky = np.mean(d/skym,axis=1)
    skystd = np.std(d/skym,axis=1)

    return np.array([sky,skystd])
        
        
def get_source(inlist):

    tmp = []
    print 'Source:'
    for frame in inlist:
        print '\t', frame
        wave, data = openfits(frame)
        idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        d = np.median(data[:,idx],axis=1)
        tmp.append(d)

    d = np.vstack(tmp)

    m = np.mean(d,axis=0)
    source = np.mean(d/m,axis=1)
    std = np.std(d/m,axis=1)

    return np.array([source,std])

def get_rms(skylist):

    skidx = np.array([1,2,18,19,20,31,32,43,44,53,54,62,63,71,79,87,88,95,102,109]) - 1
    print 'RMS:'
    bigout = []
    for frame in skylist:
        print '\t', frame
        wave, data = openfits(frame)
        idx = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        d = data[skidx]
        d = d[:,idx]
        sizeout = []
        for i in range(5):
            tmp = d[i*4:(i+1)*4,:]
            med = np.median(tmp,axis=1)
            mean = np.mean(med)
            rms = np.sqrt(np.mean((med - mean)**2))
            sizeout.append(rms)
        
        bigout.append(np.array([np.mean(sizeout),np.std(sizeout)]))

    return np.array(bigout).T

def main():

    cwd = os.getcwd()
    try:
        p = sys.argv[1]
        split = cwd.split('/')
        plot_name = "{}-P{}-{}_stats.png".format(split[-2],p,split[-1])
    except IndexError:
        p = '*'
        plot_name = "{}-{}_stats.png".format(*cwd.split('/')[-2:])

    f = go(pointing=p)
    f.savefig(plot_name)

    return 0

if __name__ == '__main__':
    main()
