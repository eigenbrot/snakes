import numpy as np
import pyfits
import matplotlib.pyplot as plt
plt.ioff()
from glob import glob
import re
import os

wavemin = 4500.
wavemax = 5500.

def go(prefix='NGC_891'):

    skylist = glob('{}*otdcs*ms_lin.fits'.format(prefix))
    skylist.sort()
    sslist = [s.replace('ms_lin','ms_s_lin') for s in skylist]
    print skylist
    print sslist

    p = re.compile('\d_(\d{3})_')
    fnum = np.array([int(p.search(i).group(1)) for i in skylist])
    
    mjd = np.array([float(pyfits.open(s)[0].header['MJD-OBS']) for s in skylist])
    mjd -= np.min(mjd)
    mjd *= 24.

    sky = get_sky(skylist)
    source_sky = get_source(sslist)
    rms = get_rms(skylist)

    fig = plt.figure()
    skax = fig.add_subplot(311)
    skax.set_ylabel('Sky flux')
    
    soax = fig.add_subplot(312)
    soax.set_ylabel('f - s/<f - s>')

    rmax = fig.add_subplot(313)
    rmax.set_xlabel('Hours since first exposure')
    rmax.set_ylabel('rms(sky - <sky>)')

    skax.errorbar(mjd, sky[0], yerr=sky[1],ls='',marker='o')
    soax.errorbar(mjd, source_sky[0], yerr=source_sky[1],ls='',marker='o')
    rmax.errorbar(mjd, rms[0], yerr=rms[1],ls='',marker='o')

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

    skax.set_ylim(skax.get_ylim()*np.array([0.9,1.]))
    soax.set_ylim(soax.get_ylim()*np.array([0.9,1.1]))
    rmax.set_ylim(rmax.get_ylim()*np.array([1.,1.1]))

    fig.subplots_adjust(hspace=0.0001)

    return fig

def openfits(image):

    h = pyfits.open(image)[0]
    data = h.data
    wave = np.arange(data.shape[1])*h.header['CDELT1'] + h.header['CRVAL1']

    return wave, data

def get_sky(skylist):

    skidx = np.array([1,2,18,19,20,31,32,43,44,53,54,62,63,71,79,87,88,95,102,109]) - 1
    print 'Sky:'
    output = []
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
        sky = np.mean(med)
        skystd = np.std(med)
        output.append(np.array([sky,skystd]))

    return np.array(output).T
        
        
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
    print d.shape

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
    plot_name = "{}-{}_stats.png".format(*cwd.split('/')[-2:])
    f = go()
    f.savefig(plot_name)

    return 0

if __name__ == '__main__':
    main()
