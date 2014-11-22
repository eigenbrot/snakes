#!/usr/bin/python

import sys
import numpy as np
import pyfits
import matplotlib.pyplot as plt

def plotstar(filename,catch=False):

    hdu = pyfits.open(filename)[0]
    data = hdu.data
    header = hdu.header

    wave = np.arange(data.shape[1])*header['CDELT1'] + header['CRVAL1']
    
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Normalized flux')
    ax.set_title(filename,fontsize=10)
    prevsize = ''
    for i in range(data.shape[0]):
        flux = data[i,:]
        flux /= flux[994]
        label = '{} ({}00 micron)'.format(*header['APNUM{}'.format(i+1)].split(' ')[0:2])
        if label[-12:] == prevsize:
            prevline = ax.plot(wave,flux,label=label,lw=0.5,color=prevline[0].get_color())
        else:
            prevline = ax.plot(wave,flux,label=label,lw=0.5)
        prevsize = label[-12:]
            
    ax.legend(title='aperture',loc=0,numpoints=1,scatterpoints=1,frameon=False,fontsize=9)
    ax.set_ylim(0,6)
    if catch:
        plt.show()
        return
    else:
        ax.figure.show()
        return ax
        
if __name__ == '__main__':

    try:
        sys.exit(plotstar(sys.argv[1],True))
    except Exception as e:
        print "The request was made but it was not good"
        sys.exit(1)
