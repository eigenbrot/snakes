#! /usr/bin/env python
#
# To be run in /usr/users/eigenbrot/research/edgeOn/sprial_profiles
#
import numpy as np
import Salty2 as salty
from glob import glob
import pyfits
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages as PDF

def galaxy_batch():
    
    for w in [0.6,0.8,0.9,1.0]:
        for pitch in np.pi/np.array([3.,4.,5.,6.]):
            for view in np.pi/np.array([2.,3.,4.,5.]):
                for N in [4,6,8]:
                    name = 'sim_{:.1f}_{:}_{:n}_{:n}.fits'.\
                        format(w,N,np.pi/pitch,np.pi/view)
                    print name
                    salty.simcurve(205,0.83,256,5.45,w=w,N=N,pitch=pitch,
                                   view_ang=view,scale=0.39,output=name)

    return

def profile_batch(radius,output):

    pp = PDF(output)

    for sim in glob('sim*.fits'):
        print sim

        header = pyfits.open(sim)[0].header
        w = header['W']
        N = header['N']
        pitch = header['PITCH']
        ang = header['VIEWANG']
        pitch = int(np.pi*2/pitch)
        ang = int(np.pi*2/ang)
        v, line, _ = salty.line_profile(sim,radius,pxbin=4.,plot=False)
        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Normalized power')
        ax.set_title(sim)
        ax.text(300,0.005,
                '$w={}$\n$N={}$\n$p=\\tau/{}$\n$\\theta_{{view}}=\\tau/{}$'.\
                    format(w,N,pitch,ang))
        ax.plot(v,line)
        pp.savefig(ax.figure)
    
    pp.close()
        
def main(radius, output):

    galaxy_batch()
    profile_batch(radius,output)

    return

if __name__ == '__main__':
    sys.exit(main(sys.argv[1],sys.argv[2]))
