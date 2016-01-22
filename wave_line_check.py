#!/usr/bin/python
import sys
import os
import time
import numpy as np
import ir2py_lines as i2p
import matplotlib.pyplot as plt
plt.ioff()
from pyraf import iraf

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

llist = ['4047','4358','HK','5461']
rlist = ['4020 4075','4325 4400','3910 4010','5437 5480']
centlist = [[4047],[4358],[3933.57,3968.53],[5461]]
numlist = [1,1,2,1]

def do_fitprof(datafile):
    """Run IRAF fitprof routine to measure line centers and widths.

    For each line specified in the module header generate the necessary inputs
    to fitprof and execture the routine. The fitting regions are also
    specified in the module header.

    Parameters
    ----------
    datafile : str
        Name of a multispec fits file to pass as input to fitprof. Must have a WCS solution in the header.

    Returns
    -------
    None
        Nothing is returned. Instead, IRAF writes the results to a file.

    """
    
    for l, cl, r in zip(llist,centlist,rlist):
        with open('{}.lines'.format(l),'w') as f:
            print l
            for c in cl:
                f.write('{} INDEF g\n'.format(c))

        iraf.fitprofs(datafile,
                      region=r,
                      positio='{}.lines'.format(l),
                      logfile='{}.fitp'.format(l))

    return

def get_results(output, threshold=3.):
    """Parse fitprof output and display results

    The line centers are taken from the output of fitprof. For each line specified in the module header the average offset and stddev across all fibers in the IFU is computed. Output is a textfile and a plot of accuracy and stochasticity as a function of wavelength.

    Parameters
    ----------
    output : str
        Name of the output text file. This file will contain the mean, offset, stddev, and number of rejected apertures.
    threshold : float, optional
        Threshold value for iterative sigma clipping in mean across IFU. The total number of rejected fibers will be recorded in the output file.

    Returns
    -------
    None :
       The result is a text file containing the results and a plot containing the accuracy as a function of wavelength.
    
    Notes
    -----
    Right now the plot is hardcoded to be writting to WLC.png

    """

    fig = plt.figure()
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('Mismatch [AA]')
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Wavelength')
    stax.set_ylabel('IFU std')

    with open(output,'a') as f:
        f.write('# {}\n'.format(time.asctime()))

        for l, n, c in zip(llist,numlist,centlist):
            proffile = '{}.fitp'.format(l)
            d = i2p.parse_fitprofs(proffile,n)[1]
            mean = np.mean(d,axis=0)
            std = np.std(d,axis=0)
            rejidx = np.where(np.abs(d - mean) > std*threshold)
            i = 0
            while rejidx[0].size != 0:
                d = np.delete(d,rejidx[0],axis=0)
                mean = np.mean(d,axis=0)
                std = np.std(d,axis=0)
                rejidx = np.where(np.abs(d - mean) > std*threshold)
                i += rejidx[0].size
                
            diff = mean - c
            outstr = '{:} ({:}):\n\t{:>7}: {:}\n\t{:>7}: {:}\n\t{:>7}: {:}\n\t{:>7}: {:}\n'.\
                     format(l,c,'numrej',i,'mean',mean,'diff',diff,'std',std)
            prtstr = ''
            for j in range(len(c)):
                prtstr += '{} {} {}\n'.format(c[j],diff[j],std[j])

            print prtstr
            f.write(outstr)
            acax.plot(c,diff,'.k')
            stax.plot(c,std,'.k')
    
        f.write('\n\n')

    fig.subplots_adjust(hspace=0.0001)
    acax.set_xlim(*stax.get_xlim())
    plt.savefig('WLC.png')

    return

if __name__ == '__main__':
    
    for l in llist:
        proffile = '{}.fitp'.format(l)
        if os.path.exists(proffile):
            print "Removing", proffile
            os.system('rm '+ proffile)

    do_fitprof(sys.argv[1])
    get_results('WLC.dat')

