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

