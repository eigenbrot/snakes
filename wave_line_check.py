#!/usr/bin/python
import sys
import os
import numpy as np
import ir2py_lines as i2p
from pyraf import iraf

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

llist = ['4047','4358','HK','5461']
rlist = ['4020 4075','4325 4400','3910 4010','5440 5480']
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

def get_results(output, threshold=2.):

    for l, n, c in zip(llist,numlist,centlist):
        proffile = '{}.fitp'.format(l)
        print proffile
        d = i2p.parse_fitprofs(proffile,n)[1]
        median = np.median(d,axis=0)
        rms = np.sqrt(np.mean((d - median)**2,axis=0))
        rejidx = np.where(np.abs(d - median) > rms*threshold)
        i = 0
        while rejidx[0].size != 0:
            d = np.delete(d,rejidx[0],axis=0)
            rms = np.sqrt(np.mean((d - median)**2,axis=0))
            rejidx = np.where(np.abs(d - median) > rms*threshold)
            i += rejidx[0].size
             
        outstr = '{:} ({:}):\n\t{:>7}: {:}\n\t{:>7}: {:}\n\t{:>7}: {:}\n'.\
            format(l,c,'numrej',i,'median',median,'rms',rms)

        print outstr
        with open(output,'a') as f:
            f.write(outstr)

    return

if __name__ == '__main__':
    
    for l in llist:
        proffile = '{}.fitp'.format(l)
        if os.path.exists(proffile):
            print "Removing", proffile
            os.system('rm '+ proffile)

    do_fitprof(sys.argv[1])
    get_results('WLC.dat')

