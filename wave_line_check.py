#!/usr/bin/python
import sys
import os
import numpy as np
import ir2py_lines as i2p
from pyraf import iraf

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)

llist = ['4047','4358','HK']
rlist = ['4020 4075','4325 4400','3910 4010']
centlist = [[4047],[4358],[3933.57,3968.53]]
numlist = [1,1,2]

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

def get_results(output):

    for l, n, c in zip(llist,numlist,centlist):
        proffile = '{}.fitp'.format(l)
        print proffile
        d = i2p.parse_fitprofs(proffile,n)
        median = np.median(d[1],axis=0)
        rms = np.sqrt(np.mean((d[1] - median)**2,axis=0))
    
        print '{:} ({:}):\n\t{:>7}: {:}\n\t{:>7}: {:}\n'.\
            format(l,c,'median',median,'rms',rms)

        with open(output,'a') as f:
            f.write('{:} ({:}):\n\t{:>7}: {:}\n\t{:>7}: {:}\n'.\
                    format(l,c,'median',median,'rms',rms))

    return

if __name__ == '__main__':
    
    for l in llist:
        proffile = '{}.fitp'.format(l)
        if os.path.exists(proffile):
            print "Removing", proffile
            os.system('rm '+ proffile)

    do_fitprof(sys.argv[1])
    get_results('WLC.dat')

