import glob
import time
import numpy as np
import pyfits
import GradPak_bin as GPB
import GradPak_plot as GPP

glob = glob.glob

def step1(SN, wavemin=3800., wavemax=6800.):

    exclude_master = [[],[],[92,93],[107,105],[92],[106,65]]
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    f = open('bin{}.pro'.format(SN),'w')

    for p in range(6):
        infile = glob('../P{}/*ms_rfsz_*fits'.format(p+1))[0]
        errfile = glob('../P{}/*me_rfz_*fits'.format(p+1))[0]
        binname = 'NGC_891_P{}_bin{}'.format(p+1,SN)
        
        print 'data:  ',infile
        print 'error: ',errfile
        print 'output:',binname

        GPB.bin(infile, errfile, SN, binname, waverange=[wavemin,wavemax],
                exclude=exclude_master[p])

        idlin = '{}.ms.fits'.format(binname)
        idlerr = '{}.me.fits'.format(binname)

        for Z in range(fraclist.size):
            idlout = '{}_Z{:04}.dat'.format(binname,int(fraclist[Z]*1000))
            f.write("do_simple, '{}', '{}', '{}', lightmin={}, lightmax={}, wavemin=3750., wavemax=6800., model='{}'\n".\
                    format(idlin, idlerr, idlout, wavemin, wavemax,modellist[Z]))

    f.close()
    return

def step2():

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
   
    for p in range(6):
    
        basename = 'NGC_891_P{}'.format(p+1)
       
        for z in range(fraclist.size):
            name = glob('{}*_Z{:04}.dat'.format(basename,int(fraclist[z]*1000)))[0]
            print name
            try:
                tmp[z] = np.loadtxt(name)
            except UnboundLocalError:
                data = np.loadtxt(name)
                tmp = np.zeros((fraclist.size,data.shape[0],data.shape[1]))
                tmp[z] = data

        bdx = np.argmin(tmp[:,:,16],axis=0)
        h = open(name,'r')
        head = h.readlines()[4]
        outname = '{}_fit.dat'.format(name.split('_Z')[0])
        f = open(outname,'w')
        f.write('# Generated on {}\n'.format(time.asctime()))
        f.write(head)
        for i in range(tmp.shape[1]):
            f.write(str('{:11n}'+12*'{:13.3e}'+'{:7.2f}{:12.3f}'+2*'{:12.3e}').format(*tmp[bdx[i],i,:-1]))
            f.write('{:10.3f}\n'.format(fraclist[bdx[i]]))
            
        f.close()
        h.close()
        del tmp

    return
