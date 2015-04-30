import numpy as np
import tau_model as tm
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os
import time

def make_monte(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60],
               N = 100):
    
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    ssp = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_008Z_ChabIMF.fits'
    for SN in SNlist:
        direc = 'SN{}'.format(SN)
        if not os.path.exists(direc):
            os.makedirs(direc)
        f = open('{}/run.pro'.format(direc),'w')
        for tau in taulist:
            for i in range(N):
                name = '{}/SN{:02}_t{:03}_N{:03}'.format(direc,SN,tau,i+1)
                print name
                tm.make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau)
                for z in range(fraclist.size):
                    f.write("do_simple, '{0:}.ms_lin.fits', '{0:}.me_lin.fits', '{0:}_Z{1:04}_fit.dat', wavemin=3750., wavemax=6800., model='{2:}', /plot\n".format(os.path.basename(name),int(fraclist[z]*1000),modellist[z]))
        f.close()
    return

def compare_SN(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60], N = 10):

    ratios = np.zeros((len(SNlist),len(taulist)))
    errs = np.zeros(ratios.shape)

    for s, SN in enumerate(SNlist):
        direc = 'SN{}'.format(SN)
        for t, tau in enumerate(taulist):
            tmp = np.zeros(N)
            for i in range(N):
                model_name = '{}/SN{:02}_t{:03}_N{:03}_model.dat'.\
                             format(direc,SN,tau,i+1)
                fit_name = '{}/SN{:02}_t{:03}_N{:03}_fit.dat'.\
                           format(direc,SN,tau,i+1)
                print model_name
                mMMWA, mMLWA = np.loadtxt(model_name,
                                          usecols=(11,12),unpack=True)
                fMMWA, fMLWA = np.loadtxt(fit_name,usecols=(11,12),unpack=True)
                tmp[i] = 1 - fMMWA/mMMWA
                ratios[s,t] = np.mean(tmp)
                errs[s,t] = np.std(tmp)

    ax = plt.figure().add_subplot(111)
    for i in range(len(taulist)):
        ax.errorbar(SNlist, ratios[:,i], yerr=errs[:,i], label=taulist[i])
#        ax.plot(SNlist,errs[:,i], label=taulist[i])

    ax.set_xlim(0,70)
    ax.set_xlabel('SNR')
#    ax.set_ylabel(r'$\sigma$(1 - MLWA$_{fit}$/MLWA$_{true}$)')
    ax.set_ylabel('1 - MMWA$_{fit}$/MMWA$_{true}$')
    ax.legend(loc=0,numpoints=1,scatterpoints=1,frameon=False,title=r'$\tau_{sf}$')

    return ratios, errs, SNlist, taulist, ax

def get_metal(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60],
               N = 100):

    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    
    for SN in SNlist:
        direc = 'SN{}'.format(SN)
        for tau in taulist:
            for i in range(N):
                tmp = np.zeros((fraclist.size,18))
                for z in range(fraclist.size):
                    name = '{}/SN{:02}_t{:03}_N{:03}_Z{:04}_fit.dat'.format(direc,SN,tau,i+1,int(fraclist[z]*1000))
                    print name
                    tmp[z] = np.loadtxt(name)

                bdx = np.argmin(tmp[:,16])
                h = open(name,'r')
                head = h.readlines()[4]
                f = open('{}/SN{:02}_t{:03}_N{:03}_fit.dat'.format(direc,SN,tau,i+1),'w')
                f.write('# Generated on {}\n'.format(time.asctime()))
                f.write(head)
                f.write(str('{:11n}'+12*'{:13.3e}'+'{:7.2f}{:12.3f}'+2*'{:12.3e}').format(*tmp[bdx][:-1]))
                f.write('{:10.3f}\n'.format(fraclist[bdx]))
                f.close()
                h.close()
    return
