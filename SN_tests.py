import numpy as np
import tau_model as tm
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import os
import time

def make_monte(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60],
               N = 10, lightmin = 5450., lightmax = 5550.):
    
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    ssp = '/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_008Z_ChabIMF.fits'
    if type(N) is not list:
        N = [N]
    for SN in SNlist:
        direc = 'SN{}'.format(SN)
        if not os.path.exists(direc):
            os.makedirs(direc)
        filename = '{}/run.pro'.format(direc)
        if os.path.exists(filename):
            n=2
            filename = '{}_{}.pro'.format(filename.split('.pro')[0],n)
            while os.path.exists(filename):
                n += 1
                filename = '{}_{}.pro'.format(filename.split('_')[0],n)

        f = open(filename,'w')
        for tau in taulist:
            for i in range(*N):
                name = '{}/SN{:02}_t{:03}_N{:03}'.format(direc,SN,tau,i+1)
                print name
                if i == N[0] - 1:
                    tm.make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau,lightmin = lightmin, lightmax = lightmax,makeplot=True)
                else:
                    tm.make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau,lightmin = lightmin, lightmax = lightmax,makeplot=False)
                for z in range(fraclist.size):
                    f.write("do_simple, '{0:}.ms_lin.fits', '{0:}.me_lin.fits', '{0:}_Z{1:04}_fit.dat', wavemin=3750., wavemax=6800., lightmin={2:}, lightmax={3:}, model='{4:}'\n".format(os.path.basename(name),int(fraclist[z]*1000),lightmin,lightmax,modellist[z]))
        f.close()
    return

def compare_SN(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60], N = 10):

    ratios = np.zeros((len(SNlist),len(taulist)))
    errs = np.zeros(ratios.shape)
    agelist = np.zeros(len(taulist))

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
                fMMWA, fMLWA, ftauV, fZ = np.loadtxt(fit_name,usecols=(11,12,13,17),unpack=True)
                # tmp[i] = 1 - fMLWA/mMLWA
                tmp[i] = 1 - fMMWA/mMMWA
                # tmp[i] = 1 - ftauV/1.5
                # tmp[i] = 1 - fZ/0.4

            agelist[t] = mMMWA
            ratios[s,t] = np.mean(tmp)
            errs[s,t] = np.sqrt(np.mean((tmp - np.mean(tmp))**2))

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(211)
    eax = fig.add_subplot(212)
    colors = ['r','g','c','b','m']
    for i in range(len(taulist)):
        ax.errorbar(SNlist, ratios[:,i], yerr=errs[:,i], color=colors[i], 
                    label='{:} $\Rightarrow$ {:4.2f} Gyr'.format(taulist[i],agelist[i]))
        # ax.plot(SNlist,errs[:,i], label=taulist[i])

        eax.plot(SNlist, errs[:,i], color=colors[i])

    fig.subplots_adjust(hspace=0.0001)
    ax.set_xticklabels([])
    ax.set_xlim(0,70)
    ax.set_ylim(ax.get_ylim()[0]*0.9,ax.get_ylim()[1])
    eax.set_xlim(ax.get_xlim())
    eax.set_ylim(eax.get_ylim()[0],eax.get_ylim()[1]*0.9)
#    ax.set_xlabel('SNR')
#    ax.set_ylabel('1 - MLWA$_{fit}$/MLWA$_{true}$')
    ax.set_ylabel('1 - MMWA$_{fit}$/MMWA$_{true}$')
#    ax.set_ylabel(r'1 - $\tau_\mathrm{V,fit}$/$\tau_\mathrm{V,true}$')
#    ax.set_ylabel(r'1 - Z$_\mathrm{fit}$/Z$_\mathrm{true}$')
#    ax.set_ylim(-2,0.5)
    ax.legend(loc=0,numpoints=1,scatterpoints=1,frameon=False,title=r'$\tau_{sf}\Rightarrow\mathrm{MMWA}$')

    eax.set_xlabel('SNR')
    eax.set_ylabel('RMS')

    return ratios, errs, SNlist, taulist, fig

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
