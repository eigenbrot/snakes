import matplotlib
#matplotlib.use('agg')
import numpy as np
import tau_model as tm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
import os
import time
plt.ioff()

def make_monte(taulist = [0.1,1,2,4,10], SNlist = [5,10,20,40,60],
               N = 10, SNmin = 5450., SNmax = 5550., lightmin = 5450., lightmax = 5550.):
    
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
                    tm.make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau,lightmin = lightmin, lightmax = lightmax,SNmin = SNmin, SNmax = SNmax,
                                   makeplot=True)
                else:
                    tm.make_galaxy(name,SSPs=ssp,SN=SN,tau_sf=tau,lightmin = lightmin, lightmax = lightmax,SNmin = SNmin, SNmax = SNmax,
                                   makeplot=False)
                for z in range(fraclist.size):
                    f.write("do_simple, '{0:}.ms_lin.fits', '{0:}.me_lin.fits', '{0:}_Z{1:04}_fit.dat', wavemin=3750., wavemax=6800., lightmin={2:}, lightmax={3:}, model='{4:}'\n".format(os.path.basename(name),int(fraclist[z]*1000),lightmin,lightmax,modellist[z]))
        f.close()
    return

def compare_SN(taulist = [0.1,1,2,4,10], SNlist = [3,5,7,10,15,20,30,40,60], 
               dirpre = '', N = 30, output=None, quant='t'):

    ratios = np.zeros((len(SNlist),len(taulist)))
    errs = np.zeros(ratios.shape)
    agelist = np.zeros(len(taulist))

    if len(dirpre) > 0 and dirpre[-1] != '/':
        dirpre += '/'

    for s, SN in enumerate(SNlist):
        direc = '{}SN{}'.format(dirpre,SN)
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
                if quant == 't':
                    tmp[i] = 1 - fMLWA/mMLWA
                elif quant == 'a':
                    tmp[i] = 1 - ftauV/1.5
                elif quant == 'z':
                    tmp[i] = 1 - fZ/0.4
                # tmp[i] = 1 - fMMWA/mMMWA

            agelist[t] = mMLWA
            ratios[s,t] = np.mean(tmp)
            errs[s,t] = np.sqrt(np.mean((tmp - np.mean(tmp))**2))

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(211)
    eax = fig.add_subplot(212)
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
    for i in range(len(taulist)):
        ax.plot(SNlist, ratios[:,i], color=colors[i], lw=1.5,
                label='{:} $\Rightarrow$ {:4.2f} Gyr'.format(taulist[i],agelist[i]))
        # ax.errorbar(SNlist, ratios[:,i], yerr=errs[:,i], color=colors[i], 
        #             label='{:} $\Rightarrow$ {:4.2f} Gyr'.format(taulist[i],agelist[i]))
        # ax.plot(SNlist,errs[:,i], label=taulist[i])

        eax.plot(SNlist, errs[:,i], color=colors[i], lw=1.5)

    fig.subplots_adjust(hspace=0.0001)
    ax.set_xticklabels([])
    ax.set_xlim(0,70)
    if quant == 't':
        ax.set_ylim(0.04,1)
        ax.set_ylabel(r'$1 - \tau_{L,\mathrm{fit}}/\tau_{L,\mathrm{true}}$')
    elif quant == 'a':
        ax.set_ylim(-0.31,0)
        ax.set_ylabel(r'$1 - A_{V,\mathrm{fit}}/A_{V,\mathrm{true}}$')
    elif quant == 'z':
        ax.set_ylabel(r'$1 - Z_\mathrm{fit}/Z_\mathrm{true}$')

    eax.set_xlim(ax.get_xlim())
    eax.set_ylim(eax.get_ylim()[0],eax.get_ylim()[1]*0.9)
#    eax.set_ylim(-0.1,2)
#    ax.set_xlabel('SNR')
#    ax.set_ylabel('1 - MMWA$_{fit}$/MMWA$_{true}$')
#    ax.set_ylim(-2,0.5)
    ax.legend(loc=0,numpoints=1,scatterpoints=1,frameon=False,title=r'$\tau_{sf}\Rightarrow\tau_L$')

    eax.set_xlabel('SNR')
    eax.set_ylabel('RMS')

    if output is not None:
        pp = PDF(output)
        pp.savefig(fig)
        pp.close()

    return ratios, errs, SNlist, taulist, fig

def plot_covariance(SN, taulist = [0.1,1,2,4,10], N = 30, output = None):

    big_age = np.zeros((len(taulist),N))
    big_Av = np.zeros(big_age.shape)
    big_Z = np.zeros(big_age.shape)
    agelist = np.zeros(len(taulist))

    fig = plt.figure(figsize=(8,8))
    age_Av_ax = fig.add_subplot(221)
    age_Z_ax = fig.add_subplot(223)
    Av_Z_ax = fig.add_subplot(224)
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']

    direc = 'SN{}'.format(SN)
    for t, tau in enumerate(taulist):
        tmp_age = np.zeros(N)
        tmp_Av = np.zeros(N)
        tmp_Z = np.zeros(N)
        for i in range(N):
            model_name = '{}/SN{:02}_t{:03}_N{:03}_model.dat'.\
                         format(direc,SN,tau,i+1)
            fit_name = '{}/SN{:02}_t{:03}_N{:03}_fit.dat'.\
                       format(direc,SN,tau,i+1)
            print model_name
            mMMWA, mMLWA = np.loadtxt(model_name,
                                      usecols=(11,12),unpack=True)
            fMMWA, fMLWA, ftauV, fZ = np.loadtxt(fit_name,usecols=(11,12,13,17),unpack=True)
            print fZ
            tmp_age[i] = 1 - fMLWA/mMLWA
            tmp_Av[i] = 1 - ftauV/1.5
            tmp_Z[i] = 1 - fZ/0.4
                
        age_Av_ax.scatter(tmp_age, tmp_Av, color=colors[t], label='{:} $\Rightarrow$ {:4.2f} Gyr'.format(tau,mMLWA))
        age_Z_ax.scatter(tmp_age, tmp_Z, color=colors[t])
        Av_Z_ax.scatter(tmp_Av, tmp_Z, color=colors[t])

    # Av_Z_ax.set_ylim(-0.1,0.1)
    # age_Z_ax.set_ylim(*Av_Z_ax.get_ylim())
    age_Av_ax.set_xticklabels([])
    Av_Z_ax.set_yticklabels([])
    age_Av_ax.set_ylabel('$A_V$')
    age_Z_ax.set_ylabel('$Z$')
    age_Z_ax.set_xlabel(r'$\tau_L$')
    Av_Z_ax.set_xlabel('$A_V$')
    fig.subplots_adjust(hspace=0.001,wspace=0.001)

    return fig
    
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
