import numpy as np
import scipy.stats as ss
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
from glob import glob
plt.ioff()

def do_allz2(stepfile, output, bluefree=652):

    MLWA, MLWZ, TauV, blueChi = np.loadtxt(stepfile,usecols=(62,64,65,69),
                                           unpack=True)
    blueChi *= bluefree
    minidx = np.where(blueChi == np.min(blueChi))
    deltaChi = blueChi - np.min(blueChi)

    ci = [0.5, 0.68, 0.9]
    dc = [ss.chi2.ppf(c,bluefree) for c in ci]
    dcidx = [np.where(deltaChi <= d) for d in dc]
    
    pidx = np.where(deltaChi < 8000)[0][::-1]
    # ax3d = plt.figure().add_subplot(111, projection='3d')
    # ax3d.set_xlabel(r'$\tau_V$')
    # ax3d.set_ylabel('MLWA')
    # ax3d.set_zlabel('MLWZ')
    # ax3d.scatter(TauV[pidx],MLWA[pidx],zs=MLWZ[pidx],c=deltaChi[pidx],
    #              alpha=0.5)
    # ax3d.figure.show()

    taumin = np.min(TauV[pidx]) * 0.95
    taumax = np.max(TauV[pidx]) * 1.05
    agemin = np.min(MLWA[pidx]) * 0.95
    agemax = np.max(MLWA[pidx]) * 1.05
    Zmin = np.min(MLWZ[pidx]) * 0.95
    Zmax = np.max(MLWZ[pidx]) * 1.05

    fig = plt.figure(figsize=(11,9))
    ax0 = fig.add_subplot(2,2,1)
    scat = ax0.scatter(TauV[pidx],MLWA[pidx],c=blueChi[pidx],alpha=0.6,s=40,
               linewidths=0,cmap=plt.cm.gnuplot2)
    ax0.set_xlim(taumin, taumax)
    ax0.set_ylim(agemin, agemax)
    ax0.set_ylabel('MLWA')
    ax0.set_xticklabels([])

    ax1 = fig.add_subplot(2,2,3)
    ax1.scatter(TauV[pidx],MLWZ[pidx],c=blueChi[pidx],alpha=0.6,s=40,
                linewidths=0,cmap=plt.cm.gnuplot2)
    ax1.set_xlim(taumin, taumax)
    ax1.set_ylim(Zmin, Zmax)
    ax1.set_xlabel(r'$\tau_V$')
    ax1.set_ylabel('MLWZ')

    ax2 = fig.add_subplot(2,2,4)
    ax2.scatter(MLWA[pidx],MLWZ[pidx],c=blueChi[pidx],alpha=0.6,s=40,
               linewidths=0,cmap=plt.cm.gnuplot2)
    ax2.set_xlim(agemin, agemax)
    ax2.set_ylim(Zmin, Zmax)
    ax2.set_xlabel('MLWA')
    ax2.set_yticklabels([])

    fig.subplots_adjust(hspace=0.0001,wspace=0.0001)
    cb = fig.colorbar(scat,ax=[ax0,ax1,ax2])
    cb.set_label(r'$\chi^2$')

    figname = '{}.png'.format(output.split('.dat')[0])
    fig.savefig(figname)
    plt.close(fig)

    bests = []
    bigCI = []
    for val, name in zip([MLWA, MLWZ, TauV], ['MLWA','MLWZ','TauV']):
        best_val = np.mean(val[minidx])
        CIs = [np.max(val[id]) - np.min(val[id]) for id in dcidx]
        bests.append(best_val)
        bigCI.append(CIs)
        print '{:} = {:5.3f} +/- {:5.3f} {:5.3f} {:5.3f}'.format(*([name, best_val] + CIs))

    return bests, np.array(bigCI)

def do_allz2_pointing(stepfolder, pointing, output, bluefree=652):
    
    aplist = glob('{}/*P{}*steps.dat'.format(stepfolder,int(pointing)))

    with open(output,'w') as f:
        
        f.write("""# Generated on {}
# Confidance intervals using data in {} for pointing {}
#\n""".format(time.asctime(), stepfolder, pointing))
        f.write('#{:>6}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n'.format('Apnum',
                                                                      'MLWA','dMLWA',
                                                                      'MLWZ','dMLWZ',
                                                                      'TauV','dTauV'))
        
        for apfile in aplist:
            num = int(apfile.split('_')[-2]) + 1
            print num, ':'
            best, CI = do_allz2(apfile, apfile, bluefree)
            f.write('{:7n}'.format(num))
            for b in zip(best, CI[:,1]):
                f.write('{:10.4f}{:10.4f}'.format(*b))
            f.write('\n')
                    
    return
