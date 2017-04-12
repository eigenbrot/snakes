import pyfits
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF


def main(model_name='DFK_allZ_vardisp.fits', output='test.pdf'):

    models = pyfits.open(model_name)[1].data[0]
    fididx = [0,0,3,3]
    clist = ['#377eb8','#4daf4a','#984ea3','#e41a1c']
    Zlist = [0.2, 0.4, 1, 2.5]
    pp = PDF(output)
    #loop over DFK ages
    for i in range(4):
        print i
        monoage = models['FLUX'][i::4,:,2]

        #take out 2 lowest metallicities
        monoage = monoage[2:,:]

        best_pars = do_fits(monoage, models['WAVE'], fididx[i]).T
        best_taus = best_pars[0]
        best_norms = best_pars[1]
        # best_taus = np.zeros(3,2) + 5
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('DFK {}'.format(i+1))
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Normalized flux')
        redd = best_norms[:,None]*np.exp(-1 * best_taus[:,None] * (models['WAVE'][None,:]/5500.)**(-0.7))
        # redd = np.zeros((4, tmpredd.shape[1]))
        # kk = 0
        # for k in range(4):
        #     if k == fididx[i]:
        #         continue
        #     else:
        #         redd[k] = tmpredd[kk]
        #         kk += 1
        print redd.shape
        print best_norms, np.sum(redd, axis=1)/np.sum(monoage[:,:], axis=1)
        for j in range(4):
            if j == fididx[i]:
                ax.plot(models['WAVE'], monoage[j], color=clist[j], label=Zlist[j], lw=2, alpha=0.5)
            else:
                ax.plot(models['WAVE'], monoage[j]*redd[j], color=clist[j], label=Zlist[j])
        # ax.plot(models['WAVE'], monoage[1]*redd[1], color=clist[1], label=Zlist[1])
        # ax.plot(models['WAVE'], monoage[2]*redd[2], color=clist[2], label=Zlist[2])
        ax.text(0.5,0.095,'{:5.3f}, {:5.3f}'.format(best_taus[0], best_norms[0]), 
                color=clist[0], transform=ax.transAxes, fontsize=9)
        ax.text(0.5,0.06,'{:5.3f}, {:5.3f}'.format(best_taus[1], best_norms[1]), 
                color=clist[1], transform=ax.transAxes, fontsize=9)
        ax.text(0.5,0.035,'{:5.3f}, {:5.3f}'.format(best_taus[2], best_norms[2]), 
                color=clist[2], transform=ax.transAxes, fontsize=9)
        ax.text(0.5,0.01,'{:5.3f}, {:5.3f}'.format(best_taus[3], best_norms[3]), 
                color=clist[3], transform=ax.transAxes, fontsize=9)
        ax.text(0.3, 0.035, r'$\tau_V$, normalization', transform=ax.transAxes, fontsize=9)
        
        ax.legend(loc=0, title=r'$Z/Z_\odot$')
        
        pp.savefig(fig)
        plt.close(fig)

    pp.close()
    return

def do_fits(monoage, wave, fididx):

    fiducial = monoage[fididx]

    results = np.zeros((4,2))
    for i in range(4):
        p0 = [1.0,1.0]
        pf = spo.fmin(red_func, p0, args=(monoage[i], fiducial, wave))
        results[i] = np.array(pf)

    return results

def red_func(X, SSP, fiducial, wave):

    redd = np.exp(-1 * X[0] * (wave/5500.)**(-0.7))
    chisq = np.sum((fiducial - X[1]*SSP*redd)**2)

    return chisq
