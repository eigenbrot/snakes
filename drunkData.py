import numpy as np
import pyfits
import matplotlib.pyplot as plt
import scipy.stats as ss
import ADESALT as sa
import ADEUtils as ADE
import time
from sklearn.mixture import GMM
from matplotlib.backends.backend_pdf import PdfPages as PDF

def get_line_distro(slayfile,radius,cent_lambda=5048.126,
                    flip=False,window=15):

    V, I, err, rwidth = sa.plot_line(slayfile,radius,wavelength=cent_lambda,
                                     velo=True,plot=False,baseline=1,
                                     window=window,flip=flip)
    factor = 1e6
    VV = V*factor
    II = np.copy(I)
    II[II < 0] = 0.0
    II /= np.sum(II)
    line_pdf = ss.rv_discrete(a=VV.min(),b=VV.max(),values=(VV,II))
    line_sample = line_pdf.rvs(size=10000)/factor

    N = np.arange(1,11)
    models = [GMM(n).fit(line_sample) for n in N]

    AIC = [m.aic(line_sample) for m in models]
    BIC = [m.bic(line_sample) for m in models]

    M_best = models[np.argmin(AIC)]
    mV = np.linspace(V.min(),V.max(),1000.)
    lp, resp = M_best.score_samples(mV)
    model_pdf = np.exp(lp)
    model_components = resp * model_pdf[:,np.newaxis]

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(mV,model_pdf,'k')
    ax1.plot(mV,model_components,'--k')
    ax1.errorbar(V,I/np.max(I)*np.max(model_pdf),
                 yerr=err/np.max(I)*np.max(model_pdf),alpha=0.4)
    ax1.hist(line_sample,50,normed=True,alpha=0.8,histtype='step')
    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Normalized counts')

    ax2 = fig.add_subplot(212)
    ax2.plot(N,BIC,'--k',label='BIC')
    ax2.plot(N,AIC,'-k',label='AIC')
    ax2.legend(loc=0)
    ax2.set_xlabel('Number of Gaussian components')
    ax2.set_ylabel('Information Criterion')

    fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))

    fig.show()

    return mV, model_pdf, M_best

def get_window(x, p, cdf_window=0.1, ax=None):

    cdf = np.cumsum(p/np.sum(p))

    xlow = np.interp(cdf_window,cdf,x)
    xhigh = np.interp(1-cdf_window,cdf,x)

    if ax:
        ax.plot(x,cdf,'k')
        ax.plot(x,p/p.max(),'--k',alpha=0.5)
        ax.axvline(x=xlow,alpha=0.7,linestyle='-.',color='k')
        ax.axvline(x=xhigh,alpha=0.7,linestyle='-.',color='k')

    return xlow, xhigh

