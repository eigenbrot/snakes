import numpy as np
import pyfits
import matplotlib
import scipy.stats as ss
import ADESALT as sa
import ADEUtils as ADE
import time
from sklearn.mixture import GMM
from matplotlib.backends.backend_pdf import PdfPages as PDF
matplotlib.pyplot = plt
matplotlib.rc('font',size=5)

def clean_model(model, tol=3.):

    strong_idx = np.argmax(model.weights_)
    strong_mean = model.means_[strong_idx]
    strong_std = np.sqrt(model.covars_[strong_idx])

    for i in range(model.n_components):
        if np.abs(model.means_[i] - strong_mean) > strong_std*tol:
            model.weights_[i] = 0.0

    return

def MGD(slayfile,radius,cent_lambda=5048.126,
        flip=False,window=15,tol=3.):

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
    clean_model(M_best,tol=tol)
    mV = np.linspace(V.min(),V.max(),1000.)
    lp, resp = M_best.score_samples(mV)
    model_pdf = np.exp(lp)
    model_components = resp * model_pdf[:,np.newaxis]

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(mV,model_pdf,'k')
    ax1.plot(mV,model_components,'--k')
    ax1.errorbar(V,I/np.max(I)*np.max(model_pdf),
                 yerr=err/np.max(I)*np.max(model_pdf),alpha=0.4)
    ax1.hist(line_sample,50,normed=True,alpha=0.8,histtype='step')
    ax1.set_xlabel('Velocity [km/s]')
    ax1.set_ylabel('Normalized counts')

    ax2 = fig.add_subplot(323)
    ax2.plot(N,BIC,'--k',label='BIC')
    ax2.plot(N,AIC,'-k',label='AIC')
    ax2.legend(loc=0)
    ax2.set_xlabel('Number of Gaussian components')
    ax2.set_ylabel('Information Criterion')

    # fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))

    # fig.show()

    return mV, model_pdf, V, I, err, rwidth, fig

def get_window(x, p, cdf_window=0.1, ax=None):

    cdf = np.cumsum(p/np.sum(p))

    xlow = np.interp(cdf_window,cdf,x)
    xhigh = np.interp(1-cdf_window,cdf,x)

    if ax:
        ax.plot(x,cdf,'k')
        ax.plot(x,p/p.max(),'--k',alpha=0.5)
        ax.axvline(x=xlow,alpha=0.7,linestyle='-.',color='k')
        ax.axvline(x=xhigh,alpha=0.7,linestyle='-.',color='k')
        ax.set_xlabel('Velocity [km/s]')

    return xlow, xhigh

def cocaine(slayfile, radius, flip=False, cdf_window=0.05, tol=3.):

    model_V, model_pdf, V, I, err, rwidth, fig = MGD(slayfile,radius,
                                                     flip=flip, tol=tol)
    
    ax = fig.add_subplot(324)
    lowV, highV = get_window(model_V, model_pdf, 
                             cdf_window=cdf_window, ax=ax)

    idx = np.where((V >= lowV) & (V <= highV))
    moments, moment_err = ADE.ADE_moments(V[idx],I[idx],err=err[idx])
    
    ax2 = fig.add_subplot(313)
    ax2.errorbar(V,I,yerr=err)
    ax2.axvline(x=lowV,alpha=0.7,linestyle='--')
    ax2.axvline(x=highV,alpha=0.7,linestyle='--')
    ax2.set_xlabel('Velocity [km/s]')
    ax2.set_ylabel('Counds [ADU]')

    # fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))
    # fig.show()

    return moments, moment_err, rwidth, fig

def get_drunk(slayfile, baseoutput, flip=False, cdf_window=0.05, tol=3.):

    fitsout = baseoutput+'.fits'
    pdfout = baseoutput+'.pdf'
    pp = PDF(pdfout)

    radii, _, _ = sa.openslay(slayfile,flip=flip)
    
    hdulist = []
    widths = np.array([])
    for i, radius in enumerate(radii):
        print 'Getting r = {:5.3} kpc'.format(radius)

        moments, moment_err, rwidth, fig = cocaine(slayfile, radius, 
                                                   flip=flip,tol=tol, 
                                                   cdf_window=cdf_window)

        widths = np.append(widths,rwidth)
        fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))
        pp.savefig(fig)
        hdu = pyfits.ImageHDU(np.transpose(
                np.dstack((moments,moment_err)),(0,2,1)))
        hdu.header.update('EXTNAME','Radius_{}'.format(i))
        hdu.header.update('RADIUS',radius,comment='Radius [kpc]')
        hdulist.append(hdu)

    radiiHDU = pyfits.PrimaryHDU(np.vstack((radii,widths)))
    radiiHDU.header.update('SLAY',slayfile,comment='Slay file')
    radiiHDU.header.update('DATE',time.asctime(),comment='Date of extraction')
    radiiHDU.header.update('FLIP',flip,comment='Flip radii around 0?')
    radiiHDU.header.update('CDF_WIN',cdf_window,comment='Limits on CDF window')
    radiiHDU.header.update('EXTNAME','Radii')
    pyfits.HDUList([radiiHDU] + hdulist).writeto(fitsout)
    pp.close()

    return
