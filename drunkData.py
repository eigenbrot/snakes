import numpy as np
import pyfits
import matplotlib
import scipy.stats as ss
import ADESALT as sa
import ADEUtils as ADE
import time
from sklearn.mixture import GMM
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt = matplotlib.pyplot

# Make the plots look nice
matplotlib.rc('axes',labelsize=9)
matplotlib.rc('xtick',labelsize=9)
matplotlib.rc('ytick',labelsize=9)
matplotlib.rc('legend',fontsize=7,frameon=False)
#matplotlib.rc('font',size=9,family='serif',serif=['Computer Modern Roman'])
#matplotlib.rc('text',usetex=True)
matplotlib.rc('axes',linewidth=0.6,labelsize=9)
matplotlib.rc('lines',linewidth=0.6)

def clean_model(model, tol=3.):

    strong_idx = np.argmax(model.weights_)
    strong_mean = model.means_[strong_idx]
    strong_std = np.sqrt(model.covars_[strong_idx])

    for i in range(model.n_components):
        if np.abs(model.means_[i] - strong_mean) > strong_std*tol:
            model.weights_[i] = 0.0
        if np.sqrt(model.covars_[i]) > 80.:
            model.weights_[i] = 0.0

    return

def MGD(slayfile,radius,cent_lambda=5048.126,
        flip=False,window=15,tol=3.,baseline=1):

    V, I, err, rwidth = sa.plot_line(slayfile,radius,wavelength=cent_lambda,
                                     velo=True,plot=False,baseline=baseline,
                                     window=window,flip=flip)
    factor = 1e6
    VV = V*factor
    II = np.copy(I)
    II[II < 0] = 0.0
    II /= np.sum(II)
    line_pdf = ss.rv_discrete(a=VV.min(),b=VV.max(),values=(VV,II))
    line_sample = line_pdf.rvs(size=10000)/factor

    N = np.arange(1,5)
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
    ax1.hist(line_sample,30,normed=True,alpha=0.8,histtype='step')
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

def cocaine(slayfile, radius, flip=False, cdf_window=0.05, tol=3.,
            window=15,baseline=1,interact=False,cent_lambda=5048.126):
    
    again = True
    argdict = {'flip': flip, 'tol': tol, 'window': window, 'baseline': baseline,
               'cent_lambda':cent_lambda}
    while again:
        model_V, model_pdf, V, I, err, rwidth, fig = MGD(slayfile,radius,
                                                         **argdict)
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
        ax2.set_ylabel('Counts [ADU]')

        # fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))
        # fig.show()
        
        if interact == True:
            print "Keys are {}".format(argdict.keys())
            print "Change any you want and then hit 'r'. 'q' will quit"
            scratch = raw_input(':')
            while scratch != 'q' or scratch != 'r':
                print 'Changing key {}'.format(scratch)
                val = input('To (current val is {}): '.\
                                format(argdict[scratch]))
                argdict[scratch] = val
                scratch = raw_input('again?:')
            if scratch == 'r':
                continue
            elif scratch == 'q':
                break
        else:
            again = False
        

    return moments, moment_err, rwidth, highV - lowV, fig

def get_drunk(slayfile, baseoutput, flip=False, cdf_window=0.05, tol=3.,
              window=15,baseline=1,skip_radii=[],cent_lambda=5048.126):

    fitsout = baseoutput+'.fits'
    pdfout = baseoutput+'.pdf'
    pp = PDF(pdfout)

    radii, _, _ = sa.openslay(slayfile,flip=flip)
    
    outradii = np.array([])
    widths = np.array([])
    vwidths = np.array([])
    m1 = np.empty((2,1))
    m2 = np.empty((2,1))
    m3 = np.empty((2,1))
    for i, radius in enumerate(radii):
        if int(np.floor(radius)) in skip_radii:
            print "user skipping radius {} kpc".format(radius)
            continue
        print 'Getting r = {:5.3} kpc'.format(radius)
        
        moments, moment_err, rwidth, vwidth, fig = cocaine(slayfile, radius, 
                                                           flip=flip,tol=tol, 
                                                           cdf_window=cdf_window,
                                                           baseline=baseline,
                                                           window=window,
                                                           cent_lambda=cent_lambda)

        if np.isnan(np.sum(moments)):
            print 'skipping radius {} kpc due to NaN condition'.format(radius)
            continue
        m1 = np.hstack((m1,np.array([[moments[0],moment_err[0]]]).T))
        m2 = np.hstack((m2,
                        np.array([[np.sqrt(moments[1]),
                                   moment_err[1]/(2*np.sqrt(moments[1]))]]).T))
        m3 = np.hstack((m3,np.array([[moments[2],moment_err[2]]]).T))
        widths = np.append(widths,rwidth)
        vwidths = np.append(vwidths,vwidth)
        outradii = np.append(outradii,radius)
        fig.tight_layout(pad=0.5)
        fig.suptitle('r = {:5.3f} kpc\n{:}'.format(radius,time.asctime()))
        pp.savefig(fig)
        
    m1 = m1[:,1:]
    m2 = m2[:,1:]
    m3 = m3[:,1:]
    radiiHDU = pyfits.PrimaryHDU(np.vstack((outradii,widths,vwidths)))
    radiiHDU.header.update('SLAY',slayfile,comment='Slay file')
    radiiHDU.header.update('DATE',time.asctime(),comment='Date of extraction')
    radiiHDU.header.update('FLIP',flip,comment='Flip radii around 0?')
    radiiHDU.header.update('CDF_WIN',cdf_window,comment='Limits on CDF window')
    radiiHDU.header.update('TOL',tol,comment='Gaussian rejection tolerance')
    radiiHDU.header.update('EXTNAME','Radii')
    
    m1HDU = pyfits.ImageHDU(m1)
    m1HDU.header.update('EXTNAME','mu1')
    m2HDU = pyfits.ImageHDU(m2)
    m2HDU.header.update('EXTNAME','sqrt(mu2)')
    m3HDU = pyfits.ImageHDU(m3)
    m3HDU.header.update('EXTNAME','mu3')

    pyfits.HDUList([radiiHDU,m1HDU,m2HDU,m3HDU]).writeto(fitsout,clobber=True)
    pp.close()
    plt.close('all')

    return

def plot_moments(moment_file):

    hdus = pyfits.open(moment_file)
    radii = hdus[0].data[0]
    
    fig = plt.figure()
    ax3 = fig.add_subplot(313)
    ax1 = fig.add_subplot(311,sharex=ax3)
    ax2 = fig.add_subplot(312,sharex=ax3)
    plt.setp(ax1.get_xticklabels(),visible=False)
    plt.setp(ax2.get_xticklabels(),visible=False)
    fig.subplots_adjust(hspace=0.0001)
    ax3.set_xlabel('radius [kpc]')
    ax1.set_ylim(-260.001,260.001)
    ax2.set_ylim(15.00001,49.99999)
    ax3.set_ylim(-1,1.001)
    matplotlib.rc('text',usetex=False)
    fig.suptitle('{}\n{}'.format(moment_file,time.asctime()))
    matplotlib.rc('text',usetex=True)

    for moment, ax, name in zip(hdus[1:],[ax1,ax2,ax3],
                                ['$\mu_1$','$\sqrt{\mu_2}$','$\mu_3$']):
        
        ax.errorbar(radii,moment.data[0],yerr=moment.data[1],fmt='.',ms=7)
        ax.set_ylabel(name)

    fig.show()

    return fig

def open_drunk(drunkfile,skip_radii=[]):

    hdus = pyfits.open(drunkfile)
    
    radii = hdus[0].data[0]
    vwidths = hdus[0].data[2]
    rwidths = hdus[0].data[1]

    m1 = hdus[1].data
    m2 = hdus[2].data
    m3 = hdus[3].data

    for i, r in enumerate(radii):
        if int(np.floor(r)) in skip_radii:
            radii = np.delete(radii,i)
            rwidths = np.delete(rwidths,i)
            vwidths = np.delete(vwidths,i)
            m1 = np.delete(m1,i,axis=1)
            m2 = np.delete(m2,i,axis=1)
            m3 = np.delete(m3,i,axis=1)

    return radii, rwidths, vwidths, m1, m2, m3
