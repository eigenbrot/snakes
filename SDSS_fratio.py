import numpy as np
import matplotlib.pyplot as plt
import pyfits
import ADEUtils as ADE
from glob import glob
from matplotlib.backends.backend_pdf import PdfPages as PDF
from datetime import datetime

def kegstand(searchstr,EEfigs,Nfig,EE=0.50):
    """Given an input list of fits files this function will try to compute the
    f-ratio of the SDSS test stand input beam.

    """
    
    R = np.array([])
    dprime = np.array([])

    pp = PDF(EEfigs)
    pp2 = PDF(Nfig)

    in_files = glob(searchstr)
    
    for image in in_files:
        print image
        HDU = pyfits.open(image)[0]
        dp = float(image.split('_')[1].split('.fits')[0])
        dprime = np.append(dprime, dp)
        radius = get_radius(HDU.data,pp,EE)
        R = np.append(R, radius)
        print "d: {:4.3f}, R: {:4.3f}".format(dp,radius)

    dprime *= -1

    # Ns = np.array([])
    # for k in range(numtrys):
    #     sampleidx = np.random.randint(dprime.size, size = dprime.size)
    #     tempdp = dprime[sampleidx]
    #     tempR = R[sampleidx]
    #     Ns = np.append(Ns,
    #                        (2.*ADE.polyclip(tempdp,tempR,1,niter=10).c[0])**-1)

    fit_coef = ADE.polyclip(dprime,R,1,niter=10).c
    slope = fit_coef[0]

    ### compute the uncertainty
    

    N = (2.*slope)**-1
#    N = np.mean(Ns)
#    N_err = np.std(Ns)
    Nfit = np.poly1d(fit_coef)
    fitd = np.linspace(dprime.min(),dprime.max(),50)
    
    ### compute the uncertainty
    see = (np.sum((R - np.polyval(fit_coef,dprime))**2)/(R.size - 2))**0.5
    slope_err = see * (1/(np.sum((dprime - np.mean(dprime))**2)))**0.5
    N_err = slope_err/(2.*slope**2)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dprime,R,marker='s',linestyle='')
    ax.plot(fitd,Nfit(fitd),'k:',label='linear fit')
    ax.set_xlabel("$-d^'$ [mm]")
    ax.set_ylabel('$R$ [mm]')
    ax.legend(loc=0)
    ax.set_title('{}\nN: {:3.2f}$\pm$ {:3.2f}'.\
                     format(datetime.now().isoformat(' '),N,N_err))

    pp2.savefig(fig)

    pp.close()
    pp2.close()

    return N, N_err

def get_radius(data,pp,EEcut):
    """Takes in a 2D numpy array and computes the radius of the beam profile.
    This function uses a parabola fit to find the true beam radius.

    """
    
    r, sb, err = ADE.fast_annulize(data,300)
    r *= 0.0044
    
    flux = np.cumsum(sb)
    EE = flux/flux.max()

    cutr = np.where(EE >= EEcut)[0][0]

    EEfit = np.poly1d(np.polyfit(r[:cutr],EE[:cutr],2))

    fitr = np.linspace(r.min(),r.max(),500)
    fitEE = EEfit(fitr)

    r1 = np.interp(1.0,fitEE,fitr)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(r,EE,marker='.',linestyle='',markersize=0.7)
    ax.plot(fitr,fitEE,'-',alpha=0.4)
    ax.set_xlabel('r [mm]')
    ax.set_ylabel('EE')
    ax.axvline(r1,linestyle='-',alpha=0.4)
    ax.axhline(1.0,linestyle=':',color='k',alpha=0.2)
    ax.set_ylim(0,1.1)
#    ax.set_xlim(0,1.5*r1)
    ax.set_title('r: {:3.2f} mm'.format(r1))

    if pp:
        pp.savefig(fig)

    return r1

def clean_VI(searchstr):
    """Takes FITS files generated by LabView and:
    1) Poor man's bkgrnd subtraction with ADE.mediclean
    2) Sets the Primary Hdu to 0 instead of LabView's 1

    """

    in_files = glob(searchstr)

    for image in in_files:

        HDU = pyfits.open(image)[1]
        cleaned_data = ADE.mediclean(HDU.data[220:870,290:960])
        clean_name = image.split('.fits')[0].split('/')[-1] + '.cleaned.fits'
        pyfits.PrimaryHDU(cleaned_data,HDU.header).writeto(clean_name)

    return
        
def radius_test(search_str):

    file_list = glob(search_str)

    zs = np.array([])
    rs = np.array([])

    for image in file_list:
        print image
        z = float(image.split('_')[1][0:4])
        data = pyfits.open(image)[0].data
        r = get_radius(data,None,0.2)/0.0044

        zs = np.append(zs,z)
        rs = np.append(rs,r)

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('z')
    ax.set_ylabel('r')
    ax.plot(zs,rs)

    ax.get_figure().show()

    return ax
