import numpy as np
import ADEUtils as ADE
import ADESALT as sa
import gc
import Salty2 as salty
import bottleneck as bn
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF

def do_line(simfile,radius,peak_scale,plot=True,Iwidth=17,rwidth=1.,ax=None,label=None):

    v, I, _ = salty.line_profile(simfile,radius,plot=False,Iwidth=Iwidth,
                                 width=rwidth,fit=False) 


#    v, I = ADE.ADE_gauss(1000,500,50)
    I *= peak_scale/I.max()
    
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    l = ax.plot(v,I,label=label)[0]

    cdf = np.cumsum(I)
    cdf /= cdf.max()
    lowidx = int(np.interp(0.1,cdf,np.arange(cdf.size)))
    highidx = int(np.interp(0.9,cdf,np.arange(cdf.size)))

    ax.axvline(x=v[lowidx],alpha=0.4,color=l.get_color())
    ax.axvline(x=v[highidx],alpha=0.4,color=l.get_color())

    moments = ADE.ADE_moments(v[lowidx:highidx+1], I[lowidx:highidx+1], threshold=np.inf)

    if plot:
        fig.show()

    return moments, v, I, ax

def noise_line(simfile,radius,peak_scale):

#    v, I, _ = salty.line_profile(simfile,radius,plot=False)
    v, I = ADE.ADE_gauss(1000,500,50)
    I *= peak_scale/I.max()
    I += 3. * np.random.randn(I.size)
#    ADE.eplot(v,I)
    moments = ADE.ADE_moments(v, I, threshold=np.inf,err=np.abs(I)**0.5)

    return moments

def monte_lines(numtrys):
    
    bigarr = np.zeros((1,3))
    for i in range(numtrys):
        v, I = ADE.ADE_gauss(1000,500,50)
        I *= 55/I.max()
        I += 3. * np.random.randn(I.size)
    #    ADE.eplot(v,I)
        moments = ADE.ADE_moments(v, I, threshold=np.inf,err=np.abs(I)**0.5)
        bigarr = np.vstack((bigarr,moments))

    bigarr = bigarr[1:]
#    print bigarr

    return bn.nanmedian(bigarr,axis=0), bn.nanstd(bigarr,axis=0)

def line_comp(simfile,slayfile,plotprefix,width=17,flip=False):
    '''compares simulation and actual data on line shapes, using a moment based
    approach'''

    pp = PDF(plotprefix+'_lines.pdf')
    ppg = PDF(plotprefix+'_linegrid.pdf')
    gridfig = plt.figure()

    dradii, dcent, derr, dm1, dm2, dm3 = sa.openslay(slayfile,moments=True,
                                                     flip=flip)
    print dradii
    rwidth = np.abs(np.mean(np.diff(dradii)))
    tm1 = np.array([])
    tm2 = np.array([])
    tm3 = np.array([])
    central_lambda=np.array([4901.416,5048.126])

    print "{:^20} {:^20} {:^20}\n{:^10}{:^10}".format('m1','m2','m3','model','data')
    print ("-"*20+'    ')*3
    for i in range(dradii.size):
        
        radius = dradii[i]
        (mm1, mm2, mm3), V, I, _ = do_line(simfile,radius,1.,
                                        plot=False,Iwidth=width,rwidth=rwidth)

        # (nm1,nm2,nm3), negV, negI, _ = do_line(simfile,-1*radius,1.,
        #                        plot=False,Iwidth=width)

        # cent_lambda = central_lambda[np.where(
        #         np.abs(central_lambda - dm1[i,1]) == np.min(
        #             np.abs(central_lambda - dm1[i,1])))]
        cent_lambda = 5048.126
        plot_center = dcent[i]/(3.e5) * cent_lambda + cent_lambda
        
        wave, spec, err = sa.plot_line(slayfile,radius,
                                       wavelength=plot_center,plot=False)
        spec -= np.mean(spec[0:10])
        dvelo = (wave - cent_lambda)/cent_lambda*3.e5
        if flip:
            dvelo *= -1.0

        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        ax0.errorbar(dvelo,spec/np.max(spec),yerr=err/np.max(spec))
        ax0.plot(V,I/np.max(I))
        #ax0.plot(negV+mm1-nm1,negI/np.max(negI))
        ax0.set_xlabel('Velocity [km/s]')
        ax0.set_ylabel('Signal')
        ax0.set_title(datetime.now().isoformat(' ')+'\nr = {:4.3f} kpc'.
                      format(dradii[i]))

        tm1 = np.append(tm1,mm1)
        tm2 = np.append(tm2,mm2)
        tm3 = np.append(tm3,mm3)
        print "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}".format(
            mm1,dm1[i,1],mm2,dm2[i,1],mm3,dm3[i,1])

        pp.savefig(fig0)
        ax0.change_geometry(4,5,i+1)
        gridfig.axes.append(ax0)
#        fig0.clf()

    pp.close()

    pp1 = PDF(plotprefix+'_moments.pdf')

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(dradii,tm1,'.',label='model')
    ax1.plot(dradii,dcent,'x',label='data centers')
    ax1.plot(dradii,dm1[:,1],'x',label='data $\mu_1$')
    ax1.set_xlabel('Radius [kpc]')
    ax1.set_ylabel('$\mu_1$')
    ax1.set_ylim(-260,260)
    ax1.legend(loc=0)
    ax1.set_title(datetime.now().isoformat(' '))
    pp1.savefig(fig1)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(dradii,tm2,'.',label='model')
    ax2.plot(dradii,dm2[:,1],'x',label='data')
    ax2.set_xlabel('Radius [kpc]')
    ax2.set_ylabel('$\mu_2$')
    ax2.set_ylim(0,2000)
    ax2.legend(loc=0)
    pp1.savefig(fig2)

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(dradii,tm3,'.',label='model')
    ax3.plot(dradii,dm3[:,1],'x',label='data')
    ax3.set_xlabel('Radius [kpc]')
    ax3.set_ylabel('$\mu_3$')
    ax3.set_ylim(-5,5)
    ax3.legend(loc=0)
    pp1.savefig(fig3)
    
    pp1.close()

    ppg.savefig(gridfig)
    ppg.close()
    
    plt.clf()
#    plt.close('all')

    return 


