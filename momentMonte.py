import numpy as np
import ADEUtils as ADE
import matplotlib.pyplot as plt
import Salty2 as salty
import bottleneck as bn
import time
from scipy import ndimage
from matplotlib.backends.backend_pdf import PdfPages as PDF

def hermnorm(N):
    # return the negatively normalized hermite polynomials up to order N-1
    #  (inclusive)
    #  using the recursive relationship
    #  p_n+1 = p_n(x)' - x*p_n(x)
    #   and p_0(x) = 1
    plist = [None]*N
    plist[0] = np.poly1d(1)
    for n in range(1,N):
        plist[n] = plist[n-1].deriv() - np.poly1d([1,0])*plist[n-1]
    return plist

def pdf_mvsk(mvsk):
    """Return the Gaussian expanded pdf function given the list of 1st, 2nd
    moment and skew and Fisher (excess) kurtosis.



    Parameters
    ----------
    mvsk : list of mu, mc2, skew, kurt
        distribution is matched to these four moments

    Returns
    -------
    pdffunc : function
        function that evaluates the pdf(x), where x is the non-standardized
        random variable.


    Notes
    -----

    Changed so it works only if four arguments are given. Uses explicit
    formula, not loop.

    This implements a Gram-Charlier expansion of the normal distribution
    where the first 2 moments coincide with those of the normal distribution
    but skew and kurtosis can deviate from it.

    In the Gram-Charlier distribution it is possible that the density
    becomes negative. This is the case when the deviation from the
    normal distribution is too large.



    References
    ----------
    http://en.wikipedia.org/wiki/Edgeworth_series
    Johnson N.L., S. Kotz, N. Balakrishnan: Continuous Univariate
    Distributions, Volume 1, 2nd ed., p.30
    """
    N = len(mvsk)
    if N < 4:
        raise ValueError("Four moments must be given to "
                         "approximate the pdf.")

    mu, mc2, skew, kurt = mvsk

    totp = np.poly1d(1)
    sig = np.sqrt(mc2)
    if N > 2:
        Dvals = hermnorm(N+1)
        C3 = skew/6.0
        C4 = kurt/24.0
        # Note: Hermite polynomial for order 3 in hermnorm is negative
        # instead of positive
        totp = totp -  C3*Dvals[3] +  C4*Dvals[4]

    def pdffunc(x):
        xn = (x-mu)/sig
        return totp(xn)*np.exp(-xn*xn/2.0)/np.sqrt(2*np.pi)/sig
    return pdffunc

def get_line(simfile, radius, observe=True):

    rwidth = 10.
    mV, mI, _ = salty.line_profile(simfile,radius,plot=False,Iwidth=17,
                                   width=rwidth,observe=observe,fit=False,
                                   verbose=False,nofits=False)

    return mV, mI

def make_line(moment_list):

    resolution = 12.06 #Pixel resolution of our data, in km/s
    x = np.arange(0,500.,12.06) - 250.
    
    func = pdf_mvsk(moment_list[1:])
    g = func(x)
    
    return x, moment_list[0]*g/np.max(g)

def get_noise(x, line, SNR):

    if SNR == np.inf:
        return np.zeros(x.shape)
    cdf = np.cumsum(line/np.sum(line))
    low, high = np.interp([0.16,0.84],cdf,x)
    idx = np.where((x > low) & (x <= high))
    sig = np.sum(line[idx]**2)
    noise = (2.*np.random.ranf(line.size) - 1.)
    noise *= np.sqrt(sig/SNR**2)
    return noise

def do_a_line(moment_list,N,line_output,monte_output):

    x, l = make_line(moment_list)
    SNRs = np.linspace(5,100,50)
    results = np.empty((SNRs.size,4,2))
    lp = PDF(line_output)

    for i, SNR in enumerate(SNRs):
        sn_res = np.empty((N,4))
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            cdf = np.cumsum(ln/np.sum(ln))
            low, high = np.interp([0.01,0.99],cdf,x)
            idx = np.where((x > low) & (x <= high))
            sn_res[j] = ADE.ADE_moments(x[idx],ln[idx])

        measured_vals = bn.nanmean(sn_res,axis=0)
        measured_stds = bn.nanstd(sn_res,axis=0)
        # print sn_res
        # print measured_stds
        # raw_input('')
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds
        
        ax = plt.figure().add_subplot(111)
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Flux')
        ax.set_title('SNR = {:5.2f}'.format(SNR))
        ax.plot(x,ln)
        lp.savefig(ax.figure)

    lp.close()
    
    mp = PDF(monte_output)
    plots = plot_results(SNRs, results, moment_list)
    for i, plot in enumerate(plots):
        if i == 2:
            plot.set_ylim(-2,2)
        if i == 3:
            plot.set_ylim(-2,5)
        mp.savefig(plot.figure)
    mp.close()
    plt.close('all')
    return SNRs, results

    
def plot_results(SNRs, results, moment_list):
    
    output = []

    for i in range(4):

        ax = plt.figure().add_subplot(111)
        y = results[:,i,0]/moment_list[i+1]
        err = results[:,i,1]/moment_list[i+1]
        ax.errorbar(SNRs,y,yerr=err)
        ax.set_xlabel('SNR')
        ax.set_ylabel('$\mu_{{{0},meas}}/\mu_{{{0},true}}$'.format(i+1))
        ax.text(0.85,0.80,
                '$\mu_1 = {}$\n$\mu_2 = {}$\n$\mu_3 = {}$\n$\mu_4 = {}$'.\
                    format(*moment_list[1:]),
                transform=ax.transAxes,fontdict={'size':14})
        output.append(ax)

    return output

def test_window(moment_list, N, SNR, output):

    x, l = make_line(moment_list)
    windows = np.linspace(0.8,0.99,50)
    results = np.empty((windows.size,4,2))
    if output: pp = PDF(output)
    ax0 = plt.figure().add_subplot(111)
    ax0.set_title('{}\nSNR = {}'.format(time.asctime(),SNR))
    ax0.set_xlabel('Velocity [km/s]')
    ax0.set_ylabel('Flux')

    for i, window in enumerate(windows):
        win_res = np.empty((N,4))
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            cdf = np.cumsum(ln/np.sum(ln))
            low, high = np.interp([1-window,window],cdf,x)
            idx = np.where((x > low) & (x <= high))
            win_res[j] = ADE.ADE_moments(x[idx],ln[idx])

        if i == 0:
            ax0.plot(x,ln)
        if i % 5 == 0 or i == windows.size - 1:
            line = ax0.plot(x,l)[0]
            ax0.axvline(x=low,label='{:5.3f}'.format(window),color=line.get_color())
            ax0.axvline(x=high,color=line.get_color())
            del ax0.lines[-3]
        measured_vals = bn.nanmean(win_res,axis=0)
        measured_stds = bn.nanstd(win_res,axis=0)
        # print win_res
        # print measured_vals
        # print measured_stds
        # raw_input('')
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds    
        
    ax0.legend(loc=0,frameon=False,fontsize=10,title='Window')
    fig = plt.figure(figsize=(8,10))
    fig.suptitle('{}\nSNR = {}'.format(time.asctime(),SNR))
    ax = fig.add_subplot(211)
    ax.set_xlabel('Window (1-X/X)')
    ax.set_ylabel('SNR')
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('Window (1-X/X)')
    ax2.set_ylabel('$\mu_{i,meas}/\mu_{i,true}$')

    for i in range(4):
        sn =  np.sqrt((results[:,i,0]/results[:,i,1])**2)
        ax.plot(windows,sn,label='$\mu_{{{}}}$'.format(i+1))
        print i
        if i == 1:
            print 'Yaya!'
            ax2.plot(windows,np.sqrt(results[:,i,0]/moment_list[i+1]),label='$\sqrt{\mu_{{{}}}}$'.format(i+1))
        else:
            ax2.plot(windows,results[:,i,0]/moment_list[i+1],label='$\mu_{{{}}}$'.format(i+1))
        
    ax.legend(loc=0)
    ax2.legend(loc=0)
    if output:
        pp.savefig(fig)
        pp.savefig(ax0.figure)
        pp.close()
    plt.close('all')
    return windows, results, [fig, ax0.figure()]

def test_window_sim(simfile, radius, N, SNR, output, observe=True):

    x, l = get_line(simfile, radius, observe=observe)
#    windows = np.linspace(0.8,0.999,50)
    windows = np.linspace(20,400,50.)
    digi_windows = np.empty(windows.shape)
    results = np.empty((windows.size,4,2))
    if output: pp = PDF(output)
    ax0 = plt.figure().add_subplot(111)
    ax0.set_title('{}\nSNR={}'.format(time.asctime(),SNR))
    ax0.set_xlabel('Velocity [km/s]')
    ax0.set_ylabel('Flux')

    true_moments = ADE.ADE_moments(x,l)

    print '{:>4}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'.\
        format('idx','window','low','high','idxsize','centerstd','mean','std','SN')
    for i, window in enumerate(windows):
        win_res = np.empty((N,4))
        window_N = np.empty((N,))
        # peak, _, _, _ = ADE.ADE_moments(x,l)
        # # cdf = np.cumsum(ln/np.sum(ln))
        # # low, high = np.interp([1-window,window],cdf,x)
        # low = peak - window/2.
        # high = peak + window/2.
        # idx = np.where((x > low) & (x <= high))[0]
        centers = np.array([])
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            peak, _, _, _ = ADE.ADE_moments(x,ln)
            centers = np.append(centers,peak)
            cdf = np.cumsum(ln/np.sum(ln))
            ##low, high = np.interp([1-window,window],cdf,x)
            low = peak - window/2.
            high = peak + window/2.
            idx = np.where((x > low) & (x <= high))[0]
            window_N[j] = x[idx[-1]] - x[idx[0]]
            if j == 0:
                print '{:4}{:10.3f}{:10.3f}{:10.3f}{:10n}'.\
                    format(i,window_N[j]/np.sqrt(true_moments[1]),low,high,idx.size),
            win_res[j] = ADE.ADE_moments(x[idx],ln[idx])
            if i == 7:
                line = ax0.plot(x,l)[0]
                ax0.axvline(x=x[idx[0]],color=line.get_color())
                ax0.axvline(x=x[idx[-1]],color=line.get_color())
                ax0.axvline(x=low,ls=':',lw=0.4,color=line.get_color())
                ax0.axvline(x=high,ls=':',lw=0.2,color=line.get_color())
                del ax0.lines[-5]

        if i == 0:
            ax0.plot(x,ln)
        #if i % 10 == 0 or i == windows.size - 1:
        if i in [88,109]:
            line = ax0.plot(x,l)[0]
            ax0.axvline(x=x[idx[0]],label='{:5.3f}'.format(window/np.sqrt(true_moments[1])),color=line.get_color())
            ax0.axvline(x=x[idx[-1]],color=line.get_color())
            ax0.axvline(x=low,ls=':',color=line.get_color())
            ax0.axvline(x=high,ls=':',color=line.get_color())
            del ax0.lines[-5]
 
        print '{:9.5f}'.format(np.std(centers)),
        measured_vals = bn.nanmean(win_res,axis=0)
        measured_stds = bn.nanstd(win_res,axis=0)
        print '{:9.3f}{:10.3f}{:10.3f}'.format(measured_vals[0],measured_stds[0],measured_vals[0]/measured_stds[0])
        digi_windows[i] = bn.nanmean(window_N)
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds    
        
    ax0.legend(loc=0,frameon=False,fontsize=10,title='Window/$\sqrt{\mu_2}$')
    fig = plt.figure(figsize=(8,10))
    fig.suptitle('{}\nSNR={}'.format(time.asctime(),SNR))
    ax = fig.add_subplot(211)
    ax.set_xlabel('Window width/$\sqrt{\mu_{2,true}}$')
    ax.set_ylabel('SNR')
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('Window width/$\sqrt{\mu_{2,true}}$')
    ax2.set_ylabel('$\mu_{i,meas}/\mu_{i,true}$')

    for i in range(4):
        sn =  np.sqrt((results[:,i,0]/results[:,i,1])**2)
#        ax.plot(digi_windows/np.sqrt(true_moments[1]),sn,label='$\mu_{{{}}}$'.format(i+1))
        ax.plot(windows/np.sqrt(true_moments[1]),sn,':')
        if i == 7:
            print 'tat'
            ax2.plot(windows/np.sqrt(true_moments[1]),np.sqrt(results[:,i,0]/true_moments[i]),label='$\sqrt{{\mu_{{{}}}}}$'.format(i+1))  
        else:
            ax2.plot(windows/np.sqrt(true_moments[1]),results[:,i,0]/true_moments[i],label='$\mu_{{{}}}$'.format(i+1))
        
    ax2.legend(loc=0)
    ax2.axhline(y=1,ls=':')
    ax2.set_ylim(-2,1.5)
    if output:
        pp.savefig(fig)
        pp.savefig(ax0.figure)
        pp.close()
    plt.close('all')
    return windows, results, [fig, ax0.figure]

def many_window(moment_list,SNRlist,output):

    pp = PDF(output)
    N = 1000
    for SNR in SNRlist:
        r = test_window(moment_list, N, SNR, False)
        pp.savefig(r[-1][0])
        pp.savefig(r[-1][1])
    
    pp.close()
    return

def test_window_int(simfile, radius, N, SNR, output, observe=True, smooth=1, resamp=False):

#    x, l = get_line(simfile, radius, observe=observe)
    
    x,g = ADE.ADE_gauss(2000,1100,120*smooth)
    x,g2 = ADE.ADE_gauss(2000,1100-80*smooth,40*smooth)
    x = (x-600)*0.25
    l = g + g2

    x = x[::48]
    l = l[::48]
    
    # if smooth:
    #     l = ndimage.gaussian_filter1d(l,smooth)
    if resamp:
        oldx = x.copy()
        resamp_x = np.linspace(x.min(),x.max(),1000.) 
        widths = np.arange(3,resamp_x.size*0.666,dtype=np.int)
    else:
        widths = np.arange(3,x.size*0.666,dtype=np.int)
    results = np.empty((widths.size,4,2))
    if output: pp = PDF(output)
    fig = plt.figure(figsize = (10,8))
    ax0 = fig.add_subplot(221)
    ax0.set_xlabel('Velocity [km/s]')
    ax0.set_ylabel('Flux')

    true_moments = ADE.ADE_moments(x,l)
    for i in range(4):
        true_moments[i] = true_moments[i]**(1./(i+1))
    print 'True moments: {}'.format(true_moments)
    print "True \\sqrt{{\\mu_2}} = {:5.2f}".format(true_moments[1])
    print "Orig px resolution = {:5.2f} km/s +/- {:5.2f}".\
        format(np.mean(np.diff(x)),np.std(np.diff(x)))
    print "Interpolated px resolution = {:5.2f} km/s +/- {:5.2f}\n".\
        format(np.mean(np.diff(resamp_x)),np.std(np.diff(resamp_x)))

    print '{:>4}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'.\
        format('idx','window','centidx','idxsize','centerstd','mean','std','SN')
    for i, width in enumerate(widths):
        win_res = np.empty((N,4))
        # peak, _, _, _ = ADE.ADE_moments(x,l)
        # cent_idx = np.argmin(np.abs(x - peak))
        # idx = np.arange(cent_idx - width/2., cent_idx + width/2., dtype=np.int)
        centers = np.array([])
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            if observe and resamp:
                ln = np.interp(resamp_x,x,ln)
                x = resamp_x
            peak, _, _, _ = ADE.ADE_moments(x,ln)
            cent_idx = np.argmin(np.abs(x - peak))
            centers = np.append(centers,cent_idx)
            idx = np.arange(cent_idx - width/2., cent_idx + width/2., dtype=np.int)
            if np.any(idx > x.size):
                idx = np.arange(idx[0],x.size - 1)
            if j == 0:
                print '{:4}{:10.3f}{:10n}{:10n}'.\
                    format(i,width,cent_idx,idx.size),
            moments = ADE.ADE_moments(x[idx],ln[idx])
            for k in range(4):
                moments[k] = moments[k]**(1./(1+k))
            win_res[j] = moments
            if resamp:
                x = oldx
            # if i == 7:
            #     line = ax0.plot(x,l)[0]
            #     ax0.axvline(x=x[idx[0]],color=line.get_color())
            #     ax0.axvline(x=x[idx[-1]],color=line.get_color())
            #     ax0.axvline(x=low,ls=':',lw=0.4,color=line.get_color())
            #     ax0.axvline(x=high,ls=':',lw=0.2,color=line.get_color())
            #     del ax0.lines[-5]

        if i == 0:
            if resamp:
                ax0.plot(resamp_x,ln)
            else:
                ax0.plot(x,ln)
        #if i % 10 == 0 or i == windows.size - 1:
        if i in [-8,-10]:
            line = ax0.plot(x,l)[0]
            ax0.axvline(x=x[idx[0]],label='{:5.3f}'.format(window/np.sqrt(true_moments[1])),color=line.get_color())
            ax0.axvline(x=x[idx[-1]],color=line.get_color())
            ax0.axvline(x=low,ls=':',color=line.get_color())
            ax0.axvline(x=high,ls=':',color=line.get_color())
            del ax0.lines[-5]
 
        print '{:9.5f}'.format(np.std(centers)),
        measured_vals = bn.nanmean(win_res,axis=0)
        measured_stds = bn.nanstd(win_res,axis=0)
        print '{:9.3f}{:10.3f}{:10.3f}'.format(measured_vals[2],measured_stds[2],measured_vals[2]/measured_stds[2])
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds    
        
    ax0.legend(loc=0,frameon=False,fontsize=10,title='Window/$\sqrt{\mu_2}$')
#    fig.subplots_adjust(hspace=0.0001)
    fig.suptitle('{}\nSNR={}, $\zeta_{{2,true}}$={:4.2f} km/s'.format(
            time.asctime(),SNR,true_moments[1]))
    ax = fig.add_subplot(222)
    ax.set_ylabel('SNR')
    ax.set_yscale('log')
    ax.set_xlabel('Window width [px]')
#    plt.setp(ax.get_xticklabels(),visible=False)
    axkm = ax.twiny()
    axkm.set_xlabel('Window width [km/s]')
    ax2 = fig.add_subplot(223)
    ax2.set_ylabel('$\zeta_{i,meas}/\zeta_{i,true}$')
    ax2.set_xlabel('Window width [km/s]')
#    plt.setp(ax2.get_xticklabels(),visible=False)
    ax3 = fig.add_subplot(224)
    ax3km = ax3.twiny()
    ax3km.set_xlabel('Window width/$\zeta_{{2,true}}$')
#    ax3.set_yscale('log')
    ax3.set_ylim(0,1.1)
    ax3.set_xlabel('Window width [px]')
    ax3.set_ylabel('$SNR \\times\, \zeta_{i,meas}/\zeta_{i,true}$ (normalized)')

    for i in range(4):
        sn =  np.sqrt((results[:,i,0]/results[:,i,1])**2)
        ax.plot(widths,sn,label='$\zeta_{{{}}}$'.format(i+1))
        if resamp:
            axkm.plot(widths*np.mean(np.diff(resamp_x)),sn)
        else:
            axkm.plot(widths*np.mean(np.diff(x)),sn)
        del axkm.lines[-1]
        # if i == 7:
        #     print 'tat'
        #     ax2.plot(windows/np.sqrt(true_moments[1]),np.sqrt(results[:,i,0]/true_moments[i]),label='$\sqrt{{\mu_{{{}}}}}$'.format(i+1))  
#        ax2.plot(widths,results[:,i,0]/true_moments[i],label='$\mu_{{{}}}$'.format(i+1))
        ax2.plot(widths,results[:,i,0]/true_moments[i],label='$\zeta_{{{}}}$'.format(i+1))
        zeta = results[:,i,0]/true_moments[i]*sn
        zeta /= zeta[50:].max()
        ax3.plot(widths,zeta)
        if resamp:
            ax3km.plot(widths*np.mean(np.diff(resamp_x))/true_moments[1],zeta)
        else:
            ax3km.plot(widths*np.mean(np.diff(x))/np.sqrt(true_moments[1]),zeta)            
        del ax3km.lines[-1]

    ax.legend(loc=0)
#    ax.set_ylim(0,1000)
    ax2.axhline(y=1,ls=':')
    ax2.set_ylim(-0.1,1.1)
    if output:
        plt.tight_layout()
        pp.savefig(fig)
        pp.close()
    plt.close('all')
    return widths, results, [fig, ax0.figure], true_moments

def test_window_line(x, l, N, SNR, output, observe=True, smooth=1, resamp=False):
    
    # if smooth:
    #     l = ndimage.gaussian_filter1d(l,smooth)
    if resamp:
        oldx = x.copy()
        resamp_x = np.linspace(x.min(),x.max(),1000.) 
        widths = np.arange(3,resamp_x.size*0.666,dtype=np.int)
    else:
        widths = np.arange(3,x.size*0.666,dtype=np.int)
    results = np.empty((widths.size,4,2))
    if output: pp = PDF(output)
    fig = plt.figure(figsize = (10,8))
    ax0 = fig.add_subplot(221)
    ax0.set_xlabel('Velocity [km/s]')
    ax0.set_ylabel('Flux')

    true_moments = ADE.ADE_moments(x,l)
    for i in range(4):
        true_moments[i] = np.sign(true_moments[i])*((np.sign(true_moments[i])*true_moments[i])**(1./(i+1)))
    ax0.text(0.7,0.8,'$\zeta_1=$ {:4.2e}\n$\zeta_2=$ {:4.2e}\n$\zeta_3=$ {:4.2e}\n$\zeta_4=$ {:4.2e}'.format(*true_moments),ha='left',va='top',fontsize=8,transform=ax0.transAxes)
    print 'True moments: {}'.format(true_moments)
    print "True \\sqrt{{\\mu_2}} = {:5.2f}".format(true_moments[1])
    print "Orig px resolution = {:5.2f} km/s +/- {:5.2f}".\
        format(np.mean(np.diff(x)),np.std(np.diff(x)))
    print "Interpolated px resolution = {:5.2f} km/s +/- {:5.2f}\n".\
        format(np.mean(np.diff(resamp_x)),np.std(np.diff(resamp_x)))

    print '{:>4}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}'.\
        format('idx','window','centidx','idxsize','centerstd','mean','std','SN')
    for i, width in enumerate(widths):
        win_res = np.empty((N,4))
        # peak, _, _, _ = ADE.ADE_moments(x,l)
        # cent_idx = np.argmin(np.abs(x - peak))
        # idx = np.arange(cent_idx - width/2., cent_idx + width/2., dtype=np.int)
        centers = np.array([])
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            if observe and resamp:
                ln = np.interp(resamp_x,x,ln)
                x = resamp_x
            peak, _, _, _ = ADE.ADE_moments(x,ln)
            cent_idx = np.argmin(np.abs(x - peak))
            centers = np.append(centers,cent_idx)
            idx = np.arange(cent_idx - width/2., cent_idx + width/2., dtype=np.int)
            if np.any(idx > x.size):
                idx = np.arange(idx[0],x.size - 1)
            if j == 0:
                print '{:4}{:10.3f}{:10n}{:10n}'.\
                    format(i,width,cent_idx,idx.size),
            moments = ADE.ADE_moments(x[idx],ln[idx])
            for k in range(4):
                moments[k] = np.sign(moments[k])*((np.sign(moments[k])*moments[k])**(1./(1+k)))
            win_res[j] = moments
            if resamp:
                x = oldx
            # if i == 7:
            #     line = ax0.plot(x,l)[0]
            #     ax0.axvline(x=x[idx[0]],color=line.get_color())
            #     ax0.axvline(x=x[idx[-1]],color=line.get_color())
            #     ax0.axvline(x=low,ls=':',lw=0.4,color=line.get_color())
            #     ax0.axvline(x=high,ls=':',lw=0.2,color=line.get_color())
            #     del ax0.lines[-5]

        if i == 0:
            if resamp:
                ax0.plot(resamp_x,ln)
            else:
                ax0.plot(x,ln)
        #if i % 10 == 0 or i == windows.size - 1:
        # if i in [-8,-10]:
        #     line = ax0.plot(x,l)[0]
        #     ax0.axvline(x=x[idx[0]],label='{:5.3f}'.format(window/np.sqrt(true_moments[1])),color=line.get_color())
        #     ax0.axvline(x=x[idx[-1]],color=line.get_color())
        #     ax0.axvline(x=low,ls=':',color=line.get_color())
        #     ax0.axvline(x=high,ls=':',color=line.get_color())
        #     del ax0.lines[-5]
 
        print '{:9.5f}'.format(np.std(centers)),
        measured_vals = bn.nanmean(win_res,axis=0)
        measured_stds = bn.nanstd(win_res,axis=0)
        print '{:9.3f}{:10.3f}{:10.3f}'.format(measured_vals[2],measured_stds[2],measured_vals[2]/measured_stds[2])
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds    
        
    ax0.legend(loc=0,frameon=False,fontsize=10,title='Window/$\sqrt{\mu_2}$')
#    fig.subplots_adjust(hspace=0.0001)
    fig.suptitle('{}\nSNR={}, $\zeta_{{2,true}}$={:4.2f} km/s'.format(
            time.asctime(),SNR,true_moments[1]))
    ax = fig.add_subplot(222)
    ax.set_ylabel('1/Noise')
    ax.set_yscale('log')
    ax.set_xlabel('Window width [px]')
#    plt.setp(ax.get_xticklabels(),visible=False)
    axkm = ax.twiny()
    axkm.set_xlabel('Window width [km/s]')
    ax2 = fig.add_subplot(223)
    ax2.set_ylabel('$\zeta_{i,meas}/\zeta_{i,true}$')
    ax2.set_xlabel('Window width [km/s]')
#    plt.setp(ax2.get_xticklabels(),visible=False)
    ax3 = fig.add_subplot(224)
    ax3km = ax3.twiny()
    ax3km.set_xlabel('Window width/$\zeta_{{2,true}}$')
#    ax3.set_yscale('log')
    ax3.set_ylim(0,1.1)
    ax3.set_xlabel('Window width [px]')
    ax3.set_ylabel('1/Noise $\\times\, \zeta_{i,meas}/\zeta_{i,true}$ (normalized)')

    for i in range(4):
        sn = 1./results[:,i,1]
        ax.plot(widths,sn/sn[100:].max(),label='$\zeta_{{{}}}$'.format(i+1))
        if resamp:
            axkm.plot(widths*np.mean(np.diff(resamp_x)),sn)
        else:
            axkm.plot(widths*np.mean(np.diff(x)),sn)
        del axkm.lines[-1]
        # if i == 7:
        #     print 'tat'
        #     ax2.plot(windows/np.sqrt(true_moments[1]),np.sqrt(results[:,i,0]/true_moments[i]),label='$\sqrt{{\mu_{{{}}}}}$'.format(i+1))  
#        ax2.plot(widths,results[:,i,0]/true_moments[i],label='$\mu_{{{}}}$'.format(i+1))
        ax2.plot(widths,results[:,i,0]/true_moments[i],label='$\zeta_{{{}}}$'.format(i+1))
        zeta = results[:,i,0]/true_moments[i]*sn
        zeta /= zeta[100:].max()
        ax3.plot(widths,zeta)
        if resamp:
            ax3km.plot(widths*np.mean(np.diff(resamp_x))/true_moments[1],zeta)
        else:
            ax3km.plot(widths*np.mean(np.diff(x))/np.sqrt(true_moments[1]),zeta)            
        del ax3km.lines[-1]

    ax.legend(loc=0)
#    ax.set_ylim(0,1000)
    ax2.axhline(y=1,ls=':')
    ax2.set_ylim(-0.1,1.1)
    if output:
        plt.tight_layout()
        pp.savefig(fig)
        pp.close()
    plt.close('all')
    return widths, results, [fig, ax0.figure], true_moments

def gogo_models():

    slist = [1,1.5,2,3]
    di = {}

    for s in slist:

        o = test_window_int('',27,1000,100,'m_s{:n}_linear.pdf'.format(s),observe=True,resamp=True,smooth=s)
        di[s] = o[-1]

    return di

def gogo_lines():
    #    x, l = get_line(simfile, radius, observe=observe)
    
    smooth = 2
    d = {}
    o = {}
    x,g = ADE.ADE_gauss(2000,1100,120*smooth)
    x,g2 = ADE.ADE_gauss(2000,1100-80*smooth,40*smooth)
    d[0] = g + g2

    x,g = ADE.ADE_gauss(2000,1000,100.*smooth)
    x,g2 = ADE.ADE_gauss(2000,1000. - 100*smooth,150.*smooth)
    x,g3 = ADE.ADE_gauss(2000,1000 + 100.*smooth,150.*smooth)
    x = (x-600)*0.25
    d[1] = g + 0.3*(g2 + g3)

    d[2] = g + g2

    x, g = ADE.ADE_gauss(2000,1000,70.*smooth)
    x, g2 = ADE.ADE_gauss(2000,1000 + 100.*smooth,250*smooth)
    d[3] = g + 0.3*g2

    x = (x-600)*0.25
    x = x[::48]
    for i in range(4):
        l = d[i]
        l = l[::48]
        o[i] = test_window_line(x, l, 1000, 100., 'l{}_s2.pdf'.format(i), observe=True, resamp=True)
        
    return o
