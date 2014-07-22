import numpy as np
import ADEUtils as ADE
import matplotlib.pyplot as plt
import bottleneck as bn
import time
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

def make_line(moment_list):

    resolution = 12.06 #Pixel resolution of our data, in km/s
    x = np.arange(0,500.,12.06) - 250.
    
    func = pdf_mvsk(moment_list[1:])
    g = func(x)
    
    return x, moment_list[0]*g/np.max(g)

def get_noise(x, line, SNR):

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
    pp = PDF(output)

    for i, window in enumerate(windows):
        win_res = np.empty((N,4))
        for j in range(N):
            noise = get_noise(x, l,SNR)
            ln = l + noise
            cdf = np.cumsum(ln/np.sum(ln))
            low, high = np.interp([1-window,window],cdf,x)
            idx = np.where((x > low) & (x <= high))
            win_res[j] = ADE.ADE_moments(x[idx],ln[idx])

        measured_vals = bn.nanmean(win_res,axis=0)
        measured_stds = bn.nanstd(win_res,axis=0)
        # print win_res
        # print measured_vals
        # print measured_stds
        # raw_input('')
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds    

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
        ax2.plot(windows,results[:,i,0]/moment_list[i+1],label='$\mu_{{{}}}$'.format(i+1))
        
    ax.legend(loc=0)
    ax2.legend(loc=0)
#    ax.figure.show()
    return windows, results, fig
