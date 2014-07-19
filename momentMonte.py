import numpy as np
import ADEUtils as ADE
import matplotlib.pyplot as plt
import bottleneck as bn

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

    x = np.arange(1000.)*1.
    
    func = pdf_mvsk(moment_list[1:])
    g = func(x)
    
    return x, moment_list[0]*g/np.max(g)

def get_noise(line, SNR):

    signal = np.sqrt(np.sum(line**2))
    noise = 2.*np.random.ranf(line.size) - 1.
    noise *= signal/(SNR*np.sqrt(np.sum(noise**2)))

    return noise

def do_a_line(moment_list,N):

    x, l = make_line(moment_list)

    SNRs = np.linspace(1,20,50)
    results = np.empty((SNRs.size,4,2))
    for i, SNR in enumerate(SNRs):
        sn_res = np.empty((N,4))
        for j in range(N):
            noise = get_noise(l,SNR)
            sn_res[j] = ADE.ADE_moments(x,l + noise)
            
        measured_vals = bn.nanmean(sn_res,axis=0)
        measured_stds = bn.nanstd(sn_res,axis=0)
        # print sn_res
        # print measured_stds
        # raw_input('')
        results[i,:,0] = measured_vals
        results[i,:,1] = measured_stds
        
    return SNRs, results

    
def plot_results(SNRs, results, moment_list):

    for i in range(4):

        ax = plt.figure().add_subplot(111)
        y = results[:,i,0]/moment_list[i+1]
        err = results[:,i,1]/moment_list[i+1]
        ax.errorbar(SNRs,y,yerr=err)
        ax.set_xlabel('SNR')
        ax.set_ylabel('$\mu_{{{0},meas}}/\mu_{{{0},true}}$'.format(i+1))
        ax.figure.show()

    return
