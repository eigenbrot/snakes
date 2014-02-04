import numpy as np
import ADEUtils as ADE
import ADESALT as sa
import Salty2 as salty
import matplotlib
import glob
import time
import pyfits
plt = matplotlib.pyplot
glob = glob.glob


def moments_notice(slayfile, simfile, plotprefix=False,
                   flip=False, cent_lambda = 5048.126):
    '''Take an extracted spectrum (slayfile) and model galaxy (simfile) and
    compute the first 3 statistical moments for both at all radii sampled by
    the spectrum. The model is binned and degraded to match the data quality
    of the actual data. The return values are the radii used, and then a tuple
    for each of the first three moments containing the moment and error on the
    data and the moment of the model. The moments are computed in the same
    window for both model and data so that we can compare the two as closely
    as possible. The exception is the first moment.
    '''

    radii, _, _ = sa.openslay(slayfile,flip=flip,moments=False)

    big_dm1 = np.array([])
    big_dm1e = np.array([])
    big_dm2 = np.array([])
    big_dm2e = np.array([])
    big_dm3 = np.array([])
    big_dm3e = np.array([])
    big_mm1 = np.array([])
    big_mm2 = np.array([])
    big_mm3 = np.array([])
    plot_radii = np.array([])
    if plotprefix:
        pp = PDF(plotprefix+'_lines.pdf')

    for radius in radii:

        dV, dI, derr, rwidth = sa.plot_line(slayfile,radius,
                                            wavelength=cent_lambda,velo=True,
                                            plot=False,baseline=1,flip=flip)

        mV, mI, _ = salty.line_profile(simfile,radius,plot=False,Iwidth=17,
                                       width=rwidth,observe=True,fit=False) 

        conv_dI = dI/np.mean(dI)
        conv_mI = mI/np.mean(mI)
        mI_pad = pad(mI,dI.size)
        corr = np.convolve(conv_dI,mI_pad[::-1],'same')
        idx = np.where(corr == corr.max())[0][0]
        corr_peak = dV[idx]
        
        fit_mV = dV + corr_peak
        cdf = np.cumsum(mI_pad/np.sum(mI_pad))
        
        # Limits defined by model data
        lowV = np.interp(0.1,cdf,fit_mV)
        highV = np.interp(0.9,cdf,fit_mV)
        dmoment_idx = np.where((dV >= lowV) & (dV <= highV))
        mmoment_idx = np.where((fit_mV >= lowV) & (fit_mV <= highV))
        print dmoment_idx

        if dmoment_idx[0].size <= 2:
            print "skipping {} kpc".format(radius)
            continue
        
        dmoments, dmerrs = ADE.ADE_moments(dV[dmoment_idx],dI[dmoment_idx],
                                         err=derr[dmoment_idx])
        mmoments = ADE.ADE_moments(fit_mV[mmoment_idx],mI_pad[mmoment_idx])

        '''We need to compute the peak of the model line separately because we
        shifted it to find the window used for the higher order moments'''
        (mpeak, _, _, _) = ADE.ADE_moments(mV,mI)

        print "moments: {}\nerrs: {}\nmmoments: {}".\
            format(dmoments,dmerrs,mmoments)
        
        big_dm1 = np.append(big_dm1,dmoments[0])
        big_mm1 = np.append(big_mm1,mpeak)
        big_dm1e = np.append(big_dm1e,dmerrs[0])

        big_dm2 = np.append(big_dm2,dmoments[1])
        big_mm2 = np.append(big_mm2,mmoments[1])
        big_dm2e = np.append(big_dm2e,dmerrs[1])

        big_dm3 = np.append(big_dm3,dmoments[2])
        big_mm3 = np.append(big_mm3,mmoments[2])
        big_dm3e = np.append(big_dm3e,dmerrs[2])

        plot_radii = np.append(plot_radii,radius)

    return plot_radii, (big_dm1, big_dm1e, big_mm1),\
        (big_dm2, big_dm2e, big_mm2),\
        (big_dm3, big_dm3e, big_mm3)

def giant_steps(slay_file, simstring, parameter_list):
    '''Take an extracted spectrum (slayfile) and a list of galaxy models and
    compute the goodness of fit for each model as a function of two model
    parameters specified in parameter_list. The goodness of fit is defined as
    the sum of the squares of the reduced chi-squared values for each of the
    first 3 statistical moments. This function is only useful if the modles
    defined by simstring vary the parameters listed in parameter_list.
    '''

    simlist = glob(simstring)
    par1_arr = np.array([]) 
    par2_arr = np.array([]) 
    value_arr = np.array([])

    for sim_file in simlist:
        header = pyfits.open(sim_file)[0].header
        par1 = float(header[parameter_list[0]])
        par2 = float(header[parameter_list[1]])
        print 'simfile: {}\n{}: {}\n{}: {}'.format(sim_file,
                                                   parameter_list[0],
                                                   par1,
                                                   parameter_list[1],
                                                   par2)

        radii, m1, m2, m3 = moments_notice(slay_file,sim_file)
        chis = np.array([])
        for moment in [m1, m2, m3]:
            red_chi = np.sum((moment[0] - moment[2])**2/moment[1])/\
                (radii.size - 1)
            chis = np.append(chis,red_chi)

        value = np.sum(chis**2)
        print 'fit goodness: {}'.format(value)
        par1_arr = np.append(par1_arr,par1)
        par2_arr = np.append(par2_arr,par2)
        value_arr = np.append(value_arr,value)
    
    return par1_arr, par2_arr, value_arr

def make_boring(vr_list, h_rot_list, z=0, size=1001):
    '''Given a list of values for v_r and h_rot, make a grid of galaxy models
    with all possible combinations of those two parameters.
    '''
    basename = 'sim_z{:n}_boring'.format(z)
    for v_r in vr_list:
        for h_rot in h_rot_list:
            name = '{:}_{:03n}_{:03n}.fits'.format(basename,v_r,h_rot*100)
            print 'building model {}:\nv_r = {} km/s\nh_rot = {} kpc'.format(
                name,v_r,h_rot)
            salty.simcurve(size,z,v_r,h_rot,output=name,scale=0.0999)

    return
    
def plot_ice(results,show=True):
    '''Given the results of moments_notice(), plot model vs. data for all
    three moments. Show the plots and return them.
    '''
    ax1 = plt.figure().add_subplot(111)
    ax2 = plt.figure().add_subplot(111)
    ax3 = plt.figure().add_subplot(111)
    ax1.set_ylabel('$\mu_1$')
    ax2.set_ylabel('$\mu_2$')
    ax3.set_ylabel('$\mu_3$')

    radius = results[0]

    for ax, data in zip([ax1, ax2, ax3],[results[1],results[2],results[3]]):
        ax.set_xlabel('radius [kpc]')
        ax.set_title(time.asctime())
        ax.errorbar(radius,data[0],yerr=data[1],fmt='.',ms=15)
        ax.plot(radius,data[2])
        if show:
            ax.get_figure().show()

    return ax1, ax2, ax3

def plot_chi(results):
    '''Given the output of giant_steps, plot a countour map of the goodness of
    fit parameter as a function of parameter 1 and parameter 2. Show this
    countour plot and return it.
    '''
    shape = (int(results[0].size**0.5),int(results[0].size**0.5))
    print shape
    P1 = results[0].reshape(shape)
    P2 = results[1].reshape(shape)
    V = results[2].reshape(shape)

    ax = plt.figure().add_subplot(111)
    ax.contourf(P1,P2,V,100,cmap=plt.cm.gray)
    ax.get_figure().show()
    
    return ax

def pad(arr,length):
    '''Take an array and pad it with zeros to the length specified. The
    original array will be in the center of the padded array.
    '''
    preroll = np.append(arr,np.zeros(length - arr.size))
    return np.roll(preroll, (length - arr.size)/2)
