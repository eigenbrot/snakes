import numpy as np
import ADEUtils as ADE
import ADESALT as sa
import Salty2 as salty
import drunkData as dd
import matplotlib
import glob
import time
import pyfits
import bottleneck as bn
import scipy.optimize as spo
import os
plt = matplotlib.pyplot
glob = glob.glob


def moments_notice(drunkdata, simfile, plotprefix=False,nofits=False,
                   flip=False, cent_lambda = 5048.126, skip_radii = []):
    '''Take an extracted spectrum (slayfile) and model galaxy (simfile) and
    compute the first 3 statistical moments for both at all radii sampled by
    the spectrum. The model is binned and degraded to match the data quality
    of the actual data. The return values are the radii used, and then a tuple
    for each of the first three moments containing the moment and error on the
    data and the moment of the model. The moments are computed in the same
    window for both model and data so that we can compare the two as closely
    as possible. The exception is the first moment.
    '''

    radii, rwidths, vwidths, m1, m2, m3 = drunkdata

    big_dm1 = m1[0]
    big_dm1e = m1[1]
    big_dm2 = m2[0]
    big_dm2e = m2[1]
    big_dm3 = m3[0]
    big_dm3e = m3[1]
    big_mm1 = np.array([])
    big_mm2 = np.array([])
    big_mm3 = np.array([])
    plot_radii = np.array([])
    if plotprefix:
        pp = PDF(plotprefix+'_lines.pdf')

    for radius, rwidth, vwidth in zip(radii,rwidths,vwidths):
        
        if int(np.floor(radius)) in skip_radii:
            print "user skipping radius {} kpc".format(radius)
            continue

        mV, mI, _ = salty.line_profile(simfile,radius,plot=False,Iwidth=17,
                                       width=rwidth,observe=True,fit=False,
                                       verbose=False,nofits=nofits) 

        '''We need to compute the peak of the model line separately because we
        shifted it to find the window used for the higher order moments'''
        (mpeak, _, _) = ADE.ADE_moments(mV,mI)

        vwidth = 135.
        mpeak_V = mV[np.argmin(np.abs(mV - mpeak))]
        lowV = mpeak_V - vwidth/2.
        highV = mpeak_V + vwidth/2.
        # cdf = np.cumsum(mI/np.sum(mI))
        # lowV, highV = np.interp([0.005,0.995],cdf,mV)
        
        mmoment_idx = np.where((mV >= lowV) & (mV <= highV))
#        print vwidth, highV - lowV, mmoment_idx[0].size
        mmoments = ADE.ADE_moments(mV[mmoment_idx],mI[mmoment_idx])

        
#        print mmoments
        big_mm1 = np.append(big_mm1,mpeak)
        big_mm2 = np.append(big_mm2,np.sqrt(mmoments[1]))
        big_mm3 = np.append(big_mm3,mmoments[2])

        plot_radii = np.append(plot_radii,radius)

    return plot_radii, [big_dm1, big_dm1e, big_mm1],\
        [big_dm2, big_dm2e, big_mm2],\
        [big_dm3, big_dm3e, big_mm3]

def giant_steps(slay_file, simstring, parameter_list,skip_radii=[]):
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

        radii, m1, m2, m3 = moments_notice(slay_file,sim_file,
                                           skip_radii=skip_radii)
        chis = np.array([])
        for moment in [m1, m2, m3]:
            red_chi = (moment[0] - moment[2])**2/moment[1]**2
            chis = np.r_[chis,red_chi]

        value = np.sum(chis)/(radii.size*3 - 1)
        print 'fit goodness: {}'.format(value)
        par1_arr = np.append(par1_arr,par1)
        par2_arr = np.append(par2_arr,par2)
        value_arr = np.append(value_arr,value)
    
    return par1_arr, par2_arr, value_arr

def cutting_session(drunkdict, output_file, 
                    pars=[239.,5.5,1.62,8.43,0.43], name='boring',
                    size=1001,fixed=[],flare=False):
    '''Use an amoeba algorithm to find a model galaxy that is the best fit (in
    a moment sense) to a set of data. Any number of heights can be fit
    simultaneously with judicious use of the slay dict, which is expected to
    have the following format:
    
        drunkdict = {hight: ['drunkfile.fits',Flip_boolean,[skip_radii]]},

    where height is a number (in scale heights) and Flip_boolean tells the
    program it needs to flip the radii around zero.

    If there are bad data at certain radii you can skip them over using
    skip_radii. You don't even need to get the exact radius correct, the check
    is done that skip_radii[i] == floor(measured_radius).

    pars and fixed are used to control which parameters get fit. Pars is a
    list with every model parameter and fixed is a list of indices for which
    parameters shoud NOT be fit (i.e., fixed). For example, the first
    parameter is the asymptotic speed so if you wanted to have all models be
    forced to have an asymptotic speed of 270 km/s then pars would be
    [270,...] and fixed=[0]

    The output of cutting session is a tuple containing the best fit
    parameters and a bar_dict (see moments_notice) for use with plotting
    functions.
    '''
    
    p0 = [pars[i] for i in range(len(pars)) if i not in fixed]

    print 'p0 is {}'.format(p0)
    #raw_input('please confirm')
    f = open(output_file,'w')
    f.write('# Generated on {}\n'.format(time.asctime())
            +'# pars: {}\n'.format(pars)
            +'# p0: {}\n'.format(p0))

    datadict = {}
    for z in drunkdict.keys():
        datadict[z] = [dd.open_drunk(drunkdict[z][0],
                                     skip_radii=drunkdict[z][2]),
                       drunkdict[z][1],
                       drunkdict[z][2]]
    print "minimizing..."
    pf = spo.fmin(solo,p0,args=(datadict,name,f,size,pars,fixed,flare))
    
    pfl = list(pf)

    parsf = [pars[j] if j in fixed else pfl.pop(0) for j in range(len(pars))]

    f.write('# final_pars: {}\n'.format(parsf))
    f.close()

    bar_dict = {}
    chis = np.array([])
    if flare:
        flarepars = {'h_zR': parsf[-1], 'ftype':flare}
    else:
        flarepars = None

    for z in datadict.keys():
        simfile = make_boring([parsf[0]],[parsf[1]],name='final_'+name,
                              size=size,z=z,
                              h_dust=parsf[3],kappa_0=parsf[2],z_d=parsf[4],
                              flarepars=flarepars)[0]
        
        bar = moments_notice(datadict[z][0],simfile,
                             skip_radii=datadict[z][2],
                             flip=datadict[z][1])
        bar_dict[z] = bar

        for moment in bar[1:]:
            red_chi = (moment[0] - moment[2])/moment[1]
            chis = np.r_[chis,red_chi]

    value = np.sum(chis**2)/(chis.size - pf.size - 1)
        
    return parsf, bar_dict, value

def solo(p,datadict,name,output_file,size,par0,fixed,flare):
    '''This is the minimizing function that is called by cutting_session. It
    constructs the necessary model galaxies, computes the moments of the
    lines, and returns a reduce chi squared value quantifying the goodness of
    fit.

    The reduced chi-squared is computed across all moments and all heights. So
    if you are fitting 3 heights that each have N radial bins (as they should,
    if using correctly binned data) then the chi-squared will be computed with
    3x3xN points (three from moments, three from heights).
    '''

    pl = list(p)
    pars = [par0[j] if j in fixed else pl.pop(0) for j in range(len(par0))]

    chis = np.array([])
    if flare:
        flarepars = {'h_zR':pars[-1],'ftype':flare}
    else:
        flarepars = None

    output_file.write(str('{:11.4f} '*len(pars)).format(*pars))

    for z in datadict.keys():
        
        simfile = make_boring([pars[0]],[pars[1]],name=name,size=size,z=z,
                              h_dust=pars[3],kappa_0=pars[2],z_d=pars[4],
                              flarepars=flarepars,nofits=True)[0]
        
        radii, m1, m2, m3 = moments_notice(datadict[z][0],simfile,
                                           skip_radii=datadict[z][2],
                                           flip=datadict[z][1],nofits=True)

        for moment in [m1,m2,m3]:#m1,m2,m3]:
            red_chi = (moment[0] - moment[2])/moment[1]
            output_file.write('{:11.4f} '.format(bn.nansum(red_chi**2)))
            chis = np.r_[chis,red_chi]

    value = bn.nansum(chis**2)/(chis.size - p.size - 1)
    output_file.write('{:11.4f}\n'.format(value))
    print 'v_r: {}\nh_rot: {}\nkappa_0: {}\nh_dust: {}\nz_dust: {}\nvalue: {}\n'.\
        format(pars[0],pars[1],pars[2],pars[3],pars[4],value)
    return value

def make_boring(vr_list, h_rot_list, h_dust=8.43, kappa_0=1.62,nofits=False,
                z=0, size=1001, z_d=0.23, name='boring',flarepars=None):
    '''Given a list of values for v_r and h_rot, make a grid of galaxy models
    with all possible combinations of those two parameters.
    '''
    basename = 'sim_z{:n}_{}'.format(z,name)
    outlist = []
    for v_r in vr_list:
        for h_rot in h_rot_list:
            name = '{}_{}.fits'.format(basename,int(round(time.time(), 3)*100))
            sim = salty.simcurve(size,z,v_r,h_rot,output=name,scale=0.0999,
                                 h_dust=h_dust,kappa_0=kappa_0,z_d=z_d,
                                 flarepars=flarepars,full=True,
                                 verbose=False,nofits=nofits)
            if nofits:
                outlist.append(sim)
            else:
                outlist.append(name)

    return outlist
    
def plot_ice(results,title=' ',show=True, rot=False, ax=False, label=''):
    '''Given the results of moments_notice(), plot model vs. data for all
    three moments. Show the plots and return them.
    '''

    if ax:
        ax1, ax2, ax3 = ax
    else:
        fig = plt.figure()
        ax3 = fig.add_subplot(313)
        ax1 = fig.add_subplot(311,sharex=ax3)
        ax2 = fig.add_subplot(312,sharex=ax3)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        fig.subplots_adjust(hspace=0.0001)
        ax1.set_ylabel('$\mu_1$')
        ax2.set_ylabel('$\sqrt{\mu_2}$')
        ax3.set_ylabel('$\mu_3$')
        ax3.set_xlabel('radius [kpc]')
        ax1.set_ylim(0.001,290.001)
        ax2.set_ylim(0.001,49.999)
        ax3.set_ylim(-1,0.999)
        fig.suptitle('{}\n{}'.format(title,time.asctime()))

    radius = results[0]
    sortidx = np.argsort(np.abs(radius))
    
    for ax, data, flip in zip([ax1, ax2, ax3],[results[1],results[2],results[3]],[True,False,True]):
        pdata = data[0]
        perr = data[1]
        pmod = data[2]
        if flip:
            pdata = np.abs(pdata)
            pmod = np.abs(pmod)
        ax.errorbar(np.abs(radius),pdata,yerr=perr,fmt='.b',ms=7)
        ax.plot(np.abs(radius)[sortidx],pmod[sortidx],label=label)

    if rot:
        rotr = np.linspace(0,np.abs(radius).max())
        rc = rot[0]*np.tanh(rotr/rot[1])
        ax1.plot(rotr,rc,'-r')

    if show: fig.show()

    return ax1, ax2, ax3

def multi_plot(bar_dict,label=''):
    '''Takes a bar_dict from a multiple-height run of cutting session and
    plots the first three moments for each height.
    '''

    figlist = []
    for z in bar_dict.keys():
        fig = plot_ice(bar_dict[z],'{}\nz={}'.format(label,z))
        figlist.append(fig)

    return figlist

def make_fits(vrlist,hrlist,footprints,output):
    '''Generate a fits version of the goodness-of-fit map from a comparison of
    data to a grid of models.
    '''
    shape = (int(footprints[0].size**0.5),int(footprints[0].size**0.5))
    hdu = pyfits.PrimaryHDU(footprints[2].reshape(shape))
    hdu.header.update('CRPIX1',1,comment='WCS: X reference pixel')
    hdu.header.update('CRVAL1',footprints[1].min(),
                      comment='WCS: X reference coordinate value')
    hdu.header.update('CDELT1',np.mean(np.diff(hrlist)),
                      comment='WCS: X pixel size')
    hdu.header.update('CTYPE1','LINEAR',comment='X type')
    hdu.header.update('CRPIX2',1,comment='WCS: Y reference pixel')
    hdu.header.update('CRVAL2',footprints[0].min(),
                      comment='WCS: Y reference coordinate value')
    hdu.header.update('CDELT2',np.mean(np.diff(vrlist)),
                      comment='WCS: Y pixel size')
    hdu.header.update('CTYPE2','LINEAR',comment='Y type')
    hdu.writeto(output)

    return


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

def mintest():

    x = np.arange(100)
    trua = np.random.random()*100
    trub = np.random.random()
    y = x*trua + trub
    y += np.random.rand(100)*10.

    p0 = [0.,0.]
    fmin = spo.fmin(fmin_func,p0,args=(x,y))
    leastsq = spo.leastsq(least_func,p0,args=(x,y))

    return trua, trub, fmin, leastsq

def fmin_func(p,x,y):

    testy = x*p[0] + p[1]
    chi = np.sum((testy - y)**2)
    return chi

def least_func(p,x,y):

    testy = x*p[0] + p[1]
    return (testy - y)

def many_zd(num):

    # zds = np.random.random(num)*3.0
    # kaps = np.random.random(num)*3.0

    zds = np.linspace(0.15,1.5,num)
    kaps = np.linspace(0.4,1.5,num)

    print zds, kaps
    raw_input('asda')

    slaydict = {0: ['../../SALT_data/tiESO_z0_MgI_binz2.slay.fits', False], 
                1: ['../../SALT_data/tiESO_z1_MgI_binz2.slay.fits', True], 
                2: ['../../SALT_data/tiESO_z2_MgI_bin60.slay.fits', True]}
    bars = []
    values = np.array([])
    outkap = np.array([])
    outzd = np.array([])

    for zd in zds:
        for kap in kaps:
            name = 'zd_{:5.3f}'.format(zd,int(round(time.time(),3)*100))
            bar = cutting_session(slaydict,name+'.out',pars=[250,3.375,kap,8.18,zd],
                                  skip_radii=[31,-33],fixed=[0,1,3],name=name)
            os.system('rm sim_z?_[^f]*.fits')
            bars.append([zd,kap,bar])
            outkap = np.append(outkap,kap)
            outzd = np.append(outzd,zd)
            values = np.append(values,bar[2])

    return bars, zds, kaps, [outzd, outkap, values]
