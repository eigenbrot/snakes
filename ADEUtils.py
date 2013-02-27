import numpy as np
import matplotlib.pyplot as plt
import pyfits
import math
import os
import scipy.optimize as spo
import bottleneck as bn
from numba.decorators import jit, autojit
from numba import float32, int16

debug=0

def ADE_gauss(length, offset, sigma, PEAK_VAL=1, FWHM=0, NORM=False):
    '''
    Calling:
        ADE_gauss(length, offset, sigma, PEAK_VAL=1, FWHM=0, NORM=True)

    Description:
        ADE_gauss generates a gaussian curve that can be customized to the
        user's desires.

    Inputs:
        length - Int
            The length of the output array
        offset - Float
            The offset that defines where the peak of the gaussian is
        sigma  - Float
            The standard deviation of the gaussian
        PEAK_VAL-Float
            The gaussian's maximum value
        FWHM   - Float
            If FWHM != 0 then the gaussian is forced to have whatever
            FWHM is specified. Setting this parameter overrides the sigma
            parameter.
        NORM   - Boolean
            If NORM is True then the gaussian is normalized to have total
            area equal to one. Setting this parameter overrides the
            PEAK_VAL parameter.

     Output:
         The output is a tuple of numpy arrays where the first array is the
         x vector used to create the gaussian and the second array is
         the actual gaussian computed using:

               gauss = PEAK_VAL*exp(-1*(x - offset)**2/s*sigma**2)
     '''

    xvec = np.arange(length, dtype=np.float64)

    if FWHM: 
        sigma = FWHM/2.35482
    
    outgauss = PEAK_VAL*np.exp(-1*(xvec - offset)**2/(2*sigma**2))

    if NORM: 
        outgauss /= np.sum(outgauss)

    return xvec,outgauss

def ADE_gauss2d(dims, offset, sigma, PEAK_VAL=1, FWHM=0, NORM=True):
    '''
    Calling:
        ADE_gauss2d(dims, offset, sigma, PEAK_VAL=1, FWHM=0, NORM=1)

    Description:
        ADE_gauss generates a 2d-gaussian array that can be customized to 
        the user's desires.

    Inputs:
        dims   - Tuple
                 The dimensions of the output array
        offset - Tuple
                 The center of the gaussian
        sigma  - Float
                 The standard deviation of the gaussian
        PEAK_VAL-Float
                 The gaussian's maximum value
        FWHM   - Float
                 If FWHM != 0 then the gaussian is forced to have whatever
                 FWHM is specified. Setting this parameter overrides the sigma
                 parameter.
        NORM   - Bool
                 If NORM then the gaussian is normalized to have total
                 area equal to one. Setting this parameter overrides the
                 PEAK_VAL parameter.

     Output:
         The output is a 2d ndarray containing the gaussian

     '''
    ids = np.indices(dims, dtype=np.float32)

    if FWHM: 
        sigma = FWHM/2.35482
    
    outgauss = PEAK_VAL*np.exp(-1*((ids[0] - offset[0])**2/(2*sigma**2) +\
                                    (ids[1] - offset[1])**2/(2*sigma**2)))

    if NORM: 
        outgauss /= np.sum(outgauss)

    return outgauss

def ADE_lorentz(length, offset, gamma, PEAK_VAL=1, FWHM=0, NORM=0):
    '''
    Calling:
       ADE_lorentz(length, offset, gamma, PEAK_VAL=1, FWHM=0, NORM=0)

    Description:
      ADE_lorentz generates a lorentzian function that can be customized
      as the user sees fit.

    Inputs:
        length - Int
            The length of the output arrays
        offset - Float
            The offset that defines the location of the peak
        gamma -  Float
            The width parameter of the lorentzian
        PEAK_VAL-Float
            The maximum value of the lorentzian
        FWHM   - Float
            If FWHM != 0 then the output lorentzian is forced to have
            the FWHM specified. If set this parameter overrides the gamma
            parameter.
        NORM   - Int
            Set NORM = 1 if you want the output lorentzian to be normalized
            to a total area of 1. If set this parameter overrides the
            PEAL_VAL paramter.

    Outputs:
        A tuple of numpy arrays, with the first being the x vector used to
        compute the lorentzian and the second being the actual lorentzian,
        which is calculated thusly:

           lorentz = PEAK_VAL*(gamma**2/((x - offset)**2 + gamma**2))
    '''

    xvec = np.arange(length, dtype=float)

    if FWHM:
        gamma = FWHM/2

    outlorentz = PEAK_VAL*(gamma**2/((xvec - offset)**2 + gamma**2))

    if NORM:
        outlorentz /= np.sum(outlorentz)

    return xvec,outlorentz

def ADE_lorentz2d(dims, offset, gamma, PEAK_VAL=1, FWHM=0, NORM=True):
    '''
    Calling:
       ADE_lorentz(dims, offset, gamma, PEAK_VAL=1, FWHM=0, NORM=0)

    Description:
      ADE_lorentz2d generates a 2 dimensional lorentzian function that 
      can be customized as the user sees fit.

    Inputs:
        dims   - Tuple
                 The size of the output array
        offset - Tuple
                 The (x,y) coordinates of the peak of the lorentzian
        gamma -  Float
                 The width parameter of the lorentzian
        PEAK_VAL-Float
                 The maximum value of the lorentzian
        FWHM   - Float
                 If FWHM != 0 then the output lorentzian is forced to have
                 the FWHM specified. If set this parameter overrides the gamma
                 parameter.
        NORM   - Boolean
            Set NORM to True if you want the output lorentzian to be normalized
            to a total area of 1. If set this parameter overrides the
            PEAL_VAL paramter.

    Outputs:
        A 2d numpy array with size specified by dims. The values of the array
        are computed via
    
           lorentz = PEAK_VAL*(gamma**2/((x - x0)**2 + (y - y0)**2 + gamma**2))
    '''

    ids = np.indices(dims, dtype=np.float32)

    if FWHM:
        gamma = FWHM/2

    outlorentz = PEAK_VAL*(gamma**2/(\
            (ids[0] - offset[0])**2 + (ids[1] - offset[1])**2 + gamma**2\
                ))
                               
    if NORM:
        outlorentz /= np.sum(outlorentz)

    return outlorentz

def multi_where(vec1, vec2):
    '''Given two vectors, multi_where returns a tuple of indices where those
    two vectors overlap.
    ****THIS FUNCTION HAS NOT BEEN TESTED ON N-DIMENSIONAL ARRAYS*******
    Inputs:
           2 numpy vectors
    Output:
           (xy, yx) where xy is a numpy vector containing the indices of the
           elements in vector 1 that are also in vector 2. yx is a vector
           containing the indices of the elements in vector 2 that are also
           in vector 1.
    Example:
           >> x = np.array([1,2,3,4,5])
           >> y = np.array([3,4,5,6,7])
           >> (xy,yx) = multi_where(x,y)
           >> xy
           array([2,3,4])
           >> yx
           array([0,1,2])
    '''

    OneInTwo = np.array([])
    TwoInOne = np.array([])
    for i in range(vec1.shape[0]):
        if np.where(vec2 == vec1[i])[0].shape[0]:
            OneInTwo = np.append(OneInTwo,i)
            TwoInOne = np.append(TwoInOne, np.where(vec2 == vec1[i])[0][0])

    return (np.int8(OneInTwo), np.int8(TwoInOne))


def mode(data):
    '''Compute the mode for an arbitrary array.
    From http://projects.scipy.org/scipy/ticket/905
    Output is a tuple of the list of mode values and how many of them there are
    '''
    counts = {}
    for x in data.flatten():
	counts[x] = counts.get(x,0) + 1
    maxcount = max(counts.values())
    modelist = []
    for x in counts:
	if counts[x] == maxcount:
	    modelist.append(x)
    return modelist,maxcount


def centroid(image, zhi=2):
    '''takes in an array and returns a list containing the
    center of mass of that array'''
  
    if debug: print('finding center...')

    '''just to make sure we don't mess up the original data somehow'''
    data = np.copy(image)

    size = data.shape
    totalMass = np.sum(data)
    
    xcom = np.sum(np.sum(data,1) * np.arange(size[0]))/totalMass
    ycom = np.sum(np.sum(data,0) * np.arange(size[1]))/totalMass
    
    return (xcom,ycom)

def dist_gen(data, center):
    '''takes in a data array with known center and returns an
    array where the value at each location is the distance of
    that point from the center.
    center - expected to be a list of two numbers'''
    
    vecvec = np.indices(data.shape,dtype=np.float32)
    distances = ((center[0] - vecvec[0,])**2 + (center[1] - vecvec[1,])**2)**0.5

    return distances

def angle_gen(data, center):
    '''
    Description:
        angle_gen generates an array of the same size as data where
        the value at each point is the angle coordinate of that point
        in polar coords with the origin at center. The angles range
        from 0 -> tau.
    '''

    vecvec = np.indices(data.shape,dtype=np.float32)
    vecvec[1] = center[1] - vecvec[1]
    vecvec[0] = center[0] - vecvec[0]
    
    'arctan2 to the rescue!'
    angles = np.arctan2(vecvec[0],vecvec[1]) + np.pi

    return angles

def appetize(image, num_ap, OUTPUT=0, NOREAD=0, DRAW=0, EXTEN=0):
    '''Appetize takes in an image and returns a vector of apperture radii
    and a vector of corresponding fluxes inside those appertures.
    image - a string with the path to the FITS file to be used
    num_ap - an integer that is the desired number of apertures
    OUTPUT - a string containing the write out file. This file will have
             both the aperture and flux vectors
    NOREAD - set this to 1 if image is a data array (i.e. already read from
             a FITS file
    EXTEN - the FITS HDU that contains the data
    '''
   
    if debug: print('appetizing...')
 
    if NOREAD: data = np.float64(image)
    else: data = np.float64(pyfits.open(image)[EXTEN].data)

    center = centroid(data)
    dims = data.shape

    distances = dist_gen(data, center)

    rlimit = distances.max()
    ap_step = rlimit/num_ap

    'this area will be used to convert flux to surface brightness'
    area_array = np.zeros(dims) + 1

    'initialize future output arrays. I think that maybe with python'
    'these can be left empty, but w/e'
    fluxes = np.zeros(num_ap)
    ap_vec = np.zeros(num_ap)
    ap_edge = 0.0


    for i in range(num_ap):
        ap_edge += ap_step
        idxs = np.where(distances < ap_edge)
        flux = np.sum(data[idxs])
        ap_area = np.sum(area_array[idxs])
        correction = math.pi*ap_edge**2/ap_area
        if debug == 2: print str(i)+"th correction is "+str(correction)
        fluxes[i] = flux*correction
        ap_vec[i] = ap_edge

    
    if OUTPUT:
        f = open(OUTPUT, 'w')
        for i in range(fluxes.shape[0]):
            np.array([ap_vec[i],fluxes[i]]).tofile(f, sep=' ',format='%9.3e')
            f.write('\n')

    if DRAW:
        fig = plt.figure(0)
        plt.plot(ap_vec,fluxes)
        fig.show()
    

    return(ap_vec, fluxes)

def annulize_sb(image, num_an, OUTPUT=0, NOREAD=0, DRAW=0, EXTEN=0, MODE=0):
    '''takes in a FITS image and computes the flux within equal area annuli
    that are as numerous as the user desires. The flux is weighted by area
    so you really are getting surface brightness. Returns a vector of annulus
    radii and a vector of surface brightnesses
    image - a string with the path to the FITS file to be used
    num_ap - an integer that is the desired number of apertures
    OUTPUT - a string containing the write out file. This file will have
             both the aperture and flux vectors
    NOREAD - set this to 1 if image is a data array (i.e. already read from
             a FITS file
    EXTEN - the FITS HDU that contains the data
    '''

    if debug: print('annulizing...')

    if NOREAD: data = np.float64(image)
    else: data = np.float32(pyfits.open(image)[EXTEN].data)

    if debug: print("data type is", data.dtype)

    center = centroid(data)
    dims = data.shape

    if MODE:
        if debug: print("subtracting mode")
        data -= MODE
       
    distances = dist_gen(data, center)

    rlimit = distances.max()
    Alimit = math.pi*(rlimit**2)
    area = Alimit/num_an

    'this area will be used to convert flux to surface brightness'
    area_array = np.zeros(dims) + 1

    'initialize future output arrays. I think that maybe with python'
    'these can be left empty, but w/e'
    fluxes = np.zeros(num_an)
    r_vec = np.zeros(num_an)
    r1 = 0.0
    r2 = (area/math.pi)**0.5
    
    for i in range(num_an):
        idx1 = np.where(distances <= r2)
        idx2 = np.where(distances[idx1] > r1)
        flux = np.sum(data[idx1][idx2])
        an_area = np.sum(area_array[idx1][idx2])
        fluxes[i] = flux/an_area

        if debug == 2: print(i,r1,r2,flux,an_area,fluxes[i]) 

        r_mid = (0.5*(r2**2 + r1**2))**0.5
        r_vec[i] = r_mid
        r1 = r2
        r2 = ((area/math.pi) + r1**2)**0.5

    if OUTPUT:
        f = open(OUTPUT, 'w')
        for i in range(fluxes.shape[0]):
            np.array([r_vec[i],fluxes[i]]).tofile(f, sep=' ',format='%9.3e')
            f.write('\n')
            
    if DRAW:
        fig = plt.figure(DRAW)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_xlabel('r [px]')
        ax.set_ylabel('Surface Brightness [ADU px$^{-1}$]')
        ax.set_title(image)
        ax.plot(r_vec,fluxes)
        fig.show()

    return(r_vec,fluxes,center)



def annulize(data, num_an, distances=np.array([0])):
    '''
    Description:
        annulize takes in a FITS image and computes the total power 
        contained within progressivly larger annuli. The number of 
        annuli used is set by the user and, unlike annulize_sb, it 
        is the width of the annuli that remains constant, not the area.

    Inputs:
        data      - ndarray
                    The data to be annulized
        num_an    - Int
                    The number of annuli to use
        distances - ndarray
                    This is expected to be a transformed distance array

    Output:
        r_vec - ndarray
                vector of radii where each entry is the radius of the
                middle of the corresponding annulus.
        fluxes- ndarray
                each entry is the flux contained within that annulus.
    '''

    if debug: print'annulizing...'
    if debug: print"data type is "+str(data.dtype)

    dims = data.shape

    '''check to see if we got some transformed distances and generate
    a basic distance array if we didn't'''
    if not(distances.any()):
        center = centroid(data)
        if debug: print "center at "+str(center)
        distances = dist_gen(data, center)


    rlimit = distances.max()
    rstep = rlimit/num_an

    'initialize future output arrays. I think that maybe with python'
    'these can be left empty, but w/e'
    fluxes = np.zeros(num_an, dtype=np.float32)
    r_vec = np.zeros(num_an, dtype=np.float32)
    outarray = np.zeros(dims, dtype=np.float32)
    area_array = np.zeros(dims, dtype=np.float32) + 1
    errors = np.zeros(num_an, dtype=np.float32)
    r1 = 0.0
    r2 = rstep
    
    for i in range(num_an):
        idx = np.where((distances <= r2) & (distances > r1))
        
        '''The correction is designed to correct annuli that are not
        entirely contained in the data array, but I still can't decide
        if it screws up the data or not'''
        correction = math.pi*(r2**2 - r1**2)/np.sum(area_array[idx])
        fluxes[i] = np.sum(data[idx])#*correction
        errors[i] = np.std(data[idx])

        if debug == 3: print(i,r1,r2,fluxes[i],correction) 
        
        '''this is only used during debug to show where the annuli
        lie on the image'''
        outarray[idx] = np.mod(i,2) + 1

        r_mid = (r1 + r2)*0.5
        r_vec[i] = r_mid
        r1 = r2
        r2 = r1 + rstep

    if debug == 2:
        pyfits.PrimaryHDU(data+outarray*data.max()/2)\
            .writeto('an_show.fits',clobber=True)
                    
    return(r_vec,fluxes,errors)

def fast_annulize(data, numan, distances=np.array([0])):

    if debug: print'annulizing...'
    if debug: print"data type is "+str(data.dtype)

    dims = data.shape

    '''check to see if we got some transformed distances and generate
    a basic distance array if we didn't'''
    if not(distances.any()):
        center = centroid(data)
        if debug: print "center at "+str(center)
        distances = dist_gen(data, center)


    rlimit = distances.max()
    ddata = np.array(data,dtype=np.float32)

    return fast_helper(ddata, distances, rlimit, numan)

@jit(argtypes=[float32[:,:],float32[:,:],float32,int16], restype=float32[:,:])
def fast_helper(data, distances, rlimit, numan):
    
    size = data.shape
    rstep = rlimit / numan

    r1 = 0.0
    r2 = float(rstep) # we need to cast this as a float so that numba doesn't
                      # complain

    r_vec = np.zeros(numan,dtype=np.float32)
    mean_vec = np.zeros(numan,dtype=np.float32)
    error_vec = np.zeros(numan,dtype=np.float32)

    for k in range(numan):
        anlist = []
        for i in range(size[0]):
            for j in range(size[1]):
                if distances[i,j] > r1:
                    if distances[i,j] <= r2:
                        anlist.append(data[i,j])
        
        anarray = np.array(anlist,dtype=np.float32)
        mean_vec[k] = bn.nansum(anarray)
        error_vec[k] = bn.nanstd(anarray)
        r_vec[k] = (r1 + r2)*0.5
        
        r1 = r2
        r2 += rstep

    return np.array([r_vec,mean_vec,error_vec])

def mediclean(data, zhi=0, keep=False):
    '''
    Description:
        mediclean is essentially a median subtraction algorithm.
        It takes in a data array and subtracts each column's median
        from the column. It then does the same thing for the rows.
        When computing the median, all data above some fixed value 
        is ignored. This algorithm is very good if you have excellent S/N.

    Inputs:
        data - ndarray
               The data to be cleaned. Should probably be 2d.
        zhi  - Float
               The zero point for defining which data points to
               ignore when computing the median. All values greater
               than zhi are ignored. The defualt is 3 standard
               deviations above the median of the data.
        keep - Bool
               Turn this on if you want to keep the intermediat steps.
               As a result you will get fits file containing the medians
               of each row and column.

    Output:
        The cleaned data array.
    '''
    size = data.shape

    'set the default value of zhi if we have to'
    std = np.std(data)
    if not zhi: zhi = np.median(data) + 3*std; k=1

    m1 = np.empty(size)
    '''we need this for loop (rather than threading) because we have
    to reject the pixels higher than zhi'''
    for r in range(data.shape[0]):
        m1[r] = np.median(data[r][np.where(data[r] < zhi)])

    if keep: pyfits.PrimaryHDU(m1).writeto('temp1.fits', clobber=True)

    'transposition keeps the code the same'
    d1 = np.transpose(data - m1)

    std = np.std(d1)
    if not k: zhi = np.median(d1) + 3*std

    m2 = np.empty(d1.shape)
    for c in range(d1.shape[0]):
        m2[c] = np.median(d1[c][np.where(d1[c] < zhi)])

    if keep: 
        pyfits.PrimaryHDU(np.transpose(m2)).writeto('temp2.fits', clobber=True)


    return  np.transpose(d1 - m2)

def polyclip(x, y0, deg, niter=20, clip_high=3., clip_low=3.):
    '''
    polyclip is a simple wrapper for Numpy's polyfit polynomial fitting routine
    that adds iterative, sigma-based rejection of "bad" points. It returns a
    np.poly1d object that can is a callable function representing the best-fit
    polynomial.

    Inputs:
        x - ndarray
            The x data. Should be 1d
        
        y0 - ndarray
            The y data to be fit. Should be 1d and the same length as x.

        deg - int
            The degree of polynomial to be fit

        niter - int
            The number of clipping iterations to perform. This algorithm is 
            pretty fast, so don't worry about setting this number too high.
            
        clip_high - float
            Any points with residuals more than clip_high*sigma above zero will 
            be replaced with values from the fit.

        clip_low - float
            Any points with residuals less than -1*|clip_low|*sigma below zero
            will be replaced with values from the fit.

    Outputs:
        A callable np.poly1d object representing the best fit to y(x).
    '''

    '''We up-type y because poly1d will always return a float64 array
    and we don't want casting errors later on'''
    y = np.array(y0,copy=True,dtype=np.float64)

    for i in range(niter):
        coef = np.polyfit(x,y,deg)
        fit = np.poly1d(coef)

        residuals = y - fit(x)
        sigma = np.std(residuals)
        
        highidx = np.where(residuals > clip_high*sigma)
        '''the np.abs is there in case the user sets clip_low to a negative
        number'''
        lowidx = np.where(residuals < -1.*np.abs(clip_low)*sigma)
        
        '''replace the bad points with points from the fit'''
        y[highidx] = fit(x[highidx])
        y[lowidx] = fit(x[lowidx])

    return np.poly1d(coef)

def eplot(x,y):
    '''a simple script for quickly plotting some data.
    It avoids all the tedious calls to create figures and axes, etc.
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    fig.show()
    return

def fit_line(x, y, error):
    '''finds a slope and intercept such that
    y - m*x + b is minimized. The initial fit is done w/ some linear algebra, 
    and then error considerations are used to refine the fit with some chisq 
    minimization.
    '''
    A = np.vstack([x, np.ones(x.size)]).T

    m, c = np.linalg.lstsq(A,y)[0]
    
    fit = spo.leastsq(chi,np.array([m,c]), args=(x,y,error))[0]

    model = fit[0]*x + fit[1]
    rchisq = ((y - model)**2/error**2).sum()/(y.size - 3)
    rms = (np.sum((y - model)**2)/y.size)**0.5
    
    return fit

def chi(p, x, y, error):
    '''the minmizing function used by fit_line'''
    model = p[0]*x + p[1]

    return (y - model)/error

def ADE_moments(x,p,threshold=np.inf,err=None):
    '''Computes the first three moments of the probability function p(x). 
    Technically, it returns the first moment (center) and then the second and
    third central moments about the first moment.
    It can try some iterative clipping to help with noisy data, but this is
    somewhat dangerous. Set threshold to something other than np.inf to perform
    this clipping.
    '''

    if err != None: p /= err

    norm_p = p/np.sum(p)

    std_p = np.array(norm_p)
    #eplot(x,std_p)
    delta_std = np.inf
    while delta_std > threshold:
        # eplot(x,std_p)
        # _ = raw_input('bel')
        std = np.std(std_p)
        sidx = np.where(np.abs(std_p) > 2*std)
        std_p[sidx] = 0.0
        newstd = np.std(std_p)
        delta_std = np.abs(std - newstd)

#    eplot(x,std_p)
    m1 = np.sum(x*norm_p)
    m2 = np.sum(norm_p * (x - m1)**2)
    
    delta_m2 = np.inf
    while delta_m2 > threshold:
        print m2, delta_m2
        std = np.std(norm_p)
        idx = np.where(np.abs(norm_p) < 1.*std)
        norm_p[idx] = 0.0
        std_p[idx] = 0.0
        norm_p /= np.sum(norm_p)
        m1 = np.sum(x*norm_p)
        new_m2 = np.sum(norm_p * (x - m1)**2)
        delta_m2 = np.abs(new_m2 - m2)
        m2 = new_m2

    m3 = np.sum(norm_p * (x - m1)**3)/m2**1.5

    return np.array([m1, m2, m3])

