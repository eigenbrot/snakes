import numpy as np
import ADEUtils as ADE
import scipy as sp
import pyfits
from pyraf.iraf import stsdas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


def chi(r, f, radius, idx):

    fn = np.copy(f)/np.sum(f)
#    idx = np.where(r >= peak)[0][0]
    size = np.float64(r.size)

    numtest = 1000
    best = 99999999.99
    goodshift = -999
    for i in range(numtest):
        i -= numtest/2
        shift = i/10.
        testcent = idx + shift
        
        theta = np.arccos((r - r[testcent])/radius)
        pn = radius*np.nan_to_num(np.sin(theta))
        
        pn /= np.sum(pn)

        nonz = np.where(pn > 0)[0]

        chisq = np.sum(((fn[nonz] - pn[nonz])**2)/pn[nonz])

        if chisq < best: 
            best = chisq
            goodshift = shift
    
    print(goodshift)
    return goodshift

def decon(r,f,peak):

    fig0 = plt.figure(0)
    plt.clf()

#    x1 = np.where(r >= r1)[0][0]
#    x2 = np.where(r >= r2)[0][0]

    radius = (500.0/5.2)/(2)#*(r2-r1)/(x2-x1))


    idx = np.where(r >= peak)[0][0]
#    idx = (x2-x1)
    idx += chi(r,f,radius, idx)

    print(idx)

    size = f.shape[0]

    theta = np.arccos((r - r[idx])/radius)
#    theta = np.arccos((np.arange(2*int(radius+1)+1) - int(radius+1))/radius)
    pn = radius*np.nan_to_num(np.sin(theta))


#    bigsig = np.zeros(2*size-1+k)
#    bigsig[0:size] = np.copy(f)/np.sum(f)
#    bigin = np.zeros(2*size -1+k)
#    bigin[0:size] = np.copy(pn)/np.sum(pn)

    ff = np.copy(f)/np.sum(f)
    pnf = np.copy(pn)/np.sum(pn)

    plt.plot(r,ff)
    plt.plot(r,pnf)
    fig0.show()

#    numfreq = ff.size + pnf.size - 1

#    ffs = np.fft.rfft(ff, n=numfreq)
#    ffi = np.fft.rfft(pnf, n=numfreq)
#    fff = ffs/ffi

#    frd = np.real(np.fft.irfft(fff))

#    plt.figure(1)
#    plt.clf()
#    plt.plot(np.arange(frd.size),frd)
    
#    pyfits.PrimaryHDU(ff).writeto('ff.fits')
#    pyfits.PrimaryHDU(pnf).writeto('pn.fits')

#    stsdas()
#    stsdas.analysis()
#    stsdas.analysis.restore()
#    stsdas.analysis.restore.lucy('ff.fits','pn.fits','frd.fits',1,0)

 #   out = raw_input('run IDL program now and enter the name of the output\n')

#    frd = pyfits.open('idl_frd.fits')[0].data

    return (r,pn)

def test(size):

    r_vec = np.arange(size)*5
    width = np.random.uniform(size/10,size/5)
    peak = np.random.uniform(width,size-width/2)
    print(peak)
    r1 = peak - width/2
    r2 = peak + width/2    

    x1 = np.where(r_vec >= r1)[0][0]
    x2 = np.where(r_vec >= r2)[0][0]

    radius = (500.0/5.2)/(2*(r2-r1)/(x2-x1))

    theta = np.arccos((np.arange(size) - peak)/radius)
    pn = radius*np.nan_to_num(np.sin(theta))

    #Now smear it
    g_width = np.random.uniform(4,5)
    print(g_width)
    gauss = np.exp((-1*(r_vec - 2*peak)**2)/(2*g_width**2))

    f = np.convolve(gauss,pn,'same')

    plt.figure(2)
    plt.clf()
    plt.plot(r_vec,pn)
    plt.plot(r_vec,f)
    plt.figure(3)
    plt.clf()
    plt.plot(r_vec,gauss)

    scratch = decon(r_vec,f,peak,r1,r2,0)

    return

def lorentz_heur(params, f, NORM=0):

    gamma = params[0]
    peak = params[1]
    offset = params[2]

    length = f.size
    lorentz = ADE.ADE_lorentz(length,offset,gamma,PEAK_VAL=peak,NORM=NORM)[1]

    'find the ChiSq difference b/t the lorentzian and the data'
    chisq = np.sum((f - lorentz)**2/lorentz)

    return chisq

def decon1d(g,f):
    length_h = len(g) - len(f) + 1
    mtx = [[0 for x in range(length_h + 1)] for y in g]
    for h_idx in range(length_h):
        for f_idx, f_val in enumerate(f):
            g_idx = h_idx + f_idx
            mtx[g_idx][h_idx] = f_val
    
    for g_idx, g_val in enumerate(g):
        mtx[g_idx][length_h] = g_val
    RRE_form(mtx)
    return [mtx[i][length_h] for i in range(length_h)]

def RRE_form(M):
    if not M: return
    lead = 0
    numrows = len(M)
    numcols = len(M[0])
    for r in range(numrows):
        if lead >= numcols:
            return
        i = r
        while M[i][lead] == 0:
            i += 1
            if i == numrows:
                i = r
                lead += 1
                if numcols == lead:
                    return
        M[i],M[r] = M[r],M[i]
        lv = M[r][lead]
        M[r] = [mrx / lv for mrx in M[r]]
        for i in range(numrows):
            if i != r:
                lv = M[i][lead]
                M[i] = [iv - lv*rv for rv,iv in zip(M[r],M[i])]
        lead += 1

def disk_matrix(size, radius):
    matrix = np.zeros(size)
    
    center = [int(size[0]/2),int(size[1]/2)]
    
    for i in range(size[0]):
        for j in range(size[1]):
            if ((i-center[0])**2 + (j-center[1])**2)**0.5 < radius:
                matrix[i][j] = 1

    return matrix

def gaussian2d(size, sigma):
    matrix = np.zeros(size)
    
    center = [int(size[0]/2),int(size[1]/2)]

    for i in range(size[0]):
        for j in range(size[1]):
            
            matrix[i][j] = np.exp(-1.*((i-center[0])**2 +(j-center[1])**2)\
                                       /(2*sigma**2))

    return matrix

def flux_test():

    g = gaussian2d((100,100),5)
    dm = disk_matrix((100,100),20)

    (_, g1d) = ADE.ADE_gauss(100,50,5,NORM=1)
    box = np.zeros(100)
    box[40:60] = 1

    l1d = ADE.ADE_lorentz(100,50,5,NORM=1)[1]

    g /= np.sum(g)

    ring = disk_matrix((100,100),21) - disk_matrix((100,100),19)
    colors = np.empty(ring.shape, dtype=str)

    for i in range(100):
        for j in range(100):
            if ring[i][j] == 1: colors[i][j] = 'b'
            else: colors[i][j] = 'w'
    
    theta = np.arccos((np.arange(100,dtype='float') - 50)/10)
    pn = 10*np.nan_to_num(np.sin(theta))


    print dm.sum()
    print dm.max()

    c2d = np.fft.ifft2(np.fft.fft2(g)*np.fft.fft2(dm)).real
    c2d = np.roll(np.roll(c2d,50,axis=0),50,axis=1)

    c1d = np.fft.irfft(np.fft.rfft(g1d)*np.fft.rfft(box)).real
    c1d = np.roll(c1d,50)

    pnc = np.fft.irfft(np.fft.rfft(g1d)*np.fft.rfft(pn)).real
    pnc = np.roll(pnc,50)

    bl = np.fft.irfft(np.fft.rfft(l1d)*np.fft.rfft(box)).real
    bl = np.roll(bl,50)

    pl = np.fft.irfft(np.fft.rfft(l1d)*np.fft.rfft(pn)).real
    pl = np.roll(pl,50)

    print c2d.sum()
    print c2d.max()

    x,y = np.meshgrid(np.arange(100),np.arange(100))


    fig0 = plt.figure(0)
    

    ax0 = fig0.add_subplot(222, projection='3d')
    ax0.clear()
    ax0.plot_surface(x,y,c2d,facecolors=colors,rstride=2,cstride=2,alpha=0.9)

#    fig1 = plt.figure(1)
    ax1 = fig0.add_subplot(221, projection='3d')
    ax1.clear()
    ax1.plot_wireframe(x,y,dm,cstride=5,rstride=5)
    
    r = np.arange(100)

    ax2 = fig0.add_subplot(223)
    ax2.clear()
    ax2.plot(r,box,label="original")
    ax2.plot(r,bl,label="with lorentzian")
    ax2.plot(r,c1d,label='with gaussian')
    ax2.plot(r,g1d,label='gaussian kernel')
    ax2.plot(r,l1d,label='lorentzian kernel')
    ax2.legend(bbox_to_anchor=(1.1,1.05))

    ax3 = fig0.add_subplot(224)
    ax3.clear()
    ax3.plot(r,pn,label='original')
    ax3.plot(r,pl,label='with lorentzian')
    ax3.plot(r,pnc,label='with gaussian')
    ax3.plot(r,g1d,label="gaussian kernel")
    ax3.plot(r,l1d,label='lorentzian kernel')
    ax3.legend(bbox_to_anchor=(1.1,1.05))

    fig0.show()

    return
