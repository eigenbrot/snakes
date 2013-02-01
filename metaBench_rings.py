import numpy as np
import ADEUtils as ADE
import pyfits
import matplotlib.pyplot as plt

def sauron(direct_image, widths_list):
    
    wf = np.loadtxt(widths_list,usecols=(0,2))
    rang = wf[:,0] - 1.0247
    zidx = np.where(rang >= 0.0)[0][0]
    ang = rang[zidx:]
    nwid = wf[:,1][:zidx]
    pwid = wf[:,1][zidx+1:]
    wid = wf[:,1][zidx:]/2

    rings = disc2rings(direct_image,ang)
    
    dims = (rings.shape[1]*2,rings.shape[2]*2)

    kk = kernel_gen(wid,dims)
    smear = ring_smear(rings,kk)
#    smear = 0
    comb = np.sum(smear,axis=0)
    
    return (rings,kk,smear,comb,wid,ang)

def disc2rings(image,angles):

    HDU = pyfits.open(image)[0]
    data = np.float32(HDU.data)
    FL = HDU.header['FOCALLEN']
    print FL
    dims = data.shape

    cent = ADE.centroid(data)
    dist = ADE.dist_gen(data,cent)
    
    ring_stack = np.zeros((1,dims[0],dims[1]),dtype=np.float32)

    sangles = np.sort(angles)
    sangles *= np.pi/180

    for i in range(sangles.size):

        r_mid = FL*np.tan(sangles[i])/24e-3

        if i == 0: r0_mid = -1*FL*np.tan(sangles[i])/24e-3
        else: r0_mid = FL*np.tan(sangles[i-1])/24e-3

        try:r2_mid = FL*np.tan(sangles[i+1])/24e-3
        except IndexError:r2_mid = r_mid + (r_mid - r0_mid)
        
        r1 = (r_mid + r0_mid)/2
        r2 = (r_mid + r2_mid)/2

        print r1, r2, r2 - r1

        idx = np.where((dist > r1) & (dist <= r2))

        temp = np.zeros(dims,dtype=np.float32)
        temp[idx] = data[idx]
        ring_stack = np.vstack((ring_stack,np.array([temp])))

    return ring_stack[1:]

def kernel_gen(widths,dims):
    
    k_stack = np.zeros((1,dims[0],dims[1]),dtype=np.float32)
    cent = (dims[0]/2,dims[1]/2)
    print cent

    for i in range(widths.size):
        k = ADE.ADE_gauss2d(dims,cent,0,FWHM=widths[i]/24e-3)

        k_stack = np.vstack((k_stack,np.array([k])))

    return k_stack[1:]

def ring_smear(rings, kernels):

    FFT = np.fft.rfft2
    IFFT = np.fft.irfft2

    pad = (rings.shape[1]*2, rings.shape[2]*2)

    conv =FFT(rings,s=pad) * FFT(kernels,s=pad)

#    return np.fft.ifftshift(IFFT(conv).real,axes=(1,2))
    return np.roll(np.roll(IFFT(conv).real,-512,axis=1),-512,axis=2)

def plot_smear(rings,title):
    FL = 50.0

    pangles = np.zeros(300)
    powers = np.zeros(300)

    for i in range(rings.shape[0]):
        print i
        
        r,f = ADE.annulize(rings[i],300)
        powers += f

    pangles = np.arctan(r*24e-3/FL)*180/np.pi


#    for i in range(angles.size):
#        r_mid = FL*np.tan(sangles[i])/24e-3
#
 #       if i == 0: r0_mid = -1*FL*np.tan(sangles[i])/24e-3
#        else: r0_mid = FL*np.tan(sangles[i-1])/24e-3
#
#        try:r2_mid = FL*np.tan(sangles[i+1])/24e-3
#        except IndexError:r2_mid = r_mid + (r_mid - r0_mid)
#        
#        r1 = (r_mid + r0_mid)/2
#        r2 = (r_mid + r2_mid)/2
#
#        #for plotting
#        a1 = np.arctan(r1*24e-3/FL)*180/np.pi
#        a2 = np.arctan(r2*24e-3/FL)*180/np.pi
#        pangles = np.append(pangles,np.linspace(a1,a2,num=100,endpoint=False))
#        powers = np.append(powers,np.zeros(100)+np.sum(rings[i]))


    fig = plt.figure(3)
    plt.clf()
    ax = fig.add_subplot('111')
    ax.plot(pangles,powers)
    ax.set_xlabel('Angle [deg]')
    ax.set_ylabel('Power')
    ax.set_xlim(0,18)
    ax.set_title(title)
    fig.show()
