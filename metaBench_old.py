import numpy as np
import ADEUtils as ADE
import pyfits
import matplotlib.pyplot as plt
import time
import glob

num_ap = 300

def sauron(direct_image, widths_list):
    
    wf = np.loadtxt(widths_list,usecols=(0,3))
    rang = wf[:,0] - 1.0247
    zidx = np.where(rang >= 0.0)[0][0]
    ang = rang[zidx:]
    nwid = wf[:,1][:zidx]
    pwid = wf[:,1][zidx+1:]
    wid = wf[:,1][zidx:]

    a, rings = disc2rings(direct_image,ang)
    
    length = rings.shape[1]

    kk = kernel_gen(wid,ang,a)
    smear = ring_smear(rings,kk)
#    smear = 0
    comb = np.sum(smear,axis=0)
    
    return (a,rings,kk,smear,comb,wid,ang)

def disc2rings(image,angles):

    HDU = pyfits.open(image)[0]
    data = np.float32(HDU.data)
    FL = HDU.header['FOCALLEN']
    length = num_ap

    r, f = ADE.annulize(data,length)

    a = np.arctan(r*24e-3/FL)*180/np.pi

    ring_stack = np.zeros((1,length),dtype=np.float32)

    sangles = np.sort(angles)
#    sangles *= np.pi/180

    for i in range(sangles.size):

        a_mid = sangles[i]

        if i == 0: a0_mid = -1*sangles[i]
        else: a0_mid = sangles[i-1]

        try: a2_mid = sangles[i+1]
        except IndexError:a2_mid = a_mid + (a_mid - a0_mid)
        
        a1 = (a_mid + a0_mid)/2
        a2 = (a_mid + a2_mid)/2

        idx = np.where((a > a1) & (a <= a2))
        counts = np.mean(f[idx])

        fi = np.zeros(length)
        
        fi[idx] = counts

        ring_stack = np.vstack((ring_stack,np.array([fi])))
        
    return (a,ring_stack[1:])

def kernel_gen(widths,angles,a):
    FL = 50.0
    length = a.size*4

    k_stack = np.zeros((1,length),dtype=np.float32)

    for i in range(widths.size):

#        dtheta = np.arctan(widths[i]*24e-3/FL)*180/np.pi
        dtheta = widths[i]
        sigma = dtheta/2.35482
        gamma = dtheta/2

        g = np.exp(-1*(a - angles[i])**2/(2*sigma**2))
        l = (gamma**2/((a - angles[i])**2 + gamma**2))
        g /= np.sum(g)
        l /= np.sum(l)
        
        k = np.fft.irfft(np.fft.rfft(g, n=length) * \
                             np.fft.rfft(l, n=length)).real
        roll = np.where(k == np.max(k))[0][0]
        k = np.roll(k,-1*roll/2)

        k_stack = np.vstack((k_stack,np.array([k])))

    return k_stack[1:,0:a.size]

def ring_smear(rings, kernels):

    length = rings.shape[1]

    FFT = np.fft.rfft
    IFFT = np.fft.irfft

    pad = max(rings.shape[1],kernels.shape[1])*2

    conv =FFT(rings,n=pad) * FFT(kernels,n=pad)

#    return np.roll(IFFT(conv).real,-1*roll,axis=1)
    
    inv = IFFT(conv).real[:,0:length]

    for i in range(inv.shape[0]):
        idx = np.where(inv[i] == np.max(inv[i]))[0][0]
        inv[i] = np.roll(inv[i],-1*idx/2)
        
    return inv

def plot_data(im_name,title,fnum=1):
    FL = 50.0

    data = pyfits.open(im_name)[0].data

    r,f = ADE.annulize(data,num_ap)

    f /= np.sum(f)

    pangles = np.arctan(r*24e-3/FL)*180/np.pi

    fig = plt.figure(fnum)
    plt.clf()
    ax = fig.add_subplot('111')
    ax.plot(pangles,f)
    ax.set_xlabel('Angle [deg]')
    ax.set_ylabel('Normalized Power')
    ax.set_xlim(0,18)
    ax.set_title(title)
    fig.show()

def plot_model(rings,a,title,clear,fnum):
    FL = 50

    if len(rings.shape) == 1: rings = np.array([rings])

    fig = plt.figure(fnum)
    if clear: plt.clf()
    ax = fig.add_subplot('111')
    ax.set_xlabel('Angle [deg]')
    ax.set_ylabel('Power')
    ax.set_title(title)

    for i in range(rings.shape[0]):
        ax.plot(a,rings[i])

    plt.suptitle(time.asctime(time.localtime()))

    fig.show()
    return

class theModule:
    
    def __init__(self, directory, datafile, nonlinmin, 
                 nonlinmax, exclude=[], threshold=0.00001):
        
        self.data = {}
        self.dir = directory
        self.file = datafile
        self.offset = \
            self.find_offset(datafile,nonlinmin,nonlinmax,exclude,threshold)
        
        file_list = glob.glob(directory+'/full*.dat')
        
        for dat_file in file_list:
            self.load_data(dat_file)

        self.angles = self.data.keys()

    def __repr__(self):
        return 'Module created using full data from '+self.dir\
            +'\nand reduced data from '+self.file

    def find_offset(self,datafile, nonlinmin, nonlinmax, exclude, threshold):

        input_a, output_a = np.loadtxt(datafile,usecols=(0,1),unpack=True)
        
        for e in exclude:
            did = np.where(input_a == e)
            output_a = np.delete(output_a, did)
            input_a = np.delete(input_a, did)

        pidx = np.where(input_a > nonlinmax)
        nidx = np.where(input_a < nonlinmin)
        
        in_a = np.append(input_a[nidx],input_a[pidx])
        out_a = np.append(-1*output_a[nidx],output_a[pidx])
        e = np.zeros(in_a.size)+1

        b = 1000.
        offset = 0.
        while abs(b) > threshold:
            m, b = ADE.fit_line(in_a,out_a,e)
            offset += b
            in_a += b

        return offset


    def load_data(self,dat_file):
        
        f = open(dat_file,'r')
        angle = float(f.readlines()[1][15:-1]) + self.offset
        
        f = open(dat_file,'r')        
        a,p = np.loadtxt(f,unpack=True)

        self.data[abs(angle)] = (a,p)

    def get_profile(self, angle):
        
        ang1, ang2, weight1, weight2 = self.find_braket(angle)

        print ang1, ang2
        
        peak_a1, peak_a2 = self.find_peaks(ang1,ang2)

        a1, p1 = self.data[ang1]
        a2, p2 = self.data[ang2]

        a1_z = a1 - peak_a1
        a2_z = a2 - peak_a2

        interpp = np.interp(a1_z, a2_z, p2)

        outputp = weight1*p1 + weight2*interpp

        shift = weight1*(peak_a1) + weight2*(peak_a2)
    
        outputa = a1_z + shift
                  

        return outputa, outputp

    def test(self,angle):

        ang1, ang2, _,_ = self.find_braket(angle)
        fig = plt.figure(1)
#        plt.clf()
        plt.plot(self.data[ang1][0],self.data[ang1][1])
        plt.plot(self.data[ang2][0],self.data[ang2][1])
        
        a, p = self.get_profile(angle)
        plt.plot(a,p,'--')
        fig.show()

    def find_braket(self, angle):
        
        i = 0
        
        k = self.data.keys()
        k.sort()

        while k[i] < angle: 
            i += 1

        ang2 = k[i]
        if ang2 == min(k): ang1 = k[i+1]
        else: ang1 = k[i-1]

        diff = abs(ang2 - ang1)

        weight1 = 1 - abs(angle - ang1)/diff
        weight2 = 1 - abs(ang2 - angle)/diff

        return ang1, ang2, weight1, weight2

    def find_peaks(self, ang1, ang2):
        
        out = []

        for k in (ang1,ang2):
            CDF = np.cumsum(self.data[k][1])
            CDF /= np.max(CDF)

            idx2 = np.where(CDF >= 0.5)[0][0]
            idx1 = idx2 - 1

            pidx = np.interp(0.5,np.array([CDF[idx1],CDF[idx2]]),\
                                 np.array([idx1,idx2]))
            
            aout = np.interp(pidx,np.arange(CDF.size),self.data[k][0])

            out.append(aout)
            
        return out
