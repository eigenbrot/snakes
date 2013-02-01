from numba.decorators import jit
from numba import float32, int16, float64
import numpy as np
import pyfits
import ADEUtils as ADE
import time
import bottleneck as bn

@jit(argtypes=[float32[:,:],float32[:,:],float32,int16],restype=float32[:,:],backend='ast')
def numba_cent(data, distances, maxd, numan):

    size = data.shape
    rstep = maxd / numan
    
    r1 = 0.0
    r2 = rstep
    
    stdarr = np.zeros(numan,dtype=np.float32)
    rarr = np.zeros(numan,dtype=np.float32)
    outarr = np.zeros(numan,dtype=np.float32)

    for k in range(numan):
        anlist = []
        for i in range(size[0]):
            for j in range(size[1]):
                if distances[i,j] > r1:
                    if distances[i,j] <= r2:
                        anlist.append(data[i,j])
#                    outarr[k] += data[i,j]

        anarray = np.array(anlist,dtype=np.float32)
        outarr[k] = bn.nansum(anarray)
        stdarr[k] = bn.nanstd(anarray)
        rarr[k] = (r1 + r2)*0.5
        r1 = r2
        r2 += rstep
        
    return np.array([rarr,outarr,stdarr])

@jit(argtypes=[int16],restype=float32[:])
def test(num):
    
    l = np.array([],dtype=np.float32)
    for i in range(num):
        l = np.append(l,i)
    
#    x = np.array(l,dtype=np.float32)
    return l

print 'Running test...'
data = pyfits.open('scratch/test.fits')['SB'].data
ddata = np.array(data,dtype=np.float32)
cent = ADE.centroid(ddata)
dist = ADE.dist_gen(ddata,cent)

at1 = time.time()
aderes = ADE.annulize(ddata,300)
at2 = time.time()

nt1 = time.time()
nres = numba_cent(ddata,dist,dist.max(),300)
nt2 = time.time()

print np.mean(nres - aderes,axis=1)
print np.std(nres - aderes,axis=1)

print 'ADE time was: {:4.5f} s\nNumba time was: {:4.5f} s'.format(at2 - at1, nt2-nt1)


# numba_cent(24)
# from numba import *

# @jit(restype=float32, argtypes=[float32[:, :]],backend='bytecode')
# def sum1d(my_double_array):
#     nsum = 0.0
#     for i in range(my_double_array.shape[0]):
#         nsum += my_double_array[i]
#     return nsum
