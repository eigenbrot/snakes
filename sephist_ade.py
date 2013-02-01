#! /usr/bin/python
# Find the minimum distance between each detection and its neighbors,
# plot result as a histogram
#
# Written by Jenna Ryon
# 09 Nov 2012
#
# Optimized by A. Eigenbrot
# 12 Nov 2012

import numpy as np
#import scipy
#from pylab import *
#import math
import time

def do_test(num):
    #x_coords,y_coords = np.loadtxt('cats/src_336.tab',comments='#',usecols=(0,1),unpack=True) # read in phot results as arrays
    x_coords = np.array(np.random.rand(num)*500.,dtype=np.float32)
    y_coords = np.array(np.random.rand(num)*500.,dtype=np.float32)

    j_t1 = time.time()
    length = len(y_coords)
    a = range(length)
    
    j_min_distances = []
    
    for i in a:
        indiv_dist = []
        for j in a:
            if i != j:
                indiv_dist.append(np.sqrt((x_coords[i]-x_coords[j])**2+(y_coords[i]-y_coords[j])**2))
        j_min_distances.append(np.amin(indiv_dist))
                
    j_t2 = time.time()
    
    ####################
    a_t1 =time.time()
    
    '''the padding is necessary to allow us to take the transpose.
    xpad is just [x]'''
    xpad = x_coords[None,:]
    
    '''xpad - xpad.T produces a 2D array with the difference between every
    combination of points in xpad. We then wrap in in a masked array.
    It's a tad bit faster to square these now.'''
    xdiff = np.ma.array(xpad - xpad.T)**2
    '''now mask out the diagonal values because those will be zero'''
    xdiff[np.diag_indices(xdiff.shape[0])] = np.ma.masked

    ypad = y_coords[None, :]
    ydiff = np.ma.array(ypad - ypad.T)**2
    ydiff[np.diag_indices(ydiff.shape[0])] = np.ma.masked

    '''find the minima. The .data makes this a standard numpy
    array rather than a masked array'''
    a_min_distances = np.min(np.sqrt(xdiff + ydiff),axis=1).data
    
    a_t2 = time.time()

    #####################
    
    avgdiff = np.mean(np.abs(j_min_distances - a_min_distances))
    stddiff = np.std(np.abs(j_min_distances - a_min_distances))
    
    print "Loop method took {:4.4f} seconds\nNumpy method took {:4.4f} seconds".format(j_t2 - j_t1, a_t2 - a_t1)
    print "Results vary by an average of {:4.4f} with a std of {:4.4f}".format(avgdiff,stddiff)
    
    return 

if __name__ == '__main__':
    do_test(5000)
