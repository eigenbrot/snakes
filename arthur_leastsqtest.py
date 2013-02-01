#!/usr/bin/python2.7

import numpy as np
import math
from scipy.optimize import leastsq
import sys
import numpy.random

def errfunc(p,x,y,err):
    print "\t{}".format(p)

    fitfunc = p[0]*x + p[1]
    print np.sum(((y-fitfunc)/err)**2,axis=0)
    return (y-fitfunc)/err

if __name__ == '__main__':

    #Some initial parameters defining the function:
    truea = 2.58
    trueb = 9.72
    randsig = 0.3
    numpoints = 50
    xinterval = (0,10)

    #Create numpy arrays for x and y, with error:
    xvals = np.linspace(xinterval[0],xinterval[1],numpoints)
    errors = np.random.normal(loc=0.0,scale=randsig,size=len(xvals))
    yvals = truea*xvals+trueb+errors
    errorvals = xvals*0 + randsig

    #Set initial guesses:
    pinit = [5,-50]#[0] is a, [1] is b
    
    #Make the functions that return residuals:
    #fitfunc = lambda p, x: p[0]*x + p[1]
    #errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err

    #Run the LM least-squares:
    fitoutput = leastsq(errfunc,pinit,args=(xvals,yvals,errorvals),full_output=1)

    print "True function: {0:.2f}x{1:+.2f}".format(truea,trueb)
    print "Initial guess: {0:.2f}x{1:+.2f}".format(pinit[0],pinit[1])
    print "Fitted function: {0:.2f}x{1:+.2f}, Number of function calls = {2:d}".format(fitoutput[0][0],fitoutput[0][1],fitoutput[2]['nfev'])
