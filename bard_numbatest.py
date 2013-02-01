from numba import double
from numba.decorators import jit

import time
import numpy as num

@jit(argtypes=[double[:,:], double[:,:]])
def pairwise_python(X, D):
    M = X.shape[0]
    N = X.shape[1]
    for i in range(M):
        for j in range(M):
            d = 0.0
            for k in range(N):
                tmp = X[i, k] - X[j, k]
                d += tmp * tmp
            D[i, j] = num.sqrt(d)

if __name__ == '__main__':
    X = num.random.random((1000,3))
    D = num.empty((1000,1000))

    t = time.time()
    pairwise_python(X,D)
    print 'time = {}'.format(time.time()-t)
