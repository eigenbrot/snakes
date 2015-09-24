import sys
import pyfits
import numpy as np


def do_single(input_image, output_image):

    h = pyfits.open(input_image)[0]
    d = h.data
    sr = np.zeros(d.shape,dtype=d.dtype)
    sr[:,0] = np.mean(d[:,0:70], axis=1)
    sr[:,-1] = np.mean(d[:,1260:], axis=1)
    
    x = np.arange(sr.shape[1])
    xp = np.array([x[0],x[-1]])
    fp = np.array([sr[:,0], sr[:,-1]])
    
    l = [np.interp(x,xp,fp[:,i]).astype(d.dtype) for i in range(fp.shape[1])]
    scat = np.vstack(l)
    
    pyfits.PrimaryHDU(d - scat, h.header).writeto(output_image)

    return

def do_multi(inputs, outputs):

    with open(inputs,'r') as fi:
        with open(outputs,'r') as fo:
            input_list = fi.readlines()
            output_list = fo.readlines()
            for i, o in zip(input_list, output_list):
                i = i.replace('\n','')
                o = o.replace('\n','')
                print i, o
                do_single(i,o)

    return

if __name__ == '__main__':
    do_multi(sys.argv[1],sys.argv[2])
