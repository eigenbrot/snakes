#! /usr/bin/env python

import pyfits
import numpy as np
import sys

def saltarith(operand1, op, operand2, result):

    if op == '/': print "Runtime errors might happen when using divide due to zeros in your images. This doesn't appear to affect the output"


    try: hdu_list1 = pyfits.open(operand1)
    except IOError: 
        print 'For now this routine only supports operand1 being a FITs file'

    head_list1 = [hdu_list1[i].header for i in range(7)]
    data_list1 = [hdu_list1[i].data for i in range(7)]
    ready_operand1 = np.dstack(data_list1[1:])

    fits = True
    try: hdu_list2 = pyfits.open(operand2)
    except IOError:
        ready_operand2 = float(operand2)
        fits = False
    
    if fits:
        data_list2 = [hdu_list2[i].data for i in range(7)]
        ready_operand2 = np.dstack(data_list2[1:])

    result_stack = eval('ready_operand1 '+op+' ready_operand2')
    
    result_hdulist = [pyfits.ImageHDU(result_stack[:,:,i],head_list1[1:][i]) for i in range(6)]
    phdu = pyfits.PrimaryHDU(None,head_list1[0])

    output_hdulist = pyfits.HDUList([phdu]+result_hdulist)
    output_hdulist.writeto(result)


def main():
    saltarith(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

if __name__ == '__main__':
    sys.exit(main())
