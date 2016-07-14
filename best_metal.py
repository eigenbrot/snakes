import pyfits
import numpy as np

def do_pointing(pointing, folder_list, suffix='bestZ'):
    
    coeflist = []
    chi = []
    for i, ff in enumerate(folder_list):
        coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(ff,pointing)
        print coeffile
        coefs = pyfits.open(coeffile)[1].data
        coeflist.append(coefs)
        chi.append(coefs['CHISQ'])

    chiarr = np.vstack(chi)
    bestmonochi = np.min(chiarr,axis=0)
    zidx = np.where(chiarr == bestmonochi)[0]

    numap = coefs.shape[0]
    outrec = np.zeros(numap, dtype=coefs.dtype)
    print numap, zidx.shape
    for i in range(numap):
        print i, zidx[i]
        if coeflist[zidx[i]].dtype != coefs.dtype:
            tmp = force_dtype(coefs.dtype, coeflist[zidx[i]].take(i))
            outrec[i] = tmp
        else:
            outrec[i] = coeflist[zidx[i]].take(i)

    pyfits.BinTableHDU(outrec).writeto('NGC_891_P{}_bin30_allz2.{}.fits'.\
                                       format(pointing, suffix), clobber=True)

    return

def do_all_pointings(folder_list = ['0.005Z', '0.02Z', '0.2Z', '0.4Z', '1Z', '2.5Z'], suffix='bestZ'):
    
    for p in range(6):
        do_pointing(p+1,folder_list,suffix)

    return

def force_dtype(dtype, data):
    #Ignores light_weight and other things that might be different
    #shapes due to multi/mono differences. All the important quantites
    #are retained.

    output = np.zeros(1, dtype=dtype)[0]
    for n in output.dtype.names:
        if output[n].shape == data[n].shape:
            output[n] = data[n]

    return output
