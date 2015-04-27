import numpy as np
import pyfits
import GradPak_plot as GPP

def bin(datafile, errfile, SNR, outputfile, waverange=None):

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    err = pyfits.open(errfile)[0].data

    if waverange is not None:
        wave = (np.arange(data.shape[1]) + hdu.header['CRPIX1'])\
               *hdu.header['CDELT1'] + hdu.header['CRVAL1']
        waveidx = np.where((wave >= waverange[0]) & (wave <= waverange[1]))[0]
        data = data[:,waveidx]
        err = err[:,waveidx]

    data *= 1e17
    err *= 1e17

    y_values = np.array([c.center[1] for c in GPP.GradPak_patches()[:,1]])
    fibnums = np.arange(109) + 1
    row_pos = np.unique(y_values)

    binf = np.zeros(data.shape[1])
    bine = np.zeros(data.shape[1])

    fibdict = {}
    for i in range(row_pos.size):

        if row_pos[i] > 80:
            continue

        idx = np.where(y_values == row_pos[i])[0]
        b = 0
        n = 0

        while n < len(idx):
            tmpf = data[idx[n]]
            tmpe = err[idx[n]]
            tmp = compute_SN(tmpf, tmpe)
            fibers = [fibnums[idx[n]]]
        
            while tmp < SNR:
                print 'fibers: {}, SNR: {}'.format(fibers, tmp)
                if n > len(idx) - 1:
                    print "WARNING, SN threshold not met in row {}, bin {}".\
                        format(i,b)
                    break
                    
                tmpf, tmpe = add_to_bin(tmpf, tmpe, data[idx[n]], err[idx[n]])
                tmp = compute_SN(tmpf, tmpe)
                fibers.append(fibnums[idx[n]])
                n += 1

            print 'fibers: {}, SNR: {}'.format(fibers, tmp)
            binf = np.vstack((binf,tmpf))
            bine = np.vstack((bine,tmpe))
            fibdict['{}_{}'.format(i,b)] = fibers
            b += 1
            n += 1

    binf = binf[1:]/1e17
    bine = bine[1:]/1e17
    pyfits.PrimaryHDU(binf, hdu.header).\
        writeto('{}.ms.fits'.format(outputfile),clobber=True)
    pyfits.PrimaryHDU(bine, hdu.header).\
        writeto('{}.me.fits'.format(outputfile),clobber=True)
    return binf, bine, fibdict

def add_to_bin(binf, bine, data, error):

    binSNR = compute_SN(binf, bine)
    addSNR = compute_SN(data, error)
    sumSNR = binSNR**2 + addSNR**2

    newbin = (binf * binSNR**2 + data * addSNR**2)/sumSNR
    newerr = np.sqrt(((binSNR**2 * bine)/sumSNR)**2 + \
                     ((addSNR**2 * error)/sumSNR)**2)

    return newbin, newerr

def compute_SN(signal, noise):

    idx = np.where(noise != 0)

    return np.sqrt(np.mean((signal[idx]/noise[idx])**2))

def plot_test():

    import matplotlib.pyplot as plt
    for SN in [0, 5, 60]:

        flux, err, _ = bin('../NGC_891_P1_final.ms_rfsz_lin.fits',
                           '../NGC_891_P1_final.me_rfz_lin.fits',
                           SN, 'NGC_891_P1_bin{}'.format(SN), 
                           waverange=[3750,6800])

        ax = plt.figure().add_subplot(111)
        ax.plot(np.arange(flux.shape[1]), flux[3])
        ax.set_title('SN > {}'.format(SN))
        ax.figure.show()


