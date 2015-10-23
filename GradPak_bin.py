import time
import numpy as np
import pyfits
import GradPak_plot as GPP

def bin(datafile, errfile, SNR, outputfile, waverange=None, exclude=[]):

    hdu = pyfits.open(datafile)[0]
    data = hdu.data
    err = pyfits.open(errfile)[0].data

    if waverange is not None:
        wave = (np.arange(data.shape[1]) + hdu.header['CRPIX1'] - 1)\
               *hdu.header['CDELT1'] + hdu.header['CRVAL1']
        waveidx = np.where((wave >= waverange[0]) & (wave <= waverange[1]))[0]
        # data = data[:,waveidx]
        # err = err[:,waveidx]
    else:
        waveidx = None

    data *= 1e17
    err *= 1e17

    y_values = np.array([c.center[1] for c in GPP.GradPak_patches()[:,1]])
    x_values = np.array([c.center[0] for c in GPP.GradPak_patches()[:,1]])
    fibnums = np.arange(109) + 1
    row_pos = np.unique(y_values)

    binf = np.zeros(data.shape[1])
    bine = np.zeros(data.shape[1])

    fibdict = {}
    binnum = 1
    for i in range(row_pos.size):

        if row_pos[i] > 80:
            continue

        idx = np.where(y_values == row_pos[i])[0]
        b = 0
        n = 0

        while fibnums[idx[n]] in exclude:
            print 'Skipping fiber {}'.format(fibnums[idx[n]])
            n += 1

        while n < len(idx):
            try:
                while fibnums[idx[n]] in exclude:
                    print 'Skipping fiber {}'.format(fibnums[idx[n]])
                    n += 1
            except IndexError: #The rest of the fibers in the row are excluded
                break

            fstack = data[idx[n]]
            estack = err[idx[n]]
            tmp = compute_SN(data[idx[n]], err[idx[n]], waveidx)
            snstack = tmp
            fibers = [fibnums[idx[n]]]
            xpos = [x_values[idx[n]]]
            ypos = [y_values[idx[n]]]

            while tmp < SNR:
                n += 1
                print 'fibers: {}, SNR: {}'.format(fibers, tmp)
                if n > len(idx) - 1:
                    print "WARNING, SN threshold not met in row {}, bin {}".\
                        format(i,b)
                    break

                if fibnums[idx[n]] in exclude:
                    print 'Skipping fiber {}'.format(fibnums[idx[n]])
                    continue
                    
                fstack = np.vstack((fstack, data[idx[n]]))
                estack = np.vstack((estack, err[idx[n]]))
                snstack = np.vstack((snstack, compute_SN(data[idx[n]], 
                                                         err[idx[n]], 
                                                         waveidx)))
                tmpf, tmpe = create_bin(fstack, estack, snstack)
                tmp = compute_SN(tmpf, tmpe, waveidx)
                fibers.append(fibnums[idx[n]])
                xpos.append(x_values[idx[n]])
                ypos.append(y_values[idx[n]])

            print 'binned aperture {}: {}, SNR: {}'.format(binnum,fibers, tmp)
            bin_x_pos = np.mean(xpos)
            bin_y_pos = np.mean(ypos)
            fibstr = [str(i) for i in fibers]
            hdu.header.update('BIN{:03}F'.format(binnum),' '.join(fibstr))
            hdu.header.update('BIN{:03}P'.format(binnum),' '.\
                              join([str(bin_x_pos),str(bin_y_pos)]))
            binf = np.vstack((binf,tmpf))
            bine = np.vstack((bine,tmpe))
            fibdict['{}_{}'.format(i,b)] = fibers
            b += 1
            n += 1
            binnum += 1

    binf = binf[1:]/1e17
    bine = bine[1:]/1e17
    pyfits.PrimaryHDU(binf, hdu.header).\
        writeto('{}.ms.fits'.format(outputfile),clobber=True)
    pyfits.PrimaryHDU(bine, hdu.header).\
        writeto('{}.me.fits'.format(outputfile),clobber=True)
    return binf, bine, fibdict

def create_bin(fstack, estack, snstack):

    sumW = np.sum(snstack**2)

    newbin = np.sum(fstack*snstack**2,axis=0)/sumW
    newerr = np.sqrt(np.sum((estack*snstack**2/sumW)**2,axis=0))

    return newbin, newerr

def compute_SN(signal, noise, idx=None):

    zidx = np.where(noise[idx] != 0)

    return np.mean(signal[idx][zidx]/(noise[idx][zidx]))

def create_locations(binfile, galcenter=[35.637962,42.347629], 
                     ifucenter=[35.637962,42.347629], reffiber=105,
                     pa=295.787, kpc_scale=0.0485):

    hdu = pyfits.open(binfile)[0]
    numaps = hdu.data.shape[0]
    binhead = hdu.header
    
    patches = GPP.get_binned_patches(binhead)
    refpatches = GPP.GradPak_patches()
    
    patches, refpatches = GPP.transform_patches(patches,refpatches=refpatches,
                                                pa=0, center=ifucenter, 
                                                reffiber=reffiber)

    decrad = galcenter[1]*2*np.pi/360.
    parad = pa*2*np.pi/360.
    ra_diff = 3600*(galcenter[0] - ifucenter[0])*np.cos(decrad)
    dec_diff = 3600*(galcenter[1] - ifucenter[1])
    print ra_diff, dec_diff, np.sqrt(ra_diff**2 + dec_diff**2)

    refr_diff = ra_diff*np.cos(parad) - dec_diff*np.sin(parad)
    refz_diff = -1*(ra_diff*np.sin(parad) + dec_diff*np.cos(parad))

    print refz_diff, refr_diff

    reffiber_r, reffiber_z = refpatches[reffiber-1,1].center

    f = open('{}_locations.dat'.format(binfile.split('.ms.fits')[0]),'w')
    f.write("""# Generated on {}
# Inpute file: {}
#
""".format(time.asctime(),binfile))
    f.write('# {:4}{:>10}{:>10}{:>10}{:>10}{:>10}\n#\n'.format('Apnum',
                                                              'size (")',
                                                              'r (")',
                                                              'z (")',
                                                              'r (kpc)',
                                                              'z (kpc)'))
            
    for i, p in enumerate(patches[:,1]):

        fibers = binhead['BIN{:03}F'.format(i+1)]
        radius = refpatches[int(fibers.split(' ')[0]) - 1][1].get_radius()

        r_diff = 3600*(p.center[0] - reffiber_r) - refr_diff
        z_diff = 3600*(p.center[1] - reffiber_z) + refz_diff
        r_diff *= -1
        print i, p.center, r_diff, z_diff

        f.write(str('{:7n}'+5*'{:10.3f}'+'\n').format(i,
                                                             radius,
                                                             r_diff,
                                                             z_diff,
                                                             r_diff*kpc_scale,
                                                             z_diff*kpc_scale))
    f.close()
    return

def plot_test():

    import matplotlib.pyplot as plt
    for SN in [0, 20, 60]:

        flux, err, _ = bin('../NGC_891_P1_final.ms_rfsz_lin.fits',
                           '../NGC_891_P1_final.me_rfz_lin.fits',
                           SN, 'NGC_891_P1_bin{}'.format(SN), 
                           waverange=[5450,5550])

        ax = plt.figure().add_subplot(111)
        ax.plot(np.arange(flux.shape[1]), flux[3])
        ax.set_title('SN > {}'.format(SN))
        ax.set_ylim(-0.2e-15,0.8e-15)
        ax.figure.show()



