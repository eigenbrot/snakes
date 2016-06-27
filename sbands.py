import numpy as np
import pyfits
import time

def read_fits(fitsfile):
    
    hdu = pyfits.open(fitsfile)[0]
    data = hdu.data

    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CDELT1'] + hdu.header['CRVAL1']

    return wave, data

def compute_Dn4000(wave, spec, err):

    redidx = np.where((wave >= 3850.) & (wave <= 3950.))[0]
    blueidx = np.where((wave >= 4000.) & (wave <= 4100.))[0]
    
    red = np.mean(spec[redidx])
    blue = np.mean(spec[blueidx])

    red_e = np.sqrt(np.mean(err[redidx]**2))
    blue_e = np.sqrt(np.mean(err[blueidx]**2))

    Dn4000 = blue/red
    Dn4000_e = np.sqrt((blue_e/red)**2 + (red_e*blue/red**2)**2)

    return Dn4000, Dn4000_e

def compute_index(wave, spec, err, centlim, bluelim, redlim):
    
    centidx = np.where((wave > centlim[0]) & (wave < centlim[1]))[0]
    redidx = np.where((wave > redlim[0]) & (wave < redlim[1]))[0]
    blueidx = np.where((wave > bluelim[0]) & (wave < bluelim[1]))[0]

    centwave = np.r_[centlim[0], wave[centidx], centlim[1]]
    redwave = np.r_[redlim[0], wave[redidx], redlim[1]]
    bluewave = np.r_[bluelim[0], wave[blueidx], bluelim[1]]

    # cent = np.mean(spec[centidx])
    # blue = np.mean(spec[blueidx])
    # red = np.mean(spec[redidx])
    cent = np.mean(np.interp(centwave, wave, spec))
    blue = np.mean(np.interp(bluewave, wave, spec))
    red = np.mean(np.interp(redwave, wave, spec))
    
    cent_e = np.sqrt(np.mean(err[centidx]**2))
    blue_e = np.sqrt(np.mean(err[blueidx]**2))
    red_e = np.sqrt(np.mean(err[redidx]**2))

    redcent = np.mean(redwave)
    bluecent = np.mean(bluewave)

    centcent = np.mean(centlim)
    dlambda = centlim[1] - centlim[0]

    cont = np.interp(centcent, [bluecent, redcent], [blue,red])
    index = (1 - cent/cont)*dlambda
    
    cont_e = np.interp(4103.9, [bluecent, redcent], [blue_e, red_e])
    err = np.sqrt((dlambda*cent_e/cont)**2 + (dlambda*cent*cont_e/cont**2)**2)

    return index, err

def main(fitsfile, errfile, output):
    
    #[[cent], [blue], [red]]
    HdAlims = [[4084.5, 4123.3], [4041.65, 4079.75], [4128.55, 4161.05]]
    Mgblims = [[5160.15, 5192.65], [5142.6, 5161.4], [5191.4, 5206.4]]
    Fe5270lims = [[5245.7, 5285.7], [5233.2, 5248.2], [5285.65, 5318.15]]
    Fe5335lims = [[5312.1, 5352.1], [5304.65, 5315.95], [5353.4, 5363.4]]

    wave, data = read_fits(fitsfile)
    _, err = read_fits(errfile)
    Dn4000 = np.zeros(data.shape[0])
    Dn4000_e = np.zeros(data.shape[0])

    HdA = np.zeros(data.shape[0])
    HdA_e = np.zeros(data.shape[0])

    Mgb = np.zeros(data.shape[0])
    Mgb_e = np.zeros(data.shape[0])

    Fe = np.zeros(data.shape[0])
    Fe_e = np.zeros(data.shape[0])

    MgFe = np.zeros(data.shape[0])
    MgFe_e = np.zeros(data.shape[0])

    for i in range(data.shape[0]):
        D, De = compute_Dn4000(wave,data[i,:],err[i,:])
        Dn4000[i] = D
        Dn4000_e[i] = De

        H, He = compute_index(wave, data[i,:], err[i,:],
                              HdAlims[0], HdAlims[1], HdAlims[2])
        HdA[i] = H
        HdA_e[i] = He

        M, Me = compute_index(wave, data[i,:], err[i,:],
                                  Mgblims[0], Mgblims[1], Mgblims[2])
        F2, F2e = compute_index(wave, data[i,:], err[i,:],
                                  Fe5270lims[0], Fe5270lims[1], Fe5270lims[2])
        F3, F3e = compute_index(wave, data[i,:], err[i,:],
                                  Fe5335lims[0], Fe5335lims[1], Fe5335lims[2])
        
        Mgb[i] = M
        Mgb_e[i] = Me

        Fe[i] = 0.5*(F2 + F3)
        Fe_e[i] = np.sqrt((F2e/2.)**2 + (F3e/2.)**2)

        MgFe[i] = np.sqrt(M * (0.72*F2 + 0.28*F3))
        Mterm = (0.72*F2 + 0.28*F3)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        F2term = (0.72*M)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        F3term = (0.28*M)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        MgFe_e[i] = np.sqrt((Mterm*Me)**2 + (F2term*F2e)**2 + (F3term*F3e)**2)


    fmt = '{:8.3f}'*10 + '\n'
    with open(output,'w') as f:
        f.write('# Generated on {}\n'.format(time.asctime()))
        f.write('# From {} and {}\n'.format(fitsfile, errfile))
        f.write(('{:>8}'*10+'\n').format('#Dn4000',
                                        'Dn4000_e',
                                        'HdA','HdA_e',
                                        'Mgb','Mgb_e',
                                        '<Fe>','<Fe>_e',
                                        'MgFe','MgFe_e'))
        for i in range(HdA.size):
            f.write(fmt.format(Dn4000[i], Dn4000_e[i], 
                               HdA[i], HdA_e[i],
                               Mgb[i], Mgb_e[i],
                               Fe[i], Fe_e[i],
                               MgFe[i], MgFe_e[i]))


    return

def do_all():

    for p in range(6):
        output = 'NGC_891_P{}_bin30.msoz.spy.dat'.format(p+1)
        main('NGC_891_P{}_bin30.msoz.fits'.format(p+1),
             'NGC_891_P{}_bin30.meoz.fits'.format(p+1),output)

    return
             
