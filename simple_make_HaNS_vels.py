import numpy as np
import pyfits
import time
import sys
from glob import glob

def consolidate(OGcoef_file, velcoef_file, location_file, HaNS_file, output,
                offset=[-56.672, -35.855, -15.033, -14.381, -4.871]):

    OGcoefs = pyfits.open(OGcoef_file)[1].data
    velcoefs = pyfits.open(velcoef_file)[1].data
    size = np.loadtxt(location_file, usecols=(1,), unpack=True)
    HaNS, HaNS_e = np.loadtxt(HaNS_file, unpack=True)

    #First, compute and correct the stellar velocity value
    starvel = OGcoefs['VSYS'] + velcoefs['VSYS']
    for i, s in enumerate([0.937, 1.406, 1.875, 2.344, 2.812]):
        idx = np.where(size == s)
        starvel[idx] += offset[i]

    #Now take the mean w/ the HaNS velocities and compute the uncertainty
    outvel = (starvel + HaNS)/2.0
    veldiff = (starvel - HaNS)/2.0
    velerr = np.sqrt(HaNS_e**2 + veldiff**2)
    
    with open(output,'w') as f:
        f.write('# Generated on {}\n#\n'.format(time.asctime()))
        f.write('# {}\n# {}\n'.format(OGcoef_file, velcoef_file))
        f.write('# Offset (by fiber size) = {} km/s\n'.format(offset))
        f.write('# Everything in km/s\n#\n')
        f.write('#{:2}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n\n'.format('Ap','V_obs','dV_obs', 'V_*', 'V_*^c', 'V_HaNS', 'dV_HaNS'))
        for i in range(outvel.size):
            f.write('{:3n}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n'.format(i+1,outvel[i],velerr[i],
                                                                                     (OGcoefs['VSYS'] + velcoefs['VSYS'])[i],
                                                                                     starvel[i],
                                                                                     HaNS[i], HaNS_e[i]))

    return

def main(offset=[-56.672, -35.855, -15.033, -14.381, -4.871]):

    baseOG = 'NGC_891_P{}_bin30_allz2.coef.prechi.fits'
    basevel = 'NGC_891_P{}_bin30_allz2.coef.vel.fits'
    baseloc = 'NGC_891_P{}_bin30_locations.dat'
    baseHaNS = 'P{}_mab_vel.txt'
    baseout = 'NGC_891_P{}_bin30_velocities.dat'

    for i in range(6):
        print baseOG.format(i+1)
        print basevel.format(i+1)
        print baseout.format(i+1)
        print baseloc.format(i+1)
        print baseHaNS.format(i+1)

        consolidate(baseOG.format(i+1),
                    basevel.format(i+1),
                    baseloc.format(i+1),
                    baseHaNS.format(i+1),
                    baseout.format(i+1),
                    offset=offset)

    return

def shift_data_files(offset=[-56.672, -35.855, -15.033, -14.381, -4.871]):
    import numpy as np

    base = 'NGC_891_P{}_bin30'

    for p in range(6):
        for suff in ['.ms','.me']:
            hdu = pyfits.open(base.format(p+1)+suff+'.fits')[0]
            location = base.format(p+1)+'_locations.dat'
            size = np.loadtxt(location, usecols=(1,), unpack=True)
            data = hdu.data
            header = hdu.header

            cdelt = header['CDELT1']
            crpix = header['CRPIX1']
            crval = header['CRVAL1']
    
            wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval

            newstack = []
            for i in range(data.shape[0]):
                oidx = np.where(np.array([0.937, 1.406, 1.875, 2.344, 2.812]) == size[i])[0][0]
                print i, oidx
                newstack.append(np.interp(wave,wave*(1 + offset[oidx]/3e5),data[i,:]*1e17))
            
            shift = np.vstack(newstack)
            
            new = base.format(p+1)+suff+'o.fits'
            print base.format(p+1)+suff+'.fits', '-->', new
            pyfits.PrimaryHDU(shift/1e17, header).writeto(new,clobber=True)

    return

if __name__ == '__main__':
    try:
        main(offset=float(sys.argv[1]))
    except IndexError:
        main()
