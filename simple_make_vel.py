import pyfits
import time
import sys
from glob import glob

def consolidate(OGcoef_file, velcoef_file, output, offset=74):

    OGcoefs = pyfits.open(OGcoef_file)[1].data
    velcoefs = pyfits.open(velcoef_file)[1].data

    outvel = OGcoefs['VSYS'] + velcoefs['VSYS'] + offset
    
    with open(output,'w') as f:
        f.write('# Generated on {}\n#\n'.format(time.asctime()))
        f.write('# {}\n# {}\n'.format(OGcoef_file, velcoef_file))
        f.write('# Offset = {} km/s\n#\n'.format(offset))
        f.write('#{:3}{:10}\n\n'.format('Ap','V [km/s]'))
        for i, v in enumerate(outvel):
            f.write('{:3n}{:10.3f}\n'.format(i+1,v))

    return

def main(offset=74):

    baseOG = 'NGC_891_P{}_bin30_allz2.coef.fits'
    basevel = 'NGC_891_P{}_bin30_allz2.coef.vel.fits'
    baseout = 'NGC_891_P{}_bin30_velocities.dat'

    for i in range(6):
        print baseOG.format(i+1)
        print basevel.format(i+1)
        print baseout.format(i+1)

        consolidate(baseOG.format(i+1),
                    basevel.format(i+1),
                    baseout.format(i+1),
                    offset=offset)

    return

def shift_data_files(offset=74.):
    import numpy as np

    base = 'NGC_891_P{}_bin30'

    for p in range(6):
        for suff in ['.ms','.me']:
            hdu = pyfits.open(base.format(p+1)+suff+'.fits')[0]
            data = hdu.data
            header = hdu.header

            cdelt = header['CDELT1']
            crpix = header['CRPIX1']
            crval = header['CRVAL1']
    
            wave = (np.arange(data.shape[1]) + crpix-1)*cdelt + crval

            shift = np.vstack([np.interp(wave,wave*(1 + offset/3e5),data[i,:]*1e17) for i in range(data.shape[0])])
            
            new = base.format(p+1)+suff+'o.fits'
            print base.format(p+1)+suff+'.fits', '-->', new
            pyfits.PrimaryHDU(shift/1e17, header).writeto(new,clobber=True)

    return

if __name__ == '__main__':
    try:
        main(offset=float(sys.argv[1]))
    except IndexError:
        main()
