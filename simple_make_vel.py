import pyfits
import time
from glob import glob

def consolidate(OGcoef_file, velcoef_file, output):

    OGcoefs = pyfits.open(OGcoef_file)[1].data
    velcoefs = pyfits.open(velcoef_file)[1].data

    outvel = OGcoefs['VSYS'] + velcoefs['VSYS']
    
    with open(output,'w') as f:
        f.write('# Generated on {}\n#\n'.format(time.asctime()))
        f.write('# {}\n# {}\n#\n'.format(OGcoef_file, velcoef_file))
        f.write('#{:3}{:8}\n\n'.format('Ap','V [km/s]'))
        for i, v in enumerate(outvel):
            f.write('{:3n}{:8.3f}\n'.format(i+1,v))

    return

def main():

    baseOG = 'NGC_891_P{}_bin30_allz2.coef.fits'
    basevel = 'NGC_891_P{}_bin30_allz2.coef.vel.fits'
    baseout = 'NGC_891_P{}_bin30_velocities.dat'

    for i in range(6):
        print baseOG.format(i+1)
        print basevel.format(i+1)
        print baseout.format(i+1)

        consolidate(baseOG.format(i+1),
                    basevel.format(i+1),
                    baseout.format(i+1))

    return

if __name__ == '__main__':
    main()
