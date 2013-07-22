import numpy as np
import Salty2 as salty
import pyfits

def make_image(output,h_zR):

    v_r = 238.3067
    hrot = 5.45
    h_z = 0.43
    size = 501
    scale = 100.0/size
    flarepars = dict(h_zR = h_zR)

    maxz = 2.0
    numz = 250

    heights = np.linspace(0,maxz,numz)
    
    SBlist = []

    for Z in heights:

        name = 'eonview_f{:04.0f}_{:04.0f}.fits'.format(h_zR*100,Z*1000)
        print Z, name
        salty.simcurve(size,Z,v_r,hrot,scale=scale,flarepars=flarepars,
                       output=name)
        SB = pyfits.open(name)['SB'].data[size/2,:]
        SBlist = [SB] + SBlist
        if Z != 0.0:
            SBlist = SBlist + [SB]

    galaxy = np.vstack(SBlist)
    HDU = pyfits.PrimaryHDU(galaxy)
    HDU.header.update('CRPIX1',size/2,comment='WCS: X reference pixel')
    HDU.header.update('CRPIX2',numz,comment='WCS: Y reference pixel')
    HDU.header.update('CRVAL1',0.0,'WCS: X reference coordinate value')
    HDU.header.update('CRVAL2',0.0,'WCS: Y reference coordinate value')    
    HDU.header.update('CDELT1',scale,comment='WCS: X pixel size')
    HDU.header.update('CDELT2',maxz/(numz-1),comment='WCS: Y pixel size')
    HDU.header.update('CTYPE1','LINEAR',comment='X type')
    HDU.header.update('CTYPE2','LINEAR',comment='Y type')
    HDU.writeto(output,clobber=True)

    return
