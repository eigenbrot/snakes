import numpy as np
import glob
import time
from yanny import yanny
glob = glob.glob

basepath = '/d/monk/eigenbrot/WIYN/14B-0456'

def create_yanny(pointing, output):
    
    numaps = get_numaps(pointing)
    
    data = np.zeros((numaps,),
                    dtype = [('ap','i'),
                             ('apsize','f4'),
                             ('vdisp','f4'),
                             ('r','f4'),
                             ('z','f4'),
                             ('SN','f4'),
                             ('gbu','i'),
                             ('ur_age','f4'),
                             ('ur_Av','f4'),
                             ('ur_Z','f4'),
                             ('ur_chisq','f4'),
                             ('ur_redchi','f4'),
                             ('ur_bluechi','f4'),
                             ('at_age','f4'),
                             ('at_Av','f4'),
                             ('at_Z','f4'),
                             ('at_chisq','f4'),
                             ('at_redchi','f4'),
                             ('at_bluechi','f4'),
                             ('az_age','f4'),
                             ('az_Av','f4'),
                             ('az_Z','f4'),
                             ('az_chisq','f4'),
                             ('az_redchi','f4'),
                             ('az_bluechi','f4')])

    data = get_basics(pointing, data)
    data = get_unreg(pointing, data)
    data = get_alltau(pointing, data)
    data = get_allz(pointing, data)

    par = yanny(filename=None, np=True, debug=False)
    struct = par.dtype_to_struct(data.dtype, 
                                 structname='APINFO',
                                 enums=dict())
    par['symbols']['struct'] = struct['struct']
    par['symbols']['enum'] = struct['enum']
    par['symbols']['APINFO'.upper()] = struct['APINFO'.upper()]
    par['APINFO'] = data

    par['pointing'] = pointing
    par['numaps'] = numaps
                             
    with open(output,'w') as f:
        f.write("""# NGC 891 aperture file for pointing {}
# File generated on {}
# 
# Uses SDSS-style Yanny .par formatting
#
# This file contains a collection of all of the important parameters derived 
# from full spectrum fitting of our NGC 891 data. The various fitting methods
# are described below:
#
# unreg: 
#    The OG. Each age is forced to have a single metallicity (found as a
#    separate step) and there is a single extinction value for all the SSPs 
#    in an aperture.
#
# all_tau:
#     Each age is forced to have a single metallicity (found as a
#     separate step), but each age is allowed to have a different A_V value
#     (limited by 0 < A_V < 20).
#
# all_Z2:
#    Each age can have an arbitrary number of metallicities with
#    arbitrary weights. There is a single extinction value for all the SSPs 
#    in an aperture.
#
# The data structure is as follows:
#
# Header keywords:
# pointing: Which pointing's data is contained in this file
# numaps: The number of binned apertures contained in this pointing
#
# Data table entries:
#  1) ap: The aperture number.
#  2) apsize: The size of the fibers combined into the aperture, in arcsec.
              Note that this is NOT the size of the final aperture.
#  3) vdisp: The resolution (in km/s) of the aperture. Different fiber sizes
             produce different resolutions.
#  4) r: Radius (in kpc) of aperture from the center of the galaxy
#  5) z: Height (in kpc) of aperture above the midplane of the galaxy
#  6) SN: Final signal to noise ratio of the binned aperture
#  7) gbu: Quality flag. 0=good, 1=bad fit, 2=ugly for some other reason
#  8) ur_age: MLWA from unreg fit, in Gyr
#  9) ur_Av: Av from unreg fit
# 10) ur_Z: Best metallicity (in chisq sense) of unreg fit
# 11) ur_chisq: Full, reduced chisq value from best unreg fit
# 12) ur_redchi: Reduced chisq for 5250 <= lambda <= 6800 from unreg
# 13) ur_bluechi: Reduced chisq for 3750 <= lambda < 5250 from unreg
# 14) at_age: MLWA from all_tau fit, in Gyr
# 15) at_Av: MLWAv all_tau from  fit
# 16) at_Z: Best metallicity (in chisq sense) of all_tau fit
# 17) at_chisq: Full, reduced chisq value from best all_tau fit
# 18) at_redchi: Reduced chisq for 5250 <= lambda <= 6800 from all_tau fit
# 19) at_bluechi: Reduced chisq for 3750 <= lambda < 5250 from all_tau fit
# 20) az_age: MLWA from all_Z2 fit, in Gyr
# 21) az_Av: Av from all_Z2 fit
# 22) az_Z: MLWZ from all_Z2 fit
# 23) az_chisq: Full, reduced chisq value from best all_Z2 fit
# 24) az_redchi: Reduced chisq for 5250 <= lambda <= 6800 from all_Z2 fit
# 25) az_bluechi: Reduced chisq for 3750 <= lambda < 5250 from all_Z2 fit
#
""".format(pointing, time.asctime()))
        
        par.write(f)
    
    return

def get_numaps(pointing):
    try:
        datfile = glob('{}/anal/unreg/good/multi_Z/NGC*P{}*_fit.dat'.\
                       format(basepath, pointing))[0]
    except IndexError:
        print 'Could not find anything in {}/anal/unreg/good/multi_Z/NGC*P{}_fit.dat'.format(basepath, pointing)
    print 'Getting numaps from {}'.format(datfile)
    d = np.loadtxt(datfile)
    
    return d.shape[0]

def get_basics(pointing, data):

    datfile = glob('{}/anal/unreg/good/multi_Z/NGC*P{}*_fit.dat'.\
                   format(basepath, pointing))[0]
    locfile = glob('{}/bin/NGC*P{}*_locations.dat'.\
                   format(basepath, pointing))[0]
    print 'Getting basic info from {}'.format(datfile)
    print '\tand {}'.format(locfile)
    
    did, SN = np.loadtxt(datfile, usecols=(0,14), unpack=True)
    lid, size, r, z = np.loadtxt(locfile, 
                                 usecols=(0,1,4,5), unpack=True)

    vdispd = {0.937: 493/2.355,
              1.406: 589/2.355,
              1.875: 691/2.355,
              2.344: 796/2.355,
              2.812: 966/2.355}

    for i in range(did.size):
        if did[i] != lid[i]+1:
            print 'WARNING: ap IDs {} and {} do not match!'.\
                format(did[i], lid[i])
        data['ap'][i] = did[i]
        data['apsize'][i] = size[i]
        data['r'][i] = r[i]
        data['z'][i] = z[i]
        data['SN'][i] = SN[i]
        data['vdisp'][i] = vdispd[size[i]]

    return data

def get_unreg(pointing, data):

    datfile = glob('{}/anal/unreg/good/multi_Z/NGC*P{}*_fit.dat'.\
                   format(basepath, pointing))[0]
    
    print 'Getting unreg info from {}'.format(datfile)

    ap, MLWA, TauV, chisq, redchi, bluechi, Z = np.loadtxt(datfile,
                                                           usecols=(0,12,
                                                                    13,15,
                                                                    16,17,19),
                                                           unpack=True)
    
    for i in range(ap.size):
        if ap[i] != data['ap'][i]:
            print 'WARNING: unreg ap {} does not match data ap {}'.\
                format(ap[i], data['ap'][i])
            raw_input('')

        data['ur_age'][i] = MLWA[i]
        data['ur_Av'][i] = 1.086*TauV[i]
        data['ur_Z'][i] = Z[i]
        data['ur_chisq'][i] = chisq[i]
        data['ur_redchi'][i] = redchi[i]
        data['ur_bluechi'][i] = bluechi[i]
    
    return data

def get_alltau(pointing, data):

    datfile = glob('{}/anal/all_tau/multi_Z/NGC*P{}*_fit.dat'.\
                   format(basepath, pointing))[0]
    
    print 'Getting all_tau info from {}'.format(datfile)

    ap, MLWA, MLWT, chisq, redchi, bluechi, Z = np.loadtxt(datfile,
                                                           usecols=(0,22,
                                                                    24,26,
                                                                    27,28,30),
                                                           unpack=True)
    
    for i in range(ap.size):
        if ap[i] != data['ap'][i]:
            print 'WARNING: all_tau ap {} does not match data ap {}'.\
                format(ap[i], data['ap'][i])
            raw_input('')

        data['at_age'][i] = MLWA[i]
        data['at_Av'][i] = 1.086*MLWT[i]
        data['at_Z'][i] = Z[i]
        data['at_chisq'][i] = chisq[i]
        data['at_redchi'][i] = redchi[i]
        data['at_bluechi'][i] = bluechi[i]
    
    return data

def get_allz(pointing, data):

    datfile = glob('{}/anal/all_Z2/NGC*P{}*_allz2.dat'.\
                   format(basepath, pointing))[0]
    
    print 'Getting all_Z info from {}'.format(datfile)

    ap, MLWA, MLWZ, TauV, chisq, redchi, bluechi = np.loadtxt(datfile,
                                                              usecols=(0,62,
                                                                       64,65,
                                                                       67,68,
                                                                       69),
                                                              unpack=True)
    
    for i in range(ap.size):
        if ap[i] != data['ap'][i]:
            print 'WARNING: all_Z ap {} does not match data ap {}'.\
                format(ap[i], data['ap'][i])
            raw_input('')

        data['az_age'][i] = MLWA[i]
        data['az_Av'][i] = 1.086*TauV[i]
        data['az_Z'][i] = MLWZ[i]
        data['az_chisq'][i] = chisq[i]
        data['az_redchi'][i] = redchi[i]
        data['az_bluechi'][i] = bluechi[i]
    
    return data
