import moment_tests as mt
import Salty2 as salty
import numpy as np

'''This file is really just a script that allows me to make a bunch of flared
and non-flared models and compare them to data. It is designed to be run in
/usr/users/eigenbrot/research/edgeOn/flares/sims
'''

def do():

    h_r = 8.3
    hzrlist = np.arange(1,2)
    filelist = ['tiESO_z0_MgI.slay.fits',
                'tiESO_z05_MgI.slay.fits',
                'tiESO_z1_MgI.slay.fits']
    zlist = [0,0.43,0.86]
    fliplist = [False,True,False]

    vclist = [salty.find_Vc(*salty.openslay(i)) for i in filelist]
    
    '''generate some sim files'''
    for filename, z, flip, vc in zip(filelist,zlist,fliplist,vclist):
        print filename+':'
        zx = filename.split('_')[1]
        noname = 'sim_{}_noflare.fits'.format(zx)
        salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output = noname)
        noprefix = 'noflare_comp_{}'.format(zx)
        mt.line_comp(noname,filename,noprefix,flip=flip)

        for hzr in hzrlist:
            print '\t'+str(hzr)
            flarename = 'sim_{:}_flare_{:3.0f}.fits'.format(zx,hzr*100.)
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,
                           output=flarename,flarepars=dict(h_zR = hzr))

            prefix = 'flare_comp_{:}_f{:3.0f}'.format(zx,hzr*100.)
            mt.line_comp(flarename,filename,prefix,flip=flip)

    return
