import moment_tests as mt
import Salty2 as salty
import numpy as np
import gc

'''This file is really just a script that allows me to make a bunch of flared
and non-flared models and compare them to data. It is designed to be run in
/usr/users/eigenbrot/research/edgeOn/flares/sims
'''

def do():

    h_r = 8.3
    hzrlist = [12.45]#np.arange(0.25,1.75,0.25) * h_r
    filelist = [
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z0_MgI_binz2.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z1_MgI_binz2.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z2_MgI_bin60.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z4_MgI_binz2.slay.fits']
    zlist = [0,0.43,0.86,1.72]
    fliplist = [True,False,False,False]

    print "hzrlist:"
    print hzrlist

    vclist = [216]*4#[salty.find_Vc(*salty.openslay(i),back=back) 
              #for (i,back) in zip(filelist,fliplist)]
    salty.plt.close('all')
    
    '''generate some sim files'''
    for filename, z, flip, vc in zip(filelist,zlist,fliplist,vclist):
        zx = filename.split('_')[2]
        noname = 'sim_{}_noflare.fits'.format(zx)
        salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output = noname)
        noprefix = 'noflare_comp_{}'.format(zx)
#        mt.line_comp(noname,filename,noprefix,flip=flip)

        for hzr in hzrlist:
            flarename = 'sim_{:}_nflarea_{:3.0f}.fits'.format(zx,hzr*100.)
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,
                           output=flarename,flarepars=dict(h_zR = hzr,
                                                           ftype='exp'))

            prefix = 'nflare_comp_{:}_f{:04.0f}'.format(zx,hzr*100.)
#            mt.line_comp(flarename,filename,prefix,flip=flip)
            salty.plt.close('all')
            mt.plt.close('all')
            salty.gc.collect()
            mt.gc.collect()
            gc.collect()

    return

def do_linear():
    '''run in /d/monk/eigenbrot/EdgeOn/flare_test/linear
    '''

    h_r = 8.3
    hzrlist = np.arange(0.1,1.0,0.25) * h_r
    filelist = [
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z0_MgI_binz2.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z1_MgI_binz2.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z2_MgI_bin60.slay.fits',
        '/d/monk/eigenbrot/EdgeOn/SALT_data/tiESO_z4_MgI_binz2.slay.fits']
    zlist = [0,0.43,0.86,1.72]
    fliplist = [True,False,False,False]

    print "hzrlist:"
    print hzrlist

    vclist = [salty.find_Vc(*salty.openslay(i),back=back) 
              for (i,back) in zip(filelist,fliplist)]
    print "vclist:"
    print vclist
    salty.plt.close('all')
    
    '''generate some sim files'''
    for filename, z, flip, vc in zip(filelist,zlist,fliplist,vclist):
        print filename+':'
        zx = filename.split('_')[2]
        noname = 'sim_{}_noflare.fits'.format(zx)
        salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output = noname)
        noprefix = 'noflare_comp_{}'.format(zx)
#        mt.line_comp(noname,filename,noprefix,flip=flip)

        for hzr in hzrlist:
            print filename+'\t'+str(hzr)
            flarename = 'sim_{:}_lflare_{:03.0f}.fits'.format(zx,hzr*100.)
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,
                           output=flarename,
                           flarepars=dict(ftype='linear',h_zR = hzr))

            prefix = 'lflare_comp_{:}_f{:04.0f}'.format(zx,hzr*100.)
#            mt.line_comp(flarename,filename,prefix,flip=flip)
            salty.plt.close('all')
            mt.plt.close('all')
            salty.gc.collect()
            mt.gc.collect()
            gc.collect()

    return    
