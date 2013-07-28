import numpy as np
import Salty2 as salty

'''This file creates a bunch of simulation files with varying degrees of
spirality
'''

def do():

    h_r = 8.3
    vc_list = [216.78,203.7565,200.516,242.0]
    zlist = [0,0.43,0.86,1.72]

    w = 0.8
    N = 6
    view_ang = 0.    
    pitch_list = np.arange(0.1,np.pi,np.pi/3.)

    for z, vc in zip(zlist,vc_list):
        for pitch in pitch_list:
            
            name = 'sim_z{:n}_spiral_p{:04.0f}.fits'.format(z/0.43,pitch*1000)
            print name
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output=name,
                           spiralpars=dict(w=w,
                                           N=N,
                                           view_ang=view_ang,
                                           pitch=pitch))
        salty.plt.close('all')
        salty.gc.collect()

    return
    
