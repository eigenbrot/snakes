import moment_tests as mt
import Salty2 as salty
import numpy as np


'''This file is for making simulations with varrying degrees of ring-itude
'''

def do():

    h_r = 8.3
    r_R = h_r*3.
    r_sig = 1.
    vc_list = [216.78,203.7565,200.516,242.0]
    zlist = [0,0.43,0.86,1.72]
    weight_list = np.arange(0.1,1.,0.25)

    for z, vc in zip(zlist,vc_list):        
        for r_w in weight_list:

            name = 'sim_z{:n}_ring_w{:03.0f}.fits'.format(z/0.43,r_w*100)
            print name
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output=name,
                           ringpars=dict(r_R=r_R,r_sig=r_sig,r_w=r_w))
            
        salty.plt.close('all')
        salty.gc.collect()

    return
    
def do_radius():

    h_r = 8.3
    r_w = 0.7
    r_sig = 1.
    vc_list = [216.78,203.7565,200.516,242.0]
    zlist = [0,0.43,0.86,1.72]
    radius_list = np.arange(1,30,10.)

    for z, vc in zip(zlist,vc_list):        
        for r_R in radius_list:

            name = 'sim_z{:n}_ring_R{:04.0f}.fits'.format(z/0.43,r_R*100)
            print name
            salty.simcurve(1001,z,vc,5.45,scale=100/1001.,output=name,
                           ringpars=dict(r_R=r_R,r_sig=r_sig,r_w=r_w))
            
        salty.plt.close('all')
        salty.gc.collect()

    return
