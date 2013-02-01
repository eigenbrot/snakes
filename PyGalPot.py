import numpy as np
import matplotlib.pyplot as plt
import os

class PyGalPot:

    def __init__(self,galpars,h_r=5.05,h_z=0.43):
        self.pars = galpars
        self.h_r = h_r
        self.h_z = h_z

        os.system('rm tempPGP.in')
        gloop = open('tempPGP.in','w')
        gloop.write('tempPGP.gp\n0.0001 50.0\n500\nn\n0.00001 5\n500\nn\n')
        gloop.close()

        os.system('rm tempPGP.gp')
        gp = open('tempPGP.gp','w')
        gp.write('1\n{:8.3E} {} {} 0 0\n1\n{:8.3E} 1.0 1.0 3.0 {:4.3f} 10000\n'.format(
                galpars[0],h_r,h_z,galpars[1],galpars[2]))
        gp.close()

        os.system('/usr/local/Gravity/galloop < tempPGP.in > tempPGP.dat')
        
        self.R, self.Z, self.phi, self.gR, self.gZ =\
            np.loadtxt('tempPGP.dat',unpack=True,usecols=(0,1,2,3,4),skiprows=18)

    def get_TVC(self,height):
        sortidx = np.argsort(np.abs(self.Z - height))
        z = self.Z[sortidx][0]

        print "extracting TVC from galpot model at z = {} kpc".format(z)

        zidx = np.where(self.Z == z)
        r = self.R[zidx]
        gr = self.gR[zidx]
        v = (gr*r*3.08e11)**0.5

        return r, v
        
    def plot_TVC(self,height):
        
        r, v = self.get_TVC(height)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Radius [kpc]')
        ax.set_ylabel('Velocity [km/s]')
        ax.plot(r,v)
        fig.show()
        return
