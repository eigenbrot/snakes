import numpy as np
import scipy.special as sps
import matplotlib.pyplot as plt

def kappa(r, z, tau=0.7, zd=0.29, hd=8.1):

    kappa0 = tau/(2 * zd)

    return kappa0 * np.exp(-1 * (r/hd + z/zd))

def A(r, z, tau0=0.7, zd=0.29, hd=8.1):

        kappa0 = tau0/(2 * zd)
        
        return 1.086 * 2 * kappa0 * r * sps.kn(1,r/hd) * np.exp(-1 * z/zd)

def A_vec(r, Z, tau0=0.7, zd=0.29, hd=8.1):
    
    output = np.zeros(Z.size)

    for i, z in enumerate(Z):
        output[i] = A(r,z,tau0=tau0, zd=zd, hd=hd)

    return output

def plot_A(r, tau0=0.7, zd=0.29, hd=8.1):
    
    zlist = np.linspace(0,2.5,20)
    #Alist = np.zeros(20)

    # for i, z in enumerate(zlist):
    #     Alist[i] = A(r,z,tau0=tau0, zd=zd, hd=hd)

    Alist = A_vec(r,zlist,tau0=tau0, zd=zd, hd=hd)

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('z [kpc]')
    ax.set_ylabel('A$_V$')
    ax.plot(zlist,Alist)

    ax.figure.show()
                
