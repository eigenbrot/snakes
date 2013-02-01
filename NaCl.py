import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import os
import time
from matplotlib.backends.backend_pdf import PdfPages as PDF

def make_galaxy(disk_frac,rot_curve=()):
    '''Makes a two-component galaxy where the stellar disk contributes disk_frac
    of the total rotation curve at 2.2h_r'''

    '''first, find the disk that matches what we want'''
    print "Finding disk..."
    t1 = time.time()
    diskpar = find_disk(disk_frac=disk_frac,rot_curve=rot_curve)[0]
    print "Finding halo..."
    halopars = find_halo(diskpar,rot_curve=rot_curve)[0]
    t2 = time.time()
    print "total time was {} seconds".format(t2-t1)
    return np.array([diskpar,halopars[0],halopars[1]])
    

def find_disk(disk_frac=0.6,rot_curve=()):
    
    h_r = 5.05 #kpc

    if len(rot_curve) > 0:
        r = rot_curve[0]
        v = rot_curve[1]
    else:
        V_r = 256.29 #km/s
        h_rot = 5.45 #kpc

        r = np.linspace(0,7*h_r,1000)
        v = V_r*np.tanh(r/h_rot)

    v22 = np.interp(11.11,r,v)
    
    os.system('rm tempd.in')
    gloop = open('tempd.in','w')
    gloop.write('tempd.gp\n0.0001 40.0\n500\nn\n0.00001 5\n500\nn\n')

    gloop.close()

    x0 = np.array([5.705e8])
    
    xf = spo.leastsq(disk_fit,x0,args=(v22,disk_frac),epsfcn=0.0001)

    m_data = np.transpose(np.loadtxt('tempd.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_r = m_data[0][m_idx]
    m_g = m_data[3][m_idx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(m_r/h_r,m_v,'k:')
    ax.plot(r/h_r,v)
    ax.set_xlabel('$r/h_r$')
    ax.set_ylabel('$V(r)$')
    fig.show()

    return xf[0]

def disk_fit(x,v22,disk_frac):

    os.system('rm tempd.gp tempd.dat')

    gp = open('tempd.gp','w')
    gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n0\n'.format(x))
    gp.close()

    os.system('/usr/local/Gravity/galloop < tempd.in > tempd.dat')
    m_data = np.transpose(np.loadtxt('tempd.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_r = m_data[0][m_idx]
    m_g = m_data[3][m_idx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5

    disk22 = np.interp(11.11,m_r,m_v)
    
    print disk22 - disk_frac*v22
    return np.array([np.abs(disk22 - disk_frac*v22)])
    
def find_halo(diskpar=[7.73e+08],rot_curve=()):
    
    h_r = 5.05 #kpc
    if len(rot_curve) > 0:
        r = rot_curve[0]
        v = rot_curve[1]
    else:
        V_r = 256.29 #km/s
        h_rot = 5.45 #kpc
        r = np.linspace(0,7*h_r,1000)
        v = V_r*np.tanh(r/h_rot)
    
    os.system('rm temph.in')
    gloop = open('temph.in','w')
    gloop.write('temph.gp\n0.0001 40.0\n500\nn\n0.00001 5\n500\nn\n')

    gloop.close()

    x0 = [7.17e+06,22.0]
    fig0 = plt.figure()

    t1 = time.time()
    xf = spo.fmin(halo_fit,x0,args=(r,v,fig0,diskpar))#,epsfcn=0.001)
    t2 = time.time()
    print 'halo time was {} seconds'.format(t2-t1)

    m_data = np.transpose(np.loadtxt('temph.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_r = m_data[0][m_idx]
    m_g = m_data[3][m_idx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(m_r/h_r,m_v,'k:')
    ax.plot(r/h_r,v)
    ax.set_xlabel('$r/h_r$')
    ax.set_ylabel('$V(r)$')
    fig.show()

    return (xf,t2-t1)

def halo_fit(x,r,v,fig,diskpar):

    print "\t{}".format(x)
    os.system('rm temph.gp temph.dat')

    gp = open('temph.gp','w')
    gp.write('1\n{0:8.3E} 5.05 0.43 0 0\n1\n{1:8.3E} 1.0 1.0 3.0 {2:4.3f} 10000\n'.format(diskpar,x[0],x[1]))
    gp.close()

    os.system('/usr/local/Gravity/galloop < temph.in > temph.dat')
    m_data = np.transpose(np.loadtxt('temph.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_zidx = np.where(m_data[0][m_idx] > 2.2*5.05)
    m_r = m_data[0][m_idx][m_zidx]
    m_g = m_data[3][m_idx][m_zidx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5

    i_v = np.interp(m_r,r,v)

    chisq = np.sum((i_v - m_v)**2)/(i_v.size - x.size - 1)

    # fig.clf()
    # ax = fig.add_subplot(111)
    # ax.plot(r, rot_curve)
    # ax.plot(m_r,m_v,':')
    # ax.text(25,100,'{:8.3E}\n{:4.3f}'.format(x[0],x[1]))
    # fig.show()

    print np.sum(((i_v - m_v)/(i_v.size - x.size - 1))**2,axis=0)
#    return (i_v - m_v)/(i_v.size - x.size - 1)
    return chisq

def find_hrot(gpars):
    
    hrot0 = [5.45]
    h_r = 5.05 #kpc
    V_r = 254.84 #km/s
    
    os.system('rm temphr.in')
    gloop = open('temphr.in','w')
    gloop.write('temphr.gp\n0.0001 {0:4.2f}\n500\nn\n0.00001 {1:4.2f}\n500\nn\n'.format(7.*5.05,5*0.43))
    
    gloop.close()
    
    gp = open('temphr.gp','w')
    gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n1\n{0[1]:8.3E} 1.0 1.0 3.0 {0[2]:8.3E} 10000\n'.format(gpars))
    gp.close()
    os.system('/usr/local/Gravity/galloop < temphr.in > temphr.dat')
    m_data = np.transpose(np.loadtxt('temphr.dat',skiprows=18))
    
    outd = {}
    pp = PDF('hrots.pdf')

    for z in np.unique(m_data[1])[::50]:
        zidx = np.where(m_data[1] == z)
        m_r = m_data[0][zidx]
        m_g = m_data[3][zidx]
        m_v = (m_g*m_r*3.08*10**11)**0.5
        
        hrotf = spo.leastsq(hrot_fit,hrot0,args=(m_r,m_v,V_r),epsfcn=0.0001)[0]

        fv = V_r*np.tanh(m_r/hrotf[0])

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(m_r,m_v,label='Galpot model')
        ax.plot(m_r,fv,':',label='h_rot fit')
        ax.set_xlabel('r [kpc]')
        ax.set_ylabel('$v_{\theta}$')
        ax.text(20,150,'$h_{{rot}} = {:4.3f}\,\mathrm{{kpc}}$\n$z/h_z={:4.3f}$'.format(hrotf[0],z/0.43),fontsize=12)
        ax.legend(loc=0)
        pp.savefig(fig)

        outd[z] = hrotf

    pp.close()
    return outd

def hrot_fit(h_rot,m_r,m_v,V_r):

    rot_curve = V_r*np.tanh(m_r/h_rot)

    print np.sum((rot_curve - m_v)**2,axis=0)
    return (rot_curve - m_v)


def make_V(gpars=np.array([4.35e+08,9.87e+06,2.19e+01]),title='',labels=''):

    realz = np.array([0,0.965*0.43,1.93*0.43,3.86*0.43])
    realv = np.array([255.377,244.504,240.967,237.265])
    realerr = np.array([3.49,5.57,20.20,19.51])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(realz,realv,yerr=realerr,marker='s',linestyle='',label='data')
    ax.set_xlabel('z [kpc]')
    ax.set_ylabel('V(r) [km/s]')
    ax.set_title(title)
    
    if len(gpars.shape) > 1:
        for pars,lab in zip(gpars,labels):
            os.system('rm tempv.in')
            gloop = open('tempv.in','w')
            gloop.write('tempv.gp\n0.0001 {0:4.2f}\n500\nn\n0.00001 5\n500\nn\n'.format(7.*5.05))
            
            gloop.close()

            gp = open('tempv.gp','w')
            if pars.size > 1:
                gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n1\n{0[1]:8.3E} 1.0 1.0 3.0 {0[2]:8.3E} 10000\n'.format(pars))
            else:
                gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n0\n'.format(pars))
            gp.close()
    
            os.system('/usr/local/Gravity/galloop < tempv.in > tempv.dat')
            m_data = np.transpose(np.loadtxt('tempv.dat',skiprows=18))
    
            rmax = m_data[0].max()
            maxidx = np.where(m_data[0] == rmax)
            
            z = m_data[1][maxidx]
            v = (m_data[3][maxidx]*m_data[0][maxidx]*3.08*10**11)**0.5
            
            ax.plot(z,v,label=lab)

    else:
        os.system('rm tempv.in')
        gloop = open('tempv.in','w')
        gloop.write('tempv.gp\n0.0001 {0:4.2f}\n500\nn\n0.00001 5\n500\nn\n'.format(7.*5.05))
        
        gloop.close()
        
        gp = open('tempv.gp','w')
        if gpars.size > 1:
            gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n1\n{0[1]:8.3E} 1.0 1.0 3.0 {0[2]:8.3E} 10000\n'.format(gpars))
        else:
            gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n0\n'.format(gpars))
        gp.close()
        
        print 'running Galloop...',
        os.system('/usr/local/Gravity/galloop < tempv.in > tempv.dat')
        print 'all done'
        m_data = np.transpose(np.loadtxt('tempv.dat',skiprows=18))
        
        # z = m_data[1]
        # v = np.array([])
        # uz = np.unique(z)
        # for zz in uz:
        #     zidx = np.where(z == zz)
        #     mr = m_data[0][zidx]
        #     mg = m_data[3][zidx]
        #     mv = (mg*mr*3.08e11)**0.5
        #     v = np.append(v,np.interp(11.11,mr,mv))

        rmax = m_data[0].max()
        maxidx = np.where(m_data[0] == rmax)
        
        z = m_data[1][maxidx]
        v = (m_data[3][maxidx]*m_data[0][maxidx]*3.08*10**11)**0.5
        
        ax.plot(z,v,label=labels)
        
    ax.legend(title='Disk maximality',loc=0)
    fig.show()

    return z,v

def find_shape():

    h_r = 5.05 #kpc
    V_r = 254.84 #km/s
    h_rot = 5.45 #kpc
    
    realz = np.array([0,0.965*0.43,1.93*0.43,3.86*0.43])
    realv = np.array([254.389,244.635,241.189,237.407])
    
    r = np.linspace(0,7*h_r,1000)
    rot_curve = V_r*np.tanh(r/h_rot)
    
    # os.system('rm tempsh.in')
    # gloop = open('tempsh.in','w')
    # gloop.write('tempsh.gp\n0.1 {0:4.2f}\n260\nn\n0.0001 3\n260\nn\n'.format(7.*h_r))

    # gloop.close()

    x0 = np.array([4.3e08])

    xf = spo.fmin(shape_fit,x0,args=(r,rot_curve),xtol=0.1,ftol=0.1)

    m_data = np.transpose(np.loadtxt('tempsh.dat',skiprows=18))
    
    rmax = m_data[0].max()
    maxidx = np.where(m_data[0] == rmax)
    
    z = m_data[1][maxidx]
    v = (m_data[3][maxidx]*m_data[0][maxidx]*3.08*10**11)**0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(z,v)
    ax.plot(realz,realv,marker='s',linestyle='')
    ax.set_xlabel('z [kpc]')
    ax.set_ylabel('V(r) [km/s]')
    fig.show()

    return


x0s = np.array([6.0e+06,30.0])

def shape_fit(x,r,rot_curve):
    '''too slow'''
    global x0s

    print "x0s: {}".format(x0s)
    realz = np.array([0,0.965*0.43,1.93*0.43,3.86*0.43])
    realv = np.array([254.389,244.635,241.189,237.407])

    xf = spo.fmin(shape_halo_fit,x0s,args=(x,r,rot_curve),xtol=0.1,ftol=0.1)
    print xf
    m_data = np.transpose(np.loadtxt('tempsh.dat',skiprows=18))
    
    rmax = m_data[0].max()
    maxidx = np.where(m_data[0] == rmax)
    
    z = m_data[1][maxidx]
    v = (m_data[3][maxidx]*m_data[0][maxidx]*3.08*10**11)**0.5

    i_v = np.interp(realz,z,v)

    chisq = np.sum((i_v - realv)**2)/(realv.size - 3. - 1.)

    x0s = xf

    print chisq
    return chisq

def shape_halo_fit(x,x1,r,rot_curve):
    '''too slow'''
    os.system('rm tempsh.gp tempsh.dat')

    gp = open('tempsh.gp','w')
    gp.write('1\n{0[0]:8.3E} 5.05 0.43 0 0\n1\n{1[0]:8.3E} 1.0 1.0 3.0 {1[1]:4.3f} 10000\n'.format(x1,x))
    gp.close()

    rrange = np.linspace(0.0001,35,10)
    zrange = np.linspace(0.0001,5,50)

#    f = open('tempsh.dat','w')
    temp_dat = np.empty((1,6))
    for rad in rrange:
        print rad
        for z in zrange:
            to = os.popen('echo "tempsh.gp\n{} {}\n" | /usr/local/Gravity/galpot'.format(rad,z))
#            f.write(to.readlines()[-1])
            temp_dat = np.vstack((temp_dat, np.array(to.readlines()[-1].split(),dtype=np.float32)))

#    f.close()
    m_data = temp_dat[1:].T
#    m_data = np.transpose(np.loadtxt('tempsh.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_zidx = np.where(m_data[0][m_idx] > 2.2*5.05)
    m_r = m_data[0][m_idx][m_zidx]
    m_g = m_data[3][m_idx][m_zidx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5

    i_v = np.interp(m_r,r,rot_curve)
    
    chisq = np.sum((i_v - m_v)**2)/(i_v.size - x.size - 1)
    
    print "\t{}".format(chisq)
    return chisq
