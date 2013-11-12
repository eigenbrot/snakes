#! /usr/bin/env python

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import sys
import scipy.optimize as spo
import os

def rotpot(infile,scale_length,z,batch,components):
    
    fig = plt.figure(0)
    fig.clf()
    ax = fig.add_subplot(111)

    sanc = np.array([70,90,120,150,150,180,190,210,230\
                         ,230,230,225,230,230,230,230,225\
                         ,225,225,225,240,225,225,210,200,230,230])
    sancr = (1.0178)*np.array([0.25,0.25,0.5,0.5,0.75,0.75,1,1.25,\
                          1.5,1.75,2,2.25,2.5,2.75,3,3.25,\
                          3.5,3.75,4,4.25,4.25,4.5,4.75,5,5.25,5.5,5.75])

    ax.plot(sancr,sanc,'g.')
    ax.set_xlabel("$r/h_r$")
    ax.set_ylabel("$v_{rot}$ [km/s]")

    
    clist = range(len(infile))

    for i in range(len(infile)):

        data = np.transpose(np.loadtxt(infile[i],skiprows=18))
        zidx = np.where(data[1] >= z)[0][0]

        print "First occurance of z="+str(data[1][zidx])+" at index "+str(zidx)

        'We only want data from a certain z'
        idx = np.where(data[1] == data[1][zidx])
        r = data[0][idx]
        g_r = data[3][idx]
    
        v_r = (g_r*r*3.08*10**11)**0.5


        clist[i], = ax.plot(r/scale_length,v_r, label='Total')
    
#    ax.set_title(ifile+', $h_r$ = '+str(scale_length)+' kpc, '\
#                     +'$h_z$ = 0.32 kpc\nz = '\
#        ax.set_title('z = '+str(round(data[1][zidx],1))+' kpc = '\
#                         +str(round(data[1][zidx]/0.32,1))+' $h_z$')

#    ax.text(5,120,'$\Sigma_d = 6.5e8 M_{sol}/kpc^2$\n'\
#                +'$ \rho_b = 1.3e9 M_{sol}/kpc^3$\n'\
#                +'$ \rho_h = 3.4e6 M_{sol}/kpc^3$\n'\
#                +'$r_e = 36 kpc$\n')
                

    if components:
        disk_g = np.transpose(np.loadtxt(ifile+'.disk',skiprows=18))[3][idx]
        bulge_g = np.transpose(np.loadtxt(ifile+'.bulge',skiprows=18))[3][idx]
        halo_g = np.transpose(np.loadtxt(ifile+'.halo',skiprows=18))[3][idx]

        disk_v = (disk_g*r*3.08*10**11)**0.5
        bulge_v = (bulge_g*r*3.08*10**11)**0.5
        halo_v = (halo_g*r*3.08*10**11)**0.5

        ax.plot(r/scale_length,disk_v, 'k--', label='Disk')
        ax.plot(r/scale_length,bulge_v, 'k:', label='Bulge')
        ax.plot(r/scale_length,halo_v, 'k-.', label='Halo')
        
        ax.legend(loc=5)
    
    if not batch:
        ax.legend(clist,["Half-max model","$\chi^2$ model","By-hand model"]\
                      ,loc=5)
        fig.show()
        scratch = raw_input("Enter to quit")
        
    return

def batch(data_file,z_file,out_file,components):
    
    pp = PdfPages(str(out_file))
    
    z_list = np.transpose(np.loadtxt(z_file))

    for i in range(z_list.size):
        rotpot(data_file,3.93,z_list[i],True,components)
        pp.savefig()

    pp.close()
    return

def amoeba(compare, rh0, r_e):
    
    c_data = np.transpose(np.loadtxt(compare,skiprows=18))
    c_zidx = np.where(c_data[1] >= 0)[0][0]

    c_idx = np.where(c_data[1] == c_data[1][c_zidx])
    c_r = c_data[0][c_idx]
    c_g = c_data[3][c_idx]
    
    c_v = (c_g*c_r*3.08*10**11)**0.5


    os.system('rm temp.in')
    gloop = open('temp.in','w')
    gloop.write('temp.gp\n0.0001 24\n500\nn\n0.000001 5\n500\nn\n')
    
    gloop.close()

    x0 = np.array([rh0,r_e])

    xopt = spo.fmin(curve_compare, x0, args=((c_v,)))
    
    print xopt

    return

def curve_compare(p, c_v):

    m_v = gen_curve(p)
    chisq = np.sum((m_v - c_v)**2)

    print chisq

    return chisq

def gen_curve(p):

    os.system('rm temp.gp temp.dat')

    gp = open('temp.gp','w')
    gp.write('1\n'\
                 +'9.705e+08 3.93 0.32 0 0\n'\
                 +'2\n'\
                 +'1.905e+09 0.75 0.211 3.35 0.24 10000\n'\
                 +'{0[0]:8.2E} 1.0 1.0 3.0 {0[1]:2.1f} 10000'.format(p)\
                 +'\n')
    gp.close()

    os.system('/usr/local/Gravity/galloop < temp.in > temp.dat')

    m_data = np.transpose(np.loadtxt('temp.dat',skiprows=18))
    m_zidx = np.where(m_data[1] >= 0)[0][0]

    m_idx = np.where(m_data[1] == m_data[1][m_zidx])
    m_r = m_data[0][m_idx]
    m_g = m_data[3][m_idx]
    
    m_v = (m_g*m_r*3.08*10**11)**0.5
    
    return m_v


def vdiff(file1, file2, zlist):
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)

    data1 = np.transpose(np.loadtxt(file1,skiprows=18))
    zidx1 = np.where(data1[1] >= 0)[0][0]

    data2 = np.transpose(np.loadtxt(file2,skiprows=18))
    zidx2 = np.where(data2[1] >= 0)[0][0]
    
    'We only want data from a certain z'
    idx1 = np.where(data1[1] == data1[1][zidx1])
    r1 = data1[0][idx1]
    g_r1 = data1[3][idx1]

    idx2 = np.where(data2[1] == data2[1][zidx2])
    r2 = data2[0][idx2]
    g_r2 = data2[3][idx2]

    v_r1 = (g_r1*r1*3.08*10**11)**0.5
    v_r2 = (g_r2*r2*3.08*10**11)**0.5

    delta_v = v_r1 - v_r2
#    ax.plot(r1/3.93, delta_v,'k--', label='$\Delta v_{mod}$')

    ax.set_xlabel("$r/h_r$")
    ax.set_ylabel("v(r,z)/v(r,z=0) [km/s]")
    ax.set_xlim(-0.1,6)
    ax.set_ylim(0,1.1)
#    ax.set_title("$\Delta v_{mod} = $"+file1+" - "+file2+"\n"\
#                     +"$z$ compaison is max - 0.5max")

    clist = ['b','r','g','c','m']

    vav = (v_r1 + v_r2)/2

    for i in range(len(zlist)):
        zidx22 = np.where(data2[1] >= zlist[i])[0][0]
        idx22 = np.where(data2[1] == data2[1][zidx22])
        g_r22 = data2[3][idx22]
        v_r22 = (g_r22*r2*3.08*10**11)**0.5

        zidx11 = np.where(data1[1] >= zlist[i])[0][0]
        idx11 = np.where(data1[1] == data1[1][zidx11])
        g_r11 = data1[3][idx11]
        v_r11 = (g_r11*r2*3.08*10**11)**0.5
        
        ax.plot(r2/3.93, (v_r11/v_r1), clist[i]+'--')
        ax.plot(r2/3.93, (v_r22/v_r2), clist[i]+'-',label=str(zlist[i]/0.32))
    
    ax.legend(loc=5,title='$z/h_z$')

    fig.show()
    
    return

if __name__ == "__main__":
    
    if sys.argv[1].find('-') != -1:
        if sys.argv[1].find('z') != -1:
            if sys.argv[1].find('c') == -1:
                rotpot((sys.argv[3],),3.93,float(sys.argv[2]),False,False)
            
            elif sys.argv[1].find('c') != -1:
                rotpot((sys.argv[3],),3.93,float(sys.argv[2]),False,True)

        if sys.argv[1].find('b') != -1:
            if sys.argv[1].find('c') == -1:
                batch((sys.argv[3],),sys.argv[2],sys.argv[4],False)
            
            elif sys.argv[1].find('c') != -1:
                batch((sys.argv[3],),sys.argv[2],sys.argv[4],True)
    
    else:
        rotpot((sys.argv[1:]),3.93,0,False,False)
