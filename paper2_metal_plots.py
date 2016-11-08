import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib import rc, rcdefaults
#All the bad and ugly apertures from the paper
pap_ex = [[2, 5, 34], 
          [1, 2, 9, 35, 37], 
          [54,56,59], 
          [2, 8, 13, 55], 
          [1, 2, 3, 11, 13, 19, 27, 28, 29, 5], 
          [35, 36, 38]]

def tau_Z(output,color=False, vmax=20, exclude=pap_ex, 
          mass=False, mass_age=False, pivot_age=6, show_fit=True):

    basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'
    #basedir = '/Users/Arthur/Documents/School/final_results'

    tlist = []
    telist = []
    Zlist = []
    Zelist = []
    rlist = []

    for p in range(6):
        coef = '{}/NGC_891_P{}_bin30_allz2.fiterr.fits'.format(basedir,p+1)
        cc = pyfits.open(coef)[1].data
        exarr = np.array(exclude[p]) - 1
        if mass_age:
            mass_A_file = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
            dd = pyfits.open(mass_A_file)[1].data
            tlist.append(np.delete(dd['MMWA'],exarr))
            Ae = cc['dMLWA']/cc['MLWA']*dd['MMWA']
            telist.append(np.delete(Ae,exarr))
        else:
            tlist.append(np.delete(cc['MLWA'],exarr))
            telist.append(np.delete(cc['dMLWA'],exarr))
        if mass:
            mass_Z_file = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
            dd = pyfits.open(mass_Z_file)[1].data
            Zlist.append(np.delete(dd['MMWZ'],exarr))
            Ze = cc['dMLWZ']/cc['MLWZ']*dd['MMWZ']
            Zelist.append(np.delete(Ze,exarr))
        else:
            Zlist.append(np.delete(cc['MLWZ'],exarr))
            Zelist.append(np.delete(cc['dMLWZ'],exarr))

        if color:
            rc('figure.subplot', right=0.99)
            rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
            rr = np.loadtxt(rphi, usecols=(1,), unpack=True)
            rlist.append(np.delete(rr,exarr))
        
    Z = np.hstack(Zlist)
    t = np.hstack(tlist)
    Ze = np.hstack(Zelist)
    te = np.hstack(telist)

    idx = np.where(t <= pivot_age)[0]
    idx2 = np.where(t > pivot_age)[0]
    p1 = np.polyfit(t[idx],Z[idx],1)
    p2 = np.polyfit(t[idx2],Z[idx2],1)

    print p1
    print p2
    
    lfit = np.poly1d(p1)
    hfit = np.poly1d(p2)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if mass_age:        
        ax.set_xlabel(r'$\tau_M\mathrm{\ [Gyr]}$')
    else:
        ax.set_xlabel(r'$\tau_L\mathrm{\ [Gyr]}$')
    if mass:
        ax.set_ylabel(r'$Z_M\ [Z_\odot]$')
    else:
        ax.set_ylabel(r'$Z_L\ [Z_\odot]$')

    ax.errorbar(t,Z,yerr=Ze,xerr=te,fmt='none', capsize=0, 
                elinewidth=1.8, ecolor='grey',zorder=0,alpha=0.6)
    
    if color:
        r = np.hstack(rlist)
        scat = ax.scatter(t,Z,c=r,s=45,linewidth=0, cmap=plt.cm.gnuplot2, vmin=0, vmax=vmax)
        cb = fig.colorbar(scat)
        cb.set_label(r'$r$')
    else:
        ax.scatter(t,Z,c='k',s=45,linewidth=0)

    if show_fit:
        ax.axvline(pivot_age,ls=':',color='k',alpha=0.8)
        ax.plot(t[idx],lfit(t[idx]),color='grey',zorder=0,lw=3)
        ax.plot(t[idx2],hfit(t[idx2]),color='grey',zorder=0,lw=3)

        ax.text(1,2.3,r'$dZ/dt \approx {:4.2f}\ Z_\odot/\mathrm{{Gyr}}$'.format(p1[0]*-1))
        ax.text(7.5,2.3,r'$dZ/dt \approx {:4.2f}\ Z_\odot/\mathrm{{Gyr}}$'.format(p2[0]*-1))
    
    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def mass_v_light(output, quant='A', exclude=pap_ex, vmax=20):

    basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'
    #basedir = '/Users/Arthur/Documents/School/final_results'

    mass_list = []
    light_list = []
    rlist = []

    for p in range(6):
        coef = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
        cc = pyfits.open(coef)[1].data
        exarr = np.array(exclude[p]) - 1
        mass_list.append(np.delete(cc['MMW{}'.format(quant)],exarr))
        light_list.append(np.delete(cc['MLW{}'.format(quant)],exarr))

        rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
        rr = np.loadtxt(rphi, usecols=(1,), unpack=True)
        rlist.append(np.delete(rr,exarr))
        
    mass = np.hstack(mass_list)
    light = np.hstack(light_list)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('MMW{}'.format(quant))
    ax.set_ylabel('MLW{}'.format(quant))    
    r = np.hstack(rlist)
    scat = ax.scatter(mass,light,c=r,s=45,linewidth=0, cmap=plt.cm.gnuplot2, vmin=0, vmax=vmax)
    cb = fig.colorbar(scat)
    cb.set_label(r'$r$')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def agemass_diff(output, exclude=pap_ex, vmax=20):

    basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'
    #basedir = '/Users/Arthur/Documents/School/final_results'

    age_list = []
    metal_list = []
    rlist = []

    for p in range(6):
        coef = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
        cc = pyfits.open(coef)[1].data
        exarr = np.array(exclude[p]) - 1
        mdiff = cc['MLWZ'] - cc['MMWZ']
        adiff = cc['MLWA'] - cc['MMWA'] 
        age_list.append(np.delete(adiff,exarr))
        metal_list.append(np.delete(mdiff,exarr))

        rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
        rr = np.loadtxt(rphi, usecols=(1,), unpack=True)
        rlist.append(np.delete(rr,exarr))
        
    age = np.hstack(age_list)
    metal = np.hstack(metal_list)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$Z_L - Z_M$')
    ax.set_xlabel(r'$\tau_L - \tau_M$')
    r = np.hstack(rlist)
    scat = ax.scatter(age,metal,c=r,s=45,linewidth=0, cmap=plt.cm.gnuplot2, vmin=0, vmax=vmax)
    cb = fig.colorbar(scat)
    cb.set_label(r'$r$')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def mass_diff(output, exclude=pap_ex, vmax=20, age_weight='L'):

    
    basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'
    #basedir = '/Users/Arthur/Documents/School/final_results'

    age_list = []
    metal_list = []
    rlist = []

    for p in range(6):
        coef = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
        cc = pyfits.open(coef)[1].data
        exarr = np.array(exclude[p]) - 1
        mdiff = cc['MLWZ'] - cc['MMWZ']
        adiff = cc['M{}WA'.format(age_weight)]
        age_list.append(np.delete(adiff,exarr))
        metal_list.append(np.delete(mdiff,exarr))

        rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
        rr = np.loadtxt(rphi, usecols=(1,), unpack=True)
        rlist.append(np.delete(rr,exarr))
        
    age = np.hstack(age_list)
    metal = np.hstack(metal_list)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$Z_L - Z_M$')
    ax.set_xlabel('M{}WA'.format(age_weight))
    r = np.hstack(rlist)
    scat = ax.scatter(age,metal,c=r,s=45,linewidth=0, cmap=plt.cm.gnuplot2, vmin=0, vmax=vmax)
    cb = fig.colorbar(scat)
    cb.set_label(r'$r$')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return
