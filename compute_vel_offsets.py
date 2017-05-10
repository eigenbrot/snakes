import numpy as np
import pyfits
from matplotlib import rc
import matplotlib.pyplot as plt
from scipy.misc import imread
from matplotlib.backends.backend_pdf import PdfPages as PDF
import allfigs

def format_mab(filename):

    P, ap, HaNS, HaNS_e = np.loadtxt(filename, unpack=True, usecols=(0,1,39,40))

    for pp in range(6):

        Pidx = np.where(P == pp+1)
        tap = ap[Pidx]
        tHaNS = HaNS[Pidx]
        tHaNS_e = HaNS_e[Pidx]

        sidx = np.argsort(tap)

        with open('P{}_mab_vel.txt'.format(pp+1), 'w') as f:
            for i in range(tap.size):
                f.write('{:7.2f}{:7.2f}\n'.format(tHaNS[sidx][i], tHaNS_e[sidx][i]))

    return

def format_mab_HbO(filename):

    P, ap, HbO, HbO_e = np.loadtxt(filename, unpack=True, usecols=(0,1,41,42))

    for pp in range(6):

        Pidx = np.where(P == pp+1)
        tap = ap[Pidx]
        tHbO = HbO[Pidx]
        tHbO_e = HbO_e[Pidx]

        sidx = np.argsort(tap)

        with open('P{}_mab_vel_HbO.txt'.format(pp+1), 'w') as f:
            for i in range(tap.size):
                f.write('{:7.2f}{:7.2f}\n'.format(tHbO[sidx][i], tHbO_e[sidx][i]))

    return

def get_data(workingdir='.',velstr='_mab_vel', velocity_dir='/Users/Arthur/Documents/School/891_research/4166_velocity_dir'):

    #Read in velocity data
    Glist = []
    Slist = []
    sizelist = []
    rlist = []
    zlist = []
    for p in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(workingdir,p+1)
        print loc
        Sfits = '{}/NGC_891_P{}_bin30_allz2.coef.prechi.fits'.format(velocity_dir,p+1)
        print Sfits
        chivfits = '{}/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(velocity_dir,p+1)
        print chivfits
        Gdat = '{}/P{}{}.txt'.format(velocity_dir,p+1,velstr)
        print Gdat
        Gv, Ge = np.loadtxt(Gdat,usecols=(0,1),unpack=True)

        size, r, z = np.loadtxt(loc, usecols=(1,4,5), unpack=True)
        z = np.abs(z)
        
        Sv = pyfits.open(Sfits)[1].data['VSYS']
        cv = pyfits.open(chivfits)[1].data['VSYS']
        Sv += cv
        print Gv.size, Sv.size
        if Gv.size != Sv.size:
            raw_input("WHAAA?")
            Sv = Sv[:Gv.size]
            size = size[:Gv.size]
            r = r[:Gv.size]
            z = z[:Gv.size]
        
        Glist.append(Gv)
        Slist.append(Sv)
        sizelist.append(size)
        rlist.append(r)
        zlist.append(z)

    G = np.hstack(Glist)
    S = np.hstack(Slist)
    Size = np.hstack(sizelist)
    R = np.hstack(rlist)
    Z = np.hstack(zlist)

    return G, S, Size, R, Z

def compute_single_offset(workingdir='.',velstr='_mab_vel'):

    G, S, _, r, z = get_data(workingdir,velstr)
    
    diff = G - S
    print np.median(diff), np.mean(diff), np.std(diff)

    return G, S, np.median(diff), r, z

def compute_fibsize_offset(workingdir='.',velstr='_mab_vel',velocity_dir='/Users/Arthur/Documents/School/891_research/4166_velocity_dir'):

    G, S, Size, r, z = get_data(workingdir,velstr,velocity_dir)

    diff = G - S
    offset = np.zeros(diff.size)
    for s in np.unique(Size):

        idx = np.where(Size == s)
        print s
        print '\t{} {} {}'.format(np.median(diff[idx]), np.mean(diff[idx]), np.std(diff[idx]))
        offset[idx] = np.median(diff[idx])

    return G, S, offset, r, z

def compute_fibsize_pointing_offset(workingdir='.',velstr='_mab_vel',final_results='../final_results'):

    Glist = []
    Slist = []
    offlist = []
    rlist = []
    zlist = []
    for p in range(6):
        print p+1
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(workingdir,p+1)
        print loc
        Sdat = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(final_results,p+1)
        print Sdat
        chivdat = '{}/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(final_results,p+1)
        print chivdat
        Gdat = '{}/P{}{}.txt'.format(workingdir,p+1,velstr)
        print Gdat
        Gv, Ge = np.loadtxt(Gdat,usecols=(0,1),unpack=True)

        size, r, z = np.loadtxt(loc, usecols=(1,2,5), unpack=True)
        z = np.abs(z)
        
        Sv = pyfits.open(Sdat)[1].data['VSYS']
        cv = pyfits.open(chivdat)[1].data['VSYS']
        Sv += cv

        if Gv.size != Sv.size:
            raw_input('WHAA?')
            Sv = Sv[:Gv.size]
            size = size[:Gv.size]
            r = r[:Gv.size]
            z = z[:Gv.size]

        diff = Gv - Sv
        offset = np.zeros(diff.size)
        for s in np.unique(size):
            idx = np.where(size == s)
            print '\t', s
            print '\t\t{} {} {}'.format(np.median(diff[idx]), np.mean(diff[idx]), np.std(diff[idx]))
            offset[idx] = np.median(diff[idx])

        Glist.append(Gv)
        Slist.append(Sv)
        offlist.append(offset)
        rlist.append(r)
        zlist.append(z)
        
    return np.hstack(Glist), np.hstack(Slist), np.hstack(offlist), np.hstack(rlist), np.hstack(zlist)

def plot_vels(gas, stars, offset, r, z, output='test.pdf'):

    gw = 0.5
    gsize = 5
    allfigs.style()

    rc('figure', figsize=(gsize*3.4,gsize*1.6))
    rc('figure.subplot', top=0.90,right=1.0)

    swatfig = '/Users/Arthur/Documents/School/891_research/Swat.png'

    im = imread(swatfig)

    b = 0.1
    h = 0.8
    w = 0.275
    l = 0.08
    pad = 0.01

    fig = plt.figure()
    ax1 = fig.add_axes([l,b,w,h])
    ax2 = fig.add_axes([l+w,b,w,h])
    ax3 = fig.add_axes([l+w*2,b,w,h])
    cax = fig.add_axes([l+w*3+pad,b,0.02,h])
    
    ax1.set_ylabel(r'$\mathrm{Velocity\ [km/s]}$')
    ax2.set_xlabel(r'$\mathrm{Distance\ from\ center\ [arcmin]}$')
    ax2.text(0.5,1.07,r'$\mathrm{[kpc]}$',ha='center',fontsize=22,
             transform=ax2.transAxes)

    ax1.text(-7.5,810,'Gas',fontsize=20)
    ax2.text(-7.5,810,'Stellar',fontsize=20)
    ax3.text(-7.5,810,'Stellar, corrected',fontsize=20)

    axlist = [ax1,ax2,ax3]

    ymin = 215
    ymax = 860
    xmin = -8.33
    xmax = 12.47

    for ax in axlist:
        ax.set_xticks([-8, -4, 0, 4, 8])
        ax.set_yticks([400,600,800])
        ax.axvline(0,ls=':',alpha=0.7,color='k')
        ax.axhline(528,ls=':',alpha=0.7,color='k')
        ax.imshow(im,zorder=0, extent=[xmin,xmax,ymin,ymax], aspect='auto')

    ax1.scatter(r/60., gas, c=z, cmap=plt.cm.gnuplot,linewidth=0,s=25,vmin=0, vmax=2.6)        
    ax2.scatter(r/60., stars, c=z, cmap=plt.cm.gnuplot,linewidth=0,s=25,vmin=0, vmax=2.6)
    scat = ax3.scatter(r/60., stars + offset, c=z, cmap=plt.cm.gnuplot,linewidth=0,s=25,vmin=0, vmax=2.6)

    fig.subplots_adjust(wspace=0.0001)
    cb = fig.colorbar(scat, cax=cax)
    cb.set_alpha(1)
    cb.set_label(r'$\mathrm{|z| [kpc]}$')
    cb.draw_all()

    for ax in axlist:
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        kax = ax.twiny()
        kax.set_xlim(*ax.get_xlim()*np.array([2.91,2.91]))
        
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def plot_diffs(gas, stars, offset, r, z, size=None, output='test.pdf', title=''):

    allfigs.style()

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('|z| [kpc]')
    ax.set_ylabel('dV')
    ax.set_title(title)

    colorlist = ['cyan','green','orange','purple','pink']
    for i, s in enumerate(np.unique(size)):
        idx = np.where(size == s)
        ax.scatter(z[idx], (gas - (stars + offset))[idx], c=colorlist[i], linewidth=0)

    ax.set_xscale('log')

    ax.set_xlim(0.01, 3)
    ax.set_ylim(-300, 300)
    
    pp = PDF(output)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)
    
    return

def plot_diffs_r(gas, stars, offset, r, z, size=None, output='test.pdf', title=''):

    allfigs.style()

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(r'$r_{proj}$ [kpc]')
    ax.set_ylabel('dV')
    ax.set_title(title)

    colorlist = ['cyan','green','orange','purple','pink']
    for i, s in enumerate(np.unique(size)):
        idx = np.where(size == s)
        ax.scatter(r[idx], (gas - (stars + offset))[idx], c=colorlist[i], linewidth=0, label="{:2.1f}".format(s*s))

    ax.legend(loc=0,numpoints=1,scatterpoints=1)
    # ax.set_xscale('log')

    # ax.set_xlim(0.01, 3)
    ax.set_ylim(-300, 300)
    
    pp = PDF(output)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)
    
    return

def gogo():

    _, _, size, _, _ = get_data()
    sing = compute_single_offset()
    fib = compute_fibsize_offset()
    point = compute_fibsize_pointing_offset()
    
    plot_diffs(*sing, size=size, output='diff_single.pdf', title='single offset')
    plot_diffs(*fib, size=size, output='diff_fibersize.pdf', title='fibersize offset')
    plot_diffs(*point, size=size, output='diff_pointing.pdf', title='fibersize and pointing offset')

    plot_diffs_r(*sing, size=size, output='diff_r_single.pdf', title='single offset')
    plot_diffs_r(*fib, size=size, output='diff_r_fibersize.pdf', title='fibersize offset')
    plot_diffs_r(*point, size=size, output='diff_r_pointing.pdf', title='fibersize and pointing offset')

    extinction_diff(*sing[:3],size=size,outpre='Avcomp_single',title='single offset')
    extinction_diff(*fib[:3],size=size,outpre='Avcomp_fibersize',title='fibersize offset')
    extinction_diff(*point[:3],size=size,outpre='Avcomp_pointing',title='fibersize and pointing offset')
    extinction_diff(sing[0], sing[1], 0, size=size,outpre='Avcomp_none',title='no offset')

    return

def gogo_HbO():

    G, _, hbsize, _, _ = get_data(velstr='_mab_vel_HbO')
    singB = compute_single_offset()
    fibB = compute_fibsize_offset()
    pointB = compute_fibsize_pointing_offset()

    singB = list(singB)
    fibB = list(fibB)
    pointB = list(pointB)
    
    singB[0] = G
    fibB[0] = G
    pointB[0] = G
    
    plot_diffs(*singB, size=hbsize, output='diff_HbO_single.pdf', title='HbO, single offset')
    plot_diffs(*fibB, size=hbsize, output='diff_HbO_fibersize.pdf', title='HbO, fibersize offset')
    plot_diffs(*pointB, size=hbsize, output='diff_HbO_pointing.pdf', title='HbO, fibersize and pointing offset')

    plot_diffs_r(*singB, size=hbsize, output='diff_HbO_r_single.pdf', title='HbO, single offset')
    plot_diffs_r(*fibB, size=hbsize, output='diff_HbO_r_fibersize.pdf', title='HbO, fibersize offset')
    plot_diffs_r(*pointB, size=hbsize, output='diff_HbO_r_pointing.pdf', title='HbO, fibersize and pointing offset')

    extinction_diff(*singB[:3],size=hbsize,outpre='Avcomp_HbO_single',title='HbO, single offset')
    extinction_diff(*fibB[:3],size=hbsize,outpre='Avcomp_HbO_fibersize',title='HbO, fibersize offset')
    extinction_diff(*pointB[:3],size=hbsize,outpre='Avcomp_HbO_pointing',title='HbO, fibersize and pointing offset')
    extinction_diff(singB[0], singB[1], 0, size=hbsize,outpre='Avcomp_HbO_none',title='HbO, no offset')

    return

def extinction_diff(gas, stars, offset, size=None, outpre='test', title=''):

    Av_g, Av_s = np.loadtxt('fpp2e.out2', usecols=(8,17), unpack=True)

    axg = plt.figure().add_subplot(111)
    axg.set_xlabel('dV')
    axg.set_ylabel('$A_{V,g}$')
    axg.set_title(title)
    
    axs = plt.figure().add_subplot(111)
    axs.set_xlabel('dV')
    axs.set_ylabel('$A_{V,*}$')
    axs.set_title(title)

    colorlist = ['cyan','green','orange','purple','pink']
    for i, s in enumerate(np.unique(size)):
        idx = np.where(size == s)
        axg.scatter((gas - (stars + offset))[idx], Av_g[idx], c=colorlist[i], linewidth=0, label="{:2.1f}".format(s*2))
        axs.scatter((gas - (stars + offset))[idx], Av_s[idx], c=colorlist[i], linewidth=0, label="{:2.1f}".format(s*2))

    axg.legend(loc=0,numpoints=1,scatterpoints=1)
    axs.legend(loc=0,numpoints=1,scatterpoints=1)

    axg.set_xlim(-300,300)
    axg.set_ylim(-2,10)
    axs.set_xlim(-300,300)
    axs.set_ylim(-2,10)
    
    pp = PDF(outpre+'_gas.pdf')
    pp.savefig(axg.figure)
    pp.close()
    plt.close(axg.figure)

    pp = PDF(outpre+'_stars.pdf')
    pp.savefig(axs.figure)
    pp.close()
    plt.close(axs.figure)

    return
