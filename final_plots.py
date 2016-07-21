import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF

exclude = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]

def combine_err():

    for p in range(6):

        basename = 'NGC_891_P{}_bin30_allz2'.format(p+1)
        print basename
        coefs = pyfits.open(basename+'.coef.fits')[1].data
        fiterr = pyfits.open(basename+'.fiterr.fits')[1].data
        syserr = pyfits.open(basename+'.syserr.fits')[1].data

        numap = coefs.shape[0]
        output = np.zeros(numap, dtype=[('MLWA','>f8'), ('dMLWA', '>f8'),
                                        ('MLWZ','>f8'), ('dMLWZ', '>f8'),
                                        ('TAUV','>f8'), ('dTAUV', '>f8')])

        for i in range(numap):
            output['MLWA'][i] = coefs['MLWA'][i]
            output['MLWZ'][i] = coefs['MLWZ'][i]
            output['TAUV'][i] = coefs['TAUV'][i]
            output['dMLWA'][i] = np.sqrt(fiterr['dMLWA'][i]**2 + syserr['dMLWA'][i]**2)
            output['dMLWZ'][i] = fiterr['dMLWZ'][i]
            output['dTAUV'][i] = fiterr['dTAUV'][i]

        outname = basename + '.final.fits'
        pyfits.BinTableHDU(output).writeto(outname, clobber=True)

    return

def val_pos_3panel(output, field='MLWA', colorfield='TAUV', basedir='.', 
                   exclude=exclude, rtrue=True):

    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)

    for p in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
        r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
        z = np.abs(z)
        r = np.abs(r)
        if rtrue:
            rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
            r = np.loadtxt(rphi, usecols=(1,), unpack=True)
        coef = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
        c = pyfits.open(coef)[1].data
        
        exarr = np.array(exclude[p]) - 1
        r = np.delete(r, exarr)
        z = np.delete(z, exarr)
        c = np.delete(c, exarr)

        ax1.scatter(r, z, c=c['TAUV']*1.086, vmin=0,vmax=5, cmap=plt.cm.gnuplot,
                    linewidth=0, alpha=0.6)
        ax2.scatter(r, c[field], c=c['TAUV']*1.086, vmin=0,vmax=5,
                    cmap=plt.cm.gnuplot,linewidth=0, alpha=0.6)
        scat = ax3.scatter(z, c[field], c=c['TAUV']*1.086, vmin=0,vmax=5,
                           cmap=plt.cm.gnuplot,linewidth=0, alpha=0.6)

    ax2.set_xlabel(r'$r_\mathrm{proj}$ [kpc]')
    if rtrue:
        ax2.set_xlabel(r'$r$ [kpc]')
        ax1.set_xlim(-0.5,20)

    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$|z|$ [kpc]')

    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylabel(r'$\tau_L$ [Gyr]')

    ax3.set_ylim(*ax2.get_ylim())
    ax3.set_yticklabels([])
    ax3.set_xlabel(r'$|z|$ [kpc]')

    
    fig.subplots_adjust(hspace=0.001,wspace=0.001)

    cax = fig.add_axes([0.66,0.55,0.04,0.32])
    cb = fig.colorbar(scat, cax=cax)
    cb.set_alpha(1)
    cb.set_label(r'$A_V$')
    cb.draw_all()

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def val_pos_2panel(output, field, errfield=None, basedir='.', exclude=exclude, rtrue=True,
                   label=None, plotfid=True, suffix='final'):

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    rmax = 20

    for p in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
        r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
        rfull = r
        z = np.abs(z)
        r = np.abs(r)
        if rtrue:
            rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
            r = np.loadtxt(rphi, usecols=(1,), unpack=True)
        coef = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,p+1,suffix)
        print coef
        c = pyfits.open(coef)[1].data
        
        exarr = np.array(exclude[p]) - 1
        r = np.delete(r, exarr)
        z = np.delete(z, exarr)
        c = np.delete(c, exarr)
        rfull = np.delete(rfull, exarr)

        posidx = np.where(rfull >= 0)[0]
        negidx = np.where(rfull < 0)[0]

        norm1 = plt.Normalize(0,2.6)
        norm2 = plt.Normalize(0,rmax)

        scat1 = ax1.scatter(r[posidx], c[field][posidx], linewidth=0, alpha=1, c=z[posidx], 
                            vmin=0, vmax=2.6, cmap=plt.cm.gnuplot)
        p1 = ax1.scatter(r[negidx], c[field][negidx], linewidth=1.1, facecolor='w',alpha=0.7)

        scat2 = ax2.scatter(z[posidx], c[field][posidx], linewidth=0, alpha=1, c=r[posidx], 
                            vmin=0, vmax=rmax, cmap=plt.cm.gnuplot)
        p2 = ax2.scatter(z[negidx], c[field][negidx], linewidth=1, facecolor='w',alpha=0.7)

        p1.set_edgecolors(plt.cm.gnuplot(norm1(z[negidx])))
        p2.set_edgecolors(plt.cm.gnuplot(norm2(r[negidx])))

        if errfield is not None:
            _,_, c1 = ax1.errorbar(r, c[field], yerr=c[errfield], capsize=0, fmt='none', alpha=1, zorder=0)
            _,_, c2 = ax2.errorbar(z, c[field], yerr=c[errfield], capsize=0, fmt='none', alpha=1, zorder=0)
            c1[0].set_color(plt.cm.gnuplot(norm1(z)))
            c2[0].set_color(plt.cm.gnuplot(norm2(r)))
        
    if label is None:
        label = field

    ax1.set_xlabel(r'$r_\mathrm{proj}$ [kpc]')
    if rtrue:
        ax1.set_xlabel(r'$r$ [kpc]')
        ax1.set_xlim(-0.5,20)

    ax1.set_ylabel(label)

    ax2.set_ylim(*ax1.get_ylim())
    ax2.set_yticklabels([])
    ax2.set_xlabel(r'$|z|$ [kpc]')
    
    if plotfid:
        ax1.axvline(3, color='k', ls=':', alpha=0.7)
        ax1.axvline(8, color='k', ls=':', alpha=0.7)
        ax2.axvline(0.4, color='k', ls=':', alpha=0.7)
        ax2.axvline(1, color='k', ls=':', alpha=0.7)

    fig.subplots_adjust(hspace=0.001,wspace=0.001)
    
    cax1 = fig.add_axes([0.17,0.95,0.3,0.03])
    cax2 = fig.add_axes([0.55,0.95,0.3,0.03])
    cb1 = fig.colorbar(scat1, cax=cax1, orientation='horizontal')
    cb2 = fig.colorbar(scat2, cax=cax2, orientation='horizontal')
    
    # if errfield is not None:
    #     for i in range(6):
    #         elist1[i][0].set_color(cb1.to_rgba(z))
    #         elist2[i][0].set_color(cb2.to_rgba(r))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def compare_weights(field, errfield, output, bins=10):

    bigf0_list = []
    bige0_list = []
    bigf1_list = []
    bige1_list = []
    for p in range(6):

        basename = 'NGC_891_P{}_bin30_allz2.fiterr.fits'.format(p+1)
        p0 = pyfits.open(basename)[1].data
        p1 = pyfits.open('IRAF_weight/'+basename)[1].data

        bigf0_list.append(p0[field])
        bige0_list.append(p0[errfield])
        bigf1_list.append(p1[field])
        bige1_list.append(p1[errfield])

    bigf0 = np.hstack(bigf0_list)
    bige0 = np.hstack(bige0_list)
    bigf1 = np.hstack(bigf1_list)
    bige1 = np.hstack(bige1_list)

    ax1 = plt.figure().add_subplot(111)
    ax1.set_xlabel('p=0')
    ax1.set_ylabel('p=1')
    ax1.set_title(field)

    ax2 = plt.figure().add_subplot(111)
    ax2.set_xlabel(errfield)
    ax2.set_ylabel('N')
        
    pe0 = bige0[np.isfinite(bige0)]
    pe1 = bige1[np.isfinite(bige1)]
        
    ax1.scatter(bigf0,bigf1,c='k')
    ax1.plot(bigf0,bigf0,color='k',ls=':',alpha=0.5)
    ax2.hist(pe0,bins=bins,color='b',alpha=0.5,histtype='stepfilled',label='p=0')
    ax2.hist(pe1,bins=bins,color='r',alpha=0.5,histtype='stepfilled',label='p=1')
    ax2.legend(loc=0)

    pp = PDF(output)
    pp.savefig(ax1.figure)
    pp.savefig(ax2.figure)
    pp.close()
    plt.close(ax1.figure)
    plt.close(ax2.figure)

    return
