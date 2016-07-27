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
                   label=None, plotfid=True, suffix='final', scale=1., componentfile=None,
                   boxcomponent=True):

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    rmax = 20

    if componentfile and not boxcomponent:
        Psub, Asub = np.loadtxt(componentfile, usecols=(1,2), unpack=True)

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
        
        if componentfile and not boxcomponent:
            Pidx = np.where(Psub == p+1)[0]
            subidx = Asub[Pidx].astype(np.int) - 1
            print subidx
            r = r[subidx]
            z = z[subidx]
            rfull = rfull[subidx]
            c = c[subidx]
        else:
            exarr = np.array(exclude[p]) - 1
            r = np.delete(r, exarr)
            z = np.delete(z, exarr)
            c = np.delete(c, exarr)
            rfull = np.delete(rfull, exarr)

        posidx = np.where(rfull >= 0)[0]
        negidx = np.where(rfull < 0)[0]

        norm1 = plt.Normalize(0,2.6)
        norm2 = plt.Normalize(0,rmax)

        scat1 = ax1.scatter(r[posidx], c[field][posidx]*scale, linewidth=0, alpha=1, c=z[posidx], 
                            vmin=0, vmax=2.6, cmap=plt.cm.gnuplot)
        p1 = ax1.scatter(r[negidx], c[field][negidx]*scale, linewidth=1.1, facecolor='w',alpha=0.7)

        scat2 = ax2.scatter(z[posidx], c[field][posidx]*scale, linewidth=0, alpha=1, c=r[posidx], 
                            vmin=0, vmax=rmax, cmap=plt.cm.gnuplot)
        p2 = ax2.scatter(z[negidx], c[field][negidx]*scale, linewidth=1, facecolor='w',alpha=0.7)

        p1.set_edgecolors(plt.cm.gnuplot(norm1(z[negidx])))
        p2.set_edgecolors(plt.cm.gnuplot(norm2(r[negidx])))

        if errfield is not None and c.size > 0:
            _,_, c1 = ax1.errorbar(r, c[field]*scale, yerr=c[errfield]*scale,
                                   capsize=0, fmt='none', alpha=1, zorder=0)
            _,_, c2 = ax2.errorbar(z, c[field]*scale, yerr=c[errfield]*scale,
                                   capsize=0, fmt='none', alpha=1, zorder=0)
            c1[0].set_color(plt.cm.gnuplot(norm1(z)))
            c2[0].set_color(plt.cm.gnuplot(norm2(r)))
        
    if componentfile and boxcomponent:
        add_IRAF_component(componentfile, ax1, ax2, field)

    if label is None:
        label = field

    ax1.set_xlabel(r'$r_\mathrm{proj}\mathrm{\ [kpc]}$')
    if rtrue:
        ax1.set_xlabel(r'$r\mathrm{\ [kpc]}$')
        ax1.set_xlim(-0.5,20)

    ax1.set_ylabel(label)
    ax1.set_xlim(0,22)
    ax2.set_xlim(-0.2,2.6)
    ax2.set_ylim(*ax1.get_ylim())
    ax2.set_yticklabels([])
    ax2.set_xlabel(r'$|z|\mathrm{\ [kpc]}$')
    
    if plotfid:
        ax1.axvline(3, color='k', ls=':', alpha=0.7)
        ax1.axvline(8, color='k', ls=':', alpha=0.7)
        ax2.axvline(0.4, color='k', ls=':', alpha=0.7)
        ax2.axvline(1, color='k', ls=':', alpha=0.7)

    fig.subplots_adjust(wspace=0.0001)
    
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    cax1 = fig.add_axes([pos1.x0 + (pos1.width-0.3)/2,0.88,0.3,0.03])
    cax2 = fig.add_axes([pos2.x0 + (pos2.width-0.3)/2,0.88,0.3,0.03])
    cb1 = fig.colorbar(scat1, cax=cax1, orientation='horizontal')
    cb2 = fig.colorbar(scat2, cax=cax2, orientation='horizontal') 
    cb1.set_ticks([0,0.5,1,1.5,2,2.5])
    cb2.set_ticks([0,4,8,12,16,20])
    cax1.text(0.5,1.3,r'$|z|\mathrm{\ [kpc]}$',va='bottom', ha='center', 
              transform=cax1.transAxes, fontsize=20)
    cax2.text(0.5,1.3,r'$r\mathrm{\ [kpc]}$', va='bottom', ha='center',
              transform=cax2.transAxes, fontsize=20)

    # if errfield is not None:
    #     for i in range(6):
    #         elist1[i][0].set_color(cb1.to_rgba(z))
    #         elist2[i][0].set_color(cb2.to_rgba(r))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def add_IRAF_component(componentfile, rax, zax, field):

    cold = {'MLWA': 7, 'MLWZ': 8, 'TAUV': 9}
    
    z, r, data = np.loadtxt(componentfile, usecols=(4,5,cold[field]), unpack=True)

    rax.scatter(r, data, color='k', marker='s', s=40, facecolor='none', alpha=0.7)
    zax.scatter(z, data, color='k', marker='s', s=40, facecolor='none', alpha=0.7)

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
    bige0 = np.hstack(bige0_list)/bigf0
    bigf1 = np.hstack(bigf1_list)
    bige1 = np.hstack(bige1_list)/bigf1

    ax1 = plt.figure().add_subplot(111)
    ax1.set_xlabel('p=0')
    ax1.set_ylabel('p=1')
    ax1.set_title(field)

    ax2 = plt.figure().add_subplot(111)
    ax2.set_xlabel(errfield)
    ax2.set_ylabel('N')
        
    pe0 = bige0[bige0 < 1]
    pe1 = bige1[bige1 < 1]
        
    ax1.errorbar(bigf0,bigf1,yerr=bige1*bigf1,xerr=bige0*bigf0,c='k',capsize=0,ls='',marker='o')
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

def weight_comp_plot(outpre, bins=20, basedir='.'):

    fields = ['TAUV','MLWA','MLWZ']
    bigD = {}
    colors = ['#1b9e77','#d95f02']

    for f in fields:
        bigD[f] = {'d0':[], 'e0': [], 'd1': [], 'e1': []}

    for p in range(6):
        
        basename = 'NGC_891_P{}_bin30_allz2.fiterr.fits'.format(p+1)
        p0 = pyfits.open(basedir+'/'+basename)[1].data
        p1 = pyfits.open(basedir+'/IRAF_weight/'+basename)[1].data

        for f in fields:
            bigD[f]['d0'].append(p0[f])
            bigD[f]['e0'].append(p0['d'+f])
            bigD[f]['d1'].append(p1[f])
            bigD[f]['e1'].append(p1['d'+f])

    for f in fields:
        bigD[f]['d0'] = np.hstack(bigD[f]['d0'])
        bigD[f]['e0'] = np.hstack(bigD[f]['e0'])
        bigD[f]['d1'] = np.hstack(bigD[f]['d1'])
        bigD[f]['e1'] = np.hstack(bigD[f]['e1'])

    sfig = plt.figure()
    sax1 = sfig.add_subplot(131)
    sax2 = sfig.add_subplot(132)
    sax3 = sfig.add_subplot(133)
    hfig = plt.figure()
    hax1 = hfig.add_subplot(131)
    hax2 = hfig.add_subplot(132)
    hax3 = hfig.add_subplot(133)

    for f, ax in zip(fields, [sax1, sax2, sax3]):
        xval = bigD[f]['d0']
        yval = bigD[f]['d1']
        if f == 'TAUV':
            xval *= 1.086
            yval *= 1.086
        lx = np.linspace(xval.min(), xval.max(), 10)
        ax.scatter(xval, yval, c='k', alpha=0.5, edgecolor='none', s=55)
        ax.plot(lx,lx,ls='--',color='k',alpha=0.9, lw=2, zorder=0)

    for f, ax in zip(fields, [hax1, hax2, hax3]):
        e0 = bigD[f]['e0']
        e1 = bigD[f]['e1']
        if f != 'TAUV':
            e0 /= bigD[f]['d0']
            e1 /= bigD[f]['d1']
        e0 = e0[e0 == e0]
        e1 = e1[e1 == e1]

        ax.hist(e0, bins=bins, color=colors[0], label='No weight', histtype='stepfilled', alpha=0.5)
        ax.hist(e1, bins=bins, color=colors[1], label='Weighted', histtype='stepfilled', alpha=0.5)
        ax.set_xlim(0,0.68)
        ax.set_ylim(0,70)

    sax2.set_xlabel('No weight')
    sax1.set_ylabel('Weighted')
    sax1.text(0.15,0.87,'$A_V$', ha='center', va='bottom', transform=sax1.transAxes, fontsize=30)
    sax2.text(0.15,0.87,r'$\tau_L$', ha='center', va='bottom', transform=sax2.transAxes, fontsize=30)
    sax3.text(0.15,0.87,r'$Z_L$', ha='center', va='bottom', transform=sax3.transAxes, fontsize=30)

    hax1.set_ylabel('$N$')
    hax1.set_xlabel(r'$\delta A_V$')
    hax2.set_xlabel(r'$\delta\tau_L/\tau_L$')
    hax3.set_xlabel(r'$\delta Z_L/Z_L$')
    hax2.legend(loc=0, frameon=False)
    hax2.set_yticklabels([])
    hax3.set_yticklabels([])
    hfig.subplots_adjust(wspace=0.001)

    spp = PDF(outpre+'_sys.pdf')
    spp.savefig(sfig)
    spp.close()
    hpp = PDF(outpre+'_hist.pdf')
    hpp.savefig(hfig)
    hpp.close()
    plt.close(sfig)
    plt.close(hfig)

    return

def SFH_cuts(output, basedir='.', exclude=exclude, rtrue=False,
             rcuts=[3,8], zcuts=[0.4,1], numZ=6, numAge=4, componentfile=None):
    import matplotlib.colorbar as mplcb

    DFK_borders = np.array([[0.9e-3,5.2e-3],
                            [5.5e-3,0.404],
                            [0.4535,5.75],
                            [6,13.5]]) * 1e3
    Zvals = np.array([0.005,0.02,0.2,0.4,1,2.5])

    fig = plt.figure()
    lax = fig.add_subplot(111, label='bigax')
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_ylabel(r'$\mathrm{Fractional }\ f_i(\lambda)\int\ \psi(t) dt$')
    lax.set_xlabel('Lookback time [Myr]')
    lax.tick_params(axis='both',pad=20,length=0)

    bigz = [0] + zcuts + [2.6]
    if rtrue:
        bigr = [0] + rcuts + [22]
    else:
        bigr = [0] + rcuts + [11]
    
    if componentfile:
        Psub, Asub = np.loadtxt(componentfile, usecols=(1,2), unpack=True)

    i = 1
    axlist = [lax]
    for z in range(len(zcuts) + 1):
        zc = [bigz[-z-2], bigz[-z-1]]
        for r in range(len(rcuts) + 1):
            rc = [bigr[r], bigr[r+1]]
            print zc, rc
            bigD_list = []
            bigE_list = []
            bigZ_list = []
            for p in range(6):
                coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(basedir,p+1)
                print coeffile
                coefs = pyfits.open(coeffile)[1].data
                loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
                print loc
                r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
                if rtrue:
                    rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, p+1)
                    r = np.loadtxt(rphi,usecols=(1,),unpack=True)
                
                r = np.abs(r)
                z = np.abs(z)

                if componentfile:
                    Pidx = np.where(Psub == p+1)[0]
                    subidx = Asub[Pidx].astype(np.int) - 1
                    print subidx
                    r = r[subidx]
                    z = z[subidx]
                    lw = coefs['LIGHT_WEIGHT'][subidx,:]
                    lwe = coefs['LIGHT_WEIGHT_ERR'][subidx,:]
                else:
                    exarr = np.array(exclude[p]) - 1
                    r = np.delete(r, exarr)
                    z = np.delete(z, exarr)
                    lw = np.delete(coefs['LIGHT_WEIGHT'], exarr, axis=0)
                    lwe = np.delete(coefs['LIGHT_WEIGHT_ERR'], exarr, axis=0)

                norm = np.sum(lw,axis=1)[:,None]
                lw /= norm
                lwe /= norm

                idx = np.where((z >= zc[0]) & (z < zc[1])
                               & (r >= rc[0]) & (r < rc[1]))[0]
                lw = lw[idx,:]
                lwe = lwe[idx,:]
                print lw.shape
                lw = np.reshape(lw,(idx.size,numZ,numAge))
                lwe = np.reshape(lwe,(idx.size,numZ,numAge))
                print lw.shape
                print np.sum(lw,axis=(1,2))
                bigD_list.append(np.sum(lw,axis=1))
                bigE_list.append(np.sqrt(np.sum(lwe**2,axis=1)))
                bigZ_list.append(np.sum(lw * Zvals[None,:,None],axis=1)/(np.sum(lw,axis=1)+0.0001))

            bigDarr = np.vstack(bigD_list)
            bigEarr = np.vstack(bigE_list)
            bigZarr = np.vstack(bigZ_list)
            print bigDarr.shape, bigEarr.shape, bigZarr.shape
            print np.mean(bigDarr,axis=0)
            print np.sum(np.mean(bigDarr,axis=0))
            bigD = np.mean(bigDarr,axis=0)
            bigE = np.sqrt(np.sum(bigEarr**2,axis=0))/(bigEarr.shape[0])
            bigZ = np.mean(bigZarr,axis=0)
            print bigD.shape, bigE.shape, bigZ.shape
            print bigZ
            norm = plt.Normalize(0,2.5)
            ax = fig.add_subplot(len(zcuts)+1,len(rcuts)+1,i)
            # ax.hlines(bigD, DFK_borders[:,0], DFK_borders[:,1],color='k', lw=2)
            # ax.hlines(bigD - bigE, DFK_borders[:,0], DFK_borders[:,1],colors='k', lw=1,linestyles=':')
            # ax.hlines(bigD + bigE, DFK_borders[:,0], DFK_borders[:,1],colors='k', lw=1,linestyles=':')
            for d in range(numAge):
                ax.fill_between(DFK_borders[d],[bigD[d] + bigE[d]]*2,[bigD[d] - bigE[d]]*2, 
                                color=plt.cm.gnuplot(norm(bigZ[d])), alpha=1)
            
            ax.set_xscale('log')
            ax.set_xlim(0.3e0,4e4)
            ax.set_ylim(-0.001,0.98)

            if i <= (len(zcuts)) * (len(rcuts) + 1):
                ax.set_xticklabels([])
            if len(rcuts) > 0 and i % (len(rcuts)+1) != 1:
                ax.set_yticklabels([])

            if  len(rcuts) == 0 or i % (len(rcuts) + 1) == 0:
                ax.text(1.09,0.5,'${}\leq |z| <{}$ kpc'.format(*zc),rotation=90,
                        ha='center',va='center',transform=ax.transAxes)

            if i <= len(rcuts) + 1:
                if rtrue:
                    rlab = '${}\leq |r| <{}$ kpc'.format(*rc)
                else:
                    rlab = '${}\leq |r_\mathrm{{proj}}| <{}$ kpc'.format(*rc)
                ax.text(0.5,1.07,rlab,
                        ha='center',va='center',transform=ax.transAxes)

            axlist.append(ax)
            i += 1

    fig.subplots_adjust(hspace=0.00001,wspace=0.0001)
    cax, _ = mplcb.make_axes(axlist,pad=0.07,aspect=30)
    cb = mplcb.ColorbarBase(cax, cmap=plt.cm.gnuplot, norm=norm)
    cb.set_label(r'$Z_L$')

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
            
    return

def age_Z_cuts(output, basedir='.', exclude=exclude, rtrue=True,
               rcuts=[3,8], zcuts=[0.4,1], suffix='fiterr'):
    
    fig = plt.figure()
    lax = fig.add_subplot(111, label='bigax')
    lax.spines['top'].set_visible(False)
    lax.spines['right'].set_visible(False)
    lax.spines['bottom'].set_visible(False)
    lax.spines['left'].set_visible(False)   
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.set_ylabel(r'$Z_L$')
    lax.set_xlabel(r'$\tau_L$')
    lax.tick_params(axis='both',pad=20,length=0)

    bigz = [0] + zcuts + [2.6]
    bigr = [0] + rcuts + [11]
    
    i = 1
    for z in range(len(zcuts) + 1):
        zc = [bigz[-z-2], bigz[-z-1]]
        for r in range(len(rcuts) + 1):
            rc = [bigr[r], bigr[r+1]]
            ax = fig.add_subplot(len(zcuts)+1,len(rcuts)+1,i)
            for p in range(6):
                coeffile = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,p+1,suffix)
                print coeffile
                coefs = pyfits.open(coeffile)[1].data
                loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
                print loc
                r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)
                fullr = r
                if rtrue:
                    rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, p+1)
                    r = np.loadtxt(rphi,usecols=(1,),unpack=True)
                
                r = np.abs(r)
                z = np.abs(z)

                exarr = np.array(exclude[p]) - 1
                r = np.delete(r, exarr)
                fullr = np.delete(fullr, exarr)
                z = np.delete(z, exarr)
                MLWA = np.delete(coefs['MLWA'], exarr, axis=0)
                MLWA_e = np.delete(coefs['dMLWA'], exarr, axis=0)
                MLWZ = np.delete(coefs['MLWZ'], exarr, axis=0)
                MLWZ_e = np.delete(coefs['dMLWZ'], exarr, axis=0)

                idx = np.where((z >= zc[0]) & (z < zc[1])
                               & (r >= rc[0]) & (r < rc[1]))[0]
                fullr = fullr[idx]
                MLWA = MLWA[idx]
                MLWZ = MLWZ[idx]
                MLWA_e = MLWA_e[idx]
                MLWZ_e = MLWZ_e[idx]
                posidx = np.where(fullr >= 0)[0]
                negidx = np.where(fullr < 0)[0]

                ax.scatter(MLWA[posidx], MLWZ[posidx], c='k', linewidths=0, alpha=0.7)
                ax.scatter(MLWA[negidx], MLWZ[negidx], linewidths=1.2, alpha=0.7,
                           facecolors='w')
                ax.errorbar(MLWA, MLWZ, yerr=MLWZ_e, xerr=MLWA_e,
                            fmt='none', capsize=0, ecolor='gray', elinewidth=1,zorder=0)

            ax.set_xlim(0,12.6)
            ax.set_ylim(0,2.6)

            if i % (len(rcuts) + 1) == 0:
                tax = ax.twinx()
                tax.set_ylabel('${}\leq |z| <{}$ kpc'.format(*zc))
                #rotation='horizontal',labelpad=20)
                tax.set_ylim(*ax.get_ylim())
                tax.set_yticklabels([])

            if i <= len(rcuts) + 1:
                tax = ax.twiny()
                tax.set_xlabel('${}\leq |r| <{}$ kpc'.format(*rc))
                #rotation='horizontal',labelpad=20)
                tax.set_xlim(*ax.get_xlim())
                tax.set_xticklabels([])

            if i <= (len(zcuts)) * (len(rcuts) + 1):
                ax.set_xticklabels([])
            if len(rcuts) > 0 and i % (len(rcuts)+1) != 1:
                ax.set_yticklabels([])

            i += 1

    fig.subplots_adjust(wspace=0.0001,hspace=0.0001)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return

def val_pos_2panel_multicomp(output, field, errfield=None, componentlist=[], basedir='.', rtrue=True,
                             label=None, plotfid=True, suffix='fiterr', scale=1.):
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    rmax = 20

    symblist = ['o','^','s']

    for p in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
        rr, zz = np.loadtxt(loc, usecols=(4,5), unpack=True)
        rrfull = rr
        zz = np.abs(zz)
        rr = np.abs(rr)
        if rtrue:
            rphi = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir,p+1)
            rr = np.loadtxt(rphi, usecols=(1,), unpack=True)
        coef = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,p+1,suffix)
        print coef
        cc = pyfits.open(coef)[1].data

        for component, symb in zip(componentlist,symblist):
            Psub, Asub = np.loadtxt(component, usecols=(1,2), unpack=True)
            Pidx = np.where(Psub == p+1)[0]
            subidx = Asub[Pidx].astype(np.int) - 1
            print subidx
            r = rr[subidx]
            z = zz[subidx]
            rfull = rrfull[subidx]
            c = cc[subidx]

            posidx = np.where(rfull >= 0)[0]
            negidx = np.where(rfull < 0)[0]

            norm1 = plt.Normalize(0,2.6)
            norm2 = plt.Normalize(0,rmax)

            scat1 = ax1.scatter(r[posidx], c[field][posidx]*scale, linewidth=0, alpha=1, c=z[posidx], 
                                vmin=0, vmax=2.6, cmap=plt.cm.gnuplot, marker=symb)
            p1 = ax1.scatter(r[negidx], c[field][negidx]*scale, linewidth=1.1, facecolor='w',alpha=0.7, marker=symb)
            
            scat2 = ax2.scatter(z[posidx], c[field][posidx]*scale, linewidth=0, alpha=1, c=r[posidx], 
                                vmin=0, vmax=rmax, cmap=plt.cm.gnuplot, marker=symb)
            p2 = ax2.scatter(z[negidx], c[field][negidx]*scale, linewidth=1, facecolor='w',alpha=0.7, marker=symb)

            p1.set_edgecolors(plt.cm.gnuplot(norm1(z[negidx])))
            p2.set_edgecolors(plt.cm.gnuplot(norm2(r[negidx])))

            if errfield is not None and c.size > 0:
                _,_, c1 = ax1.errorbar(r, c[field]*scale, yerr=c[errfield]*scale,
                                       capsize=0, fmt='none', alpha=1, zorder=0)
                _,_, c2 = ax2.errorbar(z, c[field]*scale, yerr=c[errfield]*scale,
                                       capsize=0, fmt='none', alpha=1, zorder=0)
                c1[0].set_color(plt.cm.gnuplot(norm1(z)))
                c2[0].set_color(plt.cm.gnuplot(norm2(r)))
        

    if label is None:
        label = field

    ax1.set_xlabel(r'$r_\mathrm{proj}\mathrm{\ [kpc]}$')
    if rtrue:
        ax1.set_xlabel(r'$r\mathrm{\ [kpc]}$')
        ax1.set_xlim(-0.5,20)

    ax1.set_ylabel(label)
    ax1.set_xlim(0,22)
    ax2.set_xlim(-0.2,2.6)
    ax2.set_ylim(*ax1.get_ylim())
    ax2.set_yticklabels([])
    ax2.set_xlabel(r'$|z|\mathrm{\ [kpc]}$')
    
    if plotfid:
        ax1.axvline(3, color='k', ls=':', alpha=0.7)
        ax1.axvline(8, color='k', ls=':', alpha=0.7)
        ax2.axvline(0.4, color='k', ls=':', alpha=0.7)
        ax2.axvline(1, color='k', ls=':', alpha=0.7)

    fig.subplots_adjust(wspace=0.0001)
    
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    cax1 = fig.add_axes([pos1.x0 + (pos1.width-0.3)/2,0.88,0.3,0.03])
    cax2 = fig.add_axes([pos2.x0 + (pos2.width-0.3)/2,0.88,0.3,0.03])
    cb1 = fig.colorbar(scat1, cax=cax1, orientation='horizontal')
    cb2 = fig.colorbar(scat2, cax=cax2, orientation='horizontal') 
    cb1.set_ticks([0,0.5,1,1.5,2,2.5])
    cb2.set_ticks([0,4,8,12,16,20])
    cax1.text(0.5,1.3,r'$|z|\mathrm{\ [kpc]}$',va='bottom', ha='center', 
              transform=cax1.transAxes, fontsize=20)
    cax2.text(0.5,1.3,r'$r\mathrm{\ [kpc]}$', va='bottom', ha='center',
              transform=cax2.transAxes, fontsize=20)

    # if errfield is not None:
    #     for i in range(6):
    #         elist1[i][0].set_color(cb1.to_rgba(z))
    #         elist2[i][0].set_color(cb2.to_rgba(r))

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return
