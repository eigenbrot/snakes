import sys
import time
import numpy as np
import pyfits
import scipy.ndimage as spnd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def plot_chi(chifile, coeffile, datafile,
             output=None, wavemin=3800, wavemax=6800,
             plotblue = False):

    if output == None:
        if plotblue:
            output = chifile.split('.')[0]+'.chi.blue.pdf'
        else:
            output = chifile.split('.')[0]+'.chi.pdf'

    hdu = pyfits.open(datafile)[0]
    numfibers, wavesize = hdu.data.shape

    head = hdu.header
    cdelt = head['CDELT1']
    crval = head['CRVAL1']
    crpix = head['CRPIX1']

    print 'CDELT1 = ', cdelt
    print 'CRVAL1 = ', crval
    print 'CRPIX1 = ', crpix

    wave = (np.arange(wavesize) - crpix) * cdelt + crval

    idx = np.where((wave >= wavemin) & (wave <= wavemax))
    redwl = wave[idx]
    
    redchiarray = pyfits.open(chifile)[0].data

    coeff = pyfits.open(coeffile)[1].data
    meanZ = np.mean(coeff['VSYS'])/3e5 + 1

    # Shift all chi's to rest frame
    restwl = redwl/(meanZ/3e5 + 1)
    chiarray = np.zeros(redchiarray.shape)
    for i in range(redchiarray.shape[0]):
        tl = redwl/(coeff['VSYS'][i]/3e5 + 1)
        chiarray[i,:] = np.interp(restwl,tl,redchiarray[i,:])        

    pp = PDF(output)

    medchi = np.median(chiarray,axis=0)
    mschi = spnd.filters.median_filter(medchi,50)
    rms = np.sqrt(np.mean((chiarray - medchi[None,:])**2,axis=0))
    
    nidx = np.isfinite(medchi)

    medchi = medchi[nidx]
    rms = rms[nidx]
    mschi = mschi[nidx]
    restwl = restwl[nidx]

    mchi = mschi
    stdchi = rms
    # mchi = spnd.filters.gaussian_filter(medchi,5)
    # stdchi = spnd.filters.gaussian_filter(rms,5)

    sk2 = np.array([6300., 5890., 5577.])
#    em2 = np.array([6563.8,  4861., 4959., 5006.8, 6716.0, 6583.41, 6548.04])
    em2 = np.array([6563.8, 6716.0, 6583.41, 6548.04])

    try:
        #VELSTART will be constant across the coeff array, so we'll just take
        #the first one
        em2 *= (coeff['VELSTART'][0]/3e5 + 1.)
    except KeyError:
        pass
    
    dz = 1000. / 3e5
    dzsk = 1500. / 3e5
    
    quality = np.ones(restwl.size)
    for ee in em2:
        maskout = np.where((restwl > ee*(1-dz)) & (restwl < ee*(1+dz)))
        quality[maskout] = 0
        
    for ss in sk2:
        maskout = np.where((restwl > ss*(1-dzsk)) & (restwl < ss*(1+dzsk)))
        quality[maskout] = 0
        
    ok = quality == 1

    if plotblue:
        pidx = np.where(restwl < 4500.)
        bidx = np.where((restwl > 5000.) & (restwl < 5330.))
        bwave = restwl[bidx]
        bxmin = bwave.min() - 20.
        bxmax = bwave.max() + 20.
        fbox = [0.1,0.5,0.56,0.4]
        bbox = [0.67,0.5,0.28,0.4]
    else:
        pidx = np.where(restwl == restwl)
        fbox = [0.1,0.5,0.85,0.4]

    pwave = restwl[pidx]
    xmin = pwave.min() - 20.
    xmax = pwave.max() + 20.
    
    fig = plt.figure(figsize=(11,8))
    rmax = fig.add_axes(fbox)
    rmax.set_ylabel('<Chi> - Med(<Chi>)')
    rmax.set_xlim(xmin,xmax)
    rmax.set_xticklabels([])
    rmax.set_ylim(-5,5)

    medax = fig.add_axes([fbox[0],0.1,fbox[2],0.4])
    medax.set_xlabel('Wavelength [$\AA$]')
    medax.set_ylabel('Median smoothed Chi')
    medax.set_xlim(*rmax.get_xlim())
    medax.set_ylim(-5,5)
    
    if plotblue:
        brmax = fig.add_axes(bbox)
        brmax.set_xticks([5000,5100,5200,5300])
        brmax.set_xticklabels([])
        brmax.set_yticklabels([])
        brmax.set_xlim(bxmin,bxmax)
        brmax.set_ylim(*rmax.get_ylim())
        brmax.spines['left'].set_visible(False)
        brmax.yaxis.tick_right()
        
        rmax.spines['right'].set_visible(False)
        rmax.yaxis.tick_left()

        bmedax = fig.add_axes([bbox[0],0.1,bbox[2],0.4])
        bmedax.set_xticks([5000,5100,5200,5300])
        bmedax.set_yticklabels([])
        bmedax.set_xlim(bxmin,bxmax)
        bmedax.set_ylim(*medax.get_ylim())
        bmedax.spines['left'].set_visible(False)
        bmedax.yaxis.tick_right()
        
        medax.spines['right'].set_visible(False)
        medax.yaxis.tick_left()

        medax.set_xlabel('')
        medax.text(0.76, -0.3, 'Wavelength [$\AA$]', ha='center', va='center', transform=medax.transAxes)

    prms = medchi - mschi
    mrms = np.copy(prms)
    mrms[ok] = np.NAN
    prms[~ok] = np.NAN

    rmax.plot(pwave, prms[pidx], 'k')
    rmax.plot(pwave, mrms[pidx], 'c', lw=1)

    if plotblue:
        brmax.plot(bwave, prms[bidx], 'k')
        brmax.plot(bwave, mrms[bidx], 'c', lw=1)
        
    pmchi = np.copy(mchi)
    mmchi = np.copy(mchi)
    mmchi[ok] = np.NAN
    pmchi[~ok] = np.NAN
    
    medax.plot(pwave,pmchi[pidx],'k')
    medax.plot(pwave,mmchi[pidx],'c',lw=3)
    medax.fill_between(pwave, (mchi - stdchi)[pidx],
                       (mchi + stdchi)[pidx],
                       color='k', alpha=0.2, edgecolor=None)

    if plotblue:
        bmedax.plot(bwave,pmchi[bidx],'k')
        bmedax.plot(bwave,mmchi[bidx],'c',lw=3)
        bmedax.fill_between(bwave, (mchi - stdchi)[bidx],
                           (mchi + stdchi)[bidx],
                           color='k', alpha=0.2, edgecolor=None)        

    sk = np.array([6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0])
    sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
    em = np.array([6563.8,  6716.0])
    emnam = [r'H$\alpha$', 'S2']

    ab = np.array([3820.4, 3835.4,      3889.0,     3933.7, 3968.5, 3970.18,         4304.4,   4341.,       5175.3, 5894.0, 4861.,  4102., 3820.4])
    absnam = ['L',   r'H$\eta$', r'H$\zeta$', 'K',   'H'   , r'H$\epsilon$',    'G',     r'H$\gamma$',  'MgI',   'Na',   r'H$\beta$',   r'H$\delta$',  'L']

    
    #######################
    #######################
    tlim1 = xmin + 20
    tlim2 = xmax - 20
    if plotblue:
        btlim1 = bxmin + 20
        btlim2 = bxmax - 20

    ypos = 1

    for s, sn in zip(sk/meanZ, sknam):
        tidx = np.where((restwl >= s - 10) & (restwl <= s + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        if not np.isfinite(ypos):
            ypos = 9
        if s < tlim2 and s > tlim1:
            rmax.text(s, ypos, sn, fontsize=8, ha='center', va='center')
        if plotblue:
            if s < tlim2 and s > tlim1:
                rmax.axvline(s, color='k', ls=':', alpha=0.7)
            if s > btlim1 and s < btlim2:
                brmax.text(s, ypos, sn, fontsize=8, ha='center', va='center')
                brmax.axvline(s, color='k', ls=':', alpha=0.7)
        else:
            rmax.plot((s,s), (ypos - 0.5, ypos - 1), alpha=0.8, color='k')

    prevy = 99
    for a, an in zip(ab, absnam):
        tidx = np.where((restwl >= a - 10) & (restwl <= a + 10.))
        try:
            ypos = np.min(mchi[tidx]) - 2
        except ValueError:
            pass
        if (an == r'H$\gamma$' or 
            an == r'H$\eta$' or
            an == r'H$\epsilon$') and np.abs(ypos - prevy) <= 0.5:
            ypos -= 1
        prevy = ypos
        if np.isnan(ypos) or ypos < rmax.get_ylim()[0]:
            ypos = rmax.get_ylim()[0] + 0.5
        if a > tlim1 and a < tlim2:
            rmax.text(a, ypos, an, color='r', fontsize=8, ha='center', va='center')
        if plotblue:
            if a > tlim1 and a < tlim2:
                rmax.axvline(a, color='r', ls=':', alpha=0.7)
            if a > btlim1 and a < btlim2:
                brmax.text(a, ypos, an, color='r', fontsize=8, ha='center', va='center')
                brmax.axvline(a, color='r', ls=':', alpha=0.7)
        else:
            rmax.plot((a,a), (ypos + 0.5, ypos + 1), color='r', alpha=0.8) 

    for e, en in zip(em, emnam):
        tidx = np.where((restwl >= e - 10) & (restwl <= e + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        if not np.isfinite(ypos):
            ypos = 9
        if e > tlim1 and e < tlim2:
            rmax.text(e, ypos, en, color='b', fontsize=8, ha='center', va='center')
        if plotblue:
            if e > tlim1 and e < tlim2:
                rmax.axvline(e, color='b', ls=':', alpha=0.7)
            if e > btlim1 and e < btlim2:
                brmax.text(e, ypos, en, color='b', fontsize=8, ha='center', va='center')
                brmax.axvline(e, color='b', ls=':', alpha=0.7)
        else:
            rmax.plot((e,e), (ypos - 0.5, ypos - 1), color='b', alpha=0.8)


    fig.suptitle(time.asctime())

    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return 

def parse_input(inputlist):
    
    chifile = inputlist[0]
    coeffile = inputlist[1]
    datafile = inputlist[2]

    kwar = {}

    i = 3
    while i < len(inputlist):
        
        if inputlist[i] == '-o':
            kwar['output'] = inputlist[i+1]
            i += 1
        
        if inputlist[i] == '-w':
            kwar['wavemin'] = inputlist[i+1]
            kwar['wavemax'] = inputlist[i+2]
            i += 2

        if inputlist[i] == '-b':
            kwar['plotblue'] = True
            
        if inputlist[i] == '-n':
            import nice_plots
            nice_plots.format_plots(False)

        i += 1

    return chifile, coeffile, datafile, kwar

if __name__ == '__main__':

    chi,coef,dat,kw = parse_input(sys.argv[1:])
    plot_chi(chi,coef,dat, **kw)
