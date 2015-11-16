import sys
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
    restwl = wave[idx]

    chiarray = pyfits.open(chifile)[0].data

    coeff = pyfits.open(coeffile)[1].data
    Z = np.mean(coeff['VSYS'])/3e5 + 1

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
    em2 = np.array([6563.8,  4861., 4959., 5006.8, 6716.0, 6583.41, 6548.04])
    dz = 1500. / 3e5
    dzsk = 1600. / 3e5
    
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
    else:
        pidx = np.where(restwl == restwl)

    fig = plt.figure(figsize=(11,8))
    rmax = fig.add_subplot(211)
    rmax.set_ylabel('<Chi> - Med(<Chi>)')
    rmax.set_xlim(restwl[pidx].min(),restwl[pidx].max())
    rmax.set_xticklabels([])
    rmax.set_ylim(-5,5)
    medax = fig.add_subplot(212)
    medax.set_xlabel('Wavelength [$\AA$]')
    medax.set_ylabel('Median smoothed Chi')
    medax.set_xlim(restwl[pidx].min(),restwl[pidx].max())
    medax.set_ylim(-5,5)
    
    prms = medchi - mschi
    mrms = np.copy(prms)
    mrms[ok] = np.NAN
    prms[~ok] = np.NAN

    rmax.plot(restwl[pidx], prms[pidx], 'k')
    rmax.plot(restwl[pidx], mrms[pidx], 'c', lw=3)

    pmchi = np.copy(mchi)
    mmchi = np.copy(mchi)
    mmchi[ok] = np.NAN
    pmchi[~ok] = np.NAN
    
    medax.plot(restwl[pidx],pmchi[pidx],'k')
    medax.plot(restwl[pidx],mmchi[pidx],'c',lw=3)
    medax.fill_between(restwl[pidx], (mchi - stdchi)[pidx],
                       (mchi + stdchi)[pidx],
                       color='k', alpha=0.2, edgecolor=None)

    sk = np.array([6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0])
    sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
    em = np.array([6563.8,  6716.0])
    emnam = [r'H$\alpha$', 'S2']

    ab = np.array([3820.4, 3835.4,      3889.0,     3933.7, 3968.5, 3970.18,         4304.4,   4341.,       5175.3, 5894.0, 4861.,  4102., 3820.4])
    absnam = ['L',   r'H$\eta$', r'H$\zeta$', 'K',   'H'   , r'H$\epsilon$',    'G',     r'H$\gamma$',  'Mg',   'Na',   r'H$\beta$',   r'H$\delta$',  'L']

    

    ypos = 1
    for s, sn in zip(sk, sknam):
        if s > 5500. and plotblue: continue
        tidx = np.where((restwl >= s - 10) & (restwl <= s + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        if not np.isfinite(ypos):
            ypos = 9
        rmax.text(s, ypos, sn, fontsize=8, ha='center', va='center')
        if plotblue:
            rmax.axvline(s, color='k', ls=':', alpha=0.7)
        else:
            rmax.plot((s,s), (ypos - 0.5, ypos - 1), alpha=0.8, color='k')

    prevy = 99
    for a, an in zip(ab*Z, absnam):
        if a > 5500. and plotblue: continue
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
        rmax.text(a, ypos, an, color='r', fontsize=8, ha='center', va='center')
        if plotblue:
            rmax.axvline(a, color='r', ls=':', alpha=0.7)
        else:
            rmax.plot((a,a), (ypos + 0.5, ypos + 1), color='r', alpha=0.8) 

    for e, en in zip(em*Z, emnam):
        if e > 5500. and plotblue: continue
        tidx = np.where((restwl >= e - 10) & (restwl <= e + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        if not np.isfinite(ypos):
            ypos = 9
        rmax.text(e, ypos, en, color='b', fontsize=8, ha='center', va='center')
        if plotblue:
            rmax.axvline(e, color='b', ls=':', alpha=0.7)
        else:
            rmax.plot((e,e), (ypos - 0.5, ypos - 1), color='b', alpha=0.8)

    fig.subplots_adjust(hspace=0.0001)

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
