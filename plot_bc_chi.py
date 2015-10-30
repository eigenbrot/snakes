import sys
import numpy as np
import pyfits
import scipy.ndimage as spnd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def plot_chi(chifile, datafile, output=None, wavemin=3800, wavemax=6800):

    if output == None:
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


    fig = plt.figure(figsize=(11,8))
    rmax = fig.add_subplot(211)
    rmax.set_ylabel('<Chi> - Med(<Chi>)')
    rmax.set_xlim(wavemin,wavemax)
    rmax.set_xticklabels([])
    medax = fig.add_subplot(212)
    medax.set_xlabel('Wavelength [$\AA$]')
    medax.set_ylabel('Median smoothed Chi')
    medax.set_xlim(wavemin, wavemax)
    medax.set_ylim(-5,5)
    
    rmax.plot(restwl, medchi - mschi, 'k')

    medax.plot(restwl,mchi,'k')
    medax.fill_between(restwl, mchi - stdchi, mchi + stdchi,
                    color='k', alpha=0.2, edgecolor=None)

    sk =    [6300.,        5890., 5683.8, 5577.,      5461., 5199.,      4983., 4827.32, 4665.69, 4420.23, 4358., 4165.68, 4047.0]
    sknam = ['[OI] (atm)', 'NaD', 'NaI',  'OI (atm)', 'HgI', 'NI (atm)', 'NaI', 'HgI',   'NaI',   'NaI',   'HgI', 'NaI',   'HgI']
    em = [6563.8,  6716.0]
    emnam = [r'H$\alpha$', 'S2']

    ab =    [3933.7, 3968.5, 4304.4,   5175.3, 5894.0, 4861., 4341., 4102.]
    absnam = ['K',    'H',    'G band', 'Mg',   'Na',   r'H$\beta$',  r'H$\gamma$',  r'H$\delta$']

    ypos = 1
    for s, sn in zip(sk, sknam):
        tidx = np.where((restwl >= s - 10) & (restwl <= s + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        rmax.text(s, ypos, sn, fontsize=8, ha='center', va='center')

    for a, an in zip(ab, absnam):
        tidx = np.where((restwl >= a - 10) & (restwl <= a + 10.))
        try:
            ypos = np.min(mchi[tidx]) - 2
        except ValueError:
            pass
        rmax.text(a, ypos, an, color='r', fontsize=8, ha='center', va='center')

    for e, en in zip(em, emnam):
        tidx = np.where((restwl >= e - 10) & (restwl <= e + 10.))
        try:
            ypos = np.max(mchi[tidx]) + 3
        except ValueError:
            pass
        rmax.text(e, ypos, en, color='b', fontsize=8, ha='center', va='center')

    fig.subplots_adjust(hspace=0.0001)

    pp.savefig(fig)
    pp.close()
    plt.close(fig)

    return 

def parse_input(inputlist):
    
    chifile = inputlist[0]
    datafile = inputlist[1]

    kwar = {}

    i = 2
    while i < len(inputlist):
        
        if inputlist[i] == '-o':
            kwar['output'] = inputlist[i+1]
            i += 1
        
        if inputlist[i] == '-w':
            kwar['wavemin'] = inputlist[i+1]
            kwar['wavemax'] = inputlist[i+2]
            i += 2

        i += 1

    return chifile, datafile, kwar

if __name__ == '__main__':

    chi, dat, kw = parse_input(sys.argv[1:])
    plot_chi(chi,dat, **kw)
