import numpy as np
import scipy.interpolate as spi
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
#import D4000_indices as DI
excl = [[5, 34], [1, 2, 35], [59], [2, 8], [1, 2, 3, 27, 28, 29,5], [35, 36, 38]]
def get_all_data(basedir='.'):

    zlist = []
    datalist = []

    for p in range(6):
        loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,p+1)
        fitsfile = '{}/NGC_891_P{}_bin30.msoz.fits'.format(basedir,p+1)

        rr, zz = np.loadtxt(loc,usecols=(4,5), unpack=True)
        rr = np.abs(rr)
        zz = np.abs(zz)

        hdu = pyfits.open(fitsfile)[0]
        ddata = hdu.data
        wave = (np.arange(ddata.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CDELT1'] + hdu.header['CRVAL1']

        exarr = np.array(excl[p]) - 1
        rr = np.delete(rr,exarr)
        zz = np.delete(zz,exarr)
        ddata = np.delete(ddata,exarr,axis=0)

        zlist.append(zz)
        datalist.append(ddata)
        

    z = np.hstack(zlist)
    data = np.vstack(datalist)
    idx = np.argsort(z)
    z = z[idx]
    data = data[idx,:]

    return wave, z, data

def add_line_labels(ax):

    lines = [3820.4,
             3835.4,
             3889.0,
             3933.7,
             3968.5,
#             3970.18,
             4304.4,
             4341.,
             5175.3,
             5894.0,
             4861.,
             4102.,
             3820.4]
    names = ['L',
             r'H$\eta$',
             r'H$\zeta$',
             'K',
             r'H/H$\epsilon$',
#             r'H$\epsilon$',
             'G',
             r'H$\gamma$',
             'Mg',
             'Na',
             r'H$\beta$',
             r'H$\delta$',
             'L']
    bumps = np.array([0,
                      1,
                      0,
                      1,
                      0,
                      1,
                      0,
                      0,
                      0,
                      0,
                      0,
                      0])*0.03

    yrange = np.diff(ax.get_ylim())
    ypos = ax.get_ylim()[1] + yrange*0.02

    for l, n, b in zip(lines,names,bumps):
        ax.text(l,ypos + yrange*b,n,ha='center',va='bottom', fontsize=9)

    return

def make_plot(output, basedir='.'):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    wave, z, data = get_all_data(basedir)
    widx = np.where(wave < 5500)[0]
    plotz = np.linspace(z.min(),z.max(),300)

    nd = data[:,widx]/np.mean(data[:,widx],axis=1)[:,None]

    # dfunc = spi.interp2d(wave[widx], z, nd, kind='linear')
    # plotd = dfunc(wave[widx],plotz)

    ww, zz = np.meshgrid(wave[widx],z)

    print ww.shape, zz.shape, nd.ravel().shape

    plotd = spi.griddata((ww.ravel(),zz.ravel()), nd.ravel(),
                         (wave[widx][None,:],plotz[:,None]),
                         method='nearest')

    ax.imshow(plotd, cmap=plt.cm.gnuplot2, origin='lower', 
              interpolation='none', vmax=1.8,vmin=0.2, aspect='auto',
              extent=(wave.min(), wave[widx].max(), plotz.min(), plotz.max()))
    ax.axhline(0.4, color='lime', linewidth=2, ls='--')
    
    ax.set_xlabel('$\AA$')
    ax.set_ylabel(r'|$z$| [kpc]')

    add_line_labels(ax)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return nd, plotd

def model_plot(output, infile='/d/monk/eigenbrot/WIYN/14B-0456/anal/mab_plot/zmod.const.norm_hr.ospec.fits'):

    hdu = pyfits.open(infile)[0]
    data = hdu.data
    
    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CD1_1'] + hdu.header['CRVAL1']
    
    widx = np.where((wave >= 3800) & (wave < 5500))[0]
    wave = wave[widx]
    data = data[:,widx]

    data /= np.mean(data,axis=1)[:,None]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$\AA$')
    ax.set_ylabel('|$z$| [kpc]')
    
    ax.imshow(data, cmap=plt.cm.gnuplot2, origin='lower', 
              interpolation='none', vmax=1.8,vmin=0.2, aspect='auto',
              extent=(wave.min(), wave.max(), 0, 0.01*data.shape[0]))
    ax.axhline(0.4, color='lime', linewidth=2, ls='--')

    add_line_labels(ax)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    
    return 
