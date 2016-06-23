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
             3970.18,
             4304.4,
             4341.,
             5175.3,
             5894.0,
             4861.,
             4102.]
    names = ['L',
             r'H$\eta$',
             r'H$\zeta$',
             'K',
             r'H',
             r'H$\epsilon$',
             'G',
             r'H$\gamma$',
             'Mg',
             'Na',
             r'H$\beta$',
             r'H$\delta$']
    bumps = np.array([1,
                      0,
                      0,
                      1,
                      1,
                      0,
                      1,
                      0,
                      1,
                      1,
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
    widx = np.where(wave < 4500)[0]
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

def model_plot(output, infile='/d/monk/eigenbrot/WIYN/14B-0456/anal/mab_plot/zmod.const.norm_hr.ospec.prep.fits'):

    hdu = pyfits.open(infile)[0]
    data = hdu.data
    
    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CD1_1'] + hdu.header['CRVAL1']
    
    widx = np.where((wave >= 3800) & (wave < 4500))[0]
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

def model_spec_plot(output, infile='/d/monk/eigenbrot/WIYN/14B-0456/anal/mab_plot/zmod.const.norm_hr.ospec.fits'):
    
    hdu = pyfits.open(infile)[0]
    data = hdu.data
    
    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CD1_1'] + hdu.header['CRVAL1']
    
    widx = np.where((wave >= 3800) & (wave < 4500))[0]
    wave = wave[widx]
    data = data[:,widx]

    data /= np.mean(data,axis=1)[:,None]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Wavelength [$\AA$]')
    ax.set_ylabel('Normalized Flux + offset')

    ax.plot(wave,data[16,:], 'k') #0.16 kpc
    ax.plot(wave,data[32,:] + 1, 'k') #0.32 kpc
    ax.plot(wave,data[60,:] + 2, 'k') #0.6 kpc
    
    ax.text(3850,1.1,'0.16 kpc',ha='center',va='bottom',fontsize=10)
    ax.text(3850,1.9,'0.32 kpc',ha='center',va='bottom',fontsize=10)
    ax.text(3850,2.8,'0.6 kpc',ha='center',va='bottom',fontsize=10)

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

    ypos = 3.8
    for l, n in zip(lines,names):
        ax.text(l,ypos,n,ha='center',va='top', fontsize=9)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return 

def plot_stack(output, infile='/d/monk/eigenbrot/WIYN/14B-0456/anal/mab_plot/zmod.const.norm_hr.ospec.prep.fits'):

    #Img
    hdu = pyfits.open(infile)[0]
    data = hdu.data
    
    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CDELT1'] + hdu.header['CRVAL1']
    
    widx = np.where((wave >= 3800) & (wave < 4500))[0]
    wave = wave[widx]
    data = data[:,widx]

    data /= np.mean(data,axis=1)[:,None]

    fig = plt.figure()
    ax = fig.add_subplot(211)
    # ax.set_xlabel('$\AA$')
    ax.set_xticklabels([])
    ax.set_ylabel('|$z$| [kpc]')

    ax.imshow(data, cmap=plt.cm.gnuplot2, origin='lower', 
              interpolation='none', vmax=1.8,vmin=0.2, aspect='auto',
              extent=(wave.min(), wave.max(), 0, 0.01*data.shape[0]))
    ax.axhline(0.4, color='lime', linewidth=2, ls='--')

    add_line_labels(ax)

    #Spec
    ax1 = fig.add_subplot(212)
    ax1.set_xlabel('Wavelength [$\AA$]')
    ax1.set_ylabel('Normalized Flux + offset')

    ax1.plot(wave,data[16,:], 'k') #0.16 kpc
    ax1.plot(wave,data[32,:] + 1, 'k') #0.32 kpc
    ax1.plot(wave,data[60,:] + 2, 'k') #0.6 kpc
    
    ax1.text(3870,1.1,'0.16 kpc',ha='center',va='bottom',fontsize=12)
    ax1.text(3870,1.95,'0.32 kpc',ha='center',va='bottom',fontsize=12)
    ax1.text(3870,2.9,'0.6 kpc',ha='center',va='bottom',fontsize=12)

    lines = [3820.4,
             3835.4,
             3889.0,
             3933.7,
             3968.5,
             3970.18,
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
             r'H',
             r'H$\epsilon$',
             'G',
             r'H$\gamma$',
             'Mg',
             'Na',
             r'H$\beta$',
             r'H$\delta$',
             'L']

    ypos = 3.8
    for l, n in zip(lines,names):
        if n[0:2] != r'H$':
            b = 0.2
        else:
            b = 0
        ax1.text(l,ypos+b,n,ha='center',va='top', fontsize=11)

    ax1.set_ylim(0,4.6)
    ax.set_xlim(*ax1.get_xlim())
    fig.subplots_adjust(hspace=0.0001)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return

def data_stack(output,basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/mab_plot'):

    #Img
    fig = plt.figure()
    ax = fig.add_subplot(211)

    wave, z, data = get_all_data(basedir)
    widx = np.where(wave < 4500)[0]
    plotz = np.linspace(z.min(),z.max(),300)

    nd = data[:,widx]/np.mean(data[:,widx],axis=1)[:,None]

    # dfunc = spi.interp2d(wave[widx], z, nd, kind='linear')
    # plotd = dfunc(wave[widx],plotz)

    ww, zz = np.meshgrid(wave[widx],z)

    plotd = spi.griddata((ww.ravel(),zz.ravel()), nd.ravel(),
                         (wave[widx][None,:],plotz[:,None]),
                         method='nearest')

    ax.imshow(plotd, cmap=plt.cm.gnuplot2, origin='lower', 
              interpolation='none', vmax=1.8,vmin=0.2, aspect='auto',
              extent=(wave.min(), wave[widx].max(), plotz.min(), plotz.max()))
    ax.axhline(0.4, color='lime', linewidth=2, ls='--')
    
    # ax.set_xlabel('$\AA$')
    ax.set_xticklabels([])
    ax.set_ylabel(r'|$z$| [kpc]')

    add_line_labels(ax)

    #Spec
    ax1 = fig.add_subplot(212)
    ax1.set_xlabel('Wavelength [$\AA$]')
    ax1.set_ylabel('Normalized Flux + Offset')
    ax1.set_ylim(0,5.4)

    id1 = np.argmin(np.abs(z - 0.3))
    id2 = np.argmin(np.abs(z - 0.6))
    id3 = np.argmin(np.abs(z - 1.5))

    ax1.plot(wave[widx], nd[id1,:], 'k')
    ax1.plot(wave[widx], nd[id2,:] + 1.5, 'k')
    ax1.plot(wave[widx], nd[id3,:] + 2.7, 'k')
    
    ax1.text(3870,1,'0.3 kpc',ha='center',va='bottom',fontsize=12)
    ax1.text(3870,2.4,'0.6 kpc',ha='center',va='bottom',fontsize=12)
    ax1.text(3870,3.8,'1.5 kpc',ha='center',va='bottom',fontsize=12)

    ax1.set_xlim(ax.get_xlim())

    lines = [3820.4,
             3835.4,
             3889.0,
             3933.7,
             3968.5,
             3970.18,
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
             r'H',
             r'H$\epsilon$',
             'G',
             r'H$\gamma$',
             'Mg',
             'Na',
             r'H$\beta$',
             r'H$\delta$',
             'L']

    ypos = 4.8
    for l, n in zip(lines,names):
        if n[0:2] != r'H$':
            b = 0.3
        else:
            b = 0
        ax1.text(l,ypos+b,n,ha='center',va='top', fontsize=11)
        
    fig.subplots_adjust(hspace=0.001)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return
