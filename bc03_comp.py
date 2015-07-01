from glob import glob
import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def bc03_prep():

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]

    for i, model in enumerate(modellist):
        h = pyfits.open(model)[1]
        linwave = np.arange(h.data['WAVE'].min(),h.data['WAVE'].max(),2.1)
        numage = h.data['FLUX'].shape[1]
        out = np.zeros((numage,linwave.size))
        for j in range(numage):
            out[j,:] = np.interp(linwave,h.data['WAVE'][0],
                                 h.data['FLUX'][0,j])
        
        outhdu = pyfits.PrimaryHDU(out)
        outhdu.header.update('Z',fraclist[i])
        outhdu.header.update('CTYPE1','LINEAR')
        outhdu.header.update('CRPIX1',1)
        outhdu.header.update('CRVAL1',linwave.min())
        outhdu.header.update('CDELT1',np.mean(np.diff(linwave)))
        outhdu.header.update('CTYPE2','LINEAR')
        outhdu.header.update('CRPIX2',1)
        outhdu.header.update('CRVAL2',numage)
        outhdu.header.update('CDELT2',1)
        outhdu.writeto('BC03_Z{:04n}.fits'.format(fraclist[i]*1000),clobber=True)

    return

def combine_sbands(output,numaps=10):

    numbands = 7
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)

    results = np.zeros((numaps, fraclist.size, numbands))
    for f, frac in enumerate(fraclist):
        
        fracfile = 'BC03_Z{:04n}.dat'.format(frac*1000)
        data = np.loadtxt(fracfile,usecols=(0,1,5,6),
                          dtype={'names':('aps','bands','index','eqwidth'),
                                 'formats':('S26','S11','f4','f4')})

        aps = np.unique(data['aps'])
        sar = [int(s.split(',')[-1].split(']')[0]) for s in aps]
        sidx = np.argsort(sar)
        print aps[sidx]
        for a, ap in enumerate(aps[sidx]):
            idx = np.where(data['aps'] == ap)
            results[a,f,:] = eat_index(data['index'][idx])

    
    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','aperture')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','HB')
    outhdu.header.update('i1','HdA')
    outhdu.header.update('i2','HgA')
    outhdu.header.update('i3','HdF')
    outhdu.header.update('i4','HgF')
    outhdu.header.update('i5','FeAvg')
    outhdu.header.update('i6','MgFe')
    outhdu.writeto(output,clobber=True)

def bc03_compare(databands, bc03bands, outputprefix):

    pp = PDF('{}.pdf'.format(outputprefix))
    f = open('{}.dat'.format(outputprefix), 'w')
    f.write('#{:>7}{:>10}{:>10}\n'.format('Apnum','Age [Gyr]','Z/Z_sol'))

    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    
    agelist = np.array([  5.01186000e+06,   2.51188000e+07,   1.01518000e+08,
                          2.86119008e+08,   6.40542976e+08,   9.04792000e+08,
                          1.43400000e+09,   2.50000000e+09,   5.00000000e+09,
                          1.00000000e+10])/1e9
    
    bc03data = pyfits.open(bc03bands)[0].data
    numages, numZ, numindex = bc03data.shape

    galdata = np.loadtxt(databands,usecols=(0,1,5,6),
                         dtype={'names':('aps','bands','index','eqwidth'),
                         'formats':('S35','S11','f4','f4')})

    aps = np.unique(galdata['aps'])
    sar = [int(s.split('(')[-1].split(')')[0]) for s in aps]
    sidx = np.argsort(sar)

    outdata = np.zeros((aps.size,2))

    for a, ap in enumerate(aps[sidx]):
        idx = np.where(galdata['aps'] == ap)
        bands = eat_index(galdata['index'][idx])
        compdata = np.tile(bands,(numages,numZ,1))
        chigrid = np.sum((compdata - bc03data)**2/bc03data,axis=2)
        bestage, bestZ = np.unravel_index(np.argmin(chigrid),chigrid.shape)
        outdata[a,0] = agelist[bestage]
        outdata[a,1] = fraclist[bestZ]
        
        f.write('{:8n}{:10.2f}{:10.0e}\n'.format(a+1,agelist[bestage],fraclist[bestZ]))

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        I = ax.imshow(chigrid,cmap=plt.cm.gnuplot,interpolation='nearest',
                      origin='lower')
        ax.set_title('Aperture {:}\nZ = {:5.2e} Z/Z$_{{\odot}}$\nAge = {:5.3f} Gyr'.format(a,outdata[a,1],outdata[a,0]))
        ax.set_yticks(np.arange(agelist.size))
        ax.set_yticklabels(['{:3.3f}'.format(a) for a in agelist])
        ax.set_xlabel('Z/Z$_{\odot}$')
        ax.set_xticks(np.arange(fraclist.size))
        ax.set_xticklabels(fraclist)
        ax.set_ylabel('SSP Age [Gyr]')

        ax.axvline(x=bestZ,color='white')
        ax.axhline(y=bestage,color='white')
        pp.savefig(fig)
        plt.close(fig)

    f.close()
    pp.close()
    return outdata

def eat_index(index):
    #Assumes order is:
    #0 HdA
    #1 HgA
    #2 HdF
    #3 HgF
    #4 Mgb
    #5 Fe5270
    #6 Fe5335
    #7 HB
    
    FeAvg = 0.5*(index[5] + index[6])
    MgFe = np.sqrt(index[4]*(0.72*index[5] + 0.28*index[6]))
    
    return np.array([index[7], index[0], index[1], 
                     index[2], index[3], 
                     FeAvg, MgFe])