import numpy as np
import scipy.ndimage as spnd
import pyfits
import matplotlib.pyplot as plt
plt.ioff()

def edit_errors(inputfile, velocityfile,
                output = None,
                centlist = [3969., 3934., 4102., 4341., 4861.],
                width = 20, #AA
                amp = 10):
                
    hdu = pyfits.open(inputfile)[0]
    data = hdu.data
    aps, V = np.loadtxt(velocityfile, unpack=True)
    
    wave = (np.arange(data.shape[1]) - (hdu.header['CRPIX1'] - 1))*hdu.header['CDELT1'] + hdu.header['CRVAL1']

    for i in range(data.shape[0]):
        mod = np.zeros(data.shape[1])
        for c in centlist:
            tc = np.array(c)*(V[i]/3e5 + 1)
            if tc.size == 1:
                idx = np.argmin(np.abs(wave - c))
            else:
                idx = np.where((wave >= tc[0]) & (wave <= tc[1]))[0]
            mod[idx] = 1.

        # ax = plt.figure().add_subplot(111)
        # ax.plot(wave,mod+1)
        mod = spnd.filters.gaussian_filter1d(mod,width/hdu.header['CDELT1'])
        mod *= (amp-1)/mod.max()
        mod += 1
        # ax.plot(wave,mod)
        # ax.figure.show()
        # raw_input()

        data[i,:] /= mod

    pyfits.PrimaryHDU(data,hdu.header).writeto(output, clobber=True)

    return output

def batch_edit(width=20, amp=100, centlist=None, metal=False):
    
    if centlist is None:
        if metal:
            centlist = [[6000.,6300.],4780.,4930.]
            suff = 'Z'
        else:
            centlist = [4304.4, 4341., 4102., 3933.7, 3968.5, 3970.18]
            suff = 'A'

    for p in range(6):
        inputfile = 'NGC_891_P{}_bin30.meo.fits'.format(p+1)
        velocity = 'NGC_891_P{}_bin30_velocities.dat'.format(p+1)
        output = 'NGC_891_P{}_bin30.meo{}.fits'.format(p+1,suff)
        print '{} -> {}'.format(inputfile, output)
        
        edit_errors(inputfile, velocity, output=output,
                    width=width, amp=amp, centlist=centlist)

    return

def apply_rms_spec(inputfile, velocityfile, output,
                   amp=1.):

    hdu = pyfits.open(inputfile)[0]
    data = hdu.data
    aps, V = np.loadtxt(velocityfile, unpack=True)
    
    wave = (np.arange(data.shape[1]) - (hdu.header['CRPIX1'] - 1))*hdu.header['CDELT1'] + hdu.header['CRVAL1']

    rms_wave, rms, rms_Z, rms_A = compute_rms_spec()

    for i in range(data.shape[0]):
        rms_V = np.interp(wave,rms_wave*(V[i]/3e5 + 1),rms)
        rms_V /= np.mean(rms_V)
        weights = 1./(rms_V*amp)
        data[i,:] *= weights
    
    pyfits.PrimaryHDU(data,hdu.header).writeto(output, clobber=True)

    return
        
def batch_apply(suff='R'):

    for p in range(6):
        inputfile = 'NGC_891_P{}_bin30.meo.fits'.format(p+1)
        velocity = 'NGC_891_P{}_bin30_velocities.dat'.format(p+1)
        output = 'NGC_891_P{}_bin30.meo{}.fits'.format(p+1,suff)
        print '{} -> {}'.format(inputfile, output)
        
        apply_rms_spec(inputfile, velocity, output, amp=1.)

    return

def compute_rms_spec(swindow=100):

    mlist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_02Z_vardisp.fits',
             '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_04Z_vardisp.fits',
             '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_1Z_vardisp.fits',
             '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_25Z_vardisp.fits']
    Zlist = [0.2,0.4,1,2.5]

    #Get dims
    tmp = pyfits.open(mlist[0])[1].data[0]
    nage = tmp['FLUX'].shape[0]
    nmetal = len(mlist)
    wave = tmp['WAVE']
    idx = np.where((wave >= 3800) & (wave <= 6800))[0]
    wave = wave[idx]
    age = tmp['AGE'][:,0]/1e9
    nwave = wave.size

    bigD = np.zeros((nmetal, nage, nwave))

    for i, mod in enumerate(mlist):
        data = pyfits.open(mod)[1].data[0]
        for j in range(nage):
            cont = np.convolve(data['FLUX'][j,idx,0],np.ones(swindow)/swindow,'same')
            bigD[i,j,:] = data['FLUX'][j,idx,0]/cont

    # bigAvg = np.mean(bigD,axis=(0,1))
    # smoothedAvg = np.convolve(bigAvg,np.ones(swindow)/swindow,'same')
    
    # contD = bigD / smoothedAvg
    # contAvg = np.mean(contD,axis=(0,1))

    # print bigD.shape, bigAvg.shape, contD.shape

    # sax = plt.figure().add_subplot(111)
    # sax.plot(wave, bigAvg)
    # sax.plot(wave, smoothedAvg)
    # sax.plot(wave, contAvg)
    # sax.figure.show()

    #Metal axis
    mfig = plt.figure()
    ax = mfig.add_subplot(111)
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('RMS Across age')
    rms_A = np.zeros(nwave)
    ageAvg = np.mean(bigD,axis=1)
    print bigD[0,:,:].shape, ageAvg[0].shape, (bigD[0,:,:] - ageAvg[0]).shape
    for z in range(nmetal):
        rms = np.sqrt(np.mean((bigD[z,:,:] - ageAvg[z])**2,axis=0))
        rms = np.convolve(rms,np.ones(20)/20.,'same')
        ax.plot(wave, rms,label='{}$Z_{{\odot}}$'.format(Zlist[z]))
        # ax.plot(wave, np.mean(bigD[z,:,:],axis=0),':')
        # ax.plot(wave, ageAvg[z],'--')
        rms_A += rms**2

    rms_A = np.sqrt(rms_A)/nmetal
#    ax.plot(wave,rms_Z,label='total')
    ax.legend(loc=0)
    add_line_labels(ax)
    mfig.savefig('A_rms2.png',dpi=200)

    #Age axis
    afig = plt.figure()
    ax = afig.add_subplot(111)
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('RMS across Z')
    rms_Z = np.zeros(nwave)
    ZAvg = np.mean(bigD,axis=0)
    for a in range(nage):
        rms = np.sqrt(np.mean((bigD[:,a,:] - ZAvg[a])**2,axis=0))
        rms = np.convolve(rms,np.ones(20)/20.,'same')
        ax.plot(wave, rms, label='{:7.4f} Gyr'.format(age[a]))
        rms_Z += rms**2

    rms_Z = np.sqrt(rms_Z)/nage
#    ax.plot(wave, rms_A,label='total')
    ax.legend(loc=0)
    add_line_labels(ax)
    afig.savefig('Z_rms2.png',dpi=200)

    #All
    bigfig = plt.figure()
    ax = bigfig.add_subplot(111)
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('RMS')
    bigAvg = np.mean(bigD,axis=(0,1))
    rms = np.sqrt(np.mean((bigD - bigAvg)**2,axis=(0,1)))
    rms = np.convolve(rms,np.ones(20)/20.,'same')
    ax.plot(wave, rms, label='total')
    ax.plot(wave, rms_Z/np.mean(rms_Z)*np.mean(rms_A), label=r'RMS$_Z$')
    ax.plot(wave, rms_A, label=r'RMS$_{\rm age}$')
    ax.legend(loc=0)
    add_line_labels(ax)    
#    bigfig.show()
    bigfig.savefig('big_RMS3.png',dpi=200)

    return wave, rms, rms_Z, rms_A

def add_line_labels(ax):

    em = np.array([6563.8,  6716.0, 6583.41, 6548.04, 4959., 5006.8])
    emnam = [r'H$\alpha$', 'S2', 'NII', 'NII', '[OIII]', '[OIII]']
        
    ab = np.array([3820.4, 3835.4,      3889.0,     3933.7, 3968.5, 3970.18,         4304.4,   4341.,       5175.3, 5894.0, 4861.,  4102.])
    absnam = ['L',   r'H$\eta$', r'H$\zeta$', 'K',   'H'   , r'H$\epsilon$',    'G',     r'H$\gamma$',  'Mg',   'Na',   r'H$\beta$',   r'H$\delta$']

    emposfac = [2,1,3,1,2,3]
    absposfac = [1,2,3,1,3,2,1,2,1,1,1,1]
    ymin, ymax = ax.get_ylim()
    
    for e, en, ef in zip(em, emnam, emposfac):
        top_pos = ymax + 0.01*ef*2*(ymax - ymin)
        ax.axvline(x = e, color='b', alpha=0.2)
        ax.text(e, top_pos, en, ha='center', fontsize=8)

    for a, an, af in zip(ab, absnam, absposfac):
        top_pos = ymax + 0.01*af*2*(ymax - ymin)
        ax.axvline(x = a, color='r', alpha=0.2)
        ax.text(a, top_pos, an, ha='center', fontsize=8)

    return
