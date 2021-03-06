import numpy as np
import os
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
                   power=1,
                   ml = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_02Z_vardisp.fits',
                         '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_04Z_vardisp.fits',
                         '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_1Z_vardisp.fits',
                         '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_25Z_vardisp.fits']):

    hdu = pyfits.open(inputfile)[0]
    data = hdu.data
    aps, V = np.loadtxt(velocityfile, unpack=True)
    
    wave = (np.arange(data.shape[1]) - (hdu.header['CRPIX1'] - 1))*hdu.header['CDELT1'] + hdu.header['CRVAL1']

    rms_wave, rms, rms_Z, rms_A = compute_rms_spec(mlist=ml)

    for i in range(data.shape[0]):
        rms_V = np.interp(wave,rms_wave*(V[i]/3e5 + 1),rms)
        rms_V /= np.mean(rms_V)
        weights = 1./(rms_V**power)
        data[i,:] *= weights
    
    pyfits.PrimaryHDU(data,hdu.header).writeto(output, clobber=True)

    return
        
def batch_apply(suff='R',power=1,
                ml=['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_02Z_vardisp.fits',
                    '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_04Z_vardisp.fits',
                    '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_1Z_vardisp.fits',
                    '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_25Z_vardisp.fits']):

    for p in range(6):
        inputfile = 'NGC_891_P{}_bin30.meo.fits'.format(p+1)
        velocity = 'NGC_891_P{}_bin30_velocities.dat'.format(p+1)
        output = 'NGC_891_P{}_bin30.meo{}.fits'.format(p+1,suff)
        print '{} -> {}'.format(inputfile, output)
        
        apply_rms_spec(inputfile, velocity, output, power=power, ml=ml)

    return

def convert_to_IRAF(mlist=['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_02Z_vardisp.fits',
                           '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_04Z_vardisp.fits',
                           '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_1Z_vardisp.fits',
                           '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_25Z_vardisp.fits']):
    

    for mod in mlist:
        d = pyfits.open(mod)[1].data[0]
        wave = d['WAVE']
        lin_wave = np.linspace(np.min(wave),np.max(wave),wave.size)
        print np.std(np.diff(lin_wave))
        CDELT = np.mean(np.diff(lin_wave))
        tmp = np.zeros(d['FLUX'].shape[0:2])
        for i in range(tmp.shape[0]):
            tmp[i,:] = np.interp(lin_wave, wave, d['FLUX'][i,:,0])

        HDU = pyfits.PrimaryHDU(tmp)
        HDU.header.update('CRVAL1',lin_wave.min())
        HDU.header.update('CRPIX1',1)
        HDU.header.update('CDELT1',CDELT)

        outname = '{}.iraf.fits'.format(os.path.splitext(os.path.basename(mod))[0])
        print outname
        HDU.writeto(outname,clobber=True)

    return

def compute_rms_spec(mlist=['/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_02Z_vardisp.fits',
                            '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_04Z_vardisp.fits',
                            '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_1Z_vardisp.fits',
                            '/d/monk/eigenbrot/WIYN/14B-0456/anal/DFK/models/DFK_25Z_vardisp.fits'],
                     swindow=100,order=5,paperplot=False):

    if len(mlist) == 6:
        Zlist = [0.005,0.02,0.2,0.4,1,2.5]
        clist = ['orange','purple','b','g','r','cyan']
    elif len(mlist) == 5:
        Zlist = [0.02,0.2,0.4,1,2.5]        
        clist = ['purple','b','g','r','cyan']
    else:
        Zlist = [0.2,0.4,1,2.5]
        clist = ['b','g','r','cyan']

    #Get dims, IRAF style
    tmp = pyfits.open(mlist[0])[0]
    nage, nwave = tmp.data.shape
    nmetal = len(mlist)
    wave = (np.arange(nwave) - tmp.header['CRPIX1'] - 1)*tmp.header['CDELT1'] + tmp.header['CRVAL1']
    age = np.array([2.75846757e-03,   2.16358081e-01,   3.15126872e+00, 9.68775558e+00]) #MAGIC NUMBERS!
    
    # #Get dims
    # tmp = pyfits.open(mlist[0])[1].data[0]
    # nage = tmp['FLUX'].shape[0]
    # nmetal = len(mlist)
    # wave = tmp['WAVE']
    # idx = np.where((wave >= 3800) & (wave <= 6800))[0]
    # wave = wave[idx]
    # age = tmp['AGE'][:,0]/1e9
    # nwave = wave.size

    bigD = np.zeros((nmetal, nage, nwave))

    # for i, mod in enumerate(mlist):
    #     data = pyfits.open(mod)[1].data[0]
    #     for j in range(nage):
    #         cont = np.convolve(data['FLUX'][j,:,0],np.ones(swindow)/swindow,'same')
    #         # contf = np.poly1d(np.polyfit(wave,data['FLUX'][j,idx,0],order))
    #         # cont = contf(wave)
    #         bigD[i,j,:] = data['FLUX'][j,idx,0]/cont[idx]
    #         cax = plt.figure().add_subplot(111)
    #         cax.plot(wave,data['FLUX'][j,idx,0])
    #         cax.plot(wave,cont[idx])
    #         #cax.plot(wave, bigD[i,j,:])
    #         cax.figure.show()
    #         print '{}\n\t{}'.format(mod,j)
    #         raw_input('')
    #         plt.close(cax.figure)

    for i, mod in enumerate(mlist):
        data = pyfits.open(mod)[0].data
        bigD[i,:,:] = data

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
    oldAvg = np.mean(bigD[:,2:,:],axis=1)
    youngAvg = np.mean(bigD[:,:2,:],axis=1)
    print bigD[0,:,:].shape, ageAvg[0].shape, (bigD[0,:,:] - ageAvg[0]).shape
    # saverms = np.sqrt(np.mean((bigD[3,:,:] - ageAvg[3])**2,axis=0))
    for z in range(nmetal):
        if z == 0:
            rms = np.sqrt(np.mean((bigD[z,2:,:] - oldAvg[z])**2,axis=0))
        elif z >=1 and z < 4:
            rms = np.sqrt(np.mean((bigD[z,:,:] - ageAvg[z])**2,axis=0))            
        else:
            rms = np.sqrt(np.mean((bigD[z,:2,:] - youngAvg[z])**2,axis=0))
        # rms = np.sqrt(np.mean((bigD[z,:,:] - ageAvg[z])**2,axis=0))            
        # rms /= saverms
        rms = np.convolve(rms,np.ones(20)/20.,'same')
        ax.plot(wave, rms,label='{}$Z_{{\odot}}$'.format(Zlist[z]),color=clist[z])
        # ax.plot(wave, np.mean(bigD[z,:,:],axis=0),':')
        # ax.plot(wave, ageAvg[z],'--')
        rms_A += rms**2

    rms_A = np.sqrt(rms_A)/nmetal
#    ax.plot(wave,rms_Z,label='total')
    ax.set_xlim(3400,7000)
    ax.legend(loc=0)
    add_line_labels(ax)
    if not paperplot:
        mfig.savefig('A_rms_IRAF_prior2.png',dpi=200)

    #Age axis
    afig = plt.figure()
    ax = afig.add_subplot(111)
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('RMS across Z')
    rms_Z = np.zeros(nwave)
    ZAvg = np.mean(bigD,axis=0)
    highZAvg = np.mean(bigD[1:,:,:],axis=0)
    lowZAvg = np.mean(bigD[:4,:,:],axis=0)
    # saverms = np.sqrt(np.mean((bigD[:,2,:] - ZAvg[2])**2,axis=0))
    for a in range(nage):
        if a < 2:
            rms = np.sqrt(np.mean((bigD[1:,a,:] - highZAvg[a])**2,axis=0))
        else:
            rms = np.sqrt(np.mean((bigD[:4,a,:] - lowZAvg[a])**2,axis=0))
        # rms = np.sqrt(np.mean((bigD[:,a,:] - ZAvg[a])**2,axis=0))
        # rms /= saverms
        rms = np.convolve(rms,np.ones(20)/20.,'same')
        ax.plot(wave, rms, label='{:7.4f} Gyr'.format(age[a]))
        rms_Z += rms**2

    rms_Z = np.sqrt(rms_Z)/nage
#    ax.plot(wave, rms_A,label='total')
    ax.set_xlim(3500,7000)
    ax.legend(loc=0)
    add_line_labels(ax)
    if not paperplot:
        afig.savefig('Z_rms_IRAF_prior2.png',dpi=200)
    
    #All
    bigfig = plt.figure()
    ax = bigfig.add_subplot(111)
    ax.set_xlabel(r'$\AA$')
    ax.set_ylabel('RMS')
    bigAvg = np.mean(bigD,axis=(0,1))
    rms = np.sqrt(np.mean((bigD - bigAvg)**2,axis=(0,1)))
    rms = np.convolve(rms,np.ones(20)/20.,'same')
    rms = np.sqrt(rms_Z**2 + rms_A**2)
    ax.plot(wave, rms, label='total', color='#1b9e77',zorder=100)
    ax.plot(wave, rms_Z/np.mean(rms_Z)*np.mean(rms_A), label=r'RMS$_Z$',color='#d95f02')
    ax.plot(wave, rms_A, label=r'RMS$_{\rm age}$',color='#7570b3')
    ax.set_xlim(3500,7000)
    ax.legend(loc=0)
    add_line_labels(ax)    
    # bigfig.show()
    if not paperplot:
        bigfig.savefig('big_RMS_IRAF_prior2.png',dpi=200)

    if paperplot:
        return bigfig

    return wave, rms, rms_Z, rms_A

def add_line_labels(ax):

    em = np.array([6563.8,  6716.0, 6583.41, 6548.04, 4959., 5006.8])
    emnam = [r'H$\alpha$', 'S2', 'NII', 'NII', '[OIII]', '[OIII]']
        
    ab = np.array([3820.4, 3835.4,      3889.0,     3933.7, 3968.5, 3970.18,
                   4304.4,   4341.,         5175.3, 5894.0,     4861.,        4102.,
                   5266, 5332])
    absnam = ['L',   r'H$\eta$', r'H$\zeta$', 'K',   'H'   , r'H$\epsilon$',
                     'G',     r'H$\gamma$',  'Mg',   'Na',   r'H$\beta$',   r'H$\delta$',
                    'Fe','Fe']

    emposfac = [2,1,3,1,2,3,1,1]
    absposfac = [1,2,3,1,3,2,1,2,1,1,1,1,1,1]
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
