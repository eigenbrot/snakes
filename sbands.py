import numpy as np
import pyfits
import time
from glob import glob

def read_fits(fitsfile):
    
    hdu = pyfits.open(fitsfile)[0]
    data = hdu.data

    wave = (np.arange(data.shape[1]) - hdu.header['CRPIX1'] - 1)*hdu.header['CDELT1'] + hdu.header['CRVAL1']

    return wave, data

def compute_Dn4000(wave, spec, err, doerr=True):

    redidx = np.where((wave >= 3850.) & (wave <= 3950.))[0]
    blueidx = np.where((wave >= 4000.) & (wave <= 4100.))[0]
    
    red = np.mean(spec[redidx])
    blue = np.mean(spec[blueidx])

    if doerr:
        red_e = np.sqrt(np.sum(err[redidx]**2))/redidx.size
        blue_e = np.sqrt(np.sum(err[blueidx]**2))/blueidx.size

    Dn4000 = blue/red

    if doerr:
        Dn4000_e = np.sqrt((blue_e/red)**2 + (red_e*blue/red**2)**2)
        return Dn4000, Dn4000_e
    else:
        return Dn4000

def compute_index(wave, spec, err, centlim, bluelim, redlim, doerr=True):
        
    centidx = np.where((wave > centlim[0]) & (wave < centlim[1]))[0]
    redidx = np.where((wave > redlim[0]) & (wave < redlim[1]))[0]
    blueidx = np.where((wave > bluelim[0]) & (wave < bluelim[1]))[0]

    centwave = np.r_[centlim[0], wave[centidx], centlim[1]]
    redwave = np.r_[redlim[0], wave[redidx], redlim[1]]
    bluewave = np.r_[bluelim[0], wave[blueidx], bluelim[1]]

    # cent = np.mean(spec[centidx])
    # blue = np.mean(spec[blueidx])
    # red = np.mean(spec[redidx])
    cent = np.mean(np.interp(centwave, wave, spec))
    blue = np.mean(np.interp(bluewave, wave, spec))
    red = np.mean(np.interp(redwave, wave, spec))
    
    if doerr:
        cent_e = np.sqrt(np.sum(err[centidx]**2))/centidx.size
        blue_e = np.sqrt(np.sum(err[blueidx]**2))/blueidx.size
        red_e = np.sqrt(np.sum(err[redidx]**2))/redidx.size

    redcent = np.mean(redwave)
    bluecent = np.mean(bluewave)
    centcent = np.mean(centlim)
    dlambda = centlim[1] - centlim[0]

    cont = np.interp(centcent, [bluecent, redcent], [blue,red])
    index = (1 - cent/cont)*dlambda
    
    if doerr:
#        cont_e = np.interp(4103.9, [bluecent, redcent], [blue_e, red_e])
#        cont_e = np.sqrt(0.5*(blue_e**2 + red_e**2))
        me = np.sqrt((blue_e/(bluecent - redcent))**2 + (red_e/(bluecent - redcent))**2)
        cont_e = np.sqrt((me*(centcent - redcent))**2 + red_e**2)
        err = dlambda*np.sqrt((cent_e/cont)**2 + (cent*cont_e/cont**2)**2)
        return index, err

    else:
        return index

def main(fitsfile, errfile, output):
    
    #[[cent], [blue], [red]]
    HdAlims = [[4084.5, 4123.3], [4041.65, 4079.75], [4128.55, 4161.05]]
    Mgblims = [[5160.15, 5192.65], [5142.6, 5161.4], [5191.4, 5206.4]]
    Fe5270lims = [[5245.7, 5285.7], [5233.2, 5248.2], [5285.65, 5318.15]]
    Fe5335lims = [[5312.1, 5352.1], [5304.65, 5315.95], [5353.4, 5363.4]]

    wave, data = read_fits(fitsfile)
    _, err = read_fits(errfile)
    Dn4000 = np.zeros(data.shape[0])
    Dn4000_e = np.zeros(data.shape[0])

    HdA = np.zeros(data.shape[0])
    HdA_e = np.zeros(data.shape[0])

    Mgb = np.zeros(data.shape[0])
    Mgb_e = np.zeros(data.shape[0])

    Fe = np.zeros(data.shape[0])
    Fe_e = np.zeros(data.shape[0])

    MgFe = np.zeros(data.shape[0])
    MgFe_e = np.zeros(data.shape[0])

    for i in range(data.shape[0]):
        D, De = compute_Dn4000(wave,data[i,:],err[i,:])
        Dn4000[i] = D
        Dn4000_e[i] = De

        H, He = compute_index(wave, data[i,:], err[i,:],
                              HdAlims[0], HdAlims[1], HdAlims[2])
        HdA[i] = H
        HdA_e[i] = He

        M, Me = compute_index(wave, data[i,:], err[i,:],
                                  Mgblims[0], Mgblims[1], Mgblims[2])
        F2, F2e = compute_index(wave, data[i,:], err[i,:],
                                  Fe5270lims[0], Fe5270lims[1], Fe5270lims[2])
        F3, F3e = compute_index(wave, data[i,:], err[i,:],
                                  Fe5335lims[0], Fe5335lims[1], Fe5335lims[2])
        
        Mgb[i] = M
        Mgb_e[i] = Me

        Fe[i] = 0.5*(F2 + F3)
        Fe_e[i] = np.sqrt((F2e/2.)**2 + (F3e/2.)**2)

        MgFe[i] = np.sqrt(M * (0.72*F2 + 0.28*F3))
        Mterm = (0.72*F2 + 0.28*F3)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        F2term = (0.72*M)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        F3term = (0.28*M)/(2*np.sqrt(M* (0.72*F2 + 0.28*F3)))
        MgFe_e[i] = np.sqrt((Mterm*Me)**2 + (F2term*F2e)**2 + (F3term*F3e)**2)


    fmt = '{:8.3f}'*10 + '\n'
    with open(output,'w') as f:
        f.write('# Generated on {}\n'.format(time.asctime()))
        f.write('# From {} and {}\n'.format(fitsfile, errfile))
        f.write(('{:>8}'*10+'\n').format('#Dn4000',
                                        'Dn4000_e',
                                        'HdA','HdA_e',
                                        'Mgb','Mgb_e',
                                        '<Fe>','<Fe>_e',
                                        'MgFe','MgFe_e'))
        for i in range(HdA.size):
            f.write(fmt.format(Dn4000[i], Dn4000_e[i], 
                               HdA[i], HdA_e[i],
                               Mgb[i], Mgb_e[i],
                               Fe[i], Fe_e[i],
                               MgFe[i], MgFe_e[i]))


    return

def do_all():

    for p in range(6):
        output = 'NGC_891_P{}_bin30.msoz.spy.dat'.format(p+1)
        main('NGC_891_P{}_bin30.msoz.fits'.format(p+1),
             'NGC_891_P{}_bin30.meoz.fits'.format(p+1),output)

    return
             
def run_models(output):

    deftlst = [-1,-3,-9,13,5,3,1]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    print fraclist

    HdAlims = [[4084.5, 4123.3], [4041.65, 4079.75], [4128.55, 4161.05]]
    Mgblims = [[5160.15, 5192.65], [5142.6, 5161.4], [5191.4, 5206.4]]
    Fe5270lims = [[5245.7, 5285.7], [5233.2, 5248.2], [5285.65, 5318.15]]
    Fe5335lims = [[5312.1, 5352.1], [5304.65, 5315.95], [5353.4, 5363.4]]

    results = np.zeros((len(deftlst),fraclist.size,5))

    for f, frac in enumerate(fraclist):
        modelfile = 'BC03_Z{:04n}_tau.fits'.format(frac*1000)
        print modelfile
        
        wave, data = read_fits(modelfile)
        
        for i in range(data.shape[0]):
            D = compute_Dn4000(wave,data[i,:],None,doerr=False)            
            H = compute_index(wave, data[i,:], None,
                              HdAlims[0], HdAlims[1], HdAlims[2],
                              doerr=False)
            M = compute_index(wave, data[i,:], None,
                                  Mgblims[0], Mgblims[1], Mgblims[2],
                              doerr=False)
            F2 = compute_index(wave, data[i,:],None,
                                    Fe5270lims[0], Fe5270lims[1], Fe5270lims[2],
                               doerr=False)
            F3 = compute_index(wave, data[i,:],None,
                                    Fe5335lims[0], Fe5335lims[1], Fe5335lims[2],
                               doerr=False)
            Fe = 0.5*(F2 + F3)
            MgFe = np.sqrt(M * (0.72*F2 + 0.28*F3))
            
            results[i,f,:] = np.array([D,H,M,Fe,MgFe])

    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','aperture')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','Dn4000')
    outhdu.header.update('i1','HdA')
    outhdu.header.update('i2','Mgb')
    outhdu.header.update('i3','<Fe>')
    outhdu.header.update('i4','MgFe')
    outhdu.writeto(output,clobber=True)
    
    return

def run_models_multires(output,group=1):
    dirlistd = {1: ['29377','24192','19958','18738','17947'],
                2: ['39180','34356','25751','24995','24471'],
                3: ['47014','42806','31304','30267','29751']}
    deftlst = [-1,-3,-9,13,5,3,1]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])
    fraclist = np.sort(fraclist)
    print fraclist

    HdAlims = [[4084.5, 4123.3], [4041.65, 4079.75], [4128.55, 4161.05]]
    Mgblims = [[5160.15, 5192.65], [5142.6, 5161.4], [5191.4, 5206.4]]
    Fe5270lims = [[5245.7, 5285.7], [5233.2, 5248.2], [5285.65, 5318.15]]
    Fe5335lims = [[5312.1, 5352.1], [5304.65, 5315.95], [5353.4, 5363.4]]

    results = np.zeros((len(deftlst),fraclist.size,5))

    for f, frac in enumerate(fraclist):
        reslist = []
        for resdir in dirlistd[group]:
            modelfile = '{:}/BC03_Z{:04n}_tau.fits'.format(resdir,frac*1000)
            print modelfile
        
            wave, data = read_fits(modelfile)
            aplist = []
            for i in range(data.shape[0]):
                D = compute_Dn4000(wave,data[i,:],None,doerr=False)            
                H = compute_index(wave, data[i,:], None,
                                  HdAlims[0], HdAlims[1], HdAlims[2],
                                  doerr=False)
                M = compute_index(wave, data[i,:], None,
                                  Mgblims[0], Mgblims[1], Mgblims[2],
                                  doerr=False)
                F2 = compute_index(wave, data[i,:],None,
                                   Fe5270lims[0], Fe5270lims[1], Fe5270lims[2],
                                   doerr=False)
                F3 = compute_index(wave, data[i,:],None,
                                   Fe5335lims[0], Fe5335lims[1], Fe5335lims[2],
                                   doerr=False)
                
                aplist.append([D,H,M,F2,F3])
                
            reslist.append(aplist)

        #allres = [tausf (ap),index,resolution]
        allres = np.dstack(reslist)
        print allres.shape
        tmp = np.zeros(allres.shape[:-1])
        print tmp.shape
        
        #Dn4000
        tmp[:,0] = allres[:,0,0]
        #HdA
        tmp[:,1] = allres[:,1,1]
        #Mgb
        tmp[:,2] = allres[:,2,2]
        #<Fe>
        tmp[:,3] = 0.5*(allres[:,3,3] + allres[:,4,4])
        #[MgFe]
        tmp[:,4] = np.sqrt(allres[:,2,2] * (0.72*allres[:,3,3] + 0.28*allres[:,4,4]))
        
        results[:,f,:] = tmp

    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','aperture')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','Dn4000')
    outhdu.header.update('i1','HdA')
    outhdu.header.update('i2','Mgb')
    outhdu.header.update('i3','<Fe>')
    outhdu.header.update('i4','MgFe')
    outhdu.writeto(output,clobber=True)
    
    return

def run_SSPs(output):

    HdAlims = [[4084.5, 4123.3], [4041.65, 4079.75], [4128.55, 4161.05]]
    Mgblims = [[5160.15, 5192.65], [5142.6, 5161.4], [5191.4, 5206.4]]
    Fe5270lims = [[5245.7, 5285.7], [5233.2, 5248.2], [5285.65, 5318.15]]
    Fe5335lims = [[5312.1, 5352.1], [5304.65, 5315.95], [5353.4, 5363.4]]

    SSPs = pyfits.open('/d/monk/eigenbrot/WIYN/14B-0456/anal/models/allZ2_vardisp/allz2_vardisp_batch_interp.fits')[1].data[0]
    frac = np.sort(np.unique(SSPs['Z']))
    ages = np.sort(np.unique(SSPs['AGE']))
    numZ = frac.size
    numage = ages.size

    wave = SSPs['WAVE']
    
    results = np.zeros((numage,numZ,5))
    
    for zz, Z in enumerate(frac):
        for aa, A in enumerate(ages):
            idx = np.where((SSPs['AGE'] == A) & (SSPs['Z'] == Z))[0][0]
            data = SSPs['FLUX'][idx,:,2]

            D = compute_Dn4000(wave,data,None,doerr=False)            
            H = compute_index(wave, data, None,
                              HdAlims[0], HdAlims[1], HdAlims[2],
                              doerr=False)
            M = compute_index(wave, data, None,
                              Mgblims[0], Mgblims[1], Mgblims[2],
                              doerr=False)
            F2 = compute_index(wave, data,None,
                               Fe5270lims[0], Fe5270lims[1], Fe5270lims[2],
                               doerr=False)
            F3 = compute_index(wave, data,None,
                               Fe5335lims[0], Fe5335lims[1], Fe5335lims[2],
                               doerr=False)
            
            res = [D, H, M,
                   0.5*(F2 + F3),
                   np.sqrt(M * (0.72*F2 + 0.28*F3))]

            results[aa,zz,:] = res

    
    outhdu = pyfits.PrimaryHDU(results)
    outhdu.header.update('d0','Age')
    outhdu.header.update('d1','Z')
    outhdu.header.update('d2','index')
    outhdu.header.update('i0','Dn4000')
    outhdu.header.update('i1','HdA')
    outhdu.header.update('i2','Mgb')
    outhdu.header.update('i3','<Fe>')
    outhdu.header.update('i4','MgFe')
    outhdu.writeto(output,clobber=True)
    
    return
