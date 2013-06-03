import numpy as np
import matplotlib.pyplot as plt
from MANGA_bench import Noodle
import ConfigParser
import pyfits
import ADEUtils as ADE
from glob import glob
import pickle
import os
from matplotlib.backends.backend_pdf import PdfPages as PDF

def fratio(searchstr,EEfigs,apfigs,bigfig,EE=0.50):
    '''pronounced frat-eeo, not f ratio'''
    print int(EE*100)

    pp = PDF(EEfigs)
    pp2 = PDF(apfigs)
    pp3 = PDF(bigfig)
    
    inifiles = glob(searchstr)
    datadict = get_sizes(inifiles,pp,EE=EE)
    pp.close()
    ratiodict = {}
    aps = np.array([])
    ratios = np.array([])

    for apsize in datadict.keys():

        x = np.array(datadict[apsize].keys(),dtype=float)
        r = np.array(datadict[apsize].values())

        fit_coef = ADE.polyclip(x,r,1,niter=100).c
        
        slope = fit_coef[0]
        N = (2.*slope)**-1
        fit = np.poly1d(fit_coef)

        ratiodict[apsize] = {'f/#':N,
                             'fl':float(apsize)/N}

        aps = np.append(aps,float(apsize))
        ratios = np.append(ratios,N)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,r,marker='s',linestyle='')
        fitx = np.linspace(x.min(),x.max(),50)
        ax.plot(fitx,fit(fitx),'k:',label='linear fit')
        ax.set_xlabel('back distance - C [mm]')
        ax.set_ylabel('beam radius [mm]')
        ax.legend(loc=0)
        ax.set_title('Aperture: {:n} mm\nN: {:3.2f}'.format(float(apsize),N))
        
        pp2.savefig(fig)

    pp2.close()
    
    ratiofit = ADE.polyclip(aps**-1,ratios,1,niter=50,clip_high=1,clip_low=1)
#    ratiofit = np.poly1d(np.polyfit(aps**-1,ratios,1))
    ratioslope = ratiofit.c[0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(aps**-1,ratios,marker='s',linestyle='')
    fitaps = np.linspace(aps.min(),aps.max(),50)
    ax.plot(fitaps**-1,ratiofit(fitaps**-1),':',label='linear fit')
    ax.set_xlabel('1/D [$mm^{-1}$]')
    ax.set_ylabel('f-ratio')
#    ax.set_ylim(0, ratiofit(fitaps.max())*1.1)
    ax.legend(loc=0)
    ax.set_title('slope = {:n} mm'.format(ratioslope))
    
    pp3.savefig(fig)
    pp3.close()
    
    return ratiodict

def EE_test(searchstr,EElist):

    inifiles = glob(searchstr)

    FLlist = np.array([])

    for EE in EElist:

        print '************ EE'+str(int(EE*100))+' ***************'

        datadict = get_sizes(inifiles,EErad=EE)
        ratiodict = {}
        aps = np.array([])
        ratios = np.array([])
        
        for apsize in datadict.keys():
            
            x = np.array(datadict[apsize].keys(),dtype=float)
            r = np.array(datadict[apsize].values())
            
            fit_coef = np.polyfit(x,r,1)
            
            slope = fit_coef[0]
            fit = np.poly1d(fit_coef)
            
            ratiodict[apsize] = {'f/#':(2.*slope)**-1,
                                 'fl':float(apsize)/(2.*slope)}
            
            aps = np.append(aps,float(apsize))
            ratios = np.append(ratios,(2.*slope)**-1)
            
        ratiocoef = np.polyfit(aps**-1,ratios,1)
        ratioslope = ratiocoef[0]
        ratiofit = np.poly1d(ratiocoef)
        
        FLlist = np.append(FLlist,ratioslope)
    

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.array(EElist),FLlist)
    fig.show()

    return np.array(EElist),FLlist

def ratio_boot(ratiodict,numtrys):

    D = np.array([])
    N = np.array([])

    for ap in ratiodict.keys():
        D = np.append(D,float(ap))
        N = np.append(N,ratiodict[ap]['f/#'])

    slopes = np.array([])

    for i in range(numtrys):
        
        sampleidx = np.random.randint(D.size, size = D.size)

        tempD = D[sampleidx]
        tempN = N[sampleidx]

        slopes = np.append(slopes,np.polyfit(D**-1,N,1)[0])

    print 'After {:n} trys:\n\tmean: {:3.4f} mm\n\tstd: {:3.4f} mm'.\
        format(numtrys,slopes.mean(),slopes.std())

    return slopes

    

def get_sizes(inifiles,pp,EE=0.5):
    
    datadict = {}
    
    for ini in inifiles:
        
        #for L2
        # apsize = ini.split('m')[0]
        # backd = ini.split('m')[2][1:5]

        #for L3
        apsize = ini.split('m')[0].split('_')[1]
        backd = ini.split('m')[2][1:5]

        print (apsize,backd)
        
        datafile = yaki(ini)
        print 'Got {}, reducing...'.format(datafile)
        radius = omnom(datafile,pp,EEcut=EE)

        try: datadict[apsize][backd] = radius
        except KeyError: datadict[apsize] = {backd: radius}

    return datadict

def yaki(ini):

    options = ConfigParser.ConfigParser()
    options.read(ini)

    picklename = '.'.join(ini.split('.')[:-1])+'.pkl'
    
    try:
        N = pickle.load(open(picklename,'rb'))
        print "loaded pickel from {}".format(picklename)
    except IOError:
        N = Noodle(options)
        N.build_run()
        pickle.dump(N,open(picklename,'wb'))

    key = N.ratios.keys()[0]
    
    # for L3
    return N.ratios[key]['data']['V']['L3']['final']

    #for L2
    # return N.ratios[key]['data']['V']['direct']['final']

def omnom(fitsfile,pp,EEcut=0.5,fitsexten=0):

    data = pyfits.open(fitsfile)[fitsexten].data

    r, sb, err = ADE.annulize(data,300)
#    r *= 0.024
    
    flux = np.cumsum(sb)
    EE = flux/flux.max()

    cutr = np.where(EE >= EEcut)[0][0]

    EEfit = np.poly1d(np.polyfit(r[:cutr],EE[:cutr],2))

    fitr = np.linspace(r.min(),r.max(),500)
    fitEE = EEfit(fitr)

    r1 = np.interp(1.0,fitEE,fitr)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(r,EE,marker='.',linestyle='',markersize=0.5)
    ax.plot(fitr,fitEE,':')
    ax.set_xlabel('r [mm]')
    ax.set_ylabel('EE')
    ax.axvline(r1,linestyle='-',alpha=0.4)
    ax.axhline(1.0,linestyle=':',color='k',alpha=0.2)
    ax.set_ylim(0,1.1)
    ax.set_xlim(0,1.3*r1)
    ax.set_title('{}\nr: {:3.2f} mm'.format(fitsfile,r1))

    pp.savefig(fig)

    return r1


######################################



def nomlom(fitsfile,x,fitsexten=0):

    data = pyfits.open(fitsfile)[0].data
    
    r, sb, err = ADE.annulize(data,300)
#    r *= 0.024
    flux = np.cumsum(sb)
    EE = flux/flux.max()

    radii = np.array([])

    print ' Extracting radii for EE = 1/x:'
    
    for i in x:

        print '\tx = {:n}'.format(i)
        targetEE = 1./i
        radii = np.append(radii,np.interp(targetEE,EE,r))
    
    return radii
    
def size_get(inifiles):
    
    datadict = {}
    x = np.arange(1.05,5,0.5)

    for ini in inifiles:
        
        #for L2
        apsize = ini.split('m')[0]
        backd = ini.split('m')[2][1:5]

        #for L3
        # apsize = ini.split('m')[0].split('_')[1]
        # backd = ini.split('m')[2][1:5]

        print (apsize,backd)
        
        datafile = yaki(ini)
        print 'Got {}, reducing...'.format(datafile)
        radii = nomlom(datafile,x)
        
        for i,r in zip(x,radii):
            try: 
                datadict[apsize][i] = np.vstack((datadict[apsize][i],np.array([backd,r],dtype=np.float32)))
            except KeyError: 
                try: datadict[apsize][i] = np.array([backd,r],dtype=np.float32)
                except KeyError:
                    datadict[apsize] = {i: np.array([backd,r],dtype=np.float32)}

    return datadict

def sortio(searchstr,Nfigs,apfigs,bigfig):

    inifiles = glob(searchstr)
    datadict = size_get(inifiles)
    ratiodict = {}
    aps = np.array([])
    ratios = np.array([])

    pp = PDF(Nfigs)
    pp2 = PDF(apfigs)
    pp3 = PDF(bigfig)

    for apsize in datadict.keys():

        print '\nApsize = {:n} mm:\n\t{:7}{:7}\n\t'.format(int(apsize),'x','N_x')+'-'*10

        xvec = np.array([])
        Nprime = np.array([])
        
        sortedkeys = datadict[apsize].keys()[:]
        sortedkeys.sort()

        for x in sortedkeys[::-1]:
            d = datadict[apsize][x].T[0]
            r = datadict[apsize][x].T[1]
            
            Npcoef = np.polyfit(d,r,1)
            xvec = np.append(xvec,x)
            Nprime = np.append(Nprime,(2.*Npcoef[0])**-1)

            print '\t{:4.3f}  {:4.3f}'.format(x,(2.*Npcoef[0])**-1)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(d,r,marker='s',linestyle='')
            fitvec = np.linspace(d.min(),d.max(),50)
            ax.plot(fitvec,np.poly1d(Npcoef)(fitvec),':')
            ax.set_title('Aperture: {:3n} mm\n$x={:3.2f}\Rightarrow EE={:3.2f}$\n$N_x$: ${:3.2f}$'.\
                             format(int(apsize),x,1./x,(2.*Npcoef[0])**-1),
                         fontsize=10)
            ax.set_xlabel('d')
            ax.set_ylabel('r')

            pp.savefig()

        bigNcoef = np.polyfit(xvec,Nprime,1)

        print bigNcoef

        N = bigNcoef[0]
        print (2.*N)**-1
        Nfit = np.poly1d(bigNcoef)

        ratiodict[apsize] = {'f/#':N,
                             'fl':float(apsize)*N}

        aps = np.append(aps,float(apsize))
        ratios = np.append(ratios,N)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xvec,Nprime,marker='s',linestyle='')
        fitx = np.linspace(xvec.min(),xvec.max(),50)
        ax.plot(fitx,Nfit(fitx),'k:',label='linear fit')
        ax.set_xlabel('$x$',fontsize=14)
        ax.set_ylabel('$N_x$',fontsize=14)
        ax.legend(loc=0)
        ax.set_title('Aperture: {:n} mm\nN: {:3.4f}'.format(float(apsize),N))
        
        pp2.savefig(fig)

    pp.close()
    pp2.close()

    ratiofit = ADE.polyclip(aps**-1,ratios,1,niter=50,clip_high=1,clip_low=1)
#    ratiofit = np.poly1d(np.polyfit(aps**-1,ratios,1))
    ratioslope = ratiofit.c[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(aps**-1,ratios,marker='s',linestyle='')
    fitaps = np.linspace(aps.min(),aps.max(),50)
    ax.plot(fitaps**-1,ratiofit(fitaps**-1),':',label='linear fit')
    ax.set_xlabel('1/D [$mm^{-1}$]')
    ax.set_ylabel('f-ratio')
    #    ax.set_ylim(0, ratiofit(fitaps.max())*1.1)
    ax.legend(loc=0)
    ax.set_title('slope = {:n} mm'.format(ratioslope))
    
    pp3.savefig(fig)
    pp3.close()
    
    return ratiodict
