import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF

DFK_borders = [[0.9e-3,5.2e-3],
               [5.5e-3,0.404],
               [0.4535,5.75],
               [6,13.5]]

def get_weighted_age(tau_sf, agelim):

    return tau_sf*np.log(0.5*(\
                              np.exp(agelim[1]/tau_sf) \
                              + np.exp(agelim[0]/tau_sf)))

def compare1(masses,index=None):

    # tau_sf_list = [-5.,-1.,-0.1,0.1,5.,7,10,15.]
    tau_sf_list = np.linspace(-10,20,100)
    MMWA = np.zeros(len(tau_sf_list))

    for i, tsf in enumerate(tau_sf_list):
        ages = [get_weighted_age(tsf,dkb) for dkb in DFK_borders]
        if index is None:
            MMWA[i] = np.sum(np.array(ages)*masses)/np.sum(masses)
        else:
            MMWA[i] = np.sum(np.array(ages[:index])*masses[:index])\
                      /np.sum(masses[:index])
        
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(r'$\tau_{SF}$')
    ax.set_ylabel('MMWA')
    ax.plot(tau_sf_list, MMWA)
    
    ax.figure.show()

    print (np.nanmax(MMWA) - np.nanmin(MMWA))/np.nanmean(MMWA)

    return

def compare2(masses):

    tau_sf_list = np.linspace(-10,20,100)
    age_metric = np.zeros(len(tau_sf_list))

    for i, tsf in enumerate(tau_sf_list):
        ages = [get_weighted_age(tsf,dkb) for dkb in DFK_borders]
        binsize = np.array([dkb[1] - dkb[0] for dkb in DFK_borders])
        print binsize
        age_metric[i] = np.sum(np.array(ages)*masses/binsize)/np.sum(masses)
        
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(r'$\tau_{SF}$')
    ax.set_ylabel('age metric')
    ax.plot(tau_sf_list, age_metric)

    ax.figure.show()

    return

def compare3(masses, index=-1):

    tau_sf_list = np.linspace(-10,20,100)
    youngMMWA = np.zeros(tau_sf_list.size)
    ratio = np.zeros(tau_sf_list.size)
    
    for i, tsf in enumerate(tau_sf_list):
        ages = [get_weighted_age(tsf,dkb) for dkb in DFK_borders]
        youngMMWA[i] = np.sum(np.array(ages[:index])*masses[:index])\
                       /np.sum(masses[:index])
        ratio[i] = np.sum(masses[index:])/np.sum(masses[:index])

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(r'$\tau_{L,young}$')
    ax.set_ylabel(r'$L_{old}/L_{young}$')
    s = ax.scatter(youngMMWA,ratio,c=tau_sf_list, cmap=plt.cm.gnuplot)
    cb = ax.figure.colorbar(s)
    cb.set_label(r'$\tau_{SF}$')

    ax.figure.show()

    return

def compare_Dn4000(output,tausf=5):

    import D4000_indices as DI
    import pyfits
    allZborders = DFK_borders*6

    LMdir = '.'
    D4000dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/indecies/D4000/tau_grid'
    
    ratiolist = []
    younglist = []
    D4000list = []
    # index = [0,1,4,5,8,9,12,13,16,17,20,21]
    # oindex = [2,3,6,7,10,11,14,15,18,19,22,23]
    index = [0,1,2,
             4,5,6,
             8,9,10,
             12,13,14,
             16,17,18,
             20,21,22]
    oindex = [3,7,11,15,19,23]

    for p in range(6):
        coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(LMdir,p+1)
        print coeffile
        coefs = pyfits.open(coeffile)[1].data
        ages = [get_weighted_age(tausf,dkb) for dkb in allZborders]
        num = np.sum(np.array(ages)[None,index]*\
               coefs['LIGHT_WEIGHT'][:,index],axis=1)
        print num.shape
        print coefs['LIGHT_WEIGHT'].shape
        print np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1).shape
        younglist.append(num/np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1))
        print coefs['LIGHT_WEIGHT'][:,index].shape
        print coefs['LIGHT_WEIGHT'][:,oindex].shape
        ratiolist.append(np.sum(coefs['LIGHT_WEIGHT'][:,oindex],axis=1)\
                     /np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1))

        D4000file = '{}/NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(D4000dir,p+1)
        print D4000file
        res = DI.quick_eat(D4000file)
        D4000list.append(res[:,2])

    ratio = np.hstack(ratiolist)
    print ratio.shape
    young = np.hstack(younglist)
    print young.shape
    D4000res = np.hstack(D4000list)
    print D4000res.shape

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ratio, young, c=D4000res, cmap=plt.cm.gnuplot)

    fig.show()
    
    return

def compare_ratio_Dn4000(output):

    import D4000_indices as DI
    import pyfits
    allZborders = DFK_borders*6

    LMdir = '.'
    D4000dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/indecies/D4000/tau_grid'
    
    ratiolist = []
    zlist = []
    D4000list = []
    index = [0,1,4,5,8,9,12,13,16,17,20,21]
    oindex = [2,3,6,7,10,11,14,15,18,19,22,23]
    # index = [0,1,2,
    #          4,5,6,
    #          8,9,10,
    #          12,13,14,
    #          16,17,18,
    #          20,21,22]
    # oindex = [3,7,11,15,19,23]

    for p in range(6):
        coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(LMdir,p+1)
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(p+1)
        z = np.loadtxt(loc,usecols=(5,), unpack=True)
        zlist.append(np.abs(z))
        print coeffile
        coefs = pyfits.open(coeffile)[1].data
        ratiolist.append(np.sum(coefs['LIGHT_WEIGHT'][:,oindex],axis=1)\
                     /np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1))

        D4000file = '{}/NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(D4000dir,p+1)
        print D4000file
        res = DI.quick_eat(D4000file)
        D4000list.append(res[:,2])

    ratio = np.hstack(ratiolist)
    print ratio.shape
    D4000res = np.hstack(D4000list)
    print D4000res.shape
    z = np.hstack(zlist)
    print z.shape
    gid = np.isfinite(ratio)
    ratio = ratio[gid]
    D4000res = D4000res[gid]
    z = z[gid]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(z, ratio, c=D4000res, cmap=plt.cm.gnuplot, vmax=2.0, linewidth=0)
    ax.set_ylim(0,np.median(ratio)*5)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    
    return ratio, D4000res

def compare_young_Dn4000(output,tausf=5):

    import D4000_indices as DI
    import pyfits
    allZborders = DFK_borders*6

    LMdir = '.'
    D4000dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/indecies/D4000/tau_grid'
    
    younglist = []
    zlist = []
    D4000list = []
    # index = [0,1,4,5,8,9,12,13,16,17,20,21]
    # oindex = [2,3,6,7,10,11,14,15,18,19,22,23]
    index = [0,1,2,
             4,5,6,
             8,9,10,
             12,13,14,
             16,17,18,
             20,21,22]
    oindex = [3,7,11,15,19,23]

    for p in range(6):
        coeffile = '{}/NGC_891_P{}_bin30_allz2.coef.fits'.format(LMdir,p+1)
        loc = 'NGC_891_P{}_bin30_locations.dat'.format(p+1)
        z = np.loadtxt(loc,usecols=(5,), unpack=True)
        zlist.append(np.abs(z))
        print coeffile
        coefs = pyfits.open(coeffile)[1].data
        ages = [get_weighted_age(tausf,dkb) for dkb in allZborders]
        num = np.sum(np.array(ages)[None,index]*\
               coefs['LIGHT_WEIGHT'][:,index],axis=1)
        print num.shape
        print coefs['LIGHT_WEIGHT'].shape
        print np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1).shape
        younglist.append(num/np.sum(coefs['LIGHT_WEIGHT'][:,index],axis=1))
        
        D4000file = '{}/NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(D4000dir,p+1)
        print D4000file
        res = DI.quick_eat(D4000file)
        D4000list.append(res[:,2])

    young = np.hstack(younglist)
    print young.shape
    D4000res = np.hstack(D4000list)
    print D4000res.shape
    z = np.hstack(zlist)
    print z.shape
    gid = np.isfinite(young)
    young = young[gid]
    D4000res = D4000res[gid]
    z = z[gid]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(z, young, c=D4000res, cmap=plt.cm.gnuplot, vmax=2.0, linewidth=0)
    ax.set_ylim(0,np.median(young)*5)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    
    return young, D4000res

def compute_tausf_err(tausflist=None):

    import pyfits
    allZborders = DFK_borders*6
    if tausflist is None:
        tausflist = np.linspace(-10,20,100)
    agearr = np.vstack([[get_weighted_age(tsf,dkb) for dkb in allZborders] for tsf in tausflist])
    print agearr.shape

    for p in range(6):
        coeffile = 'NGC_891_P{}_bin30_allz2.coef.fits'.format(p+1)
        coefs = pyfits.open(coeffile)[1].data
        numfibs = coefs['MLWA'].size
        
        output = np.zeros(numfibs, dtype=[('MLWA', '>f8'), ('dMLWA','>f8')])
        
        for i in range(numfibs):
            lw = coefs['light_weight'][i]
            tmp = np.zeros(len(tausflist))
            for j in range(len(tausflist)):
                tmp[j] = np.sum(agearr[j]*lw)/np.sum(lw)

            output['MLWA'][i] = np.nanmean(tmp)
            output['dMLWA'][i] = np.nanstd(tmp)
            
        pyfits.BinTableHDU(output).writeto('NGC_891_P{}_bin30_allz2.aI.fits'.format(p+1),clobber=True)

    return

def SFH_plots(output):
    import pyfits

    #Get centers (just for plotting) of DFK age bins in log-age space
    lb = np.log10(DFK_borders)
    ages = np.sum(lb,axis=1)/2.
    print ages
    print lb

    for p in range(6):
        coef = 'NGC_891_P{}_bin30_allz2.coef.fits'.format(p+1)
        data = pyfits.open(coef)[1].data
        weights = data['LIGHT_WEIGHT']
        errs = data['LIGHT_WEIGHT_ERR']

        print coef

        pp = PDF('{}_P{}.pdf'.format(output,p+1))
        for i in range(data.shape[0]):
            print '\t', i+1
            ax = plt.figure().add_subplot(111)
            ax.set_xlabel('Log(Gyr)')
            ax.set_ylabel(r'$\int \psi(t) dt$')
            
            lw = np.sum(weights[i,:].reshape(6,4),axis=0)
            lwe = np.sqrt(np.sum((errs[i,:]**2).reshape(6,4),axis=0))/4.
            #ax.plot(10**ages, lw, 'k', alpha=0.6)
            ax.scatter(10**ages, lw, marker='o', s=10, color='k')
            ax.errorbar(10**ages, lw, yerr=lwe, fmt='none',
                        capsize=0, ecolor='k')
            ax.hlines(lw,10**lb[:,0],10**lb[:,1])
            ax.set_xscale('log')
            ax.set_ylim(0,ax.get_ylim()[1])
        
            pp.savefig(ax.figure)
            plt.close(ax.figure)

        pp.close()
        plt.close('all')

    return

def compute_random_err(N=100):

    import pyfits
    borders = np.array(DFK_borders)

    agearr = np.random.random(N*borders.shape[0]).reshape(borders.shape[0],N)
    print agearr.shape
    d = np.diff(borders,axis=1)
    c = borders[:,0][:,None]
    agearr = agearr*d + c
    print agearr
    print agearr.shape
    agearr = np.vstack([agearr]*6)
    print agearr.shape

    for p in range(6):
        coeffile = 'NGC_891_P{}_bin30_allz2.coef.fits'.format(p+1)
        coefs = pyfits.open(coeffile)[1].data
        numfibs = coefs['MLWA'].size
        
        output = np.zeros(numfibs, dtype=[('MLWA', '>f8'), ('dMLWA','>f8')])
        
        for i in range(numfibs):
            lw = coefs['light_weight'][i]
            tmp = np.zeros(N)
            for j in range(N):
                tmp[j] = np.sum(agearr[:,j]*lw)/np.sum(lw)

            output['MLWA'][i] = np.nanmean(tmp)
            output['dMLWA'][i] = np.nanstd(tmp)
            
        pyfits.BinTableHDU(output).writeto('NGC_891_P{}_bin30_allz2.aIR.fits'.format(p+1),clobber=True)

    return
