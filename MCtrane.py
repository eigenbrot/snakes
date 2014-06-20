import numpy as np
import Coltrane as trane
import pymc
#import Msix
import functools
import matplotlib.pyplot as plt
import drunkData as dd

class mark_VI(object):
    
    def __init__(self, drunkdict, pardict, name, size):
        
        priors = self.make_priors(pardict)
        self.funcdict = self.make_func_pars(pardict)
        function = self.make_function(drunkdict,name,size,priors)
        data = self.get_data(drunkdict)
        detvar = pymc.Normal('gal',mu=function, tau=0.01, 
                             value=data,observed=True)
        self.model = {}
        self.model['M'] = detvar
        for k in priors.keys():
            self.model[k] = priors[k]
        

    def make_priors(self, pardict):
        
        priordict = {}
        for k in pardict.keys():
            if not pardict[k][0]:
                name = k
                lowlim = pardict[k][1]
                highlim = pardict[k][2]
                priordict[k] = pymc.Uniform(name,lowlim,highlim)

        return priordict

    def make_func_pars(self,pardict):

        funcdict = {}
        for k in pardict.keys():
            if pardict[k][0]:
                funcdict[k] = pardict[k][1]
            else:
                funcdict[k] = None

        return funcdict
    
    def make_function(self,drunkdict,name,size,priordict):
        
        def modelled_galaxy_eval(**pdict):
            print pdict
            for k in pdict.keys():
                self.funcdict[k] = pdict[k]
            simfile = trane.make_boring([self.funcdict['Vr']],
                                        [self.funcdict['hrot']],
                                        h_dust=self.funcdict['h_dust'],
                                        kappa_0=self.funcdict['kappa_0'],
                                        z_d=self.funcdict['z_d'],
                                        name=name,size=size,z=0,
                                        flarepars=False)[0]

            bar = trane.moments_notice(drunkdict[0][0],simfile,
                               skip_radii=drunkdict[0][2],
                               flip=drunkdict[0][1])

            out = np.r_[bar[1][2],bar[2][2],bar[3][2]]

            return out

        base_func = functools.partial(modelled_galaxy_eval,**priordict)
        modelled_galaxy = pymc.Deterministic(eval = modelled_galaxy_eval,
                                             name = 'modelled_galaxy',
                                             parents = priordict,
                                             doc = 'Moments of sight lines',
                                             trace = True,
                                             verbose = 0,
                                             plot = False)

        return modelled_galaxy

    def get_data(self,drunkdict):

        out = np.array([])
        for z in drunkdict.keys():
            _,_,_, m1, m2, m3 = dd.open_drunk(drunkdict[z][0])
            out = np.r_[out,m1[0],m2[0],m3[0]]

        return out

def test(drunkdict,pardict,sample,burn=10):

    # pardict = {'h_dust': [True, 8.43], 
    #            'Vr': [False, 200, 300], 
    #            'hrot': [False, 4, 5.], 
    #            'kappa_0': [True, 0.652], 
    #            'z_d': [True, 0.43]}
    # drunkdict = {0:['ESO_z0_drunk.fits',False,[]]}

    sax = mark_VI(drunkdict,pardict,'test',1001)
    S = pymc.MCMC(sax.model)
    S.sample(sample,burn=burn)
    traces = {}
    bestfit = {}
    for k in pardict.keys():
        if not pardict[k][0]:
            trace = S.trace(k)[:]
            mean = np.mean(trace)
            traces[k] = trace
            bestfit[k] = mean
            pardict[k][1] = mean

    # Vr_sample = S.trace('Vr')[:]
    # hrot_sample = S.trace('hrot')[:]
    
    # Vr = np.mean(Vr_sample)
    # hrot = np.mean(hrot_sample)
    
    simfile = trane.make_boring([pardict['Vr'][1]],
                                [pardict['hrot'][1]],
                                h_dust=pardict['h_dust'][1],
                                kappa_0=pardict['kappa_0'][1],
                                z_d=pardict['z_d'][1],
                                name='test',size=1001,z=0,
                                flarepars=False)[0]
    
    # simfile = trane.make_boring([Vr],[pars[1]],name='final',
    #                             size=1001,z=0,h_dust=pars[3],kappa_0=pars[2],
    #                             z_d=pars[4],flarepars=False)[0]

    bar = trane.moments_notice(drunkdict[0][0],simfile,
                               skip_radii=drunkdict[0][2],
                               flip=drunkdict[0][1])

    # fig = plt.figure()
    # ax = fig.add_subplot(211)
    # ax.hist(Vr_sample,bins=50,histtype='step')
    # ax.set_title('V$_r$')
    # ax2 = fig.add_subplot(212)
    # ax2.hist(hrot_sample,bins=50,histtype='step')
    # ax2.set_title('h$_{rot}$')
    # fig.show()

    return bestfit, bar
