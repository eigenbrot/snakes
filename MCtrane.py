import numpy as np
import Coltrane as trane
import pymc
import functools
import matplotlib.pyplot as plt
import drunkData as dd
import copy
from matplotlib.backends.backend_pdf import PdfPages as PDF

class mark_VI(object):
    
    def __init__(self, drunkdict, pardict, name, size):
        
        priors = self.make_priors(pardict)
        self.funcdict = self.make_func_pars(pardict)
        datadict = self.make_data_dict(drunkdict)
        function = self.make_function(datadict,name,size,priors)
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
    
    def make_data_dict(self,drunkdict):

        datadict = {}
        for z in drunkdict.keys():
            datadict[z] = [dd.open_drunk(drunkdict[z][0],skip_radii=drunkdict[z][2]),
                           drunkdict[z][1],
                           drunkdict[z][2]]

        return datadict

    def make_function(self,datadict,name,size,priordict):
        
        def modelled_galaxy_eval(**pdict):
            #print pdict
            for k in pdict.keys():
                self.funcdict[k] = pdict[k]
            bigm1 = np.array([])
            bigm2 = np.array([])
            bigm3 = np.array([])
            for z in datadict.keys():
                
                simfile = trane.make_boring([self.funcdict['Vr']],
                                            [self.funcdict['hrot']],
                                            h_dust=self.funcdict['h_dust'],
                                            kappa_0=self.funcdict['kappa_0'],
                                            z_d=self.funcdict['z_d'],
                                            name=name,size=size,z=z,
                                            flarepars=False,nofits=True)[0]
                
                _, m1, m2, m3 = trane.moments_notice(datadict[z][0],simfile,
                                           skip_radii=datadict[z][2],
                                           flip=datadict[z][1],nofits=True)

                bigm1 = np.append(bigm1,m1[2])
                bigm2 = np.append(bigm2,m2[2])
                bigm3 = np.append(bigm3,m3[2])

                out = np.r_[bigm1,bigm2]#,bigm3]
            if np.isnan(np.sum(out)):
                print '!!!!!!!!NAN!!!!!!!!'
                print out
                raw_input('')
            #print out.shape
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

        bigm1 = np.array([])
        bigm2 = np.array([])
        bigm3 = np.array([])
        for z in drunkdict.keys():
            _,_,_, m1, m2, m3 = dd.open_drunk(drunkdict[z][0],
                                              skip_radii=drunkdict[z][2])
            bigm1 = np.append(bigm1,m1[0])
            bigm2 = np.append(bigm2,m2[0])
            bigm3 = np.append(bigm3,m3[0])
        
        out = np.r_[bigm1,bigm2]#,bigm3]
        return out

def test(drunkdict,pardict,output,sample,burn=10):

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
    outdict = copy.deepcopy(pardict)
    for k in outdict.keys():
        if not outdict[k][0]:
            trace = S.trace(k)[:]
            mean = np.mean(trace)
            traces[k] = trace
            bestfit[k] = mean
            outdict[k][1] = mean
    
    bars = []
    datadict = {}
    for z in drunkdict.keys():
        datadict[z] = [dd.open_drunk(drunkdict[z][0],skip_radii=drunkdict[z][2]),
                       drunkdict[z][1],
                       drunkdict[z][2]]
    for z in datadict.keys():
        simfile = trane.make_boring([outdict['Vr'][1]],
                                    [outdict['hrot'][1]],
                                    h_dust=outdict['h_dust'][1],
                                    kappa_0=outdict['kappa_0'][1],
                                    z_d=outdict['z_d'][1],
                                    name='final',size=1001,z=z,
                                    flarepars=False)[0]
    
        bar = trane.moments_notice(datadict[z][0],simfile,
                                   skip_radii=datadict[z][2],
                                   flip=datadict[z][1])

        bars.append(bar)

    pp = PDF(output)
    for k in traces.keys():
        ax = plt.figure().add_subplot(111)
        ax.hist(traces[k],bins=50,histtype='step')
        ax.set_xlabel(k)
        ax.set_ylabel('PDF')
        ax.set_title('Most likely value:\n{:9.4f}'.format(bestfit[k]))
        pp.savefig(ax.figure)

    pp.close()

    return bestfit, bars, traces
