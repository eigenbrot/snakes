import numpy as np


def compare_metals(data_list, model_list):

    try:
        comp_list = [np.loadtxt(i,usecols=(16,),unpack=True) for i in data_list]
    except IndexError:
        print i
        return

    comp_arr = np.abs(np.array(comp_list) - 1.)

    minidx = np.argmin(comp_arr,axis=0)
    best_models = np.array(model_list)[minidx].tolist()

    return minidx, best_models

def all_pointings():

    dirlist = ['../solar_Z','../0.2solar_Z','../0.02solar_Z',
               '../0.005solar_Z','../0.4solar_Z','../2.5solar_Z']
    modellist = ['/d/monk/eigenbrot/WIYN/14B-0456/anal/models/bc03_{}_ChabIMF.fits'.format(i) for i in ['solarZ','004Z','0004Z','0001Z','008Z','05Z']]
    fraclist = np.array([1,0.2,0.02,0.005,0.4,2.5])

    for i in range(6):

        outfile = 'P{}_models.dat'.format(i+1)
        
        datalist = ['{}/P{}.dat'.format(d,i+1) for d in dirlist]
        print datalist
        minidx, best_models = compare_metals(datalist,modellist)
        fracout = fraclist[minidx]

        with open(outfile,'w') as f:
            for frac, m in zip(fracout,best_models):
                f.write('{} {}\n'.format(frac, m))

            
    return
