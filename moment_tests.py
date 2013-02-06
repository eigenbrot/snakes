import numpy as np
import ADEUtils as ADE
import ADESALT as sa
import Salty2 as salty
import bottleneck as bn
import matplotlib.pyplot as plt

# Originally lived in ~/research/edgeOn/line_profiles

def do_line(simfile,radius,peak_scale):

    v, I, _ = salty.line_profile(simfile,radius,plot=False)
#    v, I = ADE.ADE_gauss(1000,500,50)
    I *= peak_scale/I.max()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(v,I)

    cdf = np.cumsum(I)
    cdf /= cdf.max()
    lowidx = int(np.interp(0.1,cdf,np.arange(cdf.size)))
    highidx = int(np.interp(0.9,cdf,np.arange(cdf.size))) + 1

    ax.axvline(x=v[lowidx])
    ax.axvline(x=v[highidx])

    moments = ADE.ADE_moments(v[lowidx:highidx], I[lowidx:highidx], threshold=np.inf)
    fig.show()

    return moments

def noise_line(simfile,radius,peak_scale):

#    v, I, _ = salty.line_profile(simfile,radius,plot=False)
    v, I = ADE.ADE_gauss(1000,500,50)
    I *= peak_scale/I.max()
    I += 3. * np.random.randn(I.size)
#    ADE.eplot(v,I)
    moments = ADE.ADE_moments(v, I, threshold=np.inf,err=np.abs(I)**0.5)

    return moments

def monte_lines(numtrys):
    
    bigarr = np.zeros((1,3))
    for i in range(numtrys):
        v, I = ADE.ADE_gauss(1000,500,50)
        I *= 55/I.max()
        I += 3. * np.random.randn(I.size)
    #    ADE.eplot(v,I)
        moments = ADE.ADE_moments(v, I, threshold=np.inf,err=np.abs(I)**0.5)
        bigarr = np.vstack((bigarr,moments))

    bigarr = bigarr[1:]
#    print bigarr

    return bn.nanmedian(bigarr,axis=0), bn.nanstd(bigarr,axis=0)

def line_comp(simfile,slayfile):
    '''compares simulation and actual data on line shapes, using a moment based
    approach'''

    dradii, dcent, derr, dm1, dm2, dm3 = sa.openslay(slayfile,moments=True)
    tm3 = np.array([])

    print "{:^20} {:^20} {:^20}\n{:^10}{:^10}".format('m1','m2','m3','model','data')
    print ("-"*20+'    ')*3
    for i in range(dradii.size):
        mm1, mm2, mm3 = do_line(simfile,dradii[i],1.)
        tm3 = np.append(tm3,mm3)
        print "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}".format(
            mm1,dm1[i,1],mm2,dm2[i,1],mm3,dm3[i,1])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(dradii,tm3,'.',label='model')
    ax.plot(dradii,dm3[:,1],'x',label='data')
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('$\mu_3$')
    ax.legend(loc=0)
    fig.show()

    return


