import numpy as np
import GradPak_plot as GPP
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF

AGES = np.array([5,25,100,286,640,904,1434,2500,5000,1000])

def plot_idl(inputfile):

    data = np.loadtxt(inputfile)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fibers = data[:,0]
    for i in np.arange(data.shape[1] - 1) + 1:
        
        ax.plot(fibers,data[:,i],'.',label='{} Myr'.format(AGES[i-1]))
        
    ax.legend(loc=0,numpoints=1)
    
    fig.show()

    return data
    
def plot_age(inputfile):

    data = np.loadtxt(inputfile)
    fibers = data[:,0]
    ageidx = np.argmax(data[:,1:],axis=1)
    
    ax = plt.figure().add_subplot(111)
    ax.set_xlabel('Fiber number')
    ax.set_ylabel('Age of highest light fraction population [Myr]')
    ax.set_yticks(range(10))
    ax.set_yticklabels(AGES)
    ax.set_ylim(-0.5,11)
    ax.plot(fibers,ageidx,'.')
    
    ax.figure.show()

    return data

def plot_age_hist(inputfile, outputfile):

    data = np.loadtxt(inputfile)
    fibers = data[:,0]

    pp = PDF(outputfile)
    LWAs = np.zeros(fibers.size)
    for i in range(fibers.size):
        print i
        ax = plt.figure().add_subplot(111)
        ax.bar(np.arange(AGES.size),data[i,1:],align='center')
        ax.set_ylabel('Light fraction')
        ax.set_xlim(-1,AGES.size)
        ax.set_xticks(np.arange(AGES.size))
        ax.set_xticklabels(AGES)
        ax.set_xlabel('Age [Myr]')
        MLWA = np.sum(AGES*data[i,1:])/np.sum(data[i,1:])
        LWAs[i] = MLWA
        ax.set_title('Fiber {}\nMLWA = {:5.0f} Myr'.format(i+1,MLWA))
        pp.savefig(ax.figure)
        plt.close(ax.figure)

    
    ax = plt.figure().add_subplot(111)
    ax.set_ylabel('MLWA [Myr]')
    ax.set_xlabel('Fiber')
    ax.plot(fibers,LWAs,'.',ms=10)
    ax.set_yscale('log')
    pp.savefig(ax.figure)
    plt.close(ax.figure)

    ax1 = GPP.plot_img(np.log10(LWAs+1),
                       clabel='Log( Mean Light Weighted Age [Myr] )',
                       method='cubic',
                       cmap='jet')
    
    pp.savefig(ax1.figure)
    plt.close(ax1.figure)
    
    pp.close()

    return LWAs

