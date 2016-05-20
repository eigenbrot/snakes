import sys
from glob import glob
import time
import numpy as np
import scipy.stats as ss
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import ImageGrid
plt.ioff()

Zlist = [0.005,0.02,0.2,0.4,1,2.5]
flist = ['{}Z/steps'.format(i) for i in Zlist]

Zlist_all = [0.005,0.02,0.2,0.4,1,2.5,'all']
flist_all = ['{}Z/steps'.format(i) for i in Zlist_all]
offset_all = [0,0,0,0,0,0,20]
numfree_all = [1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-25-1]

Zlist_all2 = [0.005,0.02,0.2,0.4,1,2.5,'all','prior','mid']
flist_all2 = ['{}Z/steps'.format(i) for i in Zlist_all2]
offset_all2 = [0,0,0,0,0,0,20,20,20]
numfree_all2 = [1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-5-1,1428-25-1,1428-25-1,1428-25-1]

def plot(stepfile, offset=0):

    MLWA, TAUV, Chi, blueChi = np.loadtxt(stepfile,usecols=(6+offset,
                                                            10+offset,
                                                            12+offset,
                                                            14+offset),unpack=True)
    ages = np.loadtxt(stepfile,usecols=range(1,5+offset))
    t1 = np.zeros(ages.shape[0])

    #SSP_ages = [0.005, 0.025, 0.102, 0.286, 0.641, 0.905, 1.434, 2.5, 5., 10.]
    SSP_ages = [0.003,0.216,3.151,9.688]
    for s in range(ages.shape[0]):
        i = 0
        try:
            while ages[s,i] == 0:
                i += 1
            t1[s] = SSP_ages[i]
        except IndexError:
            t1[s] = 0.0

    
    print Chi.size
    Chi = Chi[MLWA == MLWA]
    blueChi = blueChi[MLWA == MLWA]
    MLWA = MLWA[MLWA == MLWA]
    TAUV = TAUV[MLWA == MLWA]
    Chi = Chi[Chi < 100]
    blueChi = blueChi[Chi < 100]
    MLWA = MLWA[Chi < 100]
    TAUV = TAUV[Chi < 100]
    print Chi.size
    Chi -= np.min(Chi)
    blueChi -= np.min(blueChi)

    ax = plt.figure().add_subplot(111)
    ax.plot(MLWA,Chi,'.',color='r',label='Full $\Delta\chi^2$')
    ax.plot(MLWA,blueChi,'.',color='b',label='Blue $\Delta\chi^2$')
    ax.set_xlabel('MLWA [Gyr]')
    ax.set_ylabel('$\Delta\chi^2$')
#    ax.set_ylim(0,100)
    ax.legend(loc=0, frameon=False, numpoints=1, scatterpoints=1)
    ax.figure.show()

    # ax2 =plt.figure().add_subplot(111)
    # ax2.plot(t1,Chi,'.',color='r',label='Full $\Delta\chi^2$')
    # ax2.plot(t1,blueChi,'.',color='b',label='Blue $\Delta\chi^2$')
    # ax2.set_xlabel('"t$_{\mathrm{form}}$"')
    # ax2.set_ylabel('$\Delta\chi^2$')
    # ax2.legend(loc=0, frameon=False, numpoints=1, scatterpoints=1)
    # ax2.figure.show()

    ax3 =plt.figure().add_subplot(111)
    ax3.plot(TAUV,Chi,'.',color='r',label='Full $\Delta\chi^2$')
    ax3.plot(TAUV,blueChi,'.',color='b',label='Blue $\Delta\chi^2$')
    ax3.set_xlabel(r'$\tau_V$')
    ax3.set_ylabel('$\Delta\chi^2$')
    ax3.set_ylim(0,100)
    ax3.legend(loc=0, frameon=False, numpoints=1, scatterpoints=1)
    ax3.figure.show()
        
    ax3d = plt.figure().add_subplot(111, projection='3d')
    ax3d.set_xlabel(r'$\tau_V$')
    ax3d.set_ylabel('MLWA')
    ax3d.set_zlabel('$\Delta\chi^2$')
    ax3d.scatter(TAUV,MLWA,zs=Chi)
    ax3d.figure.show()

    if __name__ == '__main__':
        raw_input('Enter to exit')

def plot_multiZ(prefix,minval=0,maxval=200,bluefree=702):

    fraclist = np.array([0.005,0.02,0.2,0.4,1,2.5])
    rw = 0.5
    rh = 0.05

    minchi = np.inf
    for z in range(fraclist.size):
        
        stepfile = '{}_Z{:04}_steps.dat'.format(prefix,int(fraclist[z]*1000))
        blueChi = np.loadtxt(stepfile,usecols=(17,),unpack=True)
        blueChi *= bluefree
        if blueChi[-1] < minchi:
            minchi = blueChi[-1]

    fig = plt.figure()
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%',
                     axes_class = None)
    ax = grid[0]

    ax.set_xlabel('$Z/Z_{\odot}$')
    ax.set_ylabel('MLWA')
    for z in range(fraclist.size):

        stepfile = '{}_Z{:04}_steps.dat'.format(prefix,int(fraclist[z]*1000))
        MLWA, TAUV, Chi, blueChi = np.loadtxt(stepfile,usecols=(12,13,15,17),
                                              unpack=True)
        blueChi *= bluefree
        
#        good = np.where(blueChi < 80)
#        blueChi = blueChi[good]
#        MLWA = MLWA[good]
        blueChi -= minchi
        patches = []
        print np.log10(blueChi.min()), np.log10(blueChi.max())
        for i in range(MLWA.size):
#            patches.append(Circle((z,MLWA[i]),radius=0.1,linewidth=0))
            width = rw * ((minchi+5000)/(blueChi[i]+5000))
#            width = rw * (minchi/blueChi[i])
            patches.append(Rectangle((z-width/2,MLWA[i]-rh/2),width,rh))

        collection = PatchCollection(np.array(patches)[::-1],
                                     cmap=plt.get_cmap('gnuplot'),
                                     norm=matplotlib.colors.Normalize(
                                         vmin=minval,vmax=maxval))
        collection.set_alpha(0.1)
        collection.set_linewidth(0.0)
        collection.set_array(np.log10(blueChi)[::-1])
        ax.add_collection(collection)
        
#        ax.axhline(y=MLWA[-1],color='k',ls=':',alpha=0.6)
        ax.hlines(MLWA[-1],z - 0.5, z + 0.5,color='g',lw=2)
        ax.text(z + 0.5,MLWA[-1], '{:5.2f}%, {:4.1f}'.\
                format(100 - 100*ss.chi2.cdf(blueChi[-1],bluefree+5),
                       blueChi[-1]))

    collection.set_alpha(1.0)
    cbar = ax.cax.colorbar(collection)

    ci = [68.27, 50.0, 80.0, 40.0]
    dc = [4.72,  3.36, 5.99, 2.75]
    
    # for (c, d) in zip(ci, dc):
    #     ax.cax.axvline(x = np.log10(d), color='g')
    #     ax.cax.text(np.log10(d) - 0.03, 2, '{}'.format(c),
    #                 transform=ax.cax.transData,
    #                 fontsize=8)

    collection.set_alpha(0.1)
    cbar.set_label_text('Log( $\Delta\chi^2_{\mathrm{blue}}$ )')

    ax.set_xlim(-1,fraclist.size + 1)
    ax.set_ylim(0,12)
    ax.set_xticks(np.arange(fraclist.size))
    ax.set_xticklabels(fraclist)
    fig.show()
        
    return fig

def plot_multifolder_old(folder_list, Zlist, pointing, ap):
    rw = 0.5
    rh = 0.1
    numfree = 1334 - 5 - 1

    minchi = np.inf
    maxchi = 0
    for f in folder_list:
        stepfile = glob('{:}/P{:}_*{:02n}_steps.dat'.format(f,pointing,ap))[0]
        print stepfile
        chi = np.loadtxt(stepfile,usecols=(12,),unpack=True)
        chi = chi[chi == chi]
        upperlim = np.median(chi) + 1.3*np.min(chi)
        if chi[-1] < minchi:
            minchi = chi[-1]
        if upperlim > maxchi:
            maxchi = upperlim

    fig = plt.figure()
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%',
                     axes_class = None)
    ax = grid[0]

    ax.set_xlabel('$Z/Z_{\odot}$')
    ax.set_ylabel('MLWA')
    goodage = []
    bigZ = np.zeros(len(folder_list))
    for z,f in enumerate(folder_list):

        stepfile = glob('{:}/P{:}_*{:02n}_steps.dat'.format(f,pointing,ap))[0]
        print stepfile
        MLWA, chi = np.loadtxt(stepfile,usecols=(6,12),
                               unpack=True)
        idx = chi == chi
        chi = chi[idx]
        MLWA = MLWA[idx]
        idx = MLWA == MLWA
        chi = chi[idx]
        MLWA = MLWA[idx]

#        good = np.where(chi < 80)
#        chi = chi[good]
#        MLWA = MLWA[good]
        chi -= minchi
        patches = []
        for i in range(MLWA.size):
#            patches.append(Circle((z,MLWA[i]),radius=0.1,linewidth=0))
            width = rw * ((minchi+5000)/(chi[i]+5000))
#            width = rw * (minchi/chi[i])
            patches.append(Rectangle((z-width/2,MLWA[i]-rh/2),width,rh))

        collection = PatchCollection(np.array(patches)[::-1],
                                     cmap=plt.get_cmap('gnuplot'),
                                     norm=matplotlib.colors.Normalize(
                                         vmin=minchi,vmax=maxchi))
        collection.set_alpha(0.7)
        collection.set_linewidth(0.0)
        collection.set_array(chi[::-1])
        ax.add_collection(collection)
        
#        ax.axhline(y=MLWA[-1],color='k',ls=':',alpha=0.6)
        ax.hlines(MLWA[-1],z - 0.3, z + 0.3,color='g',lw=2)
        prob = 1 - ss.chi2.cdf(chi[-1]*numfree,numfree)
        if prob >= 0.68:
            goodage.append(MLWA[-1])
            bigZ[z] = 1

        ax.text(z + 0.5,MLWA[-1],
                '{:4.2f}\n{:4.2f}'.format(prob,chi[-1]),
                fontsize=10)

    collection.set_alpha(1.0)
    cbar = ax.cax.colorbar(collection)

    ci = [68.27, 50.0, 80.0, 40.0]
    dc = [4.72,  3.36, 5.99, 2.75]
    
    # for (c, d) in zip(ci, dc):
    #     ax.cax.axvline(x = np.log10(d), color='g')
    #     ax.cax.text(np.log10(d) - 0.03, 2, '{}'.format(c),
    #                 transform=ax.cax.transData,
    #                 fontsize=8)

    collection.set_alpha(0.1)
    cbar.set_label_text(r'$\Delta\chi^2_{\nu}$')

    ax.set_xlim(-1,len(Zlist) + 1)
    ax.set_ylim(0,12)
    ax.set_xticks(np.arange(len(Zlist)))
    ax.set_xticklabels(Zlist)
#    fig.show()
        
    return fig, np.mean(goodage), np.std(goodage), bigZ

def plot_multifolder(folder_list, Zlist, pointing, ap, numfree=1428-5-1, offset=0,col=6,label='MLWA'):
    clist = ['blue','green','red','orange','purple','black','maroon']
    clist = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6']

    if not isinstance(offset,list):
        offset = [offset] * len(folder_list)
    if not isinstance(numfree,list):
        numfree = [numfree] * len(folder_list)

    minchi = np.inf
    maxchi = 0
    for i, f in enumerate(folder_list):
        stepfile = glob('{:}/P{:}_*{:02n}_steps.dat'.format(f,pointing,ap))[0]
        print stepfile
        MLWA, chi = np.loadtxt(stepfile,usecols=(col+offset[i], 12+offset[i],),unpack=True)
        idx = chi == chi
        chi = chi[idx]
        MLWA = MLWA[idx]
        upperlim = np.median(chi) + 1.3*np.min(chi)
        if np.min(chi) < minchi:
            minchi = np.min(chi)
            bestnumfree = numfree[i]
            bestMLWA = MLWA[np.argmin(chi)]
            bestZidx = i
        if upperlim > maxchi:
            maxchi = upperlim

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_title('P{}.{}\n{}'.format(pointing,ap,time.asctime()))

    ax.set_xlabel(label)
    ax.set_ylabel(r'$\Delta\chi^2_{\nu}$')
    goodage = np.array([])
    bigZ = np.zeros(len(folder_list))
    avgdiff = 0
    for z,f in enumerate(folder_list):

        stepfile = glob('{:}/P{:}_*{:02n}_steps.dat'.format(f,pointing,ap))[0]
        print stepfile
        MLWA, chi = np.loadtxt(stepfile,usecols=(col+offset[z],12+offset[z]),
                               unpack=True)
        idx = chi == chi
        chi = chi[idx]
        MLWA = MLWA[idx]
        idx = MLWA == MLWA
        chi = chi[idx]
        MLWA = MLWA[idx]

        chi -= minchi
        avgdiff += np.median(chi)

        ax.scatter(MLWA, chi, label=f, color=clist[z], s=60, alpha=0.5, edgecolors='none')
        
        prob = 1 - ss.chi2.cdf(np.min(chi)*numfree[z],numfree[z])
        if prob >= 0.68:
            goodidx = np.argmin(chi)
            goodage = np.r_[goodage,MLWA[goodidx]]
            bigZ[z] = 1
            ax.scatter(MLWA[goodidx],np.min(chi),color=clist[z],s=120)
        
        if np.min(chi) == 0.0:
            ax.axvline(MLWA[goodidx],ls='--',color=clist[z])

        # ax.text(chi[-1],MLWA[-1],
        #         '{:4.2f}\n{:4.2f}'.format(prob,chi[-1]),
        #         fontsize=10)

    avgdiff /= len(folder_list)
    ax.legend(loc=0,numpoints=1,scatterpoints=1)
    ax.set_ylim(-0.1*avgdiff,avgdiff*5)
    ax.set_xlim(0,11)
    
    cutoff = ss.chi2.ppf(0.68,bestnumfree)/bestnumfree
    ax.axhline(cutoff,ls=':',color='k',alpha=0.6)
    # ax.axvline(bestMLWA,ls='-',color='k',alpha=0.6)
    # ax.text(np.mean(goodage),avgdiff*5*0.9,'{:} = {:3.1f}'.format(label,np.mean(goodage)),
    #         ha='right',va='bottom',fontsize=8)
#    fig.show()
        
    return fig, bestMLWA, bestMLWA - np.min(goodage),\
        np.max(goodage) - bestMLWA, minchi, bestZidx, bigZ


def do_pointing(folder_list, Zlist, pointing, numaps, output, numfree=1428-5-1,offset=0,
                col=6,label='MLWA'):
    #offset = 20 for multiZ
    pp = PDF(output+'.pdf')
    txtfile = output+'.dat'
    with open(txtfile,'w') as f:
        f.write('#{:>4}{:>10}{:>13}{:>13}{:>13}{:>8}'.format('ap',label,'l'+label,'h'+label,
                                                             'minchi','Zidx'))
        f.write(str('{:>7}Z'*len(Zlist)).format(*Zlist))
        f.write('\n#\n')
        for a in range(numaps):
            o = plot_multifolder(folder_list, Zlist, pointing, a+1,
                                 numfree=numfree,offset=offset,
                                 col=col,label=label)
            fig, age, low, high, minchi, bestZidx, bigZ = o
            pp.savefig(fig)
            
            f.write('{:5n}{:10.3f}{:13.4f}{:13.4f}{:13.4f}{:8n}'.format(a+1,age,low,high,
                                                                        minchi, bestZidx))
            f.write(str('{:8n}'*bigZ.size).format(*bigZ))
            f.write('\n')

    pp.close()
    plt.close('all')
    return

def do_all_pointings(suff='',folder_list=flist, Zlist=Zlist,
                     numfree=1428-5-1,offset=0,col=6,label='MLWA'):
    
    numlist = [37,38,59,60,29,38]
    if suff != '':
        suff = '_' + suff
    for i in range(6):
        output = 'NGC_891_P{}_CI{}'.format(i+1,suff)
        do_pointing(folder_list, Zlist, i+1, numlist[i], output,
                    offset=offset,col=col,label=label,numfree=numfree)

    return

def plot_all_pointings(output,suff='',label='MLWA'):
    
    import plot_allZ2 as pa2

    pp = PDF(output)
    al = pa2.plot_heights_with_err('*'+suff,label=label,col=1,errcol=2, ylims=[0,11],
                                   plotfit=False,exclude=pa2.exclude)
    [pp.savefig(a.figure) for a in al]
    pp.close()

    return

if __name__ == '__main__':
    
    plot(sys.argv[1])
