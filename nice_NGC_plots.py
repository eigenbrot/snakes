import numpy as np
import matplotlib.pyplot as plt
import pyfits
from yanny import yanny
from glob import glob
import time
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib import rc

def make_nice():
    rc('text', usetex=False)
    rc('font', family='serif')
    rc('font', size=15.0)
    rc('axes', linewidth=0.9)
    rc('lines', linewidth=1.1)
    rc('patch', linewidth=0.1)
    rc('ps', usedistiller='Xpdf')
    rc('xtick', labelsize=12.0)
    rc('ytick', labelsize=12.0)

def height_summary(output, field, label='', title='', log=False):
    # Designed to work with drpall files, not fits

    meanbase = output.split('.pdf')[0].split('_height')[0]
    
    drpdir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/drpall'
    
    ax = plt.figure().add_subplot(111)

    symblist = ["o","^","v","s","*","x"]
    colorlist = ['b','c','g','y','m','r']
    plist = [6,3,4,2,1,5]
    rlist = [-10.0, -6.4, -2.3, 0.8, 4.4, 8.1]
    
    binidxs = [[0,1,2,3,4,5],
               [2,3],
               [1,4],
               [0,5],
               [0,1,4,5]]

    bigd = np.zeros((6,3,11))
    for i, p in enumerate(plist):
        
        parfile = glob('{}/NGC_891_P{}.drp.par'.format(drpdir,p))[0]
        print parfile
        par = yanny(parfile,np=True,debug=False)
        
        zs = par['APINFO']['z']
        zs= np.abs(zs)
        try:
            value = par['APINFO'][field]
        except ValueError:
            print 'WARNING: Could not find {} field in {}'.format(field,parfile)
            return

        if log:
            value = np.log10(value)

        unqz = np.unique(zs) #also sorts. Nice!
        abcissa = np.zeros(unqz.size)
        ordinate = np.zeros(unqz.size)
        err = np.zeros(unqz.size)
        for j, z in enumerate(unqz):
            idx = np.where(zs == unqz[j])[0]
            abcissa[j] = unqz[j]
            ordinate[j] = np.mean(value[idx])
            #err[j] = ordinate[j]*np.sqrt(1./idx.size)
            err[j] = np.std(value[zs == unqz[j]])/np.sqrt(idx.size)

        ax.errorbar(abcissa,ordinate,yerr=err,
                    color=colorlist[i],marker=symblist[i],ls='',
                    elinewidth=0.4,
                    label=str(rlist[i]))

        bigd[i,0,:] = abcissa
        bigd[i,1,:] = ordinate
        bigd[i,2,:] = err

    # num, z_edge = np.histogram(bigd[:,0,:].flatten(),11,normed=False)
    # total, z_edge2 = np.histogram(bigd[:,0,:].flatten(),11,normed=False,
    #                               weights=bigd[:,1,:].flatten())
    # print z_edge - z_edge2
    # bigz = (z_edge[1:] + z_edge[:-1])*0.5

    bigz2 = np.linspace(np.min(bigd[:,0,:]),np.max(bigd[:,0,:]),13)
    bigI = np.zeros((6,2,13))
    for k in range(len(plist)):
        print np.all(np.diff(bigd[k,0,:]) > 0)
        bigI[k,0,:] = np.interp(bigz2,bigd[k,0,:],bigd[k,1,:])
        bigI[k,1,:] = np.interp(bigz2,bigd[k,0,:],bigd[k,2,:])
        
    bigval = np.mean(bigI[:,0,:],axis=0)
    bigerr = np.sqrt(
        np.nansum(bigI[:,1,:]*(bigI[:,0,:] - bigval)**2,axis=0)/
        ((5/6.)*np.nansum(bigI[:,1,:],axis=0)))
    
    for l, b in enumerate(binidxs):
        binmean = np.mean(bigI[b,0,:],axis=0)
        binerr = np.sqrt(
            np.nansum(bigI[b,1,:]*(bigI[b,0,:] - binmean)**2,axis=0)/
            ((len(b)-1.0)/len(b) * np.nansum(bigI[b,1,:],axis=0)))
        with open('{}_bin{}_means.dat'.format(meanbase,l),'w') as f:
            f.write('# Bin {}: r = {} kpc\n'.\
                    format(l,np.array(rlist)[None,b].flatten()))
            f.write(str('#{:>12}'+'{:>13}'*2+'\n').\
                    format('height',
                           field,
                           '{} err'.format(field)))
            for Z, V, E in zip(bigz2,binmean,binerr):
                f.write(str('{:13.4f}'*3+'\n').format(Z,V,E))
    

    # ax.plot(bigz,total/num,'k-')
    ax.plot(bigz2,bigval,'b-')
    ax.fill_between(bigz2,bigval-bigerr,bigval+bigerr,alpha=0.2)
    ax.set_xlabel('z [kpc]')
    ax.set_ylabel(label)
    ax.set_title(title)
    ax.legend(loc=0,numpoints=1,title='radius [kpc]',fontsize=10)

    pp = PDF(output)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)

    return

