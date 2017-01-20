#!/usr/bin/python
import sys
import os
import time
import numpy as np
import ir2py_lines as i2p
import matplotlib.pyplot as plt
from pyraf import iraf

iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)
plt.ioff()

llist = ['4047','4358','HK','5461','5577','NaD','OI','OH1','OH2','OH3','OH4']
rlist = ['4020 4075',
         '4325 4400',
         '3910 4010',
         '5437 5480',
         '5555 5597',
         '5875 5910',
         '6264 6392',
         '6440 6520',
         '6585 6620',
         '6850 6880',
         '7030 7065']
#centlist = [[4047],[4358],[3933.57,3968.53],[5461],[5577.5],[5892.9565],[6300.3, 6363.8]]
centlist = [[4046.563],
            [4358.328],
            [3933.57,3968.53],
            [5460.735],
            [5577.5],
            [5892.9565],
            [6300.3, 6363.8],
            [6465.901, 6498.729],
            [6604.135],
            [6863.955],
            [7047.846]]
numlist = [1,1,2,1,1,1,2,2,1,1,1]

def do_fitprof(datafile):
    """Run IRAF fitprof routine to measure line centers and widths.

    For each line specified in the module header generate the necessary inputs
    to fitprof and execture the routine. The fitting regions are also
    specified in the module header.

    Parameters
    ----------
    datafile : str
        Name of a multispec fits file to pass as input to fitprof. Must have a WCS solution in the header.

    Returns
    -------
    None
        Nothing is returned. Instead, IRAF writes the results to a file.

    """
    
    for l, cl, r in zip(llist,centlist,rlist):
        with open('{}.lines'.format(l),'w') as f:
            print l
            for c in cl:
                f.write('{} INDEF g\n'.format(c))

        iraf.fitprofs(datafile,
                      region=r,
                      fitgfwhm='single',
                      nerrsample=100,
                      sigma0=0.3,
                      invgain=0,
                      positio='{}.lines'.format(l),
                      logfile='{}.fitp'.format(l))

    return

def get_results(output, threshold=3., filename=''):
    """Parse fitprof output and display results

    The line centers are taken from the output of fitprof. For each line specified in the module header the average offset and stddev across all fibers in the IFU is computed. Output is a textfile and a plot of accuracy and stochasticity as a function of wavelength.

    Parameters
    ----------
    output : str
        Name of the output text file. This file will contain the mean, offset, stddev, and number of rejected apertures.
    threshold : float, optional
        Threshold value for iterative sigma clipping in mean across IFU. The total number of rejected fibers will be recorded in the output file.

    Returns
    -------
    None :
       The result is a text file containing the results and a plot containing the accuracy as a function of wavelength.
    
    Notes
    -----
    Right now the plot is hardcoded to be writting to WLC.png

    """

    fig = plt.figure()
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('Mismatch [AA]')
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Wavelength')
    stax.set_ylabel('IFU std')

    with open(output,'a') as f:
        f.write('# {}\n'.format(time.asctime()))
        f.write('# {}\n'.format(filename))

        for l, n, c in zip(llist,numlist,centlist):
            proffile = '{}.fitp'.format(l)
            d = i2p.parse_fitprofs(proffile,n)[1]
            mean = np.mean(d,axis=0)
            std = np.std(d,axis=0)
            rejidx = np.where(np.abs(d - mean) > std*threshold)
            i = 0
            while rejidx[0].size != 0:
                d = np.delete(d,rejidx[0],axis=0)
                mean = np.mean(d,axis=0)
                std = np.std(d,axis=0)
                rejidx = np.where(np.abs(d - mean) > std*threshold)
                i += rejidx[0].size
                
            diff = mean - c
            outstr = '{:} ({:}):\n\t{:>7}: {:}\n\t{:>7}: {:}\n\t{:>7}: {:}\n\t{:>7}: {:}\n'.\
                     format(l,c,'numrej',i,'mean',mean,'diff',diff,'std',std)
            prtstr = ''
            print c, diff, std
            for j in range(len(c)):
                prtstr += '{} {} {}\n'.format(c[j],diff[j],std[j])

            print prtstr
            f.write(outstr)
            acax.plot(c,diff,'.k')
            stax.plot(c,std,'.k')
    
        f.write('\n\n')

    fig.subplots_adjust(hspace=0.0001)
    acax.set_xlim(*stax.get_xlim())
    plt.savefig('WLC.png')

    return
def fiber_size_seg(inputdir='.'):
    """Parse fitprof output and display results

    The line centers are taken from the output of fitprof. For each line specified in the module header the average offset and stddev across all fibers in the IFU is computed. Output is a textfile and a plot of accuracy and stochasticity as a function of wavelength.

    Parameters
    ----------
    output : str
        Name of the output text file. This file will contain the mean, offset, stddev, and number of rejected apertures.
    threshold : float, optional
        Threshold value for iterative sigma clipping in mean across IFU. The total number of rejected fibers will be recorded in the output file.

    Returns
    -------
    None :
       The result is a text file containing the results and a plot containing the accuracy as a function of wavelength.
    
    Notes
    -----
    Right now the plot is hardcoded to be writting to WLC.png

    """

    fib2_idx = np.arange(1,20) - 1
    fib3_idx = np.arange(20,44) - 1
    fib4_idx = np.arange(44,63) - 1
    fib5_idx = np.arange(63,88) - 1
    fib6_idx = np.arange(88,110) - 1

    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    labels = [200, 300, 400, 500, 600]

    results = np.zeros((5,np.hstack(centlist).size,2))
    
    fig = plt.figure()
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('Mismatch [AA]')
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Wavelength')
    stax.set_ylabel('IFU std')

    j = 0
    for l, n, c in zip(llist,numlist,centlist):
        proffile = '{}/{}.fitp'.format(inputdir, l)
        d = i2p.parse_fitprofs(proffile,n)[1]
        for i, idx in enumerate([fib2_idx, fib3_idx, fib4_idx, fib5_idx, fib6_idx]):
            mean = np.mean(d[idx],axis=0)
            std = np.std(d[idx],axis=0)
            
            
            diff = mean - c
            
            if len(c) > 1:
                for k in range(len(c)):
                    results[i,j+k,0] = diff[k]
                    results[i,j+k,1] = std[k]
                # results[i,j+1,0] = diff[1]
                # results[i,j+1,1] = std[1]
            else:
                results[i,j,0] = diff
                results[i,j,1] = std
                
                if c == centlist[0]:
                    acax.plot(c,diff,'.', color=colors[i], label=labels[i])
                else:
                    acax.plot(c,diff,'.', color=colors[i])
                    stax.plot(c,std,'.', color=colors[i])
                    
        if len(c) > 1:
            j += len(c)
        else:
            j += 1
                
    fig.subplots_adjust(hspace=0.0001)
    acax.set_xlim(*stax.get_xlim())
    acax.legend(loc=0, frameon=False, numpoints=1)
    plt.savefig('{}/WLC_fibsize.png'.format(inputdir))

    return results

def all_nights_fibsize(output = 'Wave_err_fibsize.pdf'):
    from matplotlib import rc
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    from glob import glob as glob

    gw = 0.5
    gsize = 5
    rc('path', simplify=True)
    rc('figure', figsize=(gsize*1.718,gsize))
    rc('font', family='serif', weight=100, size=14)
    rc('mathtext', default='regular')
    rc('xtick', labelsize=16)
    rc('ytick', labelsize=16)
    rc('xtick.major', size=8, width=gw)
    rc('ytick.major', size=8, width=gw)
    rc('xtick.minor', size=4, width=gw, visible=True)
    rc('ytick.minor', size=4, width=gw, visible=True)
    rc('lines', markeredgewidth=1)
    rc('figure.subplot', top=0.95, bottom=0.15)#, right=0.95, left=0.15)
    rc('legend', numpoints=1, scatterpoints=1, frameon=False, handletextpad=0.3)
    rc('axes', linewidth=gw, labelweight=100, labelsize=24)
    
    cents = [4047,4358,3933.57,3968.53,5461]
    night_list = glob('n*')
    
    big_results = []
    for night in night_list:
        inputdir = '{}/best_rdx'.format(night)
        print inputdir

        big_results.append(fiber_size_seg(inputdir))

    big_stack = np.stack(big_results, axis=0)
    print big_stack.shape #Should be (numnights, numsize, numline, 2)

    allnight_d = np.mean(big_stack, axis=0)

    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    labels = [200, 300, 400, 500, 600]

    fig = plt.figure(figsize=(gsize*1.718,gsize*1.3))
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('Mismatch [km/s]', fontsize=12)
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Wavelength')
    stax.set_ylabel('Fibse size stddev (across all nights)', fontsize=12)

    for i in range(allnight_d.shape[0]):
        acax.plot(cents,allnight_d[i,:,0]/cents*3e5, '.', color=colors[i])
        stax.plot(cents,allnight_d[i,:,1]/cents*3e5, '.', color=colors[i], label=labels[i])

    fig.subplots_adjust(hspace=0.0001)
    acax.axhline(0,ls='--',alpha=0.7,color='k')
    acax.set_ylim(-200,400)
    #stax.set_ylim(0,241)
    acax.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
    acax.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
    acax.set_xlim(*stax.get_xlim())
    stax.legend(loc=0, frameon=False, numpoints=1)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    
    return allnight_d

def each_night_fibsize(output = 'Wave_err_fibsize_nightly.pdf'):
    from matplotlib import rc
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    from glob import glob as glob

    gw = 0.5
    gsize = 5
    rc('path', simplify=True)
    rc('figure', figsize=(gsize*1.718,gsize))
    rc('font', family='serif', weight=100, size=14)
    rc('mathtext', default='regular')
    rc('xtick', labelsize=16)
    rc('ytick', labelsize=16)
    rc('xtick.major', size=8, width=gw)
    rc('ytick.major', size=8, width=gw)
    rc('xtick.minor', size=4, width=gw, visible=True)
    rc('ytick.minor', size=4, width=gw, visible=True)
    rc('lines', markeredgewidth=1)
    rc('figure.subplot', top=0.95, bottom=0.15)#, right=0.95, left=0.15)
    rc('legend', numpoints=1, scatterpoints=1, frameon=False, handletextpad=0.3)
    rc('axes', linewidth=gw, labelweight=100, labelsize=24)
    
    #cents = [4047,4358,3933.57,3968.53,5461]
    cents = [4047,4358,3933.57,3968.53,5461,5577.5,5892.9565,6300.3, 6363.8]
    night_list = glob('n*')
    
    big_results = []
    for night in night_list:
        inputdir = '{}/best_rdx/more_lines'.format(night)
        print inputdir

        big_results.append(fiber_size_seg(inputdir))

    big_stack = np.stack(big_results, axis=0)
    print big_stack.shape #Should be (numnights, numsize, numline, 2)

    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    labels = [200, 300, 400, 500, 600]

    fig = plt.figure(figsize=(gsize*1.718,gsize*1.3))
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('Mismatch [km/s]', fontsize=12)
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Wavelength')
    stax.set_ylabel('Fibse size stddev (across IFU on 1 night)', fontsize=12)

    for i in range(big_stack.shape[1]):
        for j in range(big_stack.shape[0]):
            acax.scatter(cents,big_stack[j,i,:,0]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, s=30)
            if j == 0:
                stax.scatter(cents,big_stack[j,i,:,1]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, label=labels[i], s=30)
            else:
                stax.scatter(cents,big_stack[j,i,:,1]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, s=30)

    fig.subplots_adjust(hspace=0.0001)
    acax.axhline(0,ls='--',alpha=0.7,color='k')
    acax.set_ylim(-200,400)
    stax.set_ylim(0,241)
    acax.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
    acax.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
    acax.set_xlim(*stax.get_xlim())
    stax.legend(loc=0, frameon=False, numpoints=1)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()
    
    return big_stack

def fib_difference(allnight_d, output):
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    
    cents = [4047,4358,3933.57,3968.53,5461]
    sizes = [200, 300, 400, 500, 600]
    
    fig = plt.figure(figsize=(5*1.718,5*1.3))
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('<Mismatch> [km/s]', fontsize=12)
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Fiber size [microns]')
    stax.set_ylabel('<Stddev>', fontsize=12)
    
    means = np.mean(allnight_d/np.array(cents)[None,:,None]*3e5, axis=1)
    
    acax.plot(sizes,means[:,0], '.k')
    stax.plot(sizes,means[:,1], '.k')

    fig.subplots_adjust(hspace=0.0001)
    acax.axhline(0,ls='--',alpha=0.7,color='k')
    acax.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
    acax.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
    acax.set_xticks(sizes)
    stax.set_xticks(sizes)
    stax.set_xlim(150,650)
    acax.set_xlim(*stax.get_xlim())
    stax.legend(loc=0, frameon=False, numpoints=1)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return

def each_fib_difference(big_stack, output):
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    
    cents = [4047,4358,3933.57,3968.53,5461]
    sizes = [200, 300, 400, 500, 600]
    
    fig = plt.figure(figsize=(5*1.718,5*1.3))
    acax = fig.add_subplot(211)
    acax.set_xticklabels([])
    acax.set_ylabel('<Mismatch> [km/s]', fontsize=12)
    
    stax = fig.add_subplot(212)
    stax.set_xlabel('Fiber size [microns]')
    stax.set_ylabel('<Stddev>', fontsize=12)
    
    means = np.mean(big_stack/np.array(cents)[None,None,:,None]*3e5, axis=2)
    
    for i in range(big_stack.shape[0]):
        acax.plot(sizes,means[i,:,0], '.k')
        stax.plot(sizes,means[i,:,1], '.k')

    fig.subplots_adjust(hspace=0.0001)
    acax.axhline(0,ls='--',alpha=0.7,color='k')
    acax.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
    acax.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
    acax.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
    acax.set_xticks(sizes)
    stax.set_xticks(sizes)
    stax.set_xlim(150,650)
    acax.set_xlim(*stax.get_xlim())
    stax.legend(loc=0, frameon=False, numpoints=1)

    pp = PDF(output)
    pp.savefig(fig)
    pp.close()

    return

def mab_panels(outpre = 'Wave_err_mab_panels'):
    from matplotlib import rc
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    from glob import glob as glob

    gw = 0.5
    gsize = 5
    rc('path', simplify=True)
    rc('figure', figsize=(gsize*1.718,gsize))
    rc('font', family='serif', weight=100, size=14)
    rc('mathtext', default='regular')
    rc('xtick', labelsize=16)
    rc('ytick', labelsize=16)
    rc('xtick.major', size=8, width=gw)
    rc('ytick.major', size=8, width=gw)
    rc('xtick.minor', size=4, width=gw, visible=True)
    rc('ytick.minor', size=4, width=gw, visible=True)
    rc('lines', markeredgewidth=1)
    rc('figure.subplot', top=0.95, bottom=0.15)#, right=0.95, left=0.15)
    rc('legend', numpoints=1, scatterpoints=1, frameon=False, handletextpad=0.3)
    rc('axes', linewidth=gw, labelweight=100, labelsize=24)
    
    #cents = [4047,4358,3933.57,3968.53,5461]
    #cents = [4047,4358,3933.57,3968.53,5461,5577.5,5892.9565,6300.3, 6363.8]
    cents = np.hstack(centlist)
    night_list = glob('n*')
    
    big_results = []
    for night in night_list:
        inputdir = '{}/best_rdx/more_lines'.format(night)
        print inputdir

        big_results.append(fiber_size_seg(inputdir))

    big_stack = np.stack(big_results, axis=0)
    print big_stack.shape #Should be (numnights, numsize, numline, 2)

    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    labels = [200, 300, 400, 500, 600]

    pp1 = PDF('{}_nightly.pdf'.format(outpre))
    pp2 = PDF('{}_all.pdf'.format(outpre))

    for i in range(big_stack.shape[1]):
        print i
        fig1 = plt.figure(figsize=(gsize*1.718,gsize*1.3))
        acax1 = fig1.add_subplot(211)
        acax1.set_xticklabels([])
        acax1.set_ylabel('Mismatch [km/s]', fontsize=12)
    
        stax1 = fig1.add_subplot(212)
        stax1.set_xlabel('Wavelength')
        stax1.set_ylabel('Fib size stddev (across IFU on 1 night)', fontsize=12)

        for j in range(big_stack.shape[0]):
            acax1.plot(cents,big_stack[j,i,:,0]/cents*3e5,'.')
            if j == 0:
                stax1.plot(cents,big_stack[j,i,:,1]/cents*3e5, '.',label=labels[i])
            else:
                stax1.plot(cents,big_stack[j,i,:,1]/cents*3e5, '.')
            # acax1.scatter(cents,big_stack[j,i,:,0]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, s=30)
            # if j == 0:
            #     stax1.scatter(cents,big_stack[j,i,:,1]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, label=labels[i], s=30)
            # else:
            #     stax1.scatter(cents,big_stack[j,i,:,1]/cents*3e5, c=colors[i], edgecolors='none', alpha=0.6, s=30)

        fig1.subplots_adjust(hspace=0.0001)
        acax1.axhline(0,ls='--',alpha=0.7,color='k')
        acax1.set_ylim(-200,400)
        stax1.set_ylim(0,241)
        acax1.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
        acax1.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
        acax1.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
        acax1.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
        acax1.set_xlim(*stax1.get_xlim())
        stax1.legend(loc=0, frameon=False, numpoints=1)

        pp1.savefig(fig1)
        plt.close(fig1)

        ###########

        fig2 = plt.figure(figsize=(5*1.718,5*1.3))
        acax2 = fig2.add_subplot(211)
        acax2.set_xticklabels([])
        acax2.set_ylabel('<Mismatch>_w [km/s]', fontsize=12)
        
        stax2 = fig2.add_subplot(212)
        stax2.set_xlabel('Wavelength')
        stax2.set_ylabel('<Fib size stddev>_w', fontsize=12)

        kms = big_stack/np.array(cents)[None,None,:,None]*3e5
        weight = 1/kms[:,:,:,1][:,:,:,None]
        avg = np.sum(kms*weight, axis=0)/np.sum(weight,axis=0)
        
        acax2.scatter(cents,avg[i,:,0], c=colors[i], s=30, alpha=0.6, edgecolors='none')
        stax2.scatter(cents,avg[i,:,1], c=colors[i], s=30, alpha=0.6, edgecolors='none', label=labels[i])
            
        fig2.subplots_adjust(hspace=0.0001)
        acax2.set_ylim(-200,400)
        acax2.axhline(0,ls='--',alpha=0.7,color='k')
        acax2.text(3951,320,'Ca H&K',fontsize=15,va='bottom',ha='center')
        acax2.text(4047+30,50,'HgI',fontsize=15,va='bottom',ha='center')
        acax2.text(4358,120,'HgI',fontsize=15,va='bottom',ha='center')
        acax2.text(5461,90,'HgI',fontsize=15,va='bottom',ha='center')
        acax2.set_xlim(*stax2.get_xlim())
        stax2.legend(loc=0, frameon=False, numpoints=1)

        pp2.savefig(fig2)
        plt.close(fig2)
        
    pp1.close()
    pp2.close()
    
    return big_stack

def night_plot(big_stack, output = 'night_by_night.pdf'):
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    gsize = 5
    pp = PDF(output)

    colors = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
    labels = [200, 300, 400, 500, 600]
    cents = [4047,4358,3933.57,3968.53,5461]
    
    for i in range(big_stack.shape[2]):
        fig = plt.figure(figsize=(gsize*1.718,gsize*1.3))
        acax = fig.add_subplot(211)
        acax.set_xticklabels([])
        acax.set_ylabel('Mismatch [km/s]', fontsize=12)
        acax.set_title(cents[i])
        
        stax = fig.add_subplot(212)
        stax.set_xlabel('Night')
        stax.set_ylabel('stddev', fontsize=12)
        
        kms = big_stack/np.array(cents)[None,None,:,None]*3e5
        
        for j in range(big_stack.shape[1]):
            acax.scatter(np.arange(big_stack.shape[0])+1, kms[:,j,i,0], c=colors[j], s=30, edgecolors='none')
            stax.scatter(np.arange(big_stack.shape[0])+1, kms[:,j,i,1], c=colors[j], s=30, edgecolors='none', label=labels[j])

        stax.legend(loc=0, frameon=False)
        fig.subplots_adjust(hspace=0.0001)
        acax.set_ylim(-200,400)
        stax.set_ylim(0,241)
        pp.savefig(fig)
        plt.close(fig)

    pp.close()

if __name__ == '__main__':
    
    for l in llist:
        proffile = '{}.fitp'.format(l)
        if os.path.exists(proffile):
            print "Removing", proffile
            os.system('rm '+ proffile)

    do_fitprof(sys.argv[1])
    get_results('WLC.dat',filename=sys.argv[1])

