import numpy as np
import matplotlib.pyplot as plt
import ADEUtils as ADE
from datetime import datetime

def plot_AR(AR19_file,AR127_file,non19_file,non127_file,err=False):

    
    try:
        fibers19AR, tput19AR, _ = open_summary(AR19_file)
        fibers127AR, tput127AR, _ = open_summary(AR127_file)
        fibers19non, tput19non, std19non = open_summary(non19_file)
        fibers127non, tput127non, std127non  = open_summary(non127_file)
    except ValueError:
        fibers19AR, tput19AR = open_summary(AR19_file)
        fibers127AR, tput127AR = open_summary(AR127_file)
        fibers19non, tput19non = open_summary(non19_file)
        fibers127non, tput127non  = open_summary(non127_file)

    AR19sidx = np.argsort(fibers19AR)
    non19sidx = np.argsort(fibers19non)
    print tput19non[non19sidx]

    sidx127 = np.array([],dtype=np.int)

    for fiberpos in fibers127AR:

        idx = np.where(fibers127non == fiberpos)[0][0]
        sidx127 = np.append(sidx127,idx)

    oneline = np.linspace(np.min(np.array([tput127non.min(),
                                           tput127AR.min(),
                                           tput19non.min(),
                                           tput19AR.min()])),
                          np.max(np.array([tput19AR.max(),
                                           tput127AR.max(),
                                           tput19non.max(),
                                           tput127non.max()]).flatten()))

    # m19, b19 = ADE.fit_line(tput19non[non19sidx],
    #                         tput19AR[AR19sidx],np.ones(tput19AR.size))
    # m127, b127 = ADE.fit_line(tput127non[sidx127],
    #                           tput127AR,np.ones(tput127AR.size))
    # line19 = tput19non[non19sidx]*m19 + b19
    # line127 = tput127non[sidx127]*m127 + b127
    mean19 = np.mean(tput19AR[AR19sidx] - tput19non[non19sidx])
    mean127 = np.mean(tput127AR - tput127non[sidx127])
    # line19 = oneline + mean19
    # line127 = oneline + mean127

    R = 0.035
    Rprime = 0.0125

    tprime = oneline*((1 - Rprime)**2/(1 - R)**2)

    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if err:
        print std19non[non19sidx]
        ax.errorbar(tput19non[non19sidx],tput19AR[AR19sidx],ls='.',color='b',xerr=std19non[non19sidx])
        ax.errorbar(tput127non[sidx127],tput127AR,ls='.',color='g',xerr=std127non[sidx127])
    else:
        ax.plot(tput19non[non19sidx],tput19AR[AR19sidx],'b.',label='19 fiber intermediate')
        ax.plot(tput127non[sidx127],tput127AR,'g.',label='127 fiber intermediate')

    # ax.plot(oneline,line19,'b:',label='19 fiber intermediate')
    # ax.plot(oneline,line127,'g:',label='127 fiber intermediate')
    # ax.plot(tput19non[non19sidx],line19,'b:')
    # ax.plot(tput127non[sidx127],line127,'g:')
    ax.plot(oneline,oneline,'k:',label="$T$ = $T^{'}$")
    ax.plot(oneline,tprime,'k--',label="$T = T^{'} \\frac{(1-0.0125)^2}{(1-0.035)^2}$")
    print np.mean(tput19AR)
    # ax.axhline(y=np.mean(tput19AR),ls=':',color='b',label='19 fiber intermediate')
    # ax.axhline(y=np.mean(tput127AR),ls=':',color='g',label='127 fiber intermediate')
    ax.set_ylim(oneline.min(),0.98)
    ax.set_xlim(oneline.min(),oneline.max())
    ax.set_xlabel('Non-coated throughput ($T$)')
    ax.set_ylabel("AR coated throughput ($T^{'}$)")
    ax.legend(loc=4,prop={'size':11})
    ax.text(0.913,0.905,
            "mean increase for 19int = {:0.3f}\nmean increase for 127int = {:0.3f}".format(mean19,mean127),fontsize=11)
    fig.suptitle('A. Eigenbrot\n'+datetime.now().isoformat(' '))

    fig.show()

    print ("{} -- {}, {} -> {}\n"*fibers127AR.size).format(
        *sum([[t[0],t[1],t[2],t[3]] for t in zip(
                    fibers127AR,fibers127non[sidx127],tput127non[sidx127],tput127AR)],[]))
    print
    print ("{} -- {}, {} -> {}\n"*fibers19AR.size).format(
        *sum([[t[0],t[1],t[2],t[3]] for t in zip(
                    fibers19AR[AR19sidx],fibers19non[non19sidx],tput19non[non19sidx],tput19AR[AR19sidx])],[]))

    return 

def open_summary(summary_file):

    fibers = np.loadtxt(summary_file,usecols=(0,),unpack=True,dtype=np.str)
    try:
        tputs = np.loadtxt(summary_file,usecols=(5,),unpack=True)
        return fibers, tputs
    except IndexError:
        tputs, errs = np.loadtxt(summary_file,usecols=(1,2),unpack=True)
        return fibers, tputs, errs

def merge_summary(summary_tup,outfile):

    bigtput = []

    for summary in summary_tup:
        
        fibers = np.loadtxt(summary,usecols=(0,),unpack=True,dtype=np.str)
        try: 
            tputs = np.loadtxt(summary,usecols=(5,),unpack=True)
        except IndexError:
            tputs = np.loadtxt(summary,usecols=(3,),unpack=True)
        sortidx = np.argsort(fibers)
        bigtput.append(tputs[sortidx])
        
    merged_tput = np.mean(np.vstack(bigtput),axis=0)
    std_tput = np.std(np.vstack(bigtput),axis=0)
    
    f = open(outfile,'w')
    f.write('#{:>10}{:>9}{:>9}\n#\n'.format('Fiber_pos','sloan','std'))
    for k in range(merged_tput.size):
        f.write('{:>11}{:9.4f}{:9.4f}\n'.format(fibers[sortidx][k],merged_tput[k],std_tput[k]))

    f.close()

    return
    
