import numpy as np
import pyfits
import matplotlib.pyplot as plt
plt.ioff()

def go_ma_bc_comp():

    bcfile = 'models/DFK_1Z_vardisp.fits'
    mafile = 'MA11/MA11_models/ma11_dfk_1Z_vardisp.fits'

    bcdata = pyfits.open(bcfile)[1].data[0]
    madata = pyfits.open(mafile)[1].data[0]

    wave = bcdata['WAVE']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Angstroms')
    ax.set_ylabel('Normalized flux')

    clist = ['blue','green','orange','red']
    dfklist = ['Y','I1','I2','O']
    for i in range(4):
        ax.plot(wave, bcdata['FLUX'][i,:,2],color=clist[i],alpha=0.8)
        ax.text(wave[-1] + 100,bcdata['FLUX'][i,-10,2],
                dfklist[i],color=clist[i],fontsize=9,va='center',ha='left')
        if i > 0:
            ax.plot(wave, madata['FLUX'][i-1,:,2],':',color=clist[i])

    ax.text(0.8,0.95,'bc03',fontsize=9,transform=ax.transAxes,va='center')
    ax.text(0.8,0.9,'ma11',fontsize=9,transform=ax.transAxes,va='center')
    ax.axhline(ax.get_ylim()[1]*0.95,xmin=0.7,xmax=0.76,color='k',ls='-')
    ax.axhline(ax.get_ylim()[1]*0.9,xmin=0.7,xmax=0.76,color='k',ls=':')

    fig.savefig('ma_bc_comp.png',dpi=200)

    return

def bc_downselect():

    bcfile = 'models/DFK_1Z_vardisp.fits'
    mabcfile = 'MABC/BC_MA_models/mabc_dfk_1Z_vardisp.fits'

    bcdata = pyfits.open(bcfile)[1].data[0]
    mabcdata = pyfits.open(mabcfile)[1].data[0]

    wave = bcdata['WAVE']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Angstroms')
    ax.set_ylabel('Normalized flux')

    clist = ['blue','green','orange','red']
    dfklist = ['Y','I1','I2','O']
    for i in range(4):
        ax.plot(wave, bcdata['FLUX'][i,:,2],
                color=clist[i],alpha=0.6)
        ax.text(wave[-1] + 100,bcdata['FLUX'][i,-10,2],
                dfklist[i],color=clist[i],fontsize=9,va='center',ha='left')
        if i > 0:
            ax.plot(wave, mabcdata['FLUX'][i-1,:,2],color=clist[i],ls='--')

    ax.text(0.8,0.95,'bc03 full set',fontsize=9,
            transform=ax.transAxes,va='center')
    ax.text(0.8,0.9,'bc03 downselected',fontsize=9,
            transform=ax.transAxes,va='center')
    ax.axhline(ax.get_ylim()[1]*0.95,xmin=0.7,xmax=0.76,color='k',
               alpha=0.7)
    ax.axhline(ax.get_ylim()[1]*0.9,xmin=0.7,xmax=0.76,color='k',ls='--')

    fig.savefig('bc_downselect.png',dpi=200)

    return

def bc_dfk_full():

    bcfile = 'models/DFK_1Z_vardisp.fits'
    allfile = '../models/allZ2_vardisp/allz2_vardisp_batch_interp.fits'

    bcdata = pyfits.open(bcfile)[1].data[0]
    alldata = pyfits.open(allfile)[1].data[0]

    wave = bcdata['WAVE']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Angstroms')
    ax.set_ylabel('Normalized flux')

    clist = ['blue','green','orange','red']
    dfklist = ['Y','I1','I2','O']

    #40 - 49 are the 1Z SSPs
    for i in range(30,40):
        age = alldata['AGE'][i,2]
        if age < 5.2e6:
            color = clist[0]
        elif 5.5e6 < age < 404e6:
            color = clist[1]
        elif 450e6 < age < 5750e6:
            color = clist[2]
        else:
            color = clist[3]
        ax.plot(wave, alldata['FLUX'][i,:,2],color=color, alpha=0.3)

    for i in range(4):
        ax.plot(wave, bcdata['FLUX'][i,:,2],color=clist[i],alpha=0.8)
        ax.text(wave[-1] + 100,bcdata['FLUX'][i,-10,2],
                dfklist[i],color=clist[i],fontsize=9,va='center',ha='left')


    # ax.text(0.8,0.95,'bc03',fontsize=9,transform=ax.transAxes,va='center')
    # ax.text(0.8,0.9,'ma11',fontsize=9,transform=ax.transAxes,va='center')
    # ax.axhline(ax.get_ylim()[1]*0.95,xmin=0.7,xmax=0.76,color='k',ls='-')
    # ax.axhline(ax.get_ylim()[1]*0.9,xmin=0.7,xmax=0.76,color='k',ls=':')
    fig.savefig('bc_dfk_comp.png',dpi=200)

    return 

def go():

    bcfile = 'models/DFK_1Z_vardisp.fits'
    cbfile = 'CB08/CB08_models/DFK_cb08_04Z_vardisp.fits'
    mafile = 'MA11/MA11_models/ma11_dfk_1Z_vardisp.fits'
    mabcfile = 'MABC/BC_MA_models/mabc_dfk_04Z_vardisp.fits'
    fspsfile = 'FSPS/FSPS_models/fsps_dfk_04Z_vardisp.fits'
    allfile = '../models/allZ2_vardisp/allz2_vardisp_batch_interp.fits'

    bcdata = pyfits.open(bcfile)[1].data[0]
    cbdata = pyfits.open(cbfile)[1].data[0]
    madata = pyfits.open(mafile)[1].data[0]
    mabcdata = pyfits.open(mabcfile)[1].data[0]
    fspsdata = pyfits.open(fspsfile)[1].data[0]
    alldata = pyfits.open(allfile)[1].data[0]

    print np.mean(bcdata['WAVE'] - cbdata['WAVE'])
    wave = bcdata['WAVE']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Angstroms')
    ax.set_ylabel('Normalized flux')

    clist = ['blue','green','orange','red']

    for i in range(4):
        ax.plot(wave, bcdata['FLUX'][i,:,2],color=clist[i],alpha=0.8)
        if i > 0:
            ax.plot(wave, madata['FLUX'][i-1,:,2],':',color=clist[i])

    #30 - 40 are the 0.4Z SSPs
    # for i in range(30,40):
    #     age = alldata['AGE'][i,2]
    #     if age < 5.2e6:
    #         color = clist[0]
    #     elif 5.5e6 < age < 404e6:
    #         color = clist[1]
    #     elif 450e6 < age < 5750e6:
    #         color = clist[2]
    #     else:
    #         color = clist[3]
    #     ax.plot(wave, alldata['FLUX'][i,:,2],color=color, alpha=0.3)

    ax.text(0.8,0.95,'bc03',fontsize=9,transform=ax.transAxes,va='center')
    ax.text(0.8,0.9,'ma11',fontsize=9,transform=ax.transAxes,va='center')
    ax.axhline(ax.get_ylim()[1]*0.95,xmin=0.7,xmax=0.76,color='k',ls='-')
    ax.axhline(ax.get_ylim()[1]*0.9,xmin=0.7,xmax=0.76,color='k',ls=':')

    return fig
