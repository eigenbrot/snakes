import numpy as np
import pyfits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid as AG
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages as PDF
from datetime import datetime

# Originally lived in /d/monk/eigenbrot/MANGA/lifetime
#

def plot_lifetime(output,title,origfits,lifefits,arfits=False):

    pp = PDF(output)

    rc('text', usetex=False)
    rc('font', family='serif')
    rc('font', size=9.0)
    rc('axes', linewidth=0.4)

    orig_dict = sort_data(origfits)
    life_dict = sort_data(lifefits)
    AR = False
    if arfits:
        ar_dict = sort_data(arfits)
        AR = True

    hi = 0
    fig = plt.figure()
    grid = AG(fig,111,
              nrows_ncols = (4,5),
              axes_pad = 0.0,
              label_mode = 'L',
              aspect = False,
              share_all = False)

    for (ax,fiber) in zip(grid,life_dict.keys()):
        print fiber

        orig_data = orig_dict[fiber].data
        orig_sloan = orig_dict[fiber].header['SLOAN']
        life_data = life_dict[fiber].data
        life_sloan = life_dict[fiber].header['SLOAN']
        if AR:
            ar_data = ar_dict[fiber].data
            ar_sloan = ar_dict[fiber].header['SLOAN']

        ax.plot(life_data[1],life_data[4],'r',lw=0.4)
        ax.plot(life_data[3],life_data[5],'b',lw=0.4)
        ax.plot(orig_data[3],orig_data[5],'g',lw=0.4)
        if AR:
            ax.plot(ar_data[3],ar_data[5],'m',lw=0.4)
        ax.text(9,1.0,'({})'.format(fiber),va='top')
        ax.text(9,0.8,'{:4.3f}'.format(orig_sloan),color='g',va='top')
        if AR:
            ax.text(9,0.65,'{:4.3f}'.format(ar_sloan),color='m',va='top')
        ax.text(9,0.5,'{:4.3f}'.format(life_sloan),color='b',va='top')

        ax.set_xlim(2,30)
        ax.set_ylim(0,1.1)
        ax.set_xscale('log')
        ax.set_xticks([5,10,20])
        ax.set_xticklabels(['$5$','$10$','$20$'])

    grid[-1].set_xlim(2,30)
    grid[-1].set_ylim(0,1.1)
    grid[-1].set_xscale('log')
    grid[-1].set_xticks([])
    grid[-1].text(2.5,0.7,'(fiber input position)\nthroughput value',
                  fontsize=8)
    grid[-1].hlines([0.5,0.4,0.3,0.2],3,5,colors=['g','m','b','r'])
#                        linestyles=['solid','dashdot','solid'])
    grid[-1].text(6,0.5,'Initial measurement',fontsize=5.0)
    grid[-1].text(6,0.4,'After AR coating',fontsize=5.0)
    grid[-1].text(6,0.3,'After lifetime tests',fontsize=5.0)
    grid[-1].text(6,0.2,'Direct normalized',fontsize=5.0)
    
    fig.text(0.5,0.04,'$f$-ratio',fontsize=12.0)
    fig.text(0.06,0.65,'Normalized Encircled Energy',
             rotation='vertical',fontsize=11.0)
    fig.suptitle(title+'\n'+datetime.now().isoformat(' '))
    
    pp.savefig(fig)
    pp.close()

    return

def sort_data(fits_file):

    hdus = pyfits.open(fits_file)[1:]

    out_dict = {}

    for h in hdus:
        fiber_pos = h.header['FIBERPOS']
        out_dict[fiber_pos] = h

    return out_dict

def merge_data(fits_files,output):

    final_hdus = []

    for fits in fits_files:
        hdus = pyfits.open(fits)[1:]
        final_hdus += hdus

    pyfits.HDUList([pyfits.PrimaryHDU(None)]+final_hdus).writeto(output)

    return
