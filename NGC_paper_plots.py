import pyfits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
plt.ioff()

def plot_psi_weights(output,
                     modelfile='/d/monk/eigenbrot/WIYN/14B-0456/anal/models/allZ2_vardisp/allz2_vardisp_batch_interp.fits'):
    #Like the last page of all the fit plots, but for all pointings at once
    #cribbed from plot_bc_vardisp.py

    m = pyfits.open(modelfile)[1].data[0]
    numZ = np.unique(m['Z'][:,0]).size
    numAge = np.unique(m['AGE'][:,0]).size
    big_W = np.zeros((numZ,numAge))
    
    for p in range(6):
        coeffile = 'NGC_891_P{}_bin30_allz2.coef.fits'.format(p+1)
        print coeffile
        coef_arr = pyfits.open(coeffile)[1].data
        numap = coef_arr['VSYS'].size
        
        for i in range(numap):
            wdata = coef_arr[i]['LIGHT_FRAC'].reshape(numZ,numAge)
            big_W += wdata/np.max(wdata)

    bwax = plt.figure().add_subplot(111)
    bwax.imshow(big_W,origin='lower',cmap='Blues',interpolation='none')
    bwax.set_xlabel('SSP Age [Gyr]')
    bwax.set_xticks(range(numAge))
    bwax.set_xticklabels(m['AGE'][:numAge,0]/1e9)
    bwax.set_ylabel(r'$Z/Z_{\odot}$')
    bwax.set_yticks(range(numZ))
    bwax.set_yticklabels(m['Z'][::numAge,0])

    pp = PDF(output)
    pp.savefig(bwax.figure)
    pp.close()
    plt.close(bwax.figure)
    
    return
