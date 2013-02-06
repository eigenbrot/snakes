import numpy as np
import matplotlib.pyplot as plt
import ADEUtils as ADE
import pyfits
from glob import glob
from MANGA_bench import thePot

#
# Originally lived in MANGA/20130122

def ratio_test(t_file,data_path):

    T = thePot(t_file)
    fits_list = glob(data_path+'/*.FIT')
    volts = np.array([])
    counts = np.array([])

    for fits in fits_list:
        print fits,
        HDU = pyfits.open(fits)[0]
        exptime = HDU.header['EXPTIME']
        timestr = HDU.header['TIME-OBS']
        obstime = np.float(timestr[6:])\
            + np.float(timestr[3:5])*60.\
            + np.float(timestr[0:2])*3600.
        
        voltage = T.get_voltage(obstime,obstime+exptime,'V')
        
        r, sb, e = ADE.fast_annulize(HDU.data,300)
        r *= 0.024
        flux = np.cumsum(sb)
        rate = flux/exptime
        ADU = np.interp(54/10.,r,rate)

        print voltage, ADU

        volts = np.append(volts,voltage)
        counts = np.append(counts,ADU)

    
    fit = ADE.fit_line(volts,counts,np.ones(counts.shape))
    print fit
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(volts,counts,'.')
    ax.set_xlabel('Voltage')
    ax.set_ylabel('ADU/s')
    fig.show()
