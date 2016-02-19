import numpy as np
import plot_simple as ps
import plot_allZ2 as pa2
import sys

if len(sys.argv) == 3:
    exclude = []
    for i in range(6):
        d, e = np.loadtxt('P{}_GBU.dat'.format(i+1),unpack=True)
        idx = (np.where(e != 0)[0]+1).tolist()
        exclude.append(idx)
else:
    exclude = [[5,34],
               [1,2,35],
               [59],
               [2,8],
               [1,2,3,27,28,29],
               [35,36,38]]

print exclude
ps.all_maps('MLWA_map.pdf',col=6,inputsuffix='allz2.dat',label='Mean Light Weighted Age [Gyr]', minval=0, maxval=12.5,binned=True,exclude=exclude)
ps.all_maps('MMWA_map.pdf',col=5,inputsuffix='allz2.dat',label='Mean Mass Weighted Age [Gyr]', minval=0, maxval=12.5,binned=True,exclude=exclude)
ps.all_maps('MLWZ_map.pdf',col=8,inputsuffix='allz2.dat',label='Mean Light Weighted Metallicity [Z_sol]', minval=-1.5, maxval=0.5,binned=True,log=True,exclude=exclude)
ps.all_maps('V_map.pdf',col=9,inputsuffix='allz2.dat',label='Velocity [km/s]', minval=298, maxval=758,binned=True,cmap='bwr',exclude=exclude)
ps.all_maps('TauV_map.pdf',col=10,inputsuffix='allz2.dat',label=r'$\tau_V$', minval=0, maxval=5,binned=True,exclude=exclude)

pa2.dfk_batch(sys.argv[1])
