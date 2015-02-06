import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

plt.ioff()

def GradPak_patches():

    patch_list = [
        Circle((-48.1542, 85.5171), radius=0.9374),
        Circle((48.1542, 85.5171), radius=0.9374),
        Circle((-15.8139, 0.0000), radius=0.9374),
        Circle((-13.5548, 0.0000), radius=0.9374),
        Circle((-11.2957, 0.0000), radius=0.9374),
        Circle((-9.0365, 0.0000), radius=0.9374),
        Circle((-6.7774, 0.0000), radius=0.9374),
        Circle((-4.5183, 0.0000), radius=0.9374),
        Circle((-2.2591, 0.0000), radius=0.9374),
        Circle((-0.0000, 0.0000), radius=0.9374),
        Circle((2.2591, 0.0000), radius=0.9374),
        Circle((4.5183, 0.0000), radius=0.9374),
        Circle((6.7774, 0.0000), radius=0.9374),
        Circle((9.0365, 0.0000), radius=0.9374),
        Circle((11.2957, 0.0000), radius=0.9374),
        Circle((13.5548, 0.0000), radius=0.9374),
        Circle((15.8139, 0.0000), radius=0.9374),
        Circle((42.3236, 85.5171), radius=0.9374),
        Circle((-42.3236, 85.5171), radius=0.9374),
        Circle((-48.9201, 90.8997), radius=1.4061),
        Circle((-15.7343, 3.1159), radius=1.4061),
        Circle((-12.2378, 3.1159), radius=1.4061),
        Circle((-8.7413, 3.1159), radius=1.4061),
        Circle((-5.2448, 3.1159), radius=1.4061),
        Circle((-1.7483, 3.1159), radius=1.4061),
        Circle((1.7483, 3.1159), radius=1.4061),
        Circle((5.2448, 3.1159), radius=1.4061),
        Circle((8.7413, 3.1159), radius=1.4061),
        Circle((12.2378, 3.1159), radius=1.4061),
        Circle((15.7343, 3.1159), radius=1.4061),
        Circle((48.9201, 90.8997), radius=1.4061),
        Circle((-41.5578, 90.8997), radius=1.4061),
        Circle((-15.7343, 6.6124), radius=1.4061),
        Circle((-12.2378, 6.6124), radius=1.4061),
        Circle((-8.7413, 6.6124), radius=1.4061),
        Circle((-5.2448, 6.6124), radius=1.4061),
        Circle((-1.7483, 6.6124), radius=1.4061),
        Circle((1.7483, 6.6124), radius=1.4061),
        Circle((5.2448, 6.6124), radius=1.4061),
        Circle((8.7413, 6.6124), radius=1.4061),
        Circle((12.2378, 6.6124), radius=1.4061),
        Circle((15.7343, 6.6124), radius=1.4061),
        Circle((41.5578, 90.8997), radius=1.4061),
        Circle((-49.7300, 97.7793), radius=1.8748),
        Circle((-15.9124, 10.8720), radius=1.8748),
        Circle((-11.3660, 10.8720), radius=1.8748),
        Circle((-6.8196, 10.8720), radius=1.8748),
        Circle((-2.2732, 10.8720), radius=1.8748),
        Circle((2.2732, 10.8720), radius=1.8748),
        Circle((6.8196, 10.8720), radius=1.8748),
        Circle((11.3660, 10.8720), radius=1.8748),
        Circle((15.9124, 10.8720), radius=1.8748),
        Circle((49.7300, 97.7793), radius=1.8748),
        Circle((-40.7478, 97.7793), radius=1.8748),
        Circle((-15.9124, 15.4184), radius=1.8748),
        Circle((-11.3660, 15.4184), radius=1.8748),
        Circle((-6.8196, 15.4184), radius=1.8748),
        Circle((-2.2732, 15.4184), radius=1.8748),
        Circle((6.8196, 15.4184), radius=1.8748),
        Circle((11.3660, 15.4184), radius=1.8748),
        Circle((15.9124, 15.4184), radius=1.8748),
        Circle((40.7478, 97.7793), radius=1.8748),
        Circle((-45.2389, 82.7434), radius=2.3435),
        Circle((-16.6201, 20.6997), radius=2.3435),
        Circle((-11.0801, 20.6997), radius=2.3435),
        Circle((-5.5400, 20.6997), radius=2.3435),
        Circle((0.0000, 20.6997), radius=2.3435),
        Circle((5.5400, 20.6997), radius=2.3435),
        Circle((11.0801, 20.6997), radius=2.3435),
        Circle((16.6201, 20.6997), radius=2.3435),
        Circle((45.2389, 82.7434), radius=2.3435),
        Circle((-16.6201, 26.2397), radius=2.3435),
        Circle((-11.0801, 26.2397), radius=2.3435),
        Circle((-5.5400, 26.2397), radius=2.3435),
        Circle((0.0000, 26.2397), radius=2.3435),
        Circle((5.5400, 26.2397), radius=2.3435),
        Circle((11.0801, 26.2397), radius=2.3435),
        Circle((16.6201, 26.2397), radius=2.3435),
        Circle((-45.2389, 88.2909), radius=2.3435),
        Circle((-16.6201, 31.7797), radius=2.3435),
        Circle((-11.0801, 31.7797), radius=2.3435),
        Circle((-5.5400, 31.7797), radius=2.3435),
        Circle((0.0000, 31.7797), radius=2.3435),
        Circle((5.5400, 31.7797), radius=2.3435),
        Circle((11.0801, 31.7797), radius=2.3435),
        Circle((16.6201, 31.7797), radius=2.3435),
        Circle((45.2389, 88.2909), radius=2.3435),
        Circle((-45.2389, 94.4215), radius=2.8122),
        Circle((-16.7092, 38.1297), radius=2.8122),
        Circle((-10.0255, 38.1297), radius=2.8122),
        Circle((-3.3418, 38.1297), radius=2.8122),
        Circle((3.3418, 38.1297), radius=2.8122),
        Circle((10.0255, 38.1297), radius=2.8122),
        Circle((16.7092, 38.1297), radius=2.8122),
        Circle((45.2389, 94.4215), radius=2.8122),
        Circle((-16.7092, 44.8133), radius=2.8122),
        Circle((-10.0255, 44.8133), radius=2.8122),
        Circle((-3.3418, 44.8133), radius=2.8122),
        Circle((3.3418, 44.8133), radius=2.8122),
        Circle((10.0255, 44.8133), radius=2.8122),
        Circle((16.7092, 44.8133), radius=2.8122),
        Circle((-45.2389, 101.1361), radius=2.8122),
        Circle((-16.7092, 51.4970), radius=2.8122),
        Circle((-10.0255, 51.4970), radius=2.8122),
        Circle((16.7092, 51.4970), radius=2.8122),
        Circle((3.3418, 51.4970), radius=2.8122),
        Circle((10.0255, 51.4970), radius=2.8122),
        Circle((-3.3418, 51.4970), radius=2.8122),
        Circle((45.2389, 101.1361), radius=2.8122)]

    return np.array(patch_list)

def plot(values, clabel='', cmap='RdYlGn', nosky=True, labelfibers = True,
         exclude=[], minval=None, maxval=None):

    fig = plt.figure(figsize=(6,6))
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%')
    ax = grid[0]
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_xlim(-59,59)
    ax.set_ylim(-2,116)
    patches = GradPak_patches()

    skyidx = [0,1,17,18,19,30,31,42,43,52,53,61,62,70,78,86,87,94,101,108]

    if labelfibers:
        for i, c in enumerate(patches):
            if not nosky or (i not in skyidx and i+1 not in exclude) :
                ax.text(c.center[0],c.center[1],i+1,fontsize=6,
                        ha='center',va='center')

    if nosky:
        exclude = np.r_[skyidx,np.array(exclude)-1]
        ax.set_ylim(-2,60)
        ax.set_xlim(-31,31)
        
    exclude = np.array(exclude)
    exclude = np.unique(exclude)
    patches = np.delete(patches, exclude)
    values = np.delete(values, exclude)
    patches = patches[values == values]
    pval = values[values == values]
    
    if minval is None:
        minval = pval.min()
    if maxval is None:
        maxval = pval.max()

    collection = PatchCollection(patches,cmap=plt.get_cmap(cmap),
                                 norm=matplotlib.colors.Normalize(
            vmin=minval,vmax=maxval))
    collection.set_array(pval)
    ax.add_collection(collection)

    cbar = ax.cax.colorbar(collection)
    cbar.set_label_text(clabel)

    return ax

def plot_img(values, clabel='', cmap='jet', nosky=True,
             numpoints=500,method='nearest',exclude=[]):
    
    fibers = GradPak_patches()
    x = np.array([c.center[0] for c in fibers])
    y = np.array([c.center[1] for c in fibers])
    fig = plt.figure()
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%')
    ax = grid[0]
    ax.set_xlabel('arcsec')
    ax.set_ylabel('arcsec')
    ax.set_xlim(-59,59)
    ax.set_ylim(-2,116)

    if nosky:
        skyidx = [0,1,17,18,19,30,31,42,43,52,53,61,62,70,78,86,87,94,101,108]
        exclude = np.r_[skyidx,np.array(exclude)-1]
        ax.set_xlim(-21,21)
        ax.set_ylim(-2,40)

    exclude = np.array(exclude)
    exclude = np.unique(exclude)
    x = np.delete(x,exclude)
    y = np.delete(y,exclude)
    values = np.delete(values,exclude)

    xi = np.linspace(x.min(),x.max(),numpoints)
    yi = np.linspace(y.min(),y.max(),numpoints)

    x = x[values == values]
    y = y[values == values]
    values = values[values == values]

    vi = spi.griddata((x,y),
                      values,
                      (xi[None,:],yi[:,None]),
                      method=method,
                      fill_value=np.inf)
    
    im = ax.imshow(vi, cmap=cmap, origin='lower', 
                   extent=(xi.min(),xi.max(),yi.min(),yi.max()))
    cbar = ax.cax.colorbar(im)
    cbar.set_label_text(clabel)

    return ax

def format_tpl(tpl):

    with open(tpl,'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line[0:6] == 'circle':
            data = line.split('(')[1].split(')')[0].replace('"','')
            x,y,r = [float(i) for i in data.split(',')]
            print 'Circle(({:6.4f}, {:6.4f}), radius={:6.4f}),'.format(
                x*3600.,y*3600,r)

    return
