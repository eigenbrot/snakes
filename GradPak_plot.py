import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pywcs
import pywcsgrid2 as wcsgrid
import pyfits
import scipy.interpolate as spi
import scipy.ndimage.interpolation as spndi
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

plt.ioff()
tau = 2 * np.pi

def GradPak_patches():

    raw_patches = [
        Circle((48.1542, 85.5171), radius=0.9374),
        Circle((-48.1542, 85.5171), radius=0.9374),
        Circle((15.8139, 0.0000), radius=0.9374),
        Circle((13.5548, 0.0000), radius=0.9374),
        Circle((11.2957, 0.0000), radius=0.9374),
        Circle((9.0365, 0.0000), radius=0.9374),
        Circle((6.7774, 0.0000), radius=0.9374),
        Circle((4.5183, 0.0000), radius=0.9374),
        Circle((2.2591, 0.0000), radius=0.9374),
        Circle((0.0000, 0.0000), radius=0.9374),
        Circle((-2.2591, 0.0000), radius=0.9374),
        Circle((-4.5183, 0.0000), radius=0.9374),
        Circle((-6.7774, 0.0000), radius=0.9374),
        Circle((-9.0365, 0.0000), radius=0.9374),
        Circle((-11.2957, 0.0000), radius=0.9374),
        Circle((-13.5548, 0.0000), radius=0.9374),
        Circle((-15.8139, 0.0000), radius=0.9374),
        Circle((-42.3236, 85.5171), radius=0.9374),
        Circle((42.3236, 85.5171), radius=0.9374),
        Circle((48.9201, 90.8997), radius=1.4061),
        Circle((15.7343, 3.1159), radius=1.4061),
        Circle((12.2378, 3.1159), radius=1.4061),
        Circle((8.7413, 3.1159), radius=1.4061),
        Circle((5.2448, 3.1159), radius=1.4061),
        Circle((1.7483, 3.1159), radius=1.4061),
        Circle((-1.7483, 3.1159), radius=1.4061),
        Circle((-5.2448, 3.1159), radius=1.4061),
        Circle((-8.7413, 3.1159), radius=1.4061),
        Circle((-12.2378, 3.1159), radius=1.4061),
        Circle((-15.7343, 3.1159), radius=1.4061),
        Circle((-48.9201, 90.8997), radius=1.4061),
        Circle((41.5578, 90.8997), radius=1.4061),
        Circle((15.7343, 6.6124), radius=1.4061),
        Circle((12.2378, 6.6124), radius=1.4061),
        Circle((8.7413, 6.6124), radius=1.4061),
        Circle((5.2448, 6.6124), radius=1.4061),
        Circle((1.7483, 6.6124), radius=1.4061),
        Circle((-1.7483, 6.6124), radius=1.4061),
        Circle((-5.2448, 6.6124), radius=1.4061),
        Circle((-8.7413, 6.6124), radius=1.4061),
        Circle((-12.2378, 6.6124), radius=1.4061),
        Circle((-15.7343, 6.6124), radius=1.4061),
        Circle((-41.5578, 90.8997), radius=1.4061),
        Circle((49.7300, 97.7793), radius=1.8748),
        Circle((15.9124, 10.8720), radius=1.8748),
        Circle((11.3660, 10.8720), radius=1.8748),
        Circle((6.8196, 10.8720), radius=1.8748),
        Circle((2.2732, 10.8720), radius=1.8748),
        Circle((-2.2732, 10.8720), radius=1.8748),
        Circle((-6.8196, 10.8720), radius=1.8748),
        Circle((-11.3660, 10.8720), radius=1.8748),
        Circle((-15.9124, 10.8720), radius=1.8748),
        Circle((-49.7300, 97.7793), radius=1.8748),
        Circle((40.7478, 97.7793), radius=1.8748),
        Circle((15.9124, 15.4184), radius=1.8748),
        Circle((11.3660, 15.4184), radius=1.8748),
        Circle((6.8196, 15.4184), radius=1.8748),
        Circle((2.2732, 15.4184), radius=1.8748),
        Circle((-6.8196, 15.4184), radius=1.8748),
        Circle((-11.3660, 15.4184), radius=1.8748),
        Circle((-15.9124, 15.4184), radius=1.8748),
        Circle((-40.7478, 97.7793), radius=1.8748),
        Circle((45.2389, 82.7434), radius=2.3435),
        Circle((16.6201, 20.6997), radius=2.3435),
        Circle((11.0801, 20.6997), radius=2.3435),
        Circle((5.5400, 20.6997), radius=2.3435),
        Circle((-0.0000, 20.6997), radius=2.3435),
        Circle((-5.5400, 20.6997), radius=2.3435),
        Circle((-11.0801, 20.6997), radius=2.3435),
        Circle((-16.6201, 20.6997), radius=2.3435),
        Circle((-45.2389, 82.7434), radius=2.3435),
        Circle((16.6201, 26.2397), radius=2.3435),
        Circle((11.0801, 26.2397), radius=2.3435),
        Circle((5.5400, 26.2397), radius=2.3435),
        Circle((-0.0000, 26.2397), radius=2.3435),
        Circle((-5.5400, 26.2397), radius=2.3435),
        Circle((-11.0801, 26.2397), radius=2.3435),
        Circle((-16.6201, 26.2397), radius=2.3435),
        Circle((45.2389, 88.2909), radius=2.3435),
        Circle((16.6201, 31.7797), radius=2.3435),
        Circle((11.0801, 31.7797), radius=2.3435),
        Circle((5.5400, 31.7797), radius=2.3435),
        Circle((-0.0000, 31.7797), radius=2.3435),
        Circle((-5.5400, 31.7797), radius=2.3435),
        Circle((-11.0801, 31.7797), radius=2.3435),
        Circle((-16.6201, 31.7797), radius=2.3435),
        Circle((-45.2389, 88.2909), radius=2.3435),
        Circle((45.2389, 94.4215), radius=2.8122),
        Circle((16.7092, 38.1297), radius=2.8122),
        Circle((10.0255, 38.1297), radius=2.8122),
        Circle((3.3418, 38.1297), radius=2.8122),
        Circle((-3.3418, 38.1297), radius=2.8122),
        Circle((-10.0255, 38.1297), radius=2.8122),
        Circle((-16.7092, 38.1297), radius=2.8122),
        Circle((-45.2389, 94.4215), radius=2.8122),
        Circle((16.7092, 44.8133), radius=2.8122),
        Circle((10.0255, 44.8133), radius=2.8122),
        Circle((3.3418, 44.8133), radius=2.8122),
        Circle((-3.3418, 44.8133), radius=2.8122),
        Circle((-10.0255, 44.8133), radius=2.8122),
        Circle((-16.7092, 44.8133), radius=2.8122),
        Circle((-45.2389, 101.1361), radius=2.8122),
        Circle((16.7092, 51.4970), radius=2.8122),
        Circle((10.0255, 51.4970), radius=2.8122),
        Circle((-16.7092, 51.4970), radius=2.8122),
        Circle((-3.3418, 51.4970), radius=2.8122),
        Circle((-10.0255, 51.4970), radius=2.8122),
        Circle((3.3418, 51.4970), radius=2.8122),
        Circle((45.2389, 101.1361), radius=2.8122)]

    patch_list = [[i+1,p] for i,p in enumerate(raw_patches)]

    return np.array(patch_list)

def transform_patches(patches, pa=0, center=[0,0], reffiber=105, scale=1.):

    
    refcenter = np.array(patches[reffiber - 1,1].center) #in arcsec
    parad = -1.*pa*tau/360.
    decrad = center[1]*tau/360.
    rotrefx = refcenter[0]*np.cos(parad) - refcenter[1]*np.sin(parad)
    rotrefy = refcenter[0]*np.sin(parad) + refcenter[1]*np.cos(parad)
    for c in patches[:,1]:
        startx = c.center[0] #in arcsec
        starty = c.center[1] 
        
        #Rotate
        rotx = startx*np.cos(parad) - starty*np.sin(parad)
        roty = startx*np.sin(parad) + starty*np.cos(parad)

        #Shift
        shiftx = (rotx - rotrefx)/(np.cos(decrad)*3600.) + center[0]
        shifty = (roty - rotrefy)/3600. + center[1]

        c.center = (shiftx, shifty)
        c.radius *= scale

    return patches

def wcs2pix(patches, header):

    header_wcs = pywcs.WCS(header)

    for c in patches[:,1]:
        c.center= tuple(header_wcs.wcs_sky2pix([c.center],0)[0])

    return patches

def prep_axis(fitsfile = None, invert = True, sky = False, imrot = False,
              centpos = None, figsize = (8,8)):

    if fitsfile:
        hdu = pyfits.open(fitsfile)[0]
        imdata = hdu.data
        if imrot:
            imdata = spndi.rotate(imdata,-1*imrot,reshape=False)
            hdu.header.update('CROTA2',imrot)
        axistype = (wcsgrid.Axes, dict(header=hdu.header))
    else:
        hdu = None
        axistype = None

    fig = plt.figure(figsize=figsize)
    grid = ImageGrid(fig, 111,
                     nrows_ncols = (1,1),
                     cbar_mode = 'each',
                     cbar_location = 'top',
                     cbar_pad = '1%',
                     axes_class = axistype)
    ax = grid[0]
    
    if fitsfile:
        if invert:
            imdata = -1*(imdata - np.max(imdata))

        ax.imshow(imdata,
                  cmap = plt.get_cmap('gray'),
                  origin = 'lower', aspect = 'equal')
        ax.set_display_coord_system('fk5')
        
        hwcs = pywcs.WCS(hdu.header)
        centpx = hwcs.wcs_sky2pix([centpos],0)[0]

        labdict = {'nbins':10}
        ax.set_ticklabel_type('arcsec','arcsec',center_pixel=centpx,
                              labtyp1_kwargs=labdict,
                              labtyp2_kwargs=labdict)
        ax.add_compass(1)

    else:
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')
        if sky:
            ax.set_xlim(59,-59)
            ax.set_ylim(-2,116)
        else:
            ax.set_ylim(-2,60)
            ax.set_xlim(31,-31)
        
    return ax, hdu

def prep_patches(values,
                 hdu = None, pa = 0, center = [0,0], reffiber = 105,
                 sky = False, exclude = []):

    patches = GradPak_patches()
    skyidx = [0,1,17,18,19,30,31,42,43,52,53,61,62,70,78,86,87,94,101,108] 
    if hdu:
        scale = 2./((np.abs(hdu.header['CDELT1']) + \
                     np.abs(hdu.header['CDELT2']))*
                    3600.)
        patches = transform_patches(patches,
                                    reffiber = reffiber,
                                    pa = pa,
                                    center = center,
                                    scale = scale)
        patches = wcs2pix(patches, hdu.header) # now in px

    refcenter = patches[reffiber - 1,1].center # in px

    if not sky:
        exclude = np.r_[skyidx,np.array(exclude)-1]
    else:
        exclude = np.array(exclude) - 1

    exclude = np.array(exclude)
    exclude = np.unique(exclude)
    patches = np.delete(patches, exclude, axis=0)
    values = np.delete(values, exclude)
    patches = patches[values == values]
    pval = values[values == values]

    return patches, pval, refcenter

def plot(values,
         ax = None, figsize = (8,8),
         fitsfile = None, imrot = False, centpos = None,
         pa = 0, center = [0,0],
         reffiber = 105, invert=True,
         clabel = '', cmap = 'gnuplot2', 
         sky = False, labelfibers = True, exclude = [], 
         minval = None, maxval = None):


    tmpax, hdu = prep_axis(fitsfile = fitsfile, 
                           invert = invert, 
                           sky = sky, 
                           imrot = imrot, 
                           centpos = centpos,
                           figsize = figsize)

    if not ax:
        ax = tmpax
  
    patches, pval, refcenter = prep_patches(values,
                                            hdu = hdu, pa = pa, 
                                            center = center,
                                            reffiber = reffiber,
                                            sky = sky, exclude = exclude)

    xdelt = 2./(60. * hdu.header['CDELT1'])
    ydelt = 2./(60. * hdu.header['CDELT2'])
    ax.set_xlim(refcenter[0] + xdelt, refcenter[0] - xdelt)
    ax.set_ylim(refcenter[1] - ydelt, refcenter[1] + ydelt)
    

    if labelfibers:
        for c in patches:
            ax.text(c[1].center[0],c[1].center[1],c[0],fontsize=7,
                    ha='center',va='center')
    
    if minval is None:
        minval = pval.min()
    if maxval is None:
        maxval = pval.max()

    collection = PatchCollection(patches[:,1],
                                 cmap=plt.get_cmap(cmap),
                                 norm=matplotlib.colors.Normalize(
                                     vmin=minval,vmax=maxval))
    collection.set_array(pval)
    ax.add_collection(collection)

    cbar = ax.cax.colorbar(collection)
    cbar.set_label_text(clabel)

    return ax

def plot_img(values, 
             ax = None, figsize = (8,8),
             fitsfile = None, imrot = False,
             pa=0, center=[0,0], 
             reffiber = 105, invert = True,
             clabel='', cmap='gnuplot2', sky = False,
             numpoints=500, method='nearest',exclude=[],
             minval = None, maxval = None):
    
    if not ax:
        ax, hdu = prep_axis(fitsfile, invert, sky, imrot, figsize)
        
    patches, pval, refcenter = prep_patches(values,
                                            hdu = hdu, pa = pa, 
                                            center = center,
                                            reffiber = reffiber,
                                            sky = sky, exclude = exclude)
    xdelt = 2./(60. * hdu.header['CDELT1'])
    ydelt = 2./(60. * hdu.header['CDELT2'])
    ax.set_xlim(refcenter[0] + xdelt, refcenter[0] - xdelt)
    ax.set_ylim(refcenter[1] - ydelt, refcenter[1] + ydelt)

    x = np.array([c.center[0] for c in patches[:,1]])
    y = np.array([c.center[1] for c in patches[:,1]])

    xi = np.linspace(x.min(),x.max(),numpoints)
    yi = np.linspace(y.min(),y.max(),numpoints)

    vi = spi.griddata((x,y),
                      pval,
                      (xi[None,:],yi[:,None]),
                      method=method,
                      fill_value=np.inf)
    if minval is None:
        minval = pval.min()
    if maxval is None:
        maxval = pval.max()

    im = ax.imshow(vi, cmap=cmap, origin='lower', 
                   extent=(xi.min(),xi.max(),yi.min(),yi.max()),
                   vmin=minval,vmax=maxval)
    cbar = ax.cax.colorbar(im)
    cbar.set_label_text(clabel)

    return ax

def plot_rows(values, ylabel='', label='',
              ax = None, fullout = False,
              weights=None, kpc_scale=None, err=True,
              **plot_kwargs):

    y_values = np.array([c.center[1] for c in GradPak_patches()[:,1]])
    row_pos = np.unique(y_values)
    binned_vals = np.array([])
    binned_errs = np.array([])
    binned_stds = np.array([])
    abcissa = np.array([])
    if weights is None:
        weights = np.ones(y_values.size)

    for i, row in enumerate(row_pos):
        if row < 80:
            idx = np.where(y_values == row)[0]
            abcissa = np.append(abcissa, row)
            mean = np.sum(values[idx]*weights[idx]/np.sum(weights[idx]))
            std = np.sqrt(
                np.sum(weights[idx]*(values[idx] - mean)**2) /
                ((idx.size - 1.)/(idx.size) * np.sum(weights[idx])))
            std /= np.sqrt(idx.size)
            err = mean*np.sqrt(1./np.sum(weights[idx]**2))
            binned_vals = np.append(binned_vals,mean)
            binned_errs = np.append(binned_errs,err)
            binned_stds = np.append(binned_stds,std)

    if kpc_scale is not None:
        abcissa *= kpc_scale
        xlabel = 'Height [kpc]'
    else:
        xlabel = 'Height [arcsec]'

    if ax is None:
        ax = plt.figure(figsize=(8,8)).add_subplot(111)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    if err:
        ax.errorbar(abcissa, binned_vals, 
                    yerr = binned_errs,label = label, **plot_kwargs)
    else:
        ax.errorbar(abcissa, binned_vals, 
                    yerr = binned_stds,label = label, **plot_kwargs)
        
    if fullout:
        return ax, abcissa, binned_vals, binned_errs, binned_stds
    else:
        return ax
    

def format_tpl(tpl):

    with open(tpl,'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if line[0:6] == 'circle':
            data = line.split('(')[1].split(')')[0].replace('"','')
            x,y,r = [float(i) for i in data.split(',')]
            print 'Circle(({:6.4f}, {:6.4f}), radius={:6.4f}),'.format(
                x*-3600.,y*3600,r)

    return
