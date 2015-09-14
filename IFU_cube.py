import pyfits
import numpy as np
from PIL import Image
import RGB_wave as RGB

def get_data(fitsfile):

    h = pyfits.open(fitsfile)[0]
    data = h.data
    wave = np.arange(data.shape[0])*h.header['CDELT3'] + h.header['CRVAL3']

    return wave, data

def get_color_overlay(wave):

    c = [np.r_[RGB.rgb(w),1] for w in wave]
    return np.vstack(c)

def cut(image, thresh, overlay):

    for i in range(image.shape[1]):
        idx = np.where(image[:,i,0] < thresh)
        overlay[idx,i,3] = 0

    return overlay

def trans(image):

    frac = np.sqrt(image[:,:,0]**2 + image[:,:,1]**2 + image[:,:,2]**2)
    image[:,:,3] = np.log10(frac*10)/np.max(np.log10(frac*10))

    return image

def wave_im(fitsfile, output, axis=1, threshold=0.1):

    wave, data = get_data(fitsfile)
    
    smoosh = np.sum(data, axis=axis)
    color_smoosh = np.repeat(smoosh[:,:,None],4,axis=2)
    uppercut = np.mean(color_smoosh) + np.std(color_smoosh)*5
    cidx = np.where(color_smoosh > uppercut)
    color_smoosh -= np.min(color_smoosh)
    color_smoosh[cidx] = 0
    color_smoosh /= np.max(color_smoosh)
    #color_smoosh /= np.mean(color_smoosh) + np.std(color_smoosh)*3
    color_smoosh[:,:,3] = 1.0

    overlay = np.repeat(get_color_overlay(wave)[:,None,:],smoosh.shape[1],axis=1)
    final_im = color_smoosh * cut(color_smoosh, threshold, overlay)

    norm = final_im * 255
    image = Image.fromarray(norm.astype(np.uint8))
    image.save(output)
    return final_im

def pos_im(fitsfile, output, waverange, threshold=0.1):

    wave, data = get_data(fitsfile)
    widx = np.where((wave >= waverange[0]) & (wave >= waverange[1]))
    smoosh = np.sum(data[widx], axis=0)
    color_smoosh = np.repeat(smoosh[:,:,None],4,axis=2)
    uppercut = np.mean(color_smoosh) + np.std(color_smoosh)*6

    cidx = np.where(color_smoosh > uppercut)
    color_smoosh -= np.min(color_smoosh)
    color_smoosh[cidx] = np.max(color_smoosh)
    color_smoosh /= np.max(color_smoosh)
    color_smoosh[:,:,3] = 1.0

    tidx = np.where(color_smoosh[:,:,0] < threshold)

    color = np.repeat(np.repeat(np.r_[RGB.rgb(np.mean(waverange)),1][None,None,:],smoosh.shape[0],axis=0),
                      smoosh.shape[1],axis=1)
    #color[:,:,3][tidx] = 0

    norm = color*trans(color_smoosh)*255
    image = Image.fromarray(norm.astype(np.uint8))
    image.save(output)

    return norm
