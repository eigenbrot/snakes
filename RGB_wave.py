# From http://hsugawa.blogspot.com/2010/01/matplotlib-colormap-for-visible.html
import numpy as np

# These functions compute the transfor from wavelength (wl) to separate R, G,
# and B values
#
# Wavelength should be in Angstroms
def factor(wl):
    return np.select(
        [ wl > 7000.,
          wl < 4200.,
          True ],
        [ .3+.7*(7800.-wl)/(7800.-7000.),
          .3+.7*(wl-3800.)/(4200.-3800.),
          1.0 ] )

def raw_r(wl):
    return np.select(
        [ wl >= 5800.,
          wl >= 5100.,
          wl >= 4400.,
          wl >= 3800.,
          True ],
        [ 1.0,
          (wl-5100.)/(5800.-5100.),
          0.0,
          (wl-4400.)/(3800.-4400.),
          0.0 ] )

def raw_g(wl):
    return np.select(
        [ wl >= 6450.,
          wl >= 5800.,
          wl >= 4900.,
          wl >= 4400.,
          True ],
        [ 0.0,
          (wl-6450.)/(5800.-6450.),
          1.0,
          (wl-4400.)/(4900.-4400.),
          0.0 ] )

def raw_b(wl):
    return np.select(
        [ wl >= 5100.,
          wl >= 4900.,
          wl >= 3800.,
          True ],
        [ 0.0,
          (wl-5100.)/(4900.-5100.),
          1.0,
          0.0 ] )

gamma = 0.80
def correct_r(wl):
    return np.power(factor(wl)*raw_r(wl),gamma)
def correct_g(wl):
    return np.power(factor(wl)*raw_g(wl),gamma)
def correct_b(wl):
    return np.power(factor(wl)*raw_b(wl),gamma)

## And here's the function you actually use. Call it like this:
#
# ax.plot(x,y,color=RGB.rgb(5500))
#
# or whatever. The wavelength needs to be in Angstroms
##
def rgb(wl):
    return np.transpose([correct_r(wl),correct_g(wl),correct_b(wl)])
