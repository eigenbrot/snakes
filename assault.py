
import numpy as np
import pyfits
import glob

from pyraf import iraf
from iraf import pysalt

from saltbias import saltbias
from saltcombine import saltcombine
from saltarith import saltarith
from saltprepare import saltprepare
from saltgain import saltgain
from saltxtalk import saltxtalk
from saltcrclean import saltcrclean
from saltmosaic import saltmosaic

from iraf import twodspec
from iraf import longslit

def assault(rawdir, biasdir, outdir):
    
    bias_list = glob.glob(biasdir+'*.fits')
    raw_list = glob.glob(rawdir+'P*.fits')

    saltbias(','.join(bias_list)+','+','.join(raw_list),'','o',subover=True,
             trim=False,subbias=False,masterbias='',median=False,
             function='polynomial',order=3,rej_lo=3,rej_hi=3,niter=10,
             plotover=False,turbo=False,clobber=False,logfile='salt.log',
             verbose=True)

    saltcombine(','.join(['o'+i for i in bias_list]),BIAS.fits,method='average',
                reject='sigclip',mask=no,weight=no,blank=0.0,scale=None,
                statsec=None
    
