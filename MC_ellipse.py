import numpy as np
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages as PDF
from matplotlib.patches import Ellipse

def do_ap(pointing, ap, nstd=1, output=None):

    apfile = 'MCdir/NGC_891_P{}_bin30_allz2.MC{:03n}.fits'.format(pointing,ap)
    print apfile
    d = pyfits.open(apfile)[1].data

    f1 = 'MLWA'
    f2 = 'TAUV'
    f3 = 'MLWZ'
    
    d1 = d[f1]
    d2 = d[f2]
    d3 = d[f3]

    ddata = np.vstack((d1,d2,d3))
    cov = np.cov(ddata)
    print cov
    vals, vecs = np.linalg.eig(cov)
    print vals
    print np.sqrt(np.diagonal(cov))
    print np.sqrt(vals)
    maxidx = np.argmax(vals)
    print maxidx
    # theta1 = np.arctan2(vecs[1,maxidx],vecs[0,maxidx])
    # theta2 = np.arctan2(vecs[2,maxidx],vecs[0,maxidx])
    # theta3 = np.arctan2(vecs[1,maxidx],vecs[2,maxidx])
    theta1 = np.arctan2(vecs[1,0],vecs[0,0])
    theta2 = np.arctan2(vecs[2,0],vecs[0,0])
    theta3 = np.arctan2(vecs[2,1],vecs[1,1])
    
    print theta1 * 180/np.pi, theta2*180/np.pi, theta3*180/np.pi

    a, b, c = nstd * np.sqrt(vals)
    
    print a*2, b*2, c*2

    # alim = a*np.cos(theta1)
    # blim = b*np.cos(theta2)
    # clim = b*np.cos(theta3)

    alim = np.sqrt(a**2 * np.cos(theta1)**2 + 
                   b**2 * np.sin(theta1)**2 +
                   c**2 * np.sin(theta2)**2)
    blim = np.sqrt(b**2 * np.cos(theta1)**2 + 
                   a**2 * np.sin(theta1)**2 +
                   c**2 * np.sin(theta3)**2)
    clim = np.sqrt(c**2 * np.cos(theta3)**2 + 
                   b**2 * np.sin(theta3)**2 +
                   a**2 * np.sin(theta2)**2)
    
    m1 = np.mean(d1)
    m2 = np.mean(d2)
    m3 = np.mean(d3)

    s1 = np.std(d1)
    s2 = np.std(d2)
    s3 = np.std(d3)

    print s1, s2, s3

    print m1 - alim, m1 + alim
    print m2 - blim, m2 + blim
    print m3 - clim, m3 + clim    

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    e1 = Ellipse(xy=(m1,m3), width=a*2, height=c*2, angle=theta2 * 180./np.pi,zorder=0,alpha=0.4)
    ax1.add_artist(e1)
    ax1.scatter(d1,d3,c='k')
    ax1.axvline(m1-alim)
    ax1.axvline(m1+alim)
    ax1.axhline(m3+clim)
    ax1.axhline(m3-clim)
    ax1.axvline(m1-s1,color='k',ls=':')
    ax1.axvline(m1+s1,color='k',ls=':')
    ax1.axhline(m3-s3,color='k',ls=':')
    ax1.axhline(m3+s3,color='k',ls=':')


    ax2 = fig.add_subplot(223)
    e2 = Ellipse(xy=(m1,m2), width=a*2, height=b*2, angle=theta1 * 180./np.pi,zorder=0,alpha=0.4)
    ax2.add_artist(e2)
    ax2.scatter(d1,d2,c='k')
    ax2.axvline(m1-alim)
    ax2.axvline(m1+alim)
    ax2.axhline(m2+blim)
    ax2.axhline(m2-blim)
    ax2.axvline(m1-s1,color='k',ls=':')
    ax2.axvline(m1+s1,color='k',ls=':')
    ax2.axhline(m2-s2,color='k',ls=':')
    ax2.axhline(m2+s2,color='k',ls=':')

    ax3 = fig.add_subplot(224)
    e3 = Ellipse(xy=(m3,m2), width=c*2, height=b*2, angle=theta3*180/np.pi,zorder=0,alpha=0.4)
    ax3.add_artist(e3)
    ax3.scatter(d3,d2,c='k')
    ax3.axvline(m3-clim)
    ax3.axvline(m3+clim)
    ax3.axhline(m2+blim)
    ax3.axhline(m2-blim)
    ax3.axvline(m3-s3,color='k',ls=':')
    ax3.axvline(m3+s3,color='k',ls=':')
    ax3.axhline(m2-s2,color='k',ls=':')
    ax3.axhline(m2+s2,color='k',ls=':')

    ax2.set_xlim(*ax1.get_xlim())
    ax2.set_ylim(*ax3.get_ylim())
    ax3.set_xlim(*ax1.get_ylim())
    ax1.set_xticklabels([])
    ax3.set_yticklabels([])
    
    ax2.set_xlabel(f1)
    ax2.set_ylabel(f2)
    ax1.set_ylabel(f3)
    ax3.set_xlabel(f3)

    fig.subplots_adjust(hspace=0.001,wspace=0.001)
    fig.suptitle('{}.{}'.format(pointing,ap))

    if output is not None:
        pp = PDF(output)
        pp.savefig(fig)
        pp.close()

    return alim, blim, clim, fig

def eigsorted(cov):
        vals, vecs = np.linalg.eig(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

def do_pointing(pointing, outpre, nstd=1, basedir='.'):

    pp = PDF(outpre+'.pdf')
    ap = 1
    
    #get number of aps
    coefs = pyfits.open('{}/NGC_891_P{}_bin30_allz2.coef.fits'.\
                      format(basedir,pointing))[1].data
    numap = coefs.shape[0]
    output = np.zeros(numap, dtype=[('MLWA','>f8'), ('dMLWA', '>f8'),
                                    ('MLWZ','>f8'), ('dMLWZ', '>f8'),
                                    ('TAUV','>f8'), ('dTAUV', '>f8')])

    for i in range(numap):
        print i+1
        a, b, c, f = do_ap(pointing, i+1, nstd)
        output['MLWA'][i] = coefs['MLWA'][i]
        output['MLWZ'][i] = coefs['MLWZ'][i]
        output['TAUV'][i] = coefs['TAUV'][i]
        output['dMLWA'][i] = a
        output['dMLWZ'][i] = b
        output['dTAUV'][i] = c
        pp.savefig(f)
        plt.close(f)

    pp.close()
    pyfits.BinTableHDU(output).writeto(outpre+'.fits',clobber=True)

    return

def do_all(nstd=1,basedir='.'):

    for p in range(6):
        outpre = 'NGC_891_P{}_bin30_allz2.fiterr'.format(p+1)
        do_pointing(p+1,outpre,nstd=nstd,basedir=basedir)

    return
