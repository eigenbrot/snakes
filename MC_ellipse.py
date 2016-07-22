import numpy as np
import pyfits
import scipy.stats as ss
import scipy.optimize as spo
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
    #vals, vecs = np.linalg.eig(cov)
    vals, vecs = eigsorted(cov)
    print vals
    print vecs
    maxidx = np.argmax(vals)
    print maxidx
    theta1 = np.arctan2(vecs[1,maxidx],vecs[0,maxidx])
    theta2 = np.arctan2(vecs[2,maxidx],vecs[0,maxidx])
    theta3 = np.arctan2(vecs[1,maxidx],vecs[2,maxidx])
    # theta1 = np.arctan2(vecs[1,0],vecs[0,0])
    # theta2 = np.arctan2(vecs[2,0],vecs[0,0])
    # theta3 = np.arctan2(vecs[2,1],vecs[1,1])
    
    print theta1 * 180/np.pi, theta2*180/np.pi, theta3*180/np.pi

    a, b, c = nstd * np.sqrt(vals)
    len1 = a*np.cos(theta1)
    len2 = a*np.sin(theta1)
    len3 = a*np.sin(theta2)

    print len1*2, len2*2, len3*2

    # alim = a*np.cos(theta1)
    # blim = b*np.cos(theta2)
    # clim = b*np.cos(theta3)

    # alim = np.sqrt(len1**2 * np.cos(theta1)**2 + 
    #                len2**2 * np.sin(theta1)**2)
    # blim = np.sqrt(len2**2 * np.cos(theta1)**2 + 
    #                len1**2 * np.sin(theta1)**2)
    # clim = np.sqrt(len3**2 * np.cos(theta3)**2 + 
    #                len2**2 * np.sin(theta3)**2)

    alim = np.sqrt(len1**2 * np.cos(theta1)**2 + 
                   len2**2 * np.sin(theta1)**2 +
                   len3**2 * np.sin(theta2)**2)
    blim = np.sqrt(len2**2 * np.cos(theta1)**2 + 
                   len1**2 * np.sin(theta1)**2 +
                   len3**2 * np.sin(theta3)**2)
    clim = np.sqrt(len3**2 * np.cos(theta3)**2 + 
                   len2**2 * np.sin(theta3)**2 +
                   len1**2 * np.sin(theta2)**2)
    
    m1 = np.mean(d1)
    m2 = np.mean(d2)
    m3 = np.mean(d3)

    s1 = np.std(d1)
    s2 = np.std(d2)
    s3 = np.std(d3)

    print s1, s2, s3

    # print m1 - alim, m1 + alim
    # print m2 - blim, m2 + blim
    # print m3 - clim, m3 + clim    

    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    e1 = Ellipse(xy=(m1,m3), width=len1*2, height=len3*2, angle=theta2 * 180./np.pi,zorder=0,alpha=0.4)
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
    e2 = Ellipse(xy=(m1,m2), width=len1*2, height=len2*2, angle=theta1 * 180./np.pi,zorder=0,alpha=0.4)
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
    e3 = Ellipse(xy=(m3,m2), width=len3*2, height=len2*2, angle=theta3*180/np.pi,zorder=0,alpha=0.4)
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

def do_2d(X, Y, CI=0.68):

    cov = np.cov(X,Y)
    print cov
    vals, vecs = eigsorted(cov)
    theta = np.arctan2(*vecs[:,0][::-1])
    a, b =  np.sqrt(vals*ss.chi2.ppf(CI,2))

    xlim = np.sqrt((a*np.cos(theta))**2 + (b*np.sin(theta))**2)
    ylim = np.sqrt((a*np.sin(theta))**2 + (b*np.cos(theta))**2)

    return (xlim, ylim), (a, b, theta)
    
def do_ap_2d(pointing, ap, CI=0.68, output=None):
    
    apfile = 'MCdir/NGC_891_P{}_bin30_allz2.MC{:03n}.fits'.format(pointing,ap)
    print apfile
    d = pyfits.open(apfile)[1].data

    MLWA = d['MLWA']
    MLWZ = d['MLWZ']
    TAUV = d['TAUV']

    fig = plt.figure()
    fig.suptitle('{}.{}'.format(pointing,ap))
    ax1 = fig.add_subplot(223)
    ax2 = fig.add_subplot(221)
    ax3 = fig.add_subplot(224)
    
    lims1, e1 = do_2d(MLWA, TAUV, CI)
    lims2, e2 = do_2d(MLWA, MLWZ, CI)
    lims3, e3 = do_2d(MLWZ, TAUV, CI)

    # lims1, e1 = scale_ellipse(MLWA, TAUV, e1, CI)
    # lims2, e2 = scale_ellipse(MLWA, MLWZ, e2, CI)
    # lims3, e3 = scale_ellipse(MLWZ, TAUV, e3, CI)

    lims1, e1, lims2, e2, lims3, e3 = scale_ellipse3(MLWA, TAUV, MLWZ, e1, e2, e3, CI)

    print lims1[0], lims2[0]
    print lims2[1], lims3[0]
    print lims1[1], lims3[1]

    dMLWA = np.sqrt(lims1[0]**2 + lims2[0]**2)/np.sqrt(2)
    dMLWZ = np.sqrt(lims2[1]**2 + lims3[0]**2)/np.sqrt(2)
    dTAUV = np.sqrt(lims1[1]**2 + lims3[1]**2)/np.sqrt(2)
    
    print dMLWA
    print dMLWZ
    print dTAUV

    MLWA_m = np.mean(MLWA)
    MLWZ_m = np.mean(MLWZ)
    TAUV_m = np.mean(TAUV)

    MLWA_s = np.std(MLWA)
    MLWZ_s = np.std(MLWZ)
    TAUV_s = np.std(TAUV)

    ax1.scatter(MLWA, TAUV, c='k')
    ax2.scatter(MLWA, MLWZ, c='k')
    ax3.scatter(MLWZ, TAUV, c='k')
    
    draw_ellipse(MLWA, TAUV, e1, ax1)
    draw_ellipse(MLWA, MLWZ, e2, ax2)
    draw_ellipse(MLWZ, TAUV, e3, ax3)

    ax1.axvline(MLWA_m + dMLWA, ls='-')
    ax1.axvline(MLWA_m - dMLWA, ls='-')
    ax1.axhline(TAUV_m + dTAUV, ls='-')
    ax1.axhline(TAUV_m - dTAUV, ls='-')
    ax1.axvline(MLWA_m + MLWA_s, ls=':',color='k')
    ax1.axvline(MLWA_m - MLWA_s, ls=':',color='k')
    ax1.axhline(TAUV_m + TAUV_s, ls=':',color='k')
    ax1.axhline(TAUV_m - TAUV_s, ls=':',color='k')

    ax2.axvline(MLWA_m + dMLWA, ls='-')
    ax2.axvline(MLWA_m - dMLWA, ls='-')
    ax2.axhline(MLWZ_m + dMLWZ, ls='-')
    ax2.axhline(MLWZ_m - dMLWZ, ls='-')
    ax2.axvline(MLWA_m + MLWA_s, ls=':',color='k')
    ax2.axvline(MLWA_m - MLWA_s, ls=':',color='k')
    ax2.axhline(MLWZ_m + MLWZ_s, ls=':',color='k')
    ax2.axhline(MLWZ_m - MLWZ_s, ls=':',color='k')

    ax3.axvline(MLWZ_m + dMLWZ, ls='-')
    ax3.axvline(MLWZ_m - dMLWZ, ls='-')
    ax3.axhline(TAUV_m + dTAUV, ls='-')
    ax3.axhline(TAUV_m - dTAUV, ls='-')
    ax3.axvline(MLWZ_m + MLWZ_s, ls=':',color='k')
    ax3.axvline(MLWZ_m - MLWZ_s, ls=':',color='k')
    ax3.axhline(TAUV_m + TAUV_s, ls=':',color='k')
    ax3.axhline(TAUV_m - TAUV_s, ls=':',color='k')

    c1 = count_in_ellipse(MLWA, TAUV, e1)
    c2 = count_in_ellipse(MLWA, MLWZ, e2)
    c3 = count_in_ellipse(MLWZ, TAUV, e3)
    
    ax1.scatter(MLWA[c1], TAUV[c1], color='r', marker='s', facecolor='none', s=30, linewidths=0.5)
    ax2.scatter(MLWA[c2], MLWZ[c2], color='r', marker='s', facecolor='none', s=30, linewidths=0.5)
    ax3.scatter(MLWZ[c3], TAUV[c3], color='r', marker='s', facecolor='none', s=30, linewidths=0.5)

    ax1.text(0.9,0.9,'{:4.2f}'.format(len(c1)*1.0/MLWA.size), fontsize=9,transform=ax1.transAxes)
    ax2.text(0.9,0.9,'{:4.2f}'.format(len(c2)*1.0/MLWA.size), fontsize=9,transform=ax2.transAxes)
    ax3.text(0.9,0.9,'{:4.2f}'.format(len(c3)*1.0/MLWA.size), fontsize=9,transform=ax3.transAxes)

    ax1.set_xlim(*ax2.get_xlim())
    ax1.set_ylim(*ax3.get_ylim())
    ax3.set_xlim(*ax2.get_ylim())
    ax2.set_xticklabels([])
    ax3.set_yticklabels([])
    
    ax1.set_xlabel(r'$\tau_L$')
    ax1.set_ylabel(r'$\tau_V$')
    ax2.set_ylabel(r'$Z_L$')
    ax3.set_xlabel(r'$Z_L$')

    fig.subplots_adjust(hspace=0.001,wspace=0.001)
    fig.suptitle('{}.{}'.format(pointing,ap))

    if output:
        pp = PDF(output)
        pp.savefig(fig)
        pp.close()

    return dMLWA, dMLWZ, dTAUV, fig

def draw_ellipse(X, Y, ellipse, ax):
    
    a, b, theta = ellipse    
    pos = [np.mean(X),np.mean(Y)]
    e = Ellipse(xy=pos, width=2*a, height=2*b, angle=theta * 180./np.pi, alpha=0.4,zorder=0)
    ax.add_artist(e)

    return

def scale_ellipse(X, Y, ellipse, CI=0.68):

    p0 = [1.0]
    pf = spo.fmin(scale_func, p0, args=(X, Y, ellipse, CI),disp=False)
    a = ellipse[0]*pf[0]
    b = ellipse[1]*pf[0]
    theta = ellipse[2]
    xlim = np.sqrt((a*np.cos(theta))**2 + (b*np.sin(theta))**2)
    ylim = np.sqrt((a*np.sin(theta))**2 + (b*np.cos(theta))**2)

    print 'pf:', pf[0]

    return (xlim, ylim), (a, b, theta)

def scale_func(p, X, Y, ellipse, CI):
    
    a, b, theta = ellipse
    new_ellipse = (a*p[0], b*p[0], theta)
    idx = count_in_ellipse(X,Y,new_ellipse)
    frac = len(idx)*1.0/X.size

    return np.abs(CI - frac)

def scale_ellipse3(X, Y, Z, e1, e2, e3, CI=0.68):

    p0 = [1.0]
    pf = spo.fmin(scale_func3, p0, args=(X, Y, Z, e1, e2, e3, CI),disp=False)
    print 'pf:', pf[0]

    a1 = e1[0]*pf[0]
    b1 = e1[1]*pf[0]
    theta1 = e1[2]
    xlim1 = np.sqrt((a1*np.cos(theta1))**2 + (b1*np.sin(theta1))**2)
    ylim1 = np.sqrt((a1*np.sin(theta1))**2 + (b1*np.cos(theta1))**2)

    a2 = e2[0]*pf[0]
    b2 = e2[1]*pf[0]
    theta2 = e2[2]
    xlim2 = np.sqrt((a2*np.cos(theta2))**2 + (b2*np.sin(theta2))**2)
    ylim2 = np.sqrt((a2*np.sin(theta2))**2 + (b2*np.cos(theta2))**2)

    a3 = e3[0]*pf[0]
    b3 = e3[1]*pf[0]
    theta3 = e3[2]
    xlim3 = np.sqrt((a3*np.cos(theta3))**2 + (b3*np.sin(theta3))**2)
    ylim3 = np.sqrt((a3*np.sin(theta3))**2 + (b3*np.cos(theta3))**2)

    return (xlim1, ylim1), (a1, b1, theta1), (xlim2, ylim2), (a2, b2, theta2), (xlim3, ylim3), (a3, b3, theta3)

def scale_func3(p, X, Y, Z, e1, e2, e3, CI):
    
    a1, b1, theta1 = e1
    ne1 = (a1*p[0], b1*p[0], theta1)
    idx1 = count_in_ellipse(X,Y,ne1)
    frac1 = len(idx1)*1.0/X.size

    a2, b2, theta2 = e2
    ne2 = (a2*p[0], b2*p[0], theta2)
    idx2 = count_in_ellipse(X,Z,ne2)
    frac2 = len(idx2)*1.0/X.size

    a3, b3, theta3 = e3
    ne3 = (a3*p[0], b3*p[0], theta3)
    idx3 = count_in_ellipse(Z,Y,ne3)
    frac3 = len(idx3)*1.0/X.size

    return np.sqrt((frac1 - CI)**2 + (frac2 - CI)**2 + (frac3 - CI)**2)

def count_in_ellipse(X, Y, ellipse):

    a, b, theta = ellipse
    mx = np.mean(X)
    my = np.mean(Y)
    in_idx = []
    for i in range(X.size):
        t1 = ((X[i] - mx)*np.cos(theta) + (Y[i] - my)*np.sin(theta))**2
        t2 = ((X[i] - mx)*np.sin(theta) - (Y[i] - my)*np.cos(theta))**2
        f = t1/a**2 + t2/b**2
        if f <= 1:
            in_idx.append(i)

#    print len(in_idx)*1.0/X.size

    return in_idx

def eigsorted(cov):
        vals, vecs = np.linalg.eig(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

def do_pointing(pointing, outpre, CI=0.68, basedir='.'):

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
        dMLWA, dMLWZ, dTAUV, f = do_ap_2d(pointing, i+1, CI)
        output['MLWA'][i] = coefs['MLWA'][i]
        output['MLWZ'][i] = coefs['MLWZ'][i]
        output['TAUV'][i] = coefs['TAUV'][i]
        output['dMLWA'][i] = dMLWA
        output['dMLWZ'][i] = dMLWZ
        output['dTAUV'][i] = dTAUV
        pp.savefig(f)
        plt.close(f)

    pp.close()
    pyfits.BinTableHDU(output).writeto(outpre+'.fits',clobber=True)

    return

def do_all(CI=0.68,basedir='.'):

    for p in range(6):
        outpre = 'NGC_891_P{}_bin30_allz2.fiterr'.format(p+1)
        do_pointing(p+1,outpre,CI=CI,basedir=basedir)

    return

def d2_test(CI=0.68):
    cov = [[1,0],[0,1]]
    # cov = [[
    # cov = [[ 0.30773231, -0.02568195],
    #        [-0.02568195,  0.00344581]]
    # cov = [[ 0.0079945,  -0.00014535],
    #        [-0.00014535,  0.00025371]]
    # cov = [[ 0.0079945,  -0.00112452],
    #        [-0.00112452,  0.00109008]]
    cov = [[ 0.01432639,  0.00315234],
           [ 0.00315234,  0.00125425]]

    t = np.random.multivariate_normal(mean=(1,1),cov=cov,size=1000)
    tX = t[:,0]
    tY = t[:,1]

    ax = plt.figure().add_subplot(111)
    ax.scatter(tX, tY, c='k')
    lim, e = do_2d(tX, tY, ax=ax, CI=CI)
    c = count_in_ellipse(tX, tY, e)
    ax.scatter(tX[c], tY[c], color='r', marker='s', facecolor='none', s=30, linewidths=0.5)

    pp = PDF('test.pdf')
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)

    return
