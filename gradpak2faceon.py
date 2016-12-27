import time
import numpy as np

def bin_r(datatable, rbins, rrange, zbins, hz, zerr):
    '''Take a data table and bin it by radius. For each radius call
    do_height() to compute the face-on integrated values

    '''

    data = np.loadtxt(datatable)
    print data.shape
    h, bin_edge = np.histogram(data[:,5], bins=rbins, range=rrange)
    
    # #get rid of bins with nobody in them
    # zero_idx = np.where(h == 0)[0]
    # bin_edge = np.delete(bin_edge, zero_idx + 1)
    # numbins = rbins - zero_idx.size

    #we do this so we can keep the < bin_edge[i+1] constraint
    bin_edge[-1] *= 1.1
    
    #choose a sensible widening value
    add_r = np.mean(np.diff(bin_edge))/3.
    rlimit = bin_edge[-1]
    
    # -4 because we don't need running seq, pointing, apnum, r_proj, or phi
    results = np.zeros((rbins, data.shape[1] - 4))

    j = 0
    for i in range(rbins):
        idx = np.where((data[:,5] >= bin_edge[i]) &
                       (data[:,5] < bin_edge[i+1]))[0]
        if idx.size == 0:
            continue
        print bin_edge[i], '<= r < ', bin_edge[i+1]
        res = do_height(data[idx,:], zbins, hz)

        # In a perfect world the average, weighted height should be
        # close to the scale height, let's make sure it's close, and
        # widen the radial bin if necessary
        if np.abs(res[0] - hz) > zerr:
            print 'Widening radial bin... [{}, {})'.format(bin_edge[i], bin_edge[i+1])
            while np.abs(res[0] - hz) > zerr:
                bin_edge[i+1:] += add_r
                
                idx = np.where((data[:,5] >= bin_edge[i]) & (data[:,5] < bin_edge[i+1]))[0]
                print '\t {} [{}, {})...'.format(i, bin_edge[i], bin_edge[i+1])
                res = do_height(data[idx,:], zbins, hz)

                if bin_edge[i+1] > rlimit:
                    print "\t reached rlimit"
                    break
                
            print 'OK!'
                
        print '#', i, idx.size, res.size, res[0]

        results[j,:] = res
        j += 1
        
        if bin_edge[i+1] > rlimit:
            break
        
    return results[:j,:]

def do_height(data, zbins=10, hz=0.42):
    '''Take in a numpy array and bin into an equal number of vertical
    slices, average across these slices. Finally average everything
    together while weighting for height

    '''

    #which columns of the original array will we keep?
    keep_idx = [4,5,6,8,9,10,11,12,13,14,15]

    #which columns (of the new array) are uncertainties?
    err_idx = [4,7,10]

    # rounding to 0.1 kpc helps us find the midplane
    roundz = np.round(data[:,4],1)
    z = np.array(data[:,4])
    
    # the first row in each pointing is actually below the plane, so
    # we'll multiply that by -1 so it doesn't get grouped with row 3
    zid = np.argmin(roundz)
    z[:zid] *= -1

    h, bin_edge = np.histogram(z, bins=zbins)

    #get rid of bins with nobody in them
    zero_idx = np.where(h == 0)[0]
    bin_edge = np.delete(bin_edge, zero_idx + 1)
    numbins = zbins - zero_idx.size

    #we do this so we can keep the < bin_edge[i+1] constraint
    bin_edge[-1] *= 1.1

    tmpres = np.zeros((numbins, data.shape[1]))

    for i in range(numbins):
        idx = np.where((z >= bin_edge[i]) & (z < bin_edge[i+1]))[0]
        print '\t\t', bin_edge[i], bin_edge[i+1], idx.size
        tmpres[i,:] = np.mean(data[idx,:], axis=0)

    #now, collapse along z with a scale-height weighting
    weight = np.exp(-tmpres[:,4]/hz)
    # print tmpres[:,4]
    res = np.sum(tmpres[:,keep_idx] * weight[:,None], axis=0)/np.sum(weight)
    # print res[0]
    
    #divide errors by root N
    res[err_idx] /= np.sqrt(tmpres.shape[0])

    #Throw a log(MMWZ) at the end for the fans
    res = np.append(res, np.log10(res[8]))

    return res

def write_header(f):

    f.write('# Generated on {}\n#\n'.format(time.asctime()))
    f.write("""#  1. Running sequence
#  2. |z| (kpc)
#  3. r (kpc) - Velocity-derived cylindrical radius
#  4. dr (kpc) - Uncertainty on r
#  5. MLWA
#  6. dMLWA
#  7. MMWA
#  8. MLWZ
#  9. dMLWZ
# 10. MMWZ
# 11. A_V
# 12. dA_V
# 13. log(MMWZ)
#
""")
    f.write(('#{:4n}'+'{:10n}'*12+'\n#\n').format(*np.arange(13)+1))
    
    return

def main(datatable='paper2_table.txt', output='paper2_faceon.dat', rbins=300, zbins=10, rrange=(0,22), hz=0.42, zerr=0.1):

    fmt = '{:5n}' + '{:10.3f}'*12 + '\n'

    results = bin_r(datatable, rbins, rrange, zbins, hz, zerr)

    print results.shape
    
    with open(output,'w') as f:
        write_header(f)
        for i in range(results.shape[0]):
            f.write(fmt.format(*([i+1] + results[i,:].tolist())))


    return

#########

def M31_comp(output, faceonfile='paper2_faceon.dat', offset=0.3):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages as PDF
    
    data = np.loadtxt(faceonfile).T

    ax = plt.figure().add_subplot(111)
    ax.set_xlabel(r'$r/h_r$ (Deprojected)')
    ax.set_ylabel('[M/H]')

    model_r = np.linspace(0,20,50)
    model_dex = model_r * -0.02 + offset# From Gregerson+ '15 Figure 9
    model_dexp = model_r * -0.024 + offset
    model_dexm = model_r * -0.016 + offset

    NGC_dex = model_r * -0.0095 + 0.27

    print np.polyfit(data[2], data[12], 1)
    
    ax.plot(model_r/5.26, model_dex, '-k', label="-0.02 dex/kpc (Gregerson+ '15) + {}".format(offset))
    ax.fill_between(model_r/5.26, model_dexp, model_dexm, color='k', alpha=0.2, interpolate=True)

    ax.plot(data[2]/5.48, data[12], ls='None', color='r', marker='o', mec='none')
    ax.plot(model_r/5.48, NGC_dex, '-r', label='-0.0095 dex/kpc + 0.27')
    ax.legend(loc=0, frameon=False)
    
    pp = PDF(output)
    pp.savefig(ax.figure)
    pp.close()
    plt.close(ax.figure)

    return
    
