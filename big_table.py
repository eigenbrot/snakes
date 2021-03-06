import numpy as np
import pyfits
import time
from glob import glob

"""For making big tables with all the CI stuff
"""

# No allZ
# uw_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/CI'
# I_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/IRAF_weight'
# I15_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/IRAF15'
# IP_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/IRAF_prior'
# IP15_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/IRAF_P15'
# IA_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/single_Z/IRAF_all'

#with allZ
uw_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/no_weight'
I_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/IRAF_weight'
I15_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/IRAF15'
IP_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/IRAF_prior'
IP15_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/IRAF_P15'
IA_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/model_balmer/bc03/full_CI/IRAF_all'


FMT = '{:8n}{:8n}{:8n}{:10.3f}{:10.3f}' + str(' '+'{:10.3f}'*3)*6

def do_all(output):

    with open(output,'w') as f:
        write_header(f)
        seq=0

        for i in range(6):
            seq = do_pointing(i+1,f,seq)


    return

def do_pointing(pointing,f,seq):

    loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)
    r, z = np.loadtxt(loc,usecols=(4,5),unpack=True)

    numap = r.size

    uw_chi, uw_t, uw_dt = get_basics(pointing, uw_dir)
    _, I_t, I_dt = get_basics(pointing, I_dir)
    _, I15_t, I15_dt = get_basics(pointing, I15_dir)
    _, IP_t, IP_dt = get_basics(pointing, IP_dir)
    _, IP15_t, IP15_dt = get_basics(pointing, IP15_dir)
    _, IA_t, IA_dt = get_basics(pointing, IA_dir)
    
    I_DC = compute_DC(pointing, I_dir, uw_chi)
    I15_DC = compute_DC(pointing, I15_dir, uw_chi)
    IP_DC = compute_DC(pointing, IP_dir, uw_chi)
    IP15_DC = compute_DC(pointing, IP15_dir, uw_chi)
    IA_DC = compute_DC(pointing, IA_dir, uw_chi)
            
    for a in range(numap):
        f.write(FMT.format(seq+1,
                           pointing,
                           a+1,
                           r[a],
                           z[a],
                           uw_chi[a], uw_t[a], uw_dt[a],
                           I_DC[a], I_t[a], I_dt[a],
                           I15_DC[a], I15_t[a], I15_dt[a],
                           IP_DC[a], IP_t[a], IP_dt[a],
                           IP15_DC[a], IP15_t[a], IP15_dt[a],
                           IA_DC[a], IA_t[a], IA_dt[a]))
        f.write('\n')
        seq += 1

    return seq

def get_basics(pointing, folder):
    
    CI_file = glob('{}/*P{}*CI*.dat'.format(folder,pointing))[0]
    print CI_file
    
    t, lt, ht, chisq = np.loadtxt(CI_file, usecols=(1,2,3,4), unpack=True)
    
    return chisq, t, ht + lt

def compute_DC(pointing, folder, uw_chi):

    CI_file = glob('{}/*P{}*CI*.dat'.format(folder,pointing))[0]
    
    bestZ = np.loadtxt(CI_file, usecols=(5,), unpack=True, dtype=np.int)
    
    fzlist = ['0.005Z','0.02Z','0.2Z','0.4Z','1Z','2.5Z','allZ']

    hdu = pyfits.open('NGC_891_P{}_bin30.mso.fits'.format(pointing))[0]
    head = hdu.header
    data = hdu.data
    error = pyfits.open('NGC_891_P{}_bin30.meo.fits'.format(pointing))[0].data

    wave = (np.arange(data.shape[1]) - head['CRPIX1'] - 1)*head['CDELT1'] + head['CRVAL1']
    idx = np.where((wave >= 3800.) & (wave <= 6800.))[0]
    wave = wave[idx]
    data = data[:,idx]
    error = error[:,idx]    

    outarr = np.zeros(data.shape[0])

    for i, bz in enumerate(bestZ):
        best_file = '{}/{}/NGC_891_P{}_bin30_allz2.fit.fits'.\
                    format(folder,fzlist[bz],pointing)
        print i+1, fzlist[bz]
        models = pyfits.open(best_file)[0].data

        coef_file = '{}/{}/NGC_891_P{}_bin30_allz2.coef.fits'.\
                    format(folder,fzlist[bz],pointing)
        coefs = pyfits.open(coef_file)[1].data
        
        chisq = np.sum((data[i,:] - models[i,:])**2/error[i,:]**2)/coefs['TOTFREE'][i]
        outarr[i] = uw_chi[i] - chisq

    return outarr

def write_header(f):

    f.write('# Generated on {}\n#\n'.format(time.asctime()))
    f.write('# Delta chisq is unweighted - XX\n')
    f.write('# allZ fits ARE included\n#\n')
    f.write("""#  1. Running sequence
#  2. Pointing
#  3. Apnum
#  4. r (kpc)
#  5. z (kpz)
#  6. unweighted chisq
#  7. unweighted  MLWA
#  8. unweighted dMLWA (min/max)
#  9. IRAF weight (old) delta unweighted chisq
# 10. IRAF weight  MLWA
# 11. IRAF weight dMLWA
# 12. IRAF weight (old) p=1.5 delta unweighted chisq
# 13. IRAF weight p=1.5  MLWA
# 14. IRAF weight p=1.5 dMLWA
# 15. IRAF w/ priors delta unweighted chisq
# 16. IRAF w/ priors  MLWA
# 17. IRAF w/ priors dMLWA
# 18. IRAF w/ priors p=1.5 delta unweighted chisq
# 19. IRAF w/ priors p=1.5  MLWA
# 20. IRAF w/ priors p=1.5 dMLWA
# 21. IRAF w/ everything delta unweighted chisq
# 22. IRAF w/ everything  MLWA
# 23. IRAF w/ everything dMLWA
""")
    f.write('#')
    f.write(' '*43)
    f.write(' {:^30}'.format('unweighted'))
    f.write(' {:^30}'.format('I'))
    f.write(' {:^30}'.format('I15'))
    f.write(' {:^30}'.format('IP'))
    f.write(' {:^30}'.format('IP15'))
    f.write(' {:^30}\n'.format('IA'))
    f.write('#'+' '*43)
    f.write(str('   '+'-'*28)*6)
    f.write('\n')
    f.write(str('#{:7n}{:8n}{:8n}{:10n}{:10n}' + str(' '+'{:10n}'*3)*6).format(*np.arange(23)+1))
    f.write('\n#\n')

    return
