import numpy as np
import time
import pyfits

"""For recording fiber locations and index measurements
"""

def write_header(f):

    f.write('# Generated on {}\n#\n'.format(time.asctime()))
    f.write("""#  1. Running sequence
#  2. Pointing
#  3. Apnum
#  4. r_proj (kpc) - Projected radius
#  5. |z| (kpc)
#  6. r (kpc) - Velocity-derived cylindrical radius
#  7. dr (kpc) -Uncertainty on r
#  8. phi (deg) - Velocity-derived cylindical angle.
#                  phi = 0 is tangent of -r_proj side of galaxy (approaching)
#  9. MLWA
# 10. dMLWA
# 11. MMWA
# 12. MLWZ
# 13. dMLWZ
# 14. MMWZ
# 15. A_V
# 16. dA_V
#
""")
    f.write(('#{:4n}{:4n}{:4n}'+'{:10n}'*13+'\n#\n').format(*np.arange(16)+1))
    
    return

def get_fitdata(pointing, basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results', suffix='fiterr'):

    coefs = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,pointing,suffix)
    d = pyfits.open(coefs)[1].data

    return d['MLWA'], d['dMLWA'], d['MLWZ'], d['dMLWZ'], d['TAUV']*1.086, d['dTAUV']*1.086

def get_massdata(pointing, basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results', suffix='coef'):
    
    coefs = '{}/NGC_891_P{}_bin30_allz2.{}.fits'.format(basedir,pointing,suffix)
    d = pyfits.open(coefs)[1].data

    return d['MMWA'], d['MMWZ']

def get_rphi(pointing, basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'):
    
    rpfile = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, pointing)
    r, phi, dr = np.loadtxt(rpfile, usecols=(1,2,3), unpack=True)
    
    return r, phi, dr

def get_rhoz(pointing, basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'):

    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)

    return r, np.abs(z)

def main(output):

    fmt = '{:5n}' + '{:4n}'*2 + '{:10.3f}'*13 + '\n'
    
    with open(output,'w') as f:
        write_header(f)
        i = 1
        for p in range(6):
            rho, z = get_rhoz(p+1)
            r, phi, dr = get_rphi(p+1)
            MLWA, dMLWA, MLWZ, dMLWZ, AV, dAV = get_fitdata(p+1)
            MMWA, MMWZ = get_massdata(p+1)

            print rho.shape, r.shape, MLWA.shape, AV.shape

            data = np.vstack((rho,z,r,dr,phi,MLWA, dMLWA, MMWA, MLWZ, dMLWZ, MMWZ, AV, dAV)).T
            print data.shape

            for j in range(data.shape[0]):
                f.write(fmt.format(*([i,p+1,j+1] + data[j].tolist())))
                i += 1

    return
