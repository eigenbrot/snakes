import numpy as np
import time
import pyfits

def write_header(f):

    f.write('# Generated on {}\n#\n'.format(time.asctime()))
    f.write("""#  1. Running sequence
#  2. Pointing
#  3. Apnum
#  4. r_proj (kpc) - Projected radius
#  5. |z| (kpc)
#  6. r (kpc) - Velocity-derived cylindrical radius
#  7. dr (kpc) - Uncertainty on r
#  8. phi (deg) - Velocity-derived cylindical angle.
#                  phi = 0 is tangent of -r_proj side of galaxy (approaching)
#  9. Vstars - Velocity from chisq SSP fit + offset of 74 km/s
# 10. dV - Total uncertainty in Vstars, comes from wavelength solnn and fit
# 11. Vgas - Velocity from H_alpha centroid
#
""")
    f.write(('#{:4n}{:4n}{:4n}'+'{:10n}'*8+'\n#\n').format(*np.arange(11)+1))
    
    return

def get_starV(pointing, basedir = '/Users/Arthur/Documents/School/891_research/final_results'):

    vel_file = '{}/NGC_891_P{}_bin30_velocities.dat'.format(basedir,pointing)
    vel = np.loadtxt(vel_file, usecols=(1,), unpack=True)

    return vel

def get_starE(pointing, dV=23.0,
              chidV_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/var_disp/final_disp/chisq_vel/4166'):

    chidV_file = '{}/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(chidV_dir, pointing)
    chidV = pyfits.open(chidV_file)[1].data['VSYS_ERROR']
    
    vel_err = np.sqrt(dV**2 + chidV**2)

    return vel_err

def get_gasV(pointing, gasdir='/d/monk/eigenbrot/WIYN/14B-0456/anal/HA_lines'):
    
    gas_file = '{}/P{}_Ha_vel.txt'.format(gasdir, pointing)
    vel = np.loadtxt(gas_file, usecols=(0,), unpack=True)
    
    return vel

def get_rphi(pointing, basedir='/Users/Arthur/Documents/School/891_research/final_results'):
    
    rpfile = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, pointing)
    r, phi, dr = np.loadtxt(rpfile, usecols=(1,2,3), unpack=True)
    
    return r, phi, dr

def get_rhoz(pointing, basedir='/Users/Arthur/Documents/School/891_research/final_results'):

    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)

    return r, np.abs(z)

def main(output):

    fmt = '{:5n}' + '{:4n}'*2 + '{:10.3f}'*8 + '\n'
    
    with open(output,'w') as f:
        write_header(f)
        i = 1
        for p in range(6):
            rho, z = get_rhoz(p+1)
            r, phi, dr = get_rphi(p+1)
            stars = get_starV(p+1)
            star_e = get_starE(p+1)
            gas = get_gasV(p+1)
            
            print rho.shape, stars.shape, star_e.shape, gas.shape

            if gas.size < stars.size:
                gas = np.r_[gas,[np.nan]*(stars.size - gas.size)]

            data = np.hstack((rho[:,None],z[:,None],r[:,None],dr[:,None],phi[:,None],stars[:,None],star_e[:,None],gas[:,None]))
            print data.shape

            for j in range(data.shape[0]):
                f.write(fmt.format(*([i,p+1,j+1] + data[j].tolist())))
                i += 1

    return
