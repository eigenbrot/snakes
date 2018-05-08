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
#  9. V_* - Velocity from chisq SSP fit
# 10. V_*^c - Velocity from chisq SSP fit + fiber-size-dependent offset
# 11. V_HaNS - Velocity from fitting Ha, [NII], SII complex
# 12. dV_HaNS - Uncertainty on above
# 13. V_obs - mean(V_*^c, V_HaNS)
# 14. dV_obs - Quadrature sum of dV_HaNS and half the difference between V_HaNS and V_*^c
#
""")
    f.write(('#{:4n}{:4n}{:4n}'+'{:10n}'*11+'\n#\n').format(*np.arange(14)+1))
    
    return

# def get_starV(pointing, basedir = '/Users/Arthur/Documents/School/891_research/final_results'):

#     vel_file = '{}/NGC_891_P{}_bin30_velocities.dat'.format(basedir,pointing)
#     vel = np.loadtxt(vel_file, usecols=(1,), unpack=True)

#     return vel

# def get_starE(pointing, dV=23.0,
#               chidV_dir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/var_disp/final_disp/chisq_vel/4166'):

#     chidV_file = '{}/NGC_891_P{}_bin30_allz2.coef.vel.fits'.format(chidV_dir, pointing)
#     chidV = pyfits.open(chidV_file)[1].data['VSYS_ERROR']
    
#     vel_err = np.sqrt(dV**2 + chidV**2)

#     return vel_err

# def get_gasV(pointing, gasdir='/d/monk/eigenbrot/WIYN/14B-0456/anal/HA_lines'):
    
#     gas_file = '{}/P{}_Ha_vel.txt'.format(gasdir, pointing)
#     vel = np.loadtxt(gas_file, usecols=(0,), unpack=True)
    
#     return vel

def get_rphi(pointing, basedir='/Users/Arthur/Documents/School/891_research/final_results'):
    
    rpfile = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, pointing)
    r, phi, dr = np.loadtxt(rpfile, usecols=(1,2,3), unpack=True)
    
    return r, phi, dr

def get_rhoz(pointing, basedir='/Users/Arthur/Documents/School/891_research/final_results'):

    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)

    return r, np.abs(z)

def get_Vs(pointing, basedir='/Users/Arthur/Documents/School/891_research/final_results'):

    vel = '{}/NGC_891_P{}_bin30_velocities.dat'.format(basedir, pointing)
    Vobs, dVobs, Vs, Vsc, Vg, dVg = np.loadtxt(vel, usecols=(1,2,3,4,5,6), unpack=True)

    return Vobs, dVobs, Vs, Vsc, Vg, dVg

def main(output):

    fmt = '{:5n}' + '{:4n}'*2 + '{:10.3f}'*11 + '\n'
    
    with open(output,'w') as f:
        write_header(f)
        i = 1
        for p in range(6):
            rho, z = get_rhoz(p+1)
            r, phi, dr = get_rphi(p+1)
            Vobs, dVobs, Vs, Vsc, Vg, dVg = get_Vs(p+1)
            
            print rho.shape, r.shape, Vobs.shape

            data = np.hstack((rho[:,None],z[:,None],r[:,None],dr[:,None],phi[:,None],
                              Vs[:,None],Vsc[:,None],Vg[:,None],dVg[:,None],Vobs[:,None],dVobs[:,None]))
            print data.shape

            for j in range(data.shape[0]):
                f.write(fmt.format(*([i,p+1,j+1] + data[j].tolist())))
                i += 1

    return
