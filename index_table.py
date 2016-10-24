import numpy as np
import time
import D4000_indices as DI
import tau_indecies as ti
#import radial_projection as rp

"""For recording fiber locations and index measurements
"""

def write_header(f):

    f.write('# Generated on {}\n#\n'.format(time.asctime()))
    f.write("""#  1. Running sequence
#  2. Pointing
#  3. Apnum
#  4. r_proj (kpc) - Projected radius
#  5. z (kpc)
#  6. r (kpc) - Velocity-derived cylindrical radius
#  7. dr (kpc) - Uncertainty on r
#  8. phi (deg) - Velocity-derived cylindical angle.
#                  phi = 0 is tangent of -r_proj side of galaxy (approaching)
#  9. HdA 
# 10. Dn4000
# 11. [MgFe]
# 12. <Fe>
# 13. Mgb
#
""")
    f.write(('#{:4n}{:4n}{:4n}'+'{:10n}'*10+'\n#\n').format(*np.arange(13)+1))
    
    return

def get_indices(pointing, basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/indecies/D4000/tau_grid'):

    DIfile = '{}/NGC_891_P{}_bin30.msoz.Dn4000.dat'.format(basedir,pointing)
    tifile = '{}/NGC_891_P{}_bin30.msoz.bands.dat'.format(basedir,pointing)

    DIres = DI.quick_eat(DIfile)
    tires = ti.quick_eat(tifile)

    #Output will be numap X 5 and 2nd dim will be:
    # Hda, Dn4000, [MgFe], <Fe>, Mgb

    return np.hstack((DIres[:,[0,2]],tires[:,[6,5,7]]))


def get_rphi(pointing):

    rphi = 'NGC_891_P{}_bin30_rphi.dat'.format(pointing)
    r, phi, dr = np.loadtxt(rphi, usecols=(1,2,3), unpack=True)
    
    return r, phi, dr

    # loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)
    # vel = 'NGC_891_P{}_bin30_velocities.dat'.format(pointing)

    # r, phi = rp.compute_rphi(loc,vel)

    # return r, phi

def get_rhoz(pointing):

    loc = 'NGC_891_P{}_bin30_locations.dat'.format(pointing)
    r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)

    return r, z

def main(output):

    fmt = '{:5n}' + '{:4n}'*2 + '{:10.3f}'*10 + '\n'
    
    with open(output,'w') as f:
        write_header(f)
        i = 1
        for p in range(6):
            rho, z = get_rhoz(p+1)
            r, phi, dr = get_rphi(p+1)
            indices = get_indices(p+1)
            
            print rho.shape, z.shape, r.shape, phi.shape, indices.shape

            data = np.hstack((rho[:,None],z[:,None],r[:,None],dr[:,None],phi[:,None],indices))
            print data.shape

            for j in range(data.shape[0]):
                f.write(fmt.format(*([i,p+1,j+1] + data[j].tolist())))
                i += 1

    return
