import numpy as np
import time

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
#  9. Dn4000
# 10. Dn4000_e
# 11. HdA
# 12. HdA_e
# 13. Mgb
# 14. Mgb_e
# 15. <Fe>
# 16. <Fe>_e
# 17. MgFe
# 18. MgFe_e
#
""")
    f.write(('#{:4n}{:4n}{:4n}'+'{:10n}'*14+'\n#\n').format(*np.arange(17)+1))
    
    return

def get_indices(pointing, basedir = '/d/monk/eigenbrot/WIYN/14B-0456/anal/indecies/D4000/py_sbands'):

    spy = '{}/NGC_891_P{}_bin30.msoz.spy.dat'.format(basedir,pointing)
    res = np.loadtxt(spy)

    return res

def get_rphi(pointing, basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'):
    
    rpfile = '{}/NGC_891_P{}_bin30_rphi.dat'.format(basedir, pointing)
    r, phi, dr = np.loadtxt(rpfile, usecols=(1,2,3), unpack=True)
    
    return r, phi, dr

def get_rhoz(pointing, basedir='/d/monk/eigenbrot/WIYN/14B-0456/anal/final_results'):

    loc = '{}/NGC_891_P{}_bin30_locations.dat'.format(basedir,pointing)
    r, z = np.loadtxt(loc, usecols=(4,5), unpack=True)

    return r, np.abs(z)

def main(output):

    fmt = '{:5n}' + '{:4n}'*2 + '{:10.3f}'*15 + '\n'
    
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
