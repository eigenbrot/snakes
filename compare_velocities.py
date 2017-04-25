import numpy as np
import pyfits

def compare(file1, file2):

    if file1[-5:] == ".fits":
        vel1 = read_fits(file1)
    elif file1[-4:] == ".dat":
        vel1 = read_dat(file1)
    else:
        print "WARNING: Could not determine file type of", file1
        return

    if file2[-5:] == ".fits":
        vel2 = read_fits(file2)
    elif file2[-4:] == ".dat":
        vel2 = read_dat(file2)
    else:
        print "WARNING: Could not determine file type of", file2
        return

    diff = vel1 - vel2
    print file1, file2
    print '\t', np.mean(diff), np.std(diff)
    return

def read_fits(filename):

    chifile = filename[:-5] + '.vel.fits'
    try:
        chivel = pyfits.open(chifile)[1].data['VSYS']
    except IOError:
        print "WARNING: Could not find chivel file associated with", filename
        print "         continuing without a chivel file"
        chivel = 0
        
    vel = pyfits.open(filename)[1].data['VSYS']

    return vel + chivel

def read_dat(filename):

    vel = np.loadtxt(filename, usecols=(1,), unpack=True)

    return vel

if __name__ == "__main__":
    import sys
    compare(sys.argv[1], sys.argv[2])
