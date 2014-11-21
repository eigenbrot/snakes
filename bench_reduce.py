#!/usr/bin/python

import sys
import os
import glob
glob = glob.glob

#Load the IRAF packages we'll need
try:
    current_dir = os.getcwd()
    if os.getlogin() == 'Arthur':
            os.chdir('/Users/Arthur/Ureka/iraf/local')
    from pyraf import iraf
    os.chdir(current_dir)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
except Exception as e:
    print "Failure: could not find pyraf/iraf"
    sys.exit(1)

def make_lists(raw_dir):

    all_raw = glob('{}/*.fits'.format(raw_dir))
    all_ot = ['{}_ot.fits'.format(os.path.basename(i).\
                                  split('.fits')[0]) for i in all_raw]
    comps = [i for i in all_ot if 'comp' in i]
    darks = [i for i in all_ot if 'dark' in i]
    flats = [i for i in all_ot if 'dflat' in i]
    zero = [i for i in all_ot if 'zero' in i]
    all_otz = [i for i in all_ot if 'zero' not in i]

    outdict = {'all_raw':all_raw,
               'all_ot':all_ot,
               'all_otz':all_otz,
               'comps':comps,
               'darks':darks,
               'zeros':zero,
               'flats':flats}

    return outdict

def write_list(inlist, name):
    
    with open(name,'w') as f:
        for i in inlist:
            f.write('{}\n'.format(i))

    return

def setup_dir(raw_dir):

    list_dict = make_lists(raw_dir)
    for k in list_dict.keys():
        if k == 'flats':
            continue
        write_list(list_dict[k],'{}.lst'.format(k))

    return list_dict

def cull_bias(biaslist,listname):

    iraf.imstat('@{}'.format(listname))
    bads = raw_input('Select any biasi you want to remove, q to continue\n')
    while bads.lower() != 'q':

        badlist = bads.split(',')
        print 'badlist is {}'.format(badlist)
        newlist = [i for i in biaslist if not any([j.strip() in i for j in badlist])]
        write_list(newlist,listname)
        iraf.imstat('@{}'.format(listname))
        bads = raw_input('Select any biasi you want to remove, q to continue\n')

    return

def detect_flats(directory):

    flatlist = glob('{}/dflat*.fits'.format(directory))
    explist = [os.path.basename(i).split('_')[1] for i in flatlist]
    unique_exp = set(explist)
    
    outlist = []
    print 'detecting flats...'
    for exp in unique_exp:
        print '\tfound {}'.format(exp)
        tmplist = [i for i in flatlist if exp in i]
        tmpname = '{}/dflats_{}.lst'.format(directory,exp)
        write_list(tmplist,tmpname)
        outlist.append(tmpname)
        
    return

def main(raw_dir):

    if not os.path.exists('{}/all_raw.lst'.format(os.getcwd())):
        print 'Setting up directory'
        list_dict = setup_dir(raw_dir)
        
        print "OT'ing"
        iraf.ccdproc('@all_raw.lst',output='@all_ot.lst',overscan=True,trim=True,
                     zerocor=False,darkcor=False,flatcor=False,illumco=False,
                     biassec='image',trimsec='image',zero='',interact=True,order=100)

        cull_bias(list_dict['zeros'],'zeros.lst')
        return
        
    print 'Making master bias'
    iraf.zerocombine('@zeros.lst',output='Zero.fits',combine='average',reject='crreject')
    print 'Subtracting bias'
    iraf.ccdproc('@all_otz.lst',overscan=False,trim=False,zerocor=True,zero='Zero.fits')

    detect_flats(os.getcwd())
    flats_to_make = glob('dflat*.lst')
    for flat in flats_to_make:
        name = 'dFlat_{}.fits'.format(flat.split('_')[1])
        print 'Making flat {}'.format(name)
        iraf.flatcombine('@{}'.format(flat),output=name,combine='average',reject='crreject')
        
    print 'Making master dark'
    iraf.darkcombine('@darks.lst',output='Dark.fits',combine='average',reject='crreject')
    print 'Making master comp'
    iraf.imcombine('@comps.lst','Comp.fits',combine='average',reject='crreject')

    return

if __name__ == '__main__':

    try:
        if not os.path.exists(sys.argv[1]):
            print 'The request was made but it was not good'
            sys.exit(1)
    except Exception as e:
        print 'The request was made but it was not good'
        sys.exit(1)

    print "Let's do it {}".format(sys.argv[1])
    main(sys.argv[1])
    
