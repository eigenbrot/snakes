import numpy as np
import ADEUtils as ADE
import pyfits
import matplotlib.pyplot as plt
import time
import glob
from scipy.interpolate import interp1d
import scipy.optimize as spo

num_ap = 300
pixelsize = 24e-3 #in mm
FL = 50. #in mm
X0 = np.array([])
X1 = np.array([])
CH = np.array([])

def sauron(direct_image, full_dir, datafile, dscale, lscale, numhats=100,
           exclude=[], threshold=0.0001):
    '''
    Description:
        Sauron, as the Lord of the Rings, is responsible for controlling
        the reconstruction of full beam bench data from a set of laser
        bench data. It breaks down a direct full beam image into a bunch of
        rings and, for each ring, constructs a power profile based on
        interpolation of the laser bench data.

    Inputs:
        direct_image - String
                       The name of the full beam, direct image from the 
                       full beam bench
        full_dir     - String
                       The directory containing the power profile data files 
                       from the laser bench
        datafile     - String
                       The name of the data file that is output from the 
                       laserBench moduel
        nonlinmin    - Float
                       The minimum angle that is not in the linear range
                       of the angle vs. radius plot
        nonlinmax    - Float
                       The maxmimum angle that is not in the linear range
                       of the angle vs. radius plot. Sauron will use all angels
                       less than nonlinmin and greater than nonlinmax to
                       find the systematic offset present in the data.
        exclude      - Python List
                       A list of input angles from the laser bench data to
                       exclude from analysis for whatever reason
        threshold    - Float
                       The value accepted to be zero intercept when
                       fitting a line to angle vs. radius data
                       The default is 0.0001.
                       
    Output:
        a        - Numpy array
                   The angle vector that all the profiles are based off of.
                   It's length is equal to num_ap.
        rings    - Numpy array
                   A N x num_ap array where N is the number of rings the
                   full beam image is broken up into. Each of the N vectors
                   is one power profile of the actual full beam, direct data.
                   i.e. all the vectors are top-hats.
        profiles - Numpy array
                   A N x num_ap array. Each vector is a simulated, smeared
                   power profile corresponding to a unsmeared profile in
                   rings.
        comb     - Numpy array
                   A num_ap length vector containing the total simulated,
                   smeared profile. It is made by summing all the vectors in
                   profiles.
        angs     - Numpy array
                   A N length vector that contains the output angles 
                   corresponding to the profiles in rings and profiles.
    '''

    'contruct the datebase for future interpolation'
    tm = theModule(full_dir,datafile,lscale,exclude=exclude,threshold=threshold)

    '''we use log spacing because we think there is more interesting stuff
    happening at smaller angles'''
#    minang = np.log10(tm.min_angle)
#    maxang = np.log10(tm.max_angle)
    minang = tm.min_angle
    maxang = tm.max_angle
    angs = np.linspace(0,maxang,num=numhats,endpoint=False)

    'break up the direct image'
    a, rings = disc2rings(direct_image,angs,dscale)
    
    'generate some smeared profiles'
    profiles = profile_gen(angs, tm, a)

    '''We now force each of the smeared profiles to have the same area as
    its corresponding direct-image profile. This makes sense because what
    we are simulating is the light from the top-hat being smeared to a wider
    range of angles so the total amount of light (power) should stay the same
    '''
    normalizations = np.array([np.sum(profiles,axis=1)]).T
    areas = np.array([np.sum(rings,axis=1)]).T

    profiles = profiles*areas/normalizations
    profiles = np.nan_to_num(profiles)

    'sum up all the smeared profiles'
    comb = np.sum(profiles,axis=0)
    
    return (a,rings,profiles,comb,angs)

def disc2rings(image,angles,dscale):
    '''takes a image and splits it into rings of constant thickness. The rings
    are centered at the angles specified in angles and the output rings
    are flattened azimuthally.'''

    HDU = pyfits.open(image)[0]
    data = np.float32(HDU.data)
    FL = HDU.header['FOCALLEN']
    length = num_ap

    r, f = ADE.annulize(data,length)

    print 'sum of data: '+str(np.sum(f))

    '''convert the pixel radii to output angle'''

    a = np.arctan(r*pixelsize/FL)*180/np.pi*dscale

    ring_stack = np.zeros((1,length),dtype=np.float32)

    sangles = np.sort(angles)

    for i in range(sangles.size):

        a_mid = sangles[i]

        '''the next few lines define the range of angles on either side of the
        requested angle we should use and deals with exceptions that happen
        at the smallest and largest angles.'''
        if i == 0: a0_mid = -1*sangles[i]
        else: a0_mid = sangles[i-1]

        try: a2_mid = sangles[i+1]
        except IndexError:a2_mid = a_mid + (a_mid - a0_mid)
        
        a1 = (a_mid + a0_mid)/2
        if i == 0: amin = a1
        a2 = (a_mid + a2_mid)/2

        idx = np.where((a > a1) & (a <= a2))
        counts = np.mean(f[idx])
        '''sometimes the range of angles is so small that the mean is nan.
        we change that to 0'''
        if np.isnan(counts): 
            counts = 0

        fi = np.zeros(length)
        
        fi[idx] = counts

        ring_stack = np.vstack((ring_stack,np.array([fi])))
    
    midx = np.where(a >= amin)[0][0]
    output_rings = ring_stack[1:]
    output_rings *= np.sum(f[midx:])/np.sum(output_rings)
    
    print 'sum of model: '+str(np.sum(output_rings))

    return (a,output_rings)

def profile_gen(angles, module, a):
    '''produces a bunch of interpolated profiles for the angles specified in
    angles using the data in module. It also forces the abcissa of the
    profiles to be a.'''

    length = a.size
    profile_stack = np.zeros((1,length),dtype=np.float32)
    
    sangles = np.sort(angles)

    for ang in sangles:
        'get the interpolated profile'
        tempa, tempp = module.get_profile(ang)

        'interpolate it once more to make it correspond to the angles in a'
        profile = np.interp(a,tempa,tempp)

        profile_stack = np.vstack((profile_stack,np.array([profile])))

    return profile_stack[1:]

def plot_data(im_name,title,clear,dscale,fnum=1):
    '''plots the power profile of some image in some specified figure number'''

    data = pyfits.open(im_name)[0].data

    r,f = ADE.annulize(data,num_ap)

    f /= np.sum(f)

    pangles = np.arctan(r*pixelsize/FL)*180./np.pi*dscale

    fig = plt.figure(fnum)
    if clear: plt.clf()
    ax = fig.add_subplot('111')
    ax.plot(pangles,f)
    ax.set_xlabel('Output Angle [deg]')
    ax.set_ylabel('Normalized Power')
    ax.set_xlim(0,18)
    ax.set_title(title)
    fig.show()

def plot_model(rings,a,title,clear,fnum):
    '''plots any of the data returned by sauron. It plots the data to figure
    number fnum and has the option of clearing the figure first or not.'''


    if len(rings.shape) == 1: rings = np.array([rings])

    fig = plt.figure(fnum)
    if clear: plt.clf()
    ax = fig.add_subplot('111')
    ax.set_xlabel('Output Angle [deg]')
    ax.set_ylabel('Power')
    ax.set_title(title)

    for i in range(rings.shape[0]):
        ax.plot(a,rings[i])

    plt.suptitle(time.asctime(time.localtime()))

    fig.show()
    return

def diff_trend(images, angle_vecs, power_vecs):
    
    fiber_powers = np.empty((1,300))
    for i in range(len(images)):
        data = pyfits.open(images[i])[0].data
        r, p = ADE.annulize(data,num_ap)
        p /= np.sum(p)
        pangles = np.arctan(r*pixelsize/FL)*180./np.pi
        ip = np.interp(angle_vecs[i],pangles,p)
        fiber_powers = np.vstack((fiber_powers,ip))

    fiber_powers = fiber_powers[1:]
    for power in power_vecs:
        power /= np.sum(power)

    print fiber_powers.shape
    print power_vecs.shape

    return np.sum(np.abs(fiber_powers - power_vecs),axis=1)

def scales(initial):
    
    global X0, X1, CH
#    X0 = X1 = CH = np.array([])

    fiber_list = ['300_1P.2P_RF3f.fits','300_1P.2P_RF42f.fits','300_1P.2P_RF63f.fits','300_1P.2P_RF135f.fits']
    
    radii = np.empty((1,num_ap))
    powers = np.empty((1,num_ap))

    for fiber in fiber_list:
        data = pyfits.open(fiber)[0].data
        r,p = ADE.annulize(data,num_ap)
        radii = np.vstack((radii,r))
        powers = np.vstack((powers,p/np.sum(p)))

    radii = radii[1:]
    powers = powers[1:]

    x0 = np.array(initial)

    print "starting minimization"
    xf = spo.fmin(func,x0,args=(radii,powers),xtol=0.0001,ftol=0.00001)

    return (xf)

def chimap(x0,x1):
    
    chi = np.empty((x0.size,x1.size))

    fiber_list = ['300_1P.2P_RF3f.fits','300_1P.2P_RF42f.fits','300_1P.2P_RF63f.fits','300_1P.2P_RF135f.fits']
    
    radii = np.empty((1,num_ap))
    powers = np.empty((1,num_ap))

    for fiber in fiber_list:
        data = pyfits.open(fiber)[0].data
        r,p = ADE.annulize(data,num_ap)
        radii = np.vstack((radii,r))
        powers = np.vstack((powers,p/np.sum(p)))

    radii = radii[1:]
    powers = powers[1:]
    
    for i in range(x0.size):
        print str(i)
        for j in range(x1.size):
            print ' '+str(j)
            chi[i,j] = func((x0[i],x1[j]),radii,powers)

    return x0,x1,chi

def func(x,radii,powers):

    direct_list = ['300_RF3d.fits','300_RF42d.fits','300_RF63d.fits','300_RF135d.fits']

    print "X0: "+str(x[0])+', X1: '+str(x[1])

    ang_stack = np.empty((1,num_ap))
    power_stack = np.empty((1,num_ap))
    diffs = np.array([])

    for direct in direct_list:
        a,_,_,comb,_ = sauron(direct,'full_better8.4.11','20110624_better8.4.11.dat',np.abs(x[0]),np.abs(x[1]),exclude=[-17.0,1.5])
        ang_stack = np.vstack((ang_stack,a))
        power_stack = np.vstack((power_stack,comb/np.sum(comb)))
        
    power_stack = power_stack[1:]
    ang_stack = ang_stack[1:]

    for i in range(len(direct_list)):
        pangles = np.arctan(radii[i]*pixelsize/FL)*180./np.pi*np.abs(x[1])
        ip = np.interp(ang_stack[i],pangles,powers[i])
        diffs = np.append(diffs,(ip - power_stack[i])**2)

        fig = plt.figure(9)
        plt.clf()
        ax = fig.add_subplot(111)
        ax.plot(a,ip,a,power_stack[i])
        fig.show()
        ch = np.sum(diffs)
        print "  chi: "+str(ch)
#        scratch = raw_input('Enter to continue...')

#    global X0, X1, CH
#    X0 = np.append(X0,x[0])
#    X1 = np.append(X1,x[1])
#    CH = np.append(CH,np.sum(diffs))

    return ch
        
def plot_all(x0,x1,sl=None):
    print x0, x1

    direct_list = ['300_RF3d.fits','300_RF42d.fits','300_RF63d.fits','300_RF135d.fits']    
    fiber_list = ['300_1P.2P_RF3f.fits','300_1P.2P_RF42f.fits','300_1P.2P_RF63f.fits','300_1P.2P_RF135f.fits']

    fig = plt.figure(5)
    plt.clf()
    for i in range(len(fiber_list)):
        
        a,_,_,comb,_ = sauron(direct_list[i],'full_better8.4.11',\
                                  '20110624_better8.4.11.dat',\
                                  x0,x1,exclude=[-17.0,1.5])

        if sl:
            print sl[i][0], sl[i][1]
            a2,_,_,comb2,_ = sauron(direct_list[i],'full_better8.4.11',\
                                      '20110624_better8.4.11.dat',\
                                      sl[i][0],sl[i][1],exclude=[-17.0,1.5])
    
        data = pyfits.open(fiber_list[i])[0].data
        r,f = ADE.annulize(data,num_ap)
        f /= np.sum(f)
        pangles = np.arctan(r*pixelsize/FL)*180./np.pi*x0
        if sl: pangles2 = np.arctan(r*pixelsize/FL)*180./np.pi*sl[i][0]
        

        ax = fig.add_subplot(2,2,i+1)
        ax.plot(pangles,f,a,comb/comb.sum())
        if sl: ax.plot(pangles2,f,'b--',a2,comb2/comb2.sum(),'g--')
        ax.set_xlabel('Output Angle [deg]')
        ax.set_ylabel('Power')
        ax.set_title(direct_list[i])

    plt.suptitle(time.asctime(time.localtime()))
    fig.show()
        
    return

##########################################################
########      theModule    ###############################
##########################################################


class theModule:
    '''This is THE Module, the final and ultimate in modular perfection. It 
    has the ability to build up a database of line profiles and then return
    a profile of _any_ line you could possibly want (that's within the range
    of the database data). What's more, theModule performs these tasks with
    a speed that incomprehensible to most sentients (man and machine).
    theModule is THE standard module.
    
    Lowely humans can only handle passing contact with theModule and are 
    advised to interact with theModule only through the get_profile function.
    theModule is THE standard module.
    '''

    def __init__(self, directory, datafile, lscale, exclude=[], threshold=0.00001):
        
        print "Initializing theModule... "

        self.clean_length = 2
        self.data = {}
        self.dir = directory
        self.file = datafile
        self.offset = self.find_offset(datafile,exclude)
        print 'offset is '+str(self.offset)+' degrees'
        
        file_list = glob.glob(directory+'/full*.dat')
        
        for dat_file in file_list:
            self.load_data(dat_file,exclude,lscale)

        self.angles = self.data.keys()
        self.max_angle = max(self.angles)
        self.min_angle = min(self.angles)

        print '...initializtion complete'

    def __repr__(self):
        return 'Module created using full data from '+self.dir\
            +'\nand reduced data from '+self.file\
            +'\ntheModule is THE standard module.'

    def find_offset_old(self,datafile, nonlinmin, nonlinmax, exclude, threshold):
        '''find_offset is used to determine the systematic offset present
        in the experimental setup that causes data to not be symmetric
        about zero input angle. It reads in the output of laserBench and
        returns the offset (in degrees)'''
        
        input_a, output_a = np.loadtxt(datafile,usecols=(0,1),unpack=True)
        
        for e in exclude:
            did = np.where(input_a == e)
            output_a = np.delete(output_a, did)
            input_a = np.delete(input_a, did)

        pidx = np.where(input_a > nonlinmax)
        nidx = np.where(input_a < nonlinmin)
        
        in_a = np.append(input_a[nidx],input_a[pidx])
        out_a = np.append(-1*output_a[nidx],output_a[pidx])
        error = np.zeros(in_a.size)+1

        b = 1000.
        offset = 0.
        while abs(b) > threshold:
            m, b = ADE.fit_line(in_a,out_a,error)
            offset += b
            in_a += b

        return offset

    def find_offset(self, datafile, exclude):
        
        input_a, output_a = np.loadtxt(datafile,usecols=(0,1),unpack=True)
        
        for e in exclude:
            did = np.where(input_a == e)
            output_a = np.delete(output_a, did)
            input_a = np.delete(input_a, did)

        initial = np.array([0.0])
        shift = spo.fmin(self.func_off,initial,args=(input_a,output_a),disp=False)

        return shift[0]

    def func_off(self, offset, in_a, out_a):
        
        input_a = in_a + offset[0]
        pidx = np.where(input_a >= 0.)
        nidx = np.where(input_a < 0.)
        
        pf = interp1d(input_a[pidx],out_a[pidx],kind='cubic')
        nf = interp1d(np.abs(input_a[nidx][::-1]),out_a[nidx][::-1]
                          ,kind='cubic')

        mina = max(np.min(input_a[pidx]),np.min(np.abs(input_a[nidx])))
        maxa = min(np.max(input_a[pidx]),np.max(np.abs(input_a[nidx])))

        a = np.linspace(mina,maxa,50)

        return np.sum((pf(a) - nf(a))**2)
        

    def load_data(self,dat_file,exclude,lscale):
        '''load_data loads data (surprise!) from a single data file and
        places it into the database of theModule. The data file should have
        2 vectors corresponding to the angle and power vectors of a power
        profile'''

        f = open(dat_file,'r')
        lines = f.readlines()
        input_angle = float(lines[1][15:-1]) + self.offset
        peak_ang = float(lines[2][16:-1])

        if input_angle in exclude: return

        f = open(dat_file,'r')        
        a,p = np.loadtxt(f,unpack=True)
#        print 'Input angle is '+str(input_angle)
#        print 'Output angle is '+str(peak_ang)+'\n'

        self.data[input_angle] = (a*lscale,p)
        return

    def get_profile(self, angle):
        '''
        Description:
            get_profile uses the database of power profiles stored in theModule
            to generate a profile at any arbitrary angle. This is the only
            function a user will likely use in theModule.

        Inputs:
            angle - Float
                    The angle you want a profile at. It must be between the
                    min and max angles contained in the database.

        Output:
            A tuple of numpy arrays. The first array is an angle of vectors
            and the second is the power profile corresponding to those angels.
            Both arrays have a length equal to the length of the arrays in
            the database.
        '''

        if (angle > self.max_angle) or (angle < self.min_angle):
            print "Angle "+str(angle)+" is not within range of database angles!"
            print "That range is "+str(self.min_angle)\
                +' to '+str(self.max_angle)
            return

        ang1, ang2, weight1, weight2 = self.find_bracket(angle)

        a1, p1 = self.data[ang1]
        a2, p2 = self.data[ang2]

        p1[np.where(p1 < 0)] = 0.0
        p2[np.where(p2 < 0)] = 0.0

        a1_z = a1 - ang1
        a2_z = a2 - ang2

        interpp = np.interp(a1_z, a2_z, p2)

        avgp = weight1*p1/p1.sum() + weight2*interpp/interpp.sum()

        shift = weight1*(ang1) + weight2*(ang2)
    
        outputp = np.interp(a1,a1_z + shift,avgp,left=0.0)
        
        if outputp[0] != 0.0: self.profile_clean(outputp,a1)

        return a1, outputp

    def profile_clean(self,profile,a):
        '''forces any generated profile to go to zero power at zero angle.
        This is done by replacing the first X points of the profile with
        a line from 0 to the power at point X+1, where X is set by
        theModule parameter clean_length.
        '''

        m = profile[self.clean_length]/(a[self.clean_length] - a[0])
        profile[0:self.clean_length+1] = \
            m*a[0:self.clean_length+1] - m*a[0]
        
    def find_bracket(self, angle):
        '''given some input angle, find_bracket returns the two angles
        in the database that bracket that angle.
        '''
        i = 0
        
        k = self.angles
        k.sort()

        while k[i] <= angle:
            if k[i] == angle: return k[i], k[i+1], 1., 0.
            i += 1

        ang2 = k[i]
        if ang2 == min(k): ang1 = k[i+1]
        else: ang1 = k[i-1]

        diff = abs(ang2 - ang1)

        weight1 = 1 - abs(angle - ang1)/diff
        weight2 = 1 - abs(ang2 - angle)/diff

        return ang1, ang2, weight1, weight2

    def find_peak(self, a, p):
        '''find the angles cooresponding to the peak in two different
        power profiles. Uses interpolation to get the angle at the peak
        not just near the peak.
        '''

        CDF = np.cumsum(p)
        CDF /= np.max(CDF)
        
        pidx = np.interp(0.5,CDF,np.arange(CDF.size))

        aout = np.interp(pidx,np.arange(CDF.size),a)

        return aout

    def test(self,angle,clear):
        '''
        Description:
            test does a quick and clean interpolation of some angle
            and plots the interpolated profile and the two bracketing
            profiles.

        Inputs:
            angle - Float
                    The angle you want a profile of
            clear - Boolean
                    If True the profile plot is cleared before plotting.
                    Otherwise the three profiles are overplotted onto whatever
                    other lines are out there.
                    
        Output:
            None
        '''
        ang1, ang2, _,_ = self.find_bracket(angle)
        
        print ang1, ang2
        
        fig = plt.figure(2)
        
        if clear:
            plt.clf()
            p1 = self.data[ang1][1]
            p2 = self.data[ang2][1]
            plt.plot(self.data[ang1][0],p1/p1.sum())
            plt.plot(self.data[ang2][0],p2/p2.sum())
            plt.xlabel('Output Angle [deg]')
            plt.ylabel('Normalized Power')
            
        a, p = self.get_profile(angle)
        plt.plot(a,p,'--')
        fig.show()
        
