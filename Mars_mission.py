import numpy as np
from datetime import datetime

def get_distance(time):

    M_earth = 5.97e24 #kg
    R_earth = 6371e3 #m
    v_0 = 4870 + 7784 #m/s
    G = 6.67e-11 #m^3 kg^-1 s^-2
    
    k = 0.5*G*M_earth*time**2
    v = v_0 * time
    
    print k * (27 * k - 4 * (R_earth + v)**3)
    
    r_1 = 2 * (-2*R_earth + v) + \
      (\
       2 * (R_earth + v)**2\
       )/\
       (-13.5 * k + (R_earth + v)**3 + 2.598 * np.sqrt(k * (27 * k + 4 * (R_earth + v)**3)))**(1/3)

    r_2 = 1.5874 * (-27 * k + 2 * (R_earth + v)**3 + 5.19615 * \
                    np.sqrt(k * (27 * k + 4 * (R_earth + v)**3)))**(1/3)

    r = 0.1666 * (r_1 + r_2)

    return r

def get_distance2(time):

    g = -9.8 #m/s^2
    v_0 = 4870 #m/s

    r = 0.5 * g * time**2 + v_0 * time

    return r

def get_distance3(time, g_tup, thresh = 1.0):

    R_earth = 6371e3 #m
    v_0 = 4870 + 7784 #m/s
    g_0 = 0 #m s^-2

    r_0 = 0.5 * g_0 * time**2 + v_0 * time
    diff = np.inf

    while diff > thresh:

        g_i = np.interp(r_0 + R_earth,g_tup[0],g_tup[1])
        r_i = 0.5 * g_i * time**2 + v_0 * time

        print '{:11.4f}, {:11.4f}'.format(r_0,r_i)
        diff = np.abs(r_i - r_0)
        r_0 = r_i

    return r_0

def get_distance4(end_time, g_tup):

    R_earth = 6371e3 #m
    v_0 = 4870 + 7784 #m/s
    g_0 = -9.8 #m s^-2

    times = np.arange(0,end_time,60.)
    distances = np.zeros(times.shape)
    velocities = np.zeros(times.shape) + v_0
    r_0 = 0.0
    
    for i in range(times.size - 1):

        v = v_0 + g_0 * 30.
        velocities[i+1] = v
        print v
        distances[i+1] = distances[i] + v * 60
        g_0 = np.interp(distances[i+1],g_tup[0],g_tup[1])
        v_0 = v

    return times, distances, velocities

def angular_size(distance):

    R_earth = 6371e3 #m

    theta_rad = 2*np.arctan(R_earth/distance) #radians

    theta_arc = theta_rad * 206265.

    return theta_arc

def do_run(output):

    M_earth = 5.97e24 #kg
    R_earth = 6371e3 #m
    gdist = np.arange(0,R_earth*1000,1000)
    g = -6.67e-11 * M_earth / (gdist + R_earth)**2
    
    moon_size = 1800. #in arcsec

    des_times = np.arange(0,2400*2,12)*3600 #every 12 hours in seconds

    times, distances, velocities = get_distance4(des_times.max(),(gdist,g))

    print "here"
    des_dist = np.interp(des_times,times,distances)
    des_velo = np.interp(des_times,times,velocities)

    print "there"
    sizes = angular_size(des_dist)

    print "elsewhere"
    moon_rel = sizes/moon_size
    
    f = open(output,'w')
    f.write('# File generated on {}\n'.format(datetime.now().isoformat(' ')) + \
            '#\n# The angular size of the Moon view from Earth is 1800 arcseconds\n' + \
            '# The angular resolution of a typical human eye is 17 - 50 arcseconds\n' + \
            '#  depending on the age and viewing conditions\n' + \
            '#\n#{:>11} = Time since departure from Low Earth Orbit [days]\n'.format('t') + \
            '#{:>11} = Distance from surface of Earth [km]\n'.format('d') + \
            '#{:>11} = Velocity [m/s]\n'.format('v') + \
            '#{:>11} = Angular size of Earth [arcsecond]\n'.format('theta_E') + \
            '#{:>11} = Angular size of Earth relative to the Moon viewed from Earth\n'.format('theta_rel') + \
            '#\n#{:>13}{:>13}{:>13}{:>13}{:>13}\n'.format('t','d','v','theta_E','theta_rel') + \
            '#{:13}{:13}{:13}{:13}{:13}\n'.format(*range(5))
            )

    for t, d, v, s, rel in zip(des_times,des_dist,des_velo,sizes,moon_rel):

        f.write('{:14.1f}{:13.3e}{:13.3f}{:13.3f}{:13.5f}\n'.format(t/3600/24.,d/1000.,v,s,rel))

    f.close()

    return
