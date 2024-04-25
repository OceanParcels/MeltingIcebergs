import math


#============================== General kernels ==============================#
def MinimumDistance(particle, fieldset, time):
    '''
    Calculate the minimum distance of the iceberg to ODP Site 696.
    '''
    
    lat_dist = math.fabs(particle.lat - fieldset.ODP_lat) * (6371e3 * (math.pi/180.))   # vertical distance [m]
    lon_dist = math.fabs(particle.lon - fieldset.ODP_lon) * ((6371e3 * math.cos(
        ((particle.lat + fieldset.ODP_lat)/2.)*(math.pi/180.))) * (math.pi/180.)) # horizontal distance [m]
    
    dist = math.sqrt(lat_dist**2 + lon_dist**2)                              # total distance [m]
    if dist < particle.prev_dist:                                            # find minimum distance
        particle.distance = dist
    else:
        particle.distance = particle.prev_dist
    particle.prev_dist = particle.distance


def SetInitialPropertiesC1(particle, fieldset, time):
    '''
    Set the initial properties of the iceberg using:
        W = L/1.5
        M = L * W * T * rho_ice
        D = T * (rho_ice/rho_ocean)
    for size class C1.
    '''
    
    if particle.check < 0.:
        ### Rackow (2017) sets the minimal iceberg volume to 1000m3 (10x10x10)
        particle.L = 17.                             # length [m]
        particle.W = particle.L/1.5                  # width [m]
        particle.T = 12.2                            # thickness [m]
        
        particle.M = particle.L * particle.W * particle.T * fieldset.rho_i # mass [kg]
        
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # draft [m]

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T
        
        particle.check = 1.
    
    particle.prev_check = particle.check


def SetInitialPropertiesC2(particle, fieldset, time):
    '''
    Set the initial properties of the iceberg using:
        W = L/1.5
        M = L * W * T * rho_ice
        D = T * (rho_ice/rho_ocean)
    for size class C2.
    '''
    
    if particle.check < 0.:
        particle.L = 100.                            # length [m]
        particle.W = particle.L/1.5                  # width [m]
        particle.T = 20.                             # thickness [m]
        
        particle.M = particle.L * particle.W * particle.T * fieldset.rho_i # mass [kg]
        
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # draft [m]

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T
        
        particle.check = 1.
    
    particle.prev_check = particle.check
    

def SetInitialPropertiesC3(particle, fieldset, time):
    '''
    Set the initial properties of the iceberg using:
        W = L/1.5
        M = L * W * T * rho_ice
        D = T * (rho_ice/rho_ocean)
    for size class C3.
    '''
    
    if particle.check < 0.:
        particle.L = 1000.                           # length [m]
        particle.W = particle.L/1.5                  # width [m]
        particle.T = 200.                            # thickness [m]
        
        particle.M = particle.L * particle.W * particle.T * fieldset.rho_i # mass [kg]
        
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # draft [m]

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T
        
        particle.check = 1.
    
    particle.prev_check = particle.check


def SetInitialPropertiesC4(particle, fieldset, time):
    '''
    Set the initial properties of the iceberg using:
        W = L/1.5
        M = L * W * T * rho_ice
        D = T * (rho_ice/rho_ocean)
    for size class C4.
    '''
    
    if particle.check < 0.:
        particle.L = 10000.                          # length [m]
        particle.W = particle.L/1.5                  # width [m]
        particle.T = 500.                            # thickness [m]
        
        particle.M = particle.L * particle.W * particle.T * fieldset.rho_i # mass [kg]
        
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # draft [m]

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T
        
        particle.check = 1.
    
    particle.prev_check = particle.check


def SetInitialPropertiesC5(particle, fieldset, time):
    '''
    Set the initial properties of the iceberg using:
        W = L/1.5
        M = L * W * T * rho_ice
        D = T * (rho_ice/rho_ocean)
    for size class C5.
    '''
    
    if particle.check < 0.:
        particle.L = 100000.                         # length [m]
        particle.W = particle.L/1.5                  # width [m]
        particle.T = 500.                            # thickness [m]
        
        particle.M = particle.L * particle.W * particle.T * fieldset.rho_i # mass [kg]
        
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # draft [m]

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T
        
        particle.check = 1.
    
    particle.prev_check = particle.check


def SampleRegion(particle, fieldset, time):
    '''
    At the iceberg's location, sample the region (see Carter et al., 2017)
    '''
    
    particle.reg = fieldset.R[time+(particle.dt/2.), 0., particle.lat, particle.lon]


def NewProperties(particle, fieldset, time):
    '''
    Calculates the new properties of the iceberg from the melt terms.
    '''
    
    if particle.check > 0.:
        # Calculate old thickness of iceberg draft
        draft = particle.prev_T * (fieldset.rho_i/fieldset.rho_o)

        # Calculate working areas for ...
        Ab = particle.prev_L * particle.prev_W                          # ... basal melt [m3]
        Av = 2. * draft * (particle.prev_L + particle.prev_W)           # ... buoyant convection on all sides below SL [m3]
        Ae = particle.prev_T * (particle.prev_L + particle.prev_W)      # ... wave erosion on two sides [m3]

        # Volume changes
        dV = (particle.prev_Mbr*Ab + particle.prev_Mvr*Av + particle.prev_Mer*Ae) * particle.dt   # total volume loss [m3]
        Vol = (particle.prev_L * particle.prev_W * particle.prev_T) - dV                          # new volume [m3]
        if Vol < 0.:
            particle.delete()

        # New dimensions
        particle.T = particle.prev_T - particle.prev_Mbr*particle.dt    # calculate new thickness directly [m]
        A = Vol/particle.T                                              # horizontal area iceberg [m2]
        particle.W = math.sqrt(A/1.5)                                   # new width follows from ratio L = 1.5W [m]
        particle.L = 1.5 * particle.W                                   # new length follows from width as L = 1.5W [m]
        particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # new draft [m]

        # New mass
        if particle.L >= 0 and particle.W >= 0 and particle.T >= 0:
            particle.M = particle.T * particle.W * particle.L * fieldset.rho_i   # iceberg mass [kg]

            ep_crit = math.sqrt(6. * (fieldset.rho_i/fieldset.rho_o) * (1. - (fieldset.rho_i/fieldset.rho_o)))
            ep = particle.W/particle.T

            if ep < ep_crit:                    # iceberg tipping (WDE17)
                particle.T = particle.W         # calculate new thickness [m]
                A = Vol/particle.T              # keep volume and ratios correct
                particle.W = math.sqrt(A/1.5)   # calculate new width [m]
                particle.L = 1.5 * particle.W   # calculate new length [m]
                particle.depth = particle.T * (fieldset.rho_i/fieldset.rho_o)   # new draft [m]
        else:
            particle.delete()

        particle.prev_L = particle.L
        particle.prev_W = particle.W
        particle.prev_T = particle.T



#========================= Kernels depth-integration =========================#
def AdvectionEE(particle, fieldset, time):
    '''
    Advection of icebergs using Explicit Euler (aka Euler Forward) integration
    using the depth-integrated ocean velocity components.
    '''
    
    particle.lon += particle.uveloD/(1852*60*math.cos(particle.lat*math.pi/180.)) * particle.dt
    particle.lat += particle.vveloD/(1852*60) * particle.dt

    
def Grounding(particle, fieldset, time):
    '''
    Set the depth-integrated velocity to zero when an iceberg has run aground.
    '''
    
    if particle.depth >= math.fabs(particle.bath):
        particle.uveloD = 0.
        particle.vveloD = 0.


def SampleFields(particle, fieldset, time):
    '''
    At the iceberg's location, sample:
     1) the surface wind field,
     2) the ocean bathymetry,
     3) the surface ocean temperature and velocity,
     4) the ocean temperature and velocity at the iceberg base, and
     5) the depth-integrated ocean temperature and velocity along the iceberg's draft.
    '''
    
    ### Wind: surface wind stress components and surface wind speed
    taux, tauy = fieldset.XY[time+(particle.dt/2.), 0., particle.lat, particle.lon]       # [g/(s2 cm)]    
    vela = math.sqrt((math.sqrt(taux**2+tauy**2)*0.1)/(fieldset.rho_a*0.0015))            # [m/s]
    particle.uvela = (taux*0.1)/(fieldset.rho_a*0.0015*math.fabs(vela))                   # [m/s]
    particle.vvela = (tauy*0.1)/(fieldset.rho_a*0.0015*math.fabs(vela))                   # [m/s]
    
    ### Bathymetry
    particle.bath = fieldset.B[time+(particle.dt/2.), particle.depth, particle.lat, particle.lon] # [m]
    
    ## Determine the deepest point between the iceberg's base and the bathymetry
    z = min(particle.depth,math.fabs(particle.bath))
    
    ### Surface (ocean): ocean temperature and velocity components at the surface
    particle.tempS = fieldset.T[time+(particle.dt/2.), 0., particle.lat, particle.lon]    # [°C]
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), 0., particle.lat, particle.lon]       # [°/s]
    particle.uveloS = uvel*1852*60*math.cos(particle.lat*math.pi/180.)                    # [m/s]
    particle.vveloS = vvel*1852*60                                                        # [m/s]
    
    ### Basal (ocean): ocean temperature and velocity components at the iceberg base
    particle.tempB = fieldset.T[time+(particle.dt/2.), z, particle.lat, particle.lon]     # [°C]
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), z, particle.lat, particle.lon]        # [°/s]
    particle.uveloB = uvel*1852*60*math.cos(particle.lat*math.pi/180.)                    # [m/s]
    particle.vveloB = vvel*1852*60                                                        # [m/s]
    
    ### Integrated (ocean): depth-integrated ocean temperature and velocity components along the iceberg's draft
    t, u, v = 0, 0, 0
    check = 1

    if z >= 10.01244 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 5.00622, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 5.00622, particle.lat, particle.lon]
        t += temp * 10.01244
        u += uv * 10.01244
        v += vv * 10.01244
    elif z < 10.01244 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 5.00622, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 5.00622, particle.lat, particle.lon]
        t += temp * (z-0.0)
        u += uv * (z-0.0)
        v += vv * (z-0.0)
        check = -1
    
    if z >= 20.12502 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 15.068729, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 15.068729, particle.lat, particle.lon]
        t += temp * (20.12502-10.01244)
        u += uv * (20.12502-10.01244)
        v += vv * (20.12502-10.01244)
    elif z < 20.12502 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 15.068729, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 15.068729, particle.lat, particle.lon]
        t += temp * (z-10.01244)
        u += uv * (z-10.01244)
        v += vv * (z-10.01244)
        check = -1
    
    if z >= 30.44184 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 25.283428, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 25.283428, particle.lat, particle.lon]
        t += temp * (30.44184-20.12502)
        u += uv * (30.44184-20.12502)
        v += vv * (30.44184-20.12502)
    elif z < 30.44184 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 25.283428, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 25.283428, particle.lat, particle.lon]
        t += temp * (z-20.12502)
        u += uv * (z-20.12502)
        v += vv * (z-20.12502)
        check = -1
    
    if z >= 41.075138 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 35.758488, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 35.758488, particle.lat, particle.lon]
        t += temp * (41.075138-30.44184)
        u += uv * (41.075138-30.44184)
        v += vv * (41.075138-30.44184)
    elif z < 41.075138 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 35.758488, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 35.758488, particle.lat, particle.lon]
        t += temp * (z-30.44184)
        u += uv * (z-30.44184)
        v += vv * (z-30.44184)
        check = -1
    
    if z >= 52.150238 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 46.612686, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 46.612686, particle.lat, particle.lon]
        t += temp * (52.150238-41.075138)
        u += uv * (52.150238-41.075138)
        v += vv * (52.150238-41.075138)
    elif z < 52.150238 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 46.612686, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 46.612686, particle.lat, particle.lon]
        t += temp * (z-41.075138)
        u += uv * (z-41.075138)
        v += vv * (z-41.075138)
        check = -1
    
    if z >= 63.811737 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 57.980988, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 57.980988, particle.lat, particle.lon]
        t += temp * (63.811737-52.150238)
        u += uv * (63.811737-52.150238)
        v += vv * (63.811737-52.150238)
    elif z < 63.811737 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 57.980988, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 57.980988, particle.lat, particle.lon]
        t += temp * (z-52.150238)
        u += uv * (z-52.150238)
        v += vv * (z-52.150238)
        check = -1
    
    if z >= 76.23103 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 70.021385, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 70.021385, particle.lat, particle.lon]
        t += temp * (76.23103-63.811737)
        u += uv * (76.23103-63.811737)
        v += vv * (76.23103-63.811737)
    elif z < 76.23103 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 70.021385, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 70.021385, particle.lat, particle.lon]
        t += temp * (z-63.811737)
        u += uv * (z-63.811737)
        v += vv * (z-63.811737)
        check = -1
    
    if z >= 89.617134 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 82.92409, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 82.92409, particle.lat, particle.lon]
        t += temp * (89.617134-76.23103)
        u += uv * (89.617134-76.23103)
        v += vv * (89.617134-76.23103)
    elif z < 89.617134 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 82.92409, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 82.92409, particle.lat, particle.lon]
        t += temp * (z-76.23103)
        u += uv * (z-76.23103)
        v += vv * (z-76.23103)
        check = -1
    
    if z >= 104.23113 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 96.92413, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 96.92413, particle.lat, particle.lon]
        t += temp * (104.23113-89.617134)
        u += uv * (104.23113-89.617134)
        v += vv * (104.23113-89.617134)
    elif z < 104.23113 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 96.92413, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 96.92413, particle.lat, particle.lon]
        t += temp * (z-89.617134)
        u += uv * (z-89.617134)
        v += vv * (z-89.617134)
        check = -1
        
    if z >= 120.40673 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 112.31893, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 112.31893, particle.lat, particle.lon]
        t += temp * (120.40673-104.23113)
        u += uv * (120.40673-104.23113)
        v += vv * (120.40673-104.23113)
    elif z < 120.40673 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 112.31893, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 112.31893, particle.lat, particle.lon]
        t += temp * (z-104.23113)
        u += uv * (z-104.23113)
        v += vv * (z-104.23113)
        check = -1
    
    if z >= 138.58043 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 129.49358, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 129.49358, particle.lat, particle.lon]
        t += temp * (138.58043-120.40673)
        u += uv * (138.58043-120.40673)
        v += vv * (138.58043-120.40673)
    elif z < 138.58043 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 129.49358, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 129.49358, particle.lat, particle.lon]
        t += temp * (z-120.40673)
        u += uv * (z-120.40673)
        v += vv * (z-120.40673)
        check = -1
    
    if z >= 159.33603 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 148.95822, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 148.95822, particle.lat, particle.lon]
        t += temp * (159.33603-138.58043)
        u += uv * (159.33603-138.58043)
        v += vv * (159.33603-138.58043)
    elif z < 159.33603 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 148.95822, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 148.95822, particle.lat, particle.lon]
        t += temp * (z-138.58043)
        u += uv * (z-138.58043)
        v += vv * (z-138.58043)
        check = -1
    
    if z >= 183.47282 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 171.40442, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 171.40442, particle.lat, particle.lon]
        t += temp * (183.47282-159.33603)
        u += uv * (183.47282-159.33603)
        v += vv * (183.47282-159.33603)
    elif z < 183.47282 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 171.40442, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 171.40442, particle.lat, particle.lon]
        t += temp * (z-159.33603)
        u += uv * (z-159.33603)
        v += vv * (z-159.33603)
        check = -1
    
    if z >= 212.11102 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 197.79192, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 197.79192, particle.lat, particle.lon]
        t += temp * (212.11102-183.47282)
        u += uv * (212.11102-183.47282)
        v += vv * (212.11102-183.47282)
    elif z < 212.11102 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 197.79192, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 197.79192, particle.lat, particle.lon]
        t += temp * (z-183.47282)
        u += uv * (z-183.47282)
        v += vv * (z-183.47282)
        check = -1
    
    if z >= 246.85742 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 229.48422, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 229.48422, particle.lat, particle.lon]
        t += temp * (246.85742-212.11102)
        u += uv * (246.85742-212.11102)
        v += vv * (246.85742-212.11102)
    elif z < 246.85742 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 229.48422, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 229.48422, particle.lat, particle.lon]
        t += temp * (z-212.11102)
        u += uv * (z-212.11102)
        v += vv * (z-212.11102)
        check = -1
    
    if z >= 290.066 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 268.46173, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 268.46173, particle.lat, particle.lon]
        t += temp * (290.066-246.85742)
        u += uv * (290.066-246.85742)
        v += vv * (290.066-246.85742)
    elif z < 290.066 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 268.46173, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 268.46173, particle.lat, particle.lon]
        t += temp * (z-246.8574)
        u += uv * (z-246.8574)
        v += vv * (z-246.8574)
        check = -1
    
    if z >= 345.2341 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 317.6501, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 317.6501, particle.lat, particle.lon]
        t += temp * (345.2341-290.066)
        u += uv * (345.2341-290.066)
        v += vv * (345.2341-290.066)
    elif z < 345.2341 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 317.6501, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 317.6501, particle.lat, particle.lon]
        t += temp * (z-290.066)
        u += uv * (z-290.066)
        v += vv * (z-290.066)
        check = -1
    
    if z >= 417.5387 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 381.38644, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 381.38644, particle.lat, particle.lon]
        t += temp * (417.5387-345.2341)
        u += uv * (417.5387-345.2341)
        v += vv * (417.5387-345.2341)
    elif z < 417.5387 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 381.38644, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 381.38644, particle.lat, particle.lon]
        t += temp * (z-345.2341)
        u += uv * (z-345.2341)
        v += vv * (z-345.2341)
        check = -1
    
    if z >= 514.2877 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 465.91324, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 465.91324, particle.lat, particle.lon]
        t += temp * (514.2877-417.5387)
        u += uv * (514.2877-417.5387)
        v += vv * (514.2877-417.5387)
    elif z < 514.2877 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 465.91324, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 465.91324, particle.lat, particle.lon]
        t += temp * (z-417.5387)
        u += uv * (z-417.5387)
        v += vv * (z-417.5387)
        check = -1
    
    if z >= 644.3267 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 579.30725, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 579.30725, particle.lat, particle.lon]
        t += temp * (644.3267-514.2877)
        u += uv * (644.3267-514.2877)
        v += vv * (644.3267-514.2877)
    elif z < 644.3267 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 579.30725, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 579.30725, particle.lat, particle.lon]
        t += temp * (z-514.2877)
        u += uv * (z-514.2877)
        v += vv * (z-514.2877)
        check = -1
    
    if z >= 814.37573 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 729.35126, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 729.35126, particle.lat, particle.lon]
        t += temp * (814.37573-644.3267)
        u += uv * (814.37573-644.3267)
        v += vv * (814.37573-644.3267)
    elif z < 814.37573 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 729.35126, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 729.35126, particle.lat, particle.lon]
        t += temp * (z-644.3267)
        u += uv * (z-644.3267)
        v += vv * (z-644.3267)
        check = -1
    
    particle.tempD = t / z                                               # [°C]
    uvel = u / z                                                         # [°/s]
    vvel = v / z                                                         # [°/s]
    
    particle.uveloD = uvel*1852*60*math.cos(particle.lat*math.pi/180.)   # [m/s]
    particle.vveloD = vvel*1852*60                                       # [m/s]


def BuoyantConvection(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through buoyant convection at the icebergs sides following:
        Mv = b1 * T + b2 * T^2
    with
        b1 = 7.62E-3 [m/(d °C)]
        b2 = 1.29E-3 [m/(d °C)]
        T: depth-integrated temperature [°C]
    See also Merino et al. (2016).
    '''
    
    if particle.depth < math.fabs(particle.bath):                        ## if not grounded
        Mvr = 7.62e-3 * particle.tempD + 1.29e-3 * (particle.tempD)**2   # [m/d]
    else:                                                                ## if grounded
        # meltrate till seafloor
        Mv1 = 7.62e-3 * particle.tempD + 1.29e-3 * (particle.tempD)**2
        Mvr1 = Mv1 * math.fabs(particle.bath)
        # meltrate 'below' seafloor
        Mv2 = 7.62e-3 * particle.tempB + 1.29e-3 * (particle.tempB)**2
        Mvr2 = Mv2 * (particle.depth-math.fabs(particle.bath))
        # total meltrate
        Mvr = (Mvr1+Mvr2)/particle.depth                                 # [m/d]
    
    particle.Mvr = Mvr/fieldset.sec_to_day                               # [m/s]
    particle.prev_Mvr = particle.Mvr


def BasalMelt(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through basal melt at the iceberg base following:
        Mb = C * |v_ice - v_ocean|^0.8 * (T_ocean - T_ice) / (L^0.2) 
    with
        C = 0.58 [°C^-1 m^0.4 d^-1 s^0.8]
        v_ice: iceberg velocity (depth-integrated) [m/s]
        v_ocean: ocean velocity at the iceberg base [m/s]
        T_ocean: ocean temperature at the iceberg base [°C]
        T_ice = -1.92 [°C]; freezing temperature
        L: iceberg length [m]
    See FitzMaurice & Stern (2018) and references therein.
    '''
    
    dvelabs = math.sqrt((particle.uveloD-particle.uveloB)**2 + (particle.vveloD-particle.vveloB)**2)
    
    Mbr = 0.58 * (dvelabs**0.8) * ((particle.tempB-fieldset.Ti)/(particle.L**0.2))   # [m/d]
    particle.Mbr = Mbr/fieldset.sec_to_day                                           # [m/s]
    particle.prev_Mbr = particle.Mbr


def BasalMeltAlt(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through basal melt at the iceberg base following the bulk
    equation when its length is smaller than 1E2 m and the adapted equation when larger.
    '''
    
    if particle.L <= 15000.:
        dvelabs = math.sqrt((particle.uveloD-particle.uveloB)**2 + (particle.vveloD-particle.vveloB)**2)

        Mbr = 0.58 * (dvelabs**0.8) * ((particle.tempB-fieldset.Ti)/(particle.L**0.2))   # [m/d]
        particle.Mbr = Mbr/fieldset.sec_to_day                                           # [m/s]
        particle.prev_Mbr = particle.Mbr
    else:    
        f = math.fabs(2 * fieldset.Om * math.sin(particle.lat*math.pi/180.))
        u_st = math.sqrt(fieldset.cd) * math.sqrt(particle.uveloD**2 + particle.vveloD**2)
        eta_st = 1./math.sqrt(1 + fieldset.xiN*u_st/(f*fieldset.Lo*fieldset.Rc))
        h_nu = 5*fieldset.nu/u_st

        GammaT = fieldset.k_r*math.log(
            u_st*fieldset.xiN*eta_st**2/(f*h_nu)) + 1./(
            2*fieldset.xiN*eta_st) - fieldset.k_r                                                   # [-]
        GammaM = 12.5 * fieldset.Pr**(2./3.) - 6                                                    # [-]
        gamma_T3EM = (fieldset.rho_o*fieldset.cpo*math.sqrt(fieldset.cd)*(
            0.004*particle.tempS+0.02)) / (GammaT + GammaM)                                         # [W/(m2 K)]

        particle.Mbr = gamma_T3EM * (0.7*particle.tempS-1.2-fieldset.Ti) / (fieldset.rho_i * fieldset.Lf)   # [m/s]
        particle.prev_Mbr = particle.Mbr


def WaveErosion(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through wave erosion at the icebergs sides:
        Me = 1/12 Ss * (1+cos(pi*A_ice^3)) * (T_ocean + 2)
    with sea state:
        Ss = 3/2 (|v_atm - v_ocean|)^0.5 + 1/10 |v_atm - v_ocean|
    and
        A_ice: fractional sea-ice area [m2]; negligible in the Eocene
        T_ocean: ocean surface temperature [°C]
        v_atm: surface wind velocity [m/s]
        v_ocean: ocean surface velocity [m/s]
    '''
    
    dvelabs = math.sqrt((particle.uvela-particle.uveloS)**2 + (particle.vvela-particle.vveloS)**2)
    Ss = (3./2.) * math.sqrt(dvelabs) + (1./10.) * dvelabs
    
    Ai = 0.
    damping = 0.5 * (1 + math.cos(math.pi * Ai**3))
    
    Mer = (1./6.) * Ss * damping * (particle.tempS - - 2.)   # [m/d]
    particle.Mer = Mer/fieldset.sec_to_day                   # [m/s]
    particle.prev_Mer = particle.Mer


def SampleFieldsModern(particle, fieldset, time):
    '''
    At the iceberg's location, sample:
     1) the surface wind field,
     2) the ocean bathymetry,
     3) the surface ocean temperature and velocity,
     4) the ocean temperature and velocity at the iceberg base, and
     5) the depth-integrated ocean temperature and velocity along the iceberg's draft.
    '''
    
    ### Wind: surface wind components
    uvel, vvel = fieldset.UV10[time+(particle.dt/2.), 0., particle.lat, particle.lon]     # [m/s]
    
    particle.uvela = uvel                                                                 # [m/s]
    particle.vvela = vvel                                                                 # [m/s]
    
    ### Bathymetry
    particle.bath = fieldset.B[time+(particle.dt/2.), particle.depth, particle.lat, particle.lon] # [m]
    
    ## Determine the deepest point between the iceberg's base and the bathymetry
    z = min(particle.depth,math.fabs(particle.bath))
    
    ### Surface (ocean): ocean temperature and velocity components at the surface
    particle.tempS = fieldset.T[time+(particle.dt/2.), 0., particle.lat, particle.lon]    # [°C]
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), 0., particle.lat, particle.lon]       # [°/s]
    particle.uveloS = uvel*1852*60*math.cos(particle.lat*math.pi/180.)                    # [m/s]
    particle.vveloS = vvel*1852*60                                                        # [m/s]
    
    ### Basal (ocean): ocean temperature and velocity components at the iceberg base
    particle.tempB = fieldset.T[time+(particle.dt/2.), z, particle.lat, particle.lon]     # [°C]
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), z, particle.lat, particle.lon]        # [°/s]
    particle.uveloB = uvel*1852*60*math.cos(particle.lat*math.pi/180.)                    # [m/s]
    particle.vveloB = vvel*1852*60                                                        # [m/s]
    
    ### Integrated (ocean): depth-integrated ocean temperature and velocity components along the iceberg's draft
    t, u, v = 0, 0, 0
    check = 1

    if z >= 1.01127517 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 0.494025379, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 0.494025379, particle.lat, particle.lon]
        t += temp * 1.01127517
        u += uv * 1.01127517
        v += vv * 1.01127517
    elif z < 1.01127517 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 0.494025379, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 0.494025379, particle.lat, particle.lon]
        t += temp * (z-0.0)
        u += uv * (z-0.0)
        v += vv * (z-0.0)
        check = -1
    
    if z >= 2.08567595 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 1.54137540, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 1.54137540, particle.lat, particle.lon]
        t += temp * (2.08567595-1.01127517)
        u += uv * (2.08567595-1.01127517)
        v += vv * (2.08567595-1.01127517)
    elif z < 2.08567595 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 1.54137540, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 1.54137540, particle.lat, particle.lon]
        t += temp * (z-1.01127517)
        u += uv * (z-1.01127517)
        v += vv * (z-1.01127517)
        check = -1
    
    if z >= 3.22300124 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 2.64566851, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 2.64566851, particle.lat, particle.lon]
        t += temp * (3.22300124-2.08567595)
        u += uv * (3.22300124-2.08567595)
        v += vv * (3.22300124-2.08567595)
    elif z < 3.22300124 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 2.64566851, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 2.64566851, particle.lat, particle.lon]
        t += temp * (z-2.08567595)
        u += uv * (z-2.08567595)
        v += vv * (z-2.08567595)
        check = -1
    
    if z >= 4.43716145 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 3.81949472, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 3.81949472, particle.lat, particle.lon]
        t += temp * (4.43716145-3.22300124)
        u += uv * (4.43716145-3.22300124)
        v += vv * (4.43716145-3.22300124)
    elif z < 4.43716145 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 3.81949472, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 3.81949472, particle.lat, particle.lon]
        t += temp * (z-3.22300124)
        u += uv * (z-3.22300124)
        v += vv * (z-3.22300124)
        check = -1
    
    if z >= 5.74513674 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 5.07822371, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 5.07822371, particle.lat, particle.lon]
        t += temp * (5.74513674-4.43716145)
        u += uv * (5.74513674-4.43716145)
        v += vv * (5.74513674-4.43716145)
    elif z < 5.74513674 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 5.07822371, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 5.07822371, particle.lat, particle.lon]
        t += temp * (z-4.43716145)
        u += uv * (z-4.43716145)
        v += vv * (z-4.43716145)
        check = -1
    
    if z >= 7.16765165 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 6.44061422, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 6.44061422, particle.lat, particle.lon]
        t += temp * (7.16765165-5.74513674)
        u += uv * (7.16765165-5.74513674)
        v += vv * (7.16765165-5.74513674)
    elif z < 7.16765165 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 6.44061422, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 6.44061422, particle.lat, particle.lon]
        t += temp * (z-5.74513674)
        u += uv * (z-5.74513674)
        v += vv * (z-5.74513674)
        check = -1
    
    if z >= 8.72999573 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 7.92956018, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 7.92956018, particle.lat, particle.lon]
        t += temp * (8.72999573-7.16765165)
        u += uv * (8.72999573-7.16765165)
        v += vv * (8.72999573-7.16765165)
    elif z < 8.72999573 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 7.92956018, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 7.92956018, particle.lat, particle.lon]
        t += temp * (z-7.16765165)
        u += uv * (z-7.16765165)
        v += vv * (z-7.16765165)
        check = -1
    
    if z >= 10.46302414 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 9.57299709, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 9.57299709, particle.lat, particle.lon]
        t += temp * (10.46302414-8.72999573)
        u += uv * (10.46302414-8.72999573)
        v += vv * (10.46302414-8.72999573)
    elif z < 10.46302414 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 9.57299709, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 9.57299709, particle.lat, particle.lon]
        t += temp * (z-8.72999573)
        u += uv * (z-8.72999573)
        v += vv * (z-8.72999573)
        check = -1
    
    if z >= 12.40437222 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 11.4050026, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 11.4050026, particle.lat, particle.lon]
        t += temp * (12.40437222-10.46302414)
        u += uv * (12.40437222-10.46302414)
        v += vv * (12.40437222-10.46302414)
    elif z < 12.40437222 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 11.4050026, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 11.4050026, particle.lat, particle.lon]
        t += temp * (z-10.46302414)
        u += uv * (z-10.46302414)
        v += vv * (z-10.46302414)
        check = -1
        
    if z >= 14.59993172 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 13.4671383, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 13.4671383, particle.lat, particle.lon]
        t += temp * (14.59993172-12.40437222)
        u += uv * (14.59993172-12.40437222)
        v += vv * (14.59993172-12.40437222)
    elif z < 14.59993172 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 13.4671383, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 13.4671383, particle.lat, particle.lon]
        t += temp * (z-12.40437222)
        u += uv * (z-12.40437222)
        v += vv * (z-12.40437222)
        check = -1
    
    if z >= 17.10564232 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 15.8100729, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 15.8100729, particle.lat, particle.lon]
        t += temp * (17.10564232-14.59993172)
        u += uv * (17.10564232-14.59993172)
        v += vv * (17.10564232-14.59993172)
    elif z < 17.10564232 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 15.8100729, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 15.8100729, particle.lat, particle.lon]
        t += temp * (z-14.59993172)
        u += uv * (z-14.59993172)
        v += vv * (z-14.59993172)
        check = -1
    
    if z >= 19.98966408 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 18.4955597, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 18.4955597, particle.lat, particle.lon]
        t += temp * (19.98966408-17.10564232)
        u += uv * (19.98966408-17.10564232)
        v += vv * (19.98966408-17.10564232)
    elif z < 19.98966408 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 18.4955597, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 18.4955597, particle.lat, particle.lon]
        t += temp * (z-17.10564232)
        u += uv * (z-17.10564232)
        v += vv * (z-17.10564232)
        check = -1
    
    if z >= 23.33499336 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 21.5988159, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 21.5988159, particle.lat, particle.lon]
        t += temp * (23.33499336-19.98966408)
        u += uv * (23.33499336-19.98966408)
        v += vv * (23.33499336-19.98966408)
    elif z < 23.33499336 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 21.5988159, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 21.5988159, particle.lat, particle.lon]
        t += temp * (z-19.98966408)
        u += uv * (z-19.98966408)
        v += vv * (z-19.98966408)
        check = -1
    
    if z >= 27.24263382 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 25.2114086, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 25.2114086, particle.lat, particle.lon]
        t += temp * (27.24263382-23.33499336)
        u += uv * (27.24263382-23.33499336)
        v += vv * (27.24263382-23.33499336)
    elif z < 27.24263382 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 25.2114086, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 25.2114086, particle.lat, particle.lon]
        t += temp * (z-23.33499336)
        u += uv * (z-23.33499336)
        v += vv * (z-23.33499336)
        check = -1
    
    if z >= 31.83539963 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 29.4447289, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 29.4447289, particle.lat, particle.lon]
        t += temp * (31.83539963-27.24263382)
        u += uv * (31.83539963-27.24263382)
        v += vv * (31.83539963-27.24263382)
    elif z < 31.83539963 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 29.4447289, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 29.4447289, particle.lat, particle.lon]
        t += temp * (z-27.24263382)
        u += uv * (z-27.24263382)
        v += vv * (z-27.24263382)
        check = -1
    
    if z >= 37.2624855 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 34.4341545, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 34.4341545, particle.lat, particle.lon]
        t += temp * (37.2624855-31.83539963)
        u += uv * (37.2624855-31.83539963)
        v += vv * (37.2624855-31.83539963)
    elif z < 37.2624855 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 34.4341545, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 34.4341545, particle.lat, particle.lon]
        t += temp * (z-246.8574)
        u += uv * (z-246.8574)
        v += vv * (z-246.8574)
        check = -1
    
    if z >= 43.70489883 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 40.3440514, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 40.3440514, particle.lat, particle.lon]
        t += temp * (43.70489883-37.2624855)
        u += uv * (43.70489883-37.2624855)
        v += vv * (43.70489883-37.2624855)
    elif z < 43.70489883 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 40.3440514, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 40.3440514, particle.lat, particle.lon]
        t += temp * (z-37.2624855)
        u += uv * (z-37.2624855)
        v += vv * (z-37.2624855)
        check = -1
    
    if z >= 51.38193512 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 47.3736877, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 47.3736877, particle.lat, particle.lon]
        t += temp * (51.38193512-43.70489883)
        u += uv * (51.38193512-43.70489883)
        v += vv * (51.38193512-43.70489883)
    elif z < 51.38193512 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 47.3736877, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 47.3736877, particle.lat, particle.lon]
        t += temp * (z-43.70489883)
        u += uv * (z-43.70489883)
        v += vv * (z-43.70489883)
        check = -1
    
    if z >= 60.55880737 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 55.7642899, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 55.7642899, particle.lat, particle.lon]
        t += temp * (60.55880737-51.38193512)
        u += uv * (60.55880737-51.38193512)
        v += vv * (60.55880737-51.38193512)
    elif z < 60.55880737 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 55.7642899, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 55.7642899, particle.lat, particle.lon]
        t += temp * (z-51.38193512)
        u += uv * (z-51.38193512)
        v += vv * (z-51.38193512)
        check = -1
    
    if z >= 71.55553436 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 65.8072739, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 65.8072739, particle.lat, particle.lon]
        t += temp * (71.55553436-60.55880737)
        u += uv * (71.55553436-60.55880737)
        v += vv * (71.55553436-60.55880737)
    elif z < 71.55553436 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 65.8072739, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 65.8072739, particle.lat, particle.lon]
        t += temp * (z-60.55880737)
        u += uv * (z-60.55880737)
        v += vv * (z-60.55880737)
        check = -1
    
    if z >= 84.75727844 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 77.8538513, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 77.8538513, particle.lat, particle.lon]
        t += temp * (84.75727844-71.55553436)
        u += uv * (84.75727844-71.55553436)
        v += vv * (84.75727844-71.55553436)
    elif z < 84.75727844 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 77.8538513, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 77.8538513, particle.lat, particle.lon]
        t += temp * (z-71.55553436)
        u += uv * (z-71.55553436)
        v += vv * (z-71.55553436)
        check = -1

    if z >= 100.626091 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 92.3260727, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 92.3260727, particle.lat, particle.lon]
        t += temp * (100.626091-84.75727844)
        u += uv * (100.626091-84.75727844)
        v += vv * (100.626091-84.75727844)
    elif z < 100.626091 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 92.3260727, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 92.3260727, particle.lat, particle.lon]
        t += temp * (z-84.75727844)
        u += uv * (z-84.75727844)
        v += vv * (z-84.75727844)
        check = -1

    if z >= 119.71409607 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 109.729279, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 109.729279, particle.lat, particle.lon]
        t += temp * (119.71409607-100.626091)
        u += uv * (119.71409607-100.626091)
        v += vv * (119.71409607-100.626091)
    elif z < 119.71409607 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 109.729279, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 109.729279, particle.lat, particle.lon]
        t += temp * (z-100.626091)
        u += uv * (z-100.626091)
        v += vv * (z-100.626091)
        check = -1

    if z >= 142.67788696 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 130.665985, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 130.665985, particle.lat, particle.lon]
        t += temp * (142.67788696-119.71409607)
        u += uv * (142.67788696-119.71409607)
        v += vv * (142.67788696-119.71409607)
    elif z < 142.67788696 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 130.665985, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 130.665985, particle.lat, particle.lon]
        t += temp * (z-119.71409607)
        u += uv * (z-119.71409607)
        v += vv * (z-119.71409607)
        check = -1

    if z >= 170.29385376 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 155.850723, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 155.850723, particle.lat, particle.lon]
        t += temp * (170.29385376-142.67788696)
        u += uv * (170.29385376-142.67788696)
        v += vv * (170.29385376-142.67788696)
    elif z < 170.29385376 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 155.850723, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 155.850723, particle.lat, particle.lon]
        t += temp * (z-142.67788696)
        u += uv * (z-142.67788696)
        v += vv * (z-142.67788696)
        check = -1

    if z >= 203.473526 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 186.125565, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 186.125565, particle.lat, particle.lon]
        t += temp * (203.473526-170.29385376)
        u += uv * (203.473526-170.29385376)
        v += vv * (203.473526-170.29385376)
    elif z < 203.473526 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 186.125565, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 186.125565, particle.lat, particle.lon]
        t += temp * (z-170.29385376)
        u += uv * (z-170.29385376)
        v += vv * (z-170.29385376)
        check = -1

    if z >= 243.27807617 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 222.475174, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 222.475174, particle.lat, particle.lon]
        t += temp * (243.27807617-203.473526)
        u += uv * (243.27807617-203.473526)
        v += vv * (243.27807617-203.473526)
    elif z < 243.27807617 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 222.475174, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 222.475174, particle.lat, particle.lon]
        t += temp * (z-203.473526)
        u += uv * (z-203.473526)
        v += vv * (z-203.473526)
        check = -1

    if z >= 290.93026733 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 266.040253, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 266.040253, particle.lat, particle.lon]
        t += temp * (290.93026733-243.27807617)
        u += uv * (290.93026733-243.27807617)
        v += vv * (290.93026733-243.27807617)
    elif z < 290.93026733 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 266.040253, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 266.040253, particle.lat, particle.lon]
        t += temp * (z-243.27807617)
        u += uv * (z-243.27807617)
        v += vv * (z-243.27807617)
        check = -1

    if z >= 347.82159424 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 318.127441, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 318.127441, particle.lat, particle.lon]
        t += temp * (347.82159424-290.93026733)
        u += uv * (347.82159424-290.93026733)
        v += vv * (347.82159424-290.93026733)
    elif z < 347.82159424 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 318.127441, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 318.127441, particle.lat, particle.lon]
        t += temp * (z-290.93026733)
        u += uv * (z-290.93026733)
        v += vv * (z-290.93026733)
        check = -1

    if z >= 415.51190186 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 380.213013, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 380.213013, particle.lat, particle.lon]
        t += temp * (415.51190186-347.82159424)
        u += uv * (415.51190186-347.82159424)
        v += vv * (415.51190186-347.82159424)
    elif z < 415.51190186 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 380.213013, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 380.213013, particle.lat, particle.lon]
        t += temp * (z-347.82159424)
        u += uv * (z-347.82159424)
        v += vv * (z-347.82159424)
        check = -1

    if z >= 495.71838379 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 453.937744, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 453.937744, particle.lat, particle.lon]
        t += temp * (495.71838379-415.51190186)
        u += uv * (495.71838379-415.51190186)
        v += vv * (495.71838379-415.51190186)
    elif z < 495.71838379 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 453.937744, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 453.937744, particle.lat, particle.lon]
        t += temp * (z-415.51190186)
        u += uv * (z-415.51190186)
        v += vv * (z-415.51190186)
        check = -1

    if z >= 590.2901001 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 541.088928, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 541.088928, particle.lat, particle.lon]
        t += temp * (590.2901001-495.71838379)
        u += uv * (590.2901001-495.71838379)
        v += vv * (590.2901001-495.71838379)
    elif z < 590.2901001 and check > 0:
        temp = fieldset.T[time+(particle.dt/2.), 541.088928, particle.lat, particle.lon]
        uv, vv = fieldset.UV[time+(particle.dt/2.), 541.088928, particle.lat, particle.lon]
        t += temp * (z-495.71838379)
        u += uv * (z-495.71838379)
        v += vv * (z-495.71838379)
        check = -1
    
    particle.tempD = t / z                                               # [°C]
    uvel = u / z                                                         # [°/s]
    vvel = v / z                                                         # [°/s]
    
    particle.uveloD = uvel*1852*60*math.cos(particle.lat*math.pi/180.)   # [m/s]
    particle.vveloD = vvel*1852*60                                       # [m/s]



#============================== Kernels surface ==============================#
def SampleFieldsSurf(particle, fieldset, time):
    '''
    At the iceberg's location, sample:
     1) the surface wind field,
     2) the surface ocean temperature and velocity.
    '''
    
    ### Wind: surface wind stress components and surface wind speed
    taux, tauy = fieldset.XY[time+(particle.dt/2.), 0., particle.lat, particle.lon]       # [g/(s2 cm)]    
    vela = math.sqrt((math.sqrt(taux**2+tauy**2)*0.1)/(fieldset.rho_a*0.0015))            # [m/s]
    particle.uvela = (taux*0.1)/(fieldset.rho_a*0.0015*math.fabs(vela))                   # [m/s]
    particle.vvela = (tauy*0.1)/(fieldset.rho_a*0.0015*math.fabs(vela))                   # [m/s]
    
    ### Surface (ocean): ocean temperature and velocity components at the surface
    particle.temp = fieldset.T[time+(particle.dt/2.), 0., particle.lat, particle.lon]     # [°C]
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), 0., particle.lat, particle.lon]       # [°/s]
    particle.uvelo = uvel*1852*60*math.cos(particle.lat*math.pi/180.)                     # [m/s]
    particle.vvelo = vvel*1852*60                                                         # [m/s]


def AdvectionEESurf(particle, fieldset, time):
    '''
    Advection of icebergs using Explicit Euler (aka Euler Forward) integration
    using the surface ocean velocity components.
    '''
    
    particle.lon += particle.uvelo/(1852*60*math.cos(particle.lat*math.pi/180.)) * particle.dt
    particle.lat += particle.vvelo/(1852*60) * particle.dt


def BuoyantConvectionSurf(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through buoyant convection at the icebergs sides following:
        Mv = b1 * T + b2 * T^2
    with
        b1 = 7.62E-3 [m/(d °C)]
        b2 = 1.29E-3 [m/(d °C)]
        T: ocean temperature [°C]
    See also Merino et al. (2016).
    '''
    
    Mvr = 7.62e-3 * particle.temp + 1.29e-3 * (particle.temp)**2   # [m/d]
    
    particle.Mvr = Mvr/fieldset.sec_to_day                         # [m/s]
    particle.prev_Mvr = particle.Mvr


def BasalMeltSurf(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through basal melt at the iceberg base following:
        Mb = C * |v_ice - v_ocean|^0.8 * (T_ocean - T_ice) / (L^0.2) 
    with
        C = 0.58 [°C^-1 m^0.4 d^-1 s^0.8]
        v_ice: iceberg velocity [m/s]
        v_ocean: ocean velocity [m/s]
        T_ocean: ocean temperature [°C]
        T_ice = -1.92 [°C]; freezing temperature
        L: iceberg length [m]
    See FitzMaurice & Stern (2018) and references therein.
    '''
    
    # Let's sample the second ocean layer so we can calculate some basal melt
    # Nnormally, vIB =/= vS. Now, let's calculate uB as mean of u0 and u1.
    uvel, vvel = fieldset.UV[time+(particle.dt/2.), 10.1, particle.lat, particle.lon]     # [°/s]
    uvel = (uvel*1852*60*math.cos(particle.lat*math.pi/180.)+particle.uvelo)/2.           # [m/s]
    vvel = (vvel*1852*60+particle.vvelo)/2.                                               # [m/s]
    
    dvelabs = math.sqrt((uvel-particle.uvelo)**2 + (vvel-particle.vvelo)**2) # actually 0 as IB moves with surface flow
    
    Mbr = 0.58 * (dvelabs**0.8) * ((particle.temp-fieldset.Ti)/(particle.L**0.2))   # [m/d]
    particle.Mbr = Mbr/fieldset.sec_to_day                                          # [m/s]
    particle.prev_Mbr = particle.Mbr


def WaveErosionSurf(particle, fieldset, time):
    '''
    Calculates the mass loss rate [m/s] through wave erosion at the icebergs sides:
        Me = 1/12 Ss * (1+cos(pi*A_ice^3)) * (T_ocean + 2)
    with sea state:
        Ss = 3/2 (|v_atm - v_ocean|)^0.5 + 1/10 |v_atm - v_ocean|
    and
        A_ice: fractional sea-ice area [m2]; negligible in the Eocene
        T_ocean: ocean surface temperature [°C]
        v_atm: surface wind velocity [m/s]
        v_ocean: ocean surface velocity [m/s]
    '''
    
    dvelabs = math.sqrt((particle.uvela-particle.uvelo)**2 + (particle.vvela-particle.vvelo)**2)
    Ss = (3./2.) * math.sqrt(dvelabs) + (1./10.) * dvelabs
    
    Ai = 0.
    damping = 0.5 * (1 + math.cos(math.pi * Ai**3))
    
    Mer = (1./6.) * Ss * damping * (particle.temp - - 2.)    # [m/d]
    particle.Mer = Mer/fieldset.sec_to_day                   # [m/s]
    particle.prev_Mer = particle.Mer