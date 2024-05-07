'''
Code for passive particle simulations
---------------------------------------------
#(!)#  marks changes required every run (filename)
#(?)#  marks changes required to run backward simulation
'''

#================================== IMPORT ===================================#
from parcels import (FieldSet, ParticleSet, Variable, JITParticle, ErrorCode,
                     ParticleFile, UnitConverter, Field, AdvectionRK4)
from datetime import datetime, date
from datetime import timedelta as delta
from glob import glob

import math
import numpy as np
import xarray as xr
import pandas as pd

from operator import attrgetter



#==================== LOAD/SELECT DATA & DEFINE FIELDSET ====================#
### Grid (POP)
mesh_mask = '/nethome/5867800/grid/edited_grid_coordinates_pop_tx0.1_38ma.nc'

### Eocene data Nooteboom et al. (2022)
data_path_ocean = '/storage/shared/pop/p21a.EO38Ma.tx0.1.2pic_control/daily/'

files = sorted(glob(data_path_ocean+'vars_first600m_eocene_2pic_pop_00*.nc'))
filename = data_path_ocean+'vars_first600m_eocene_2pic_pop_00400301.nc'
for i, name in enumerate(files):
    if name == filename:
        files.pop(i)   # remove empty file

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': mesh_mask, 'data': files},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': mesh_mask, 'data': files},
             'T': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': mesh_mask, 'data': files}}#,
             # 'X': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': mesh_mask, 'data': files},
             # 'Y': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': mesh_mask, 'data': files}}

variables = {'U': 'UVEL', 'V': 'VVEL', 'T': 'TEMP'}#, 'X': 'TAUX', 'Y': 'TAUY'}

dimensions = {'U': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep', 'time': 'time'},   # Should all be equal
              'V': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep', 'time': 'time'},
              'T': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'depth': 'w_dep', 'time': 'time'}}#,   # Use U even though defined op T-grid
              # 'X': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time': 'time'},
              # 'Y': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D', 'time': 'time'}}

### Antarctica indices
indices = {'lat': np.arange(0, 752).tolist()} # ~90-50°S

### Time
dr = pd.date_range(start='01/01/2038', end='12/31/2042', freq='d')            # exact dates not relevant, only frequency !!! US format !!!
dates1 = dr[(dr.day != 29) | (dr.month != 2)]
dates = dates1[(dates1.day != 1) | (dates1.month != 3) | (dates1.year != 2040)]   # this file is empty
timestamps = np.expand_dims(np.array(dates), axis=1)

### Load fieldset
fieldset = FieldSet.from_pop(filenames, variables, dimensions, indices=indices, timestamps=timestamps)   # from_pop assumes velocities in cm/s
# XY = VectorField('XY', fieldset.X, fieldset.Y)
# fieldset.add_vector_field(XY)

### Bathymetry and regions #(?)#
# br = '/nethome/5867800/grid/bathymetry_regions.nc'

# filenames_b = {'B': {'lon': mesh_mask, 'lat': mesh_mask, 'data': br}}
# variables_b = {'B': 'bathymetry'}
# dimensions_b = {'B': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'}}

# fieldset_b = FieldSet.from_netcdf(filenames_b, variables_b, dimensions_b)
# fieldset.add_field(fieldset_b.B)

### Regions based on Carter et al. (2017) --- Add only for backward simulations #(?)#
# filenames_r = {'R': {'lon': mesh_mask, 'lat': mesh_mask, 'data': br}}
# variables_r = {'R': 'region'}
# dimensions_r = {'R': {'lon': 'U_LON_2D', 'lat': 'U_LAT_2D'}}
# fieldset_r = FieldSet.from_netcdf(filenames_r, variables_r, dimensions_r)
# fieldset.add_field(fieldset_r.R)
# fieldset.R.interp_method = ('cgrid_tracer')



#==================== DEFINE PARTICLESET ====================#
### Define particle class #(?)#
class fieldsampler(JITParticle):
    uvel = Variable('uvel', initial=0, to_write=False)
    vvel = Variable('vvel', initial=0, to_write=False)
    velo = Variable('velo', initial=0)
    temp = Variable('temp', initial=0)   
    # taux = Variable('taux', initial=0)
    # tauy = Variable('tauy', initial=0)
    vela = Variable('vela', initial=0) 
    # bath = Variable('bath', initial=0, to_write=False)    # bathymetry [m]
    
    ## Backwards only
    # reg = Variable('reg', initial=0, dtype=np.int32)      # region Carter
    
    ## Forwards only
    distance = Variable('distance', initial=1000000.)  # distance to ODP Site 696
    prev_dist = Variable('prev_dist', to_write=False, initial=attrgetter('distance'))


### Define release locations #(?)#
## Forwards
lons = np.array([-55., -50., -45., -40., -35., -30., -25., -20.,
                 -15., -10.,  -5.,   0.,   5.,  10.,  15.,  20.,
                  25.,  30.,  35.,  40.,  45.,  50.,  55.,  60.,
                 -57.82330253, -59.25, -60.08662877, -60.01036191, -58.82106809])
lats = np.array([-77.99537227, -78.45787684, -79.17977465, -79.75459498, -79.99027112, -80.07689735, -79.9997914 , -80.12097168,
                 -79.24559138, -78.39033646, -77.33178711, -76.13251966, -75.17986379, -74.68456009, -73.6896764 , -73.53134648,
                 -72.9724492 , -72.3117251 , -72.29006817, -72.16268516, -70.86573029, -68.58359528, -68.58359528, -66.59729004,
                 -77.        , -75.        , -73.        , -71.        , -69.        ])
## Backwards
# lons = np.array([-57.04998779, -57.0249939,  -57.,         -56.9750061,  -56.95001221,
#                  -57.04998779, -57.0249939,  -57.,         -56.9750061,  -56.95001221,
#                  -57.04998779, -57.0249939,  -57.,         -56.9750061,  -56.95001221,
#                  -57.04998779, -57.0249939,  -57.,         -56.9750061,  -56.95001221,
#                  -57.04998779, -57.0249939,  -57.,         -56.9750061,  -56.95001221])
# lats = np.array([-67.5270462,  -67.5270462,  -67.5270462,  -67.5270462,  -67.5270462,
#                  -67.5164814,  -67.5164814,  -67.5164814,  -67.5164814,  -67.5164814,
#                  -67.5059166,  -67.5059166,  -67.5059166,  -67.5059166,  -67.5059166,
#                  -67.49535179, -67.49535179, -67.49535179, -67.49535179, -67.49535179,
#                  -67.48478699, -67.48478699, -67.48478699, -67.48478699, -67.48478699])

lonsar = np.tile(lons, len(timestamps))
latsar = np.tile(lats, len(timestamps))
timear = np.repeat(timestamps, len(lons))

### Define global parameters
fieldset.add_constant('ODP_lat', -67.5)       # Latitude ODP 696 Eocene [°]
fieldset.add_constant('ODP_lon', -57.)        # Longitude ODP 696 Eocene [°]
fieldset.add_constant('rho_a', 1.293)         # Air density [kg/m3]
fieldset.add_constant('halo_west', min(fieldset.U.grid.lon))
fieldset.add_constant('halo_east', max(fieldset.U.grid.lon))
fieldset.add_periodic_halo(zonal=True)

### Load particleset  #(?)#
## FW and BW runs
pset = ParticleSet.from_list(fieldset = fieldset,
                             pclass   = fieldsampler,
                             lon      = lonsar,
                             lat      = latsar,
                             depth    = np.full(len(latsar), [0.]),
                             time     = timear)
## Circumpolar run
# pset = ParticleSet.from_line(fieldset = fieldset,
#                              pclass   = JITParticle,
#                              size     = 20,
#                              start    = (-57, -67),
#                              finish   = (-57, -66),
#                              repeatdt = delta(days=1))



#==================== DEFINE KERNELS ====================#
def periodicBC(particle, fieldset, time):
    '''
    Add periodic boundary conditions to the domain.
    '''
    
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west


def DeleteParticle(particle, fieldset, time):
    '''
    Remove icebergs that leave the domain.
    '''
    
    particle.delete()


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

    
### Sampling kernel #(?)#
def SampleUVT(particle, fieldset, time):
    ''' Samples the temperature, velocity, stress fields and water depth at the location of the particle. '''
    uv, vv = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    particle.uvel = uv*1852*60*math.cos(particle.lat*math.pi/180.)                 # [m/s]
    particle.vvel = vv*1852*60                                                     # [m/s]
    particle.velo = math.sqrt(particle.uvel**2 + particle.vvel**2)
    particle.temp = fieldset.T[time, particle.depth, particle.lat, particle.lon]   # [°C]
    taux = fieldset.X[time, 0., particle.lat, particle.lon]   # [g/(s2 cm)]
    tauy = fieldset.Y[time, 0., particle.lat, particle.lon]   # [g/(s2 cm)]
    particle.vela = math.sqrt((math.sqrt(taux**2+tauy**2)*0.1)/(fieldset.rho_a*0.0015)) # [m/s]
    # particle.bathy = fieldset.B[time, particle.depth, particle.lat, particle.lon]
    
    ## Backwards
    # particle.reg = fieldset.R[time, 0., particle.lat, particle.lon]



#==================== EXECUTE ====================#
### Define output file #(!)#
out = pset.ParticleFile(name='29p_5y_pass_30d_1hdtar.zarr', outputdt=delta(days=1))

### Execute #(?)#
## Record initial values of field
pset.execute(pset.Kernel(SampleUVT), dt=0)

## Forward
pset.execute(pset.Kernel(SampleUVT) + pset.Kernel(MinimumDistance) + AdvectionRK4 + pset.Kernel(periodicBC),# + pset.Kernel(TotalDistance),
             runtime=delta(days=(len(dr)-2)),
             dt=delta(hours=1),
             output_file=out,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

# ## Backwards
# pset.execute(pset.Kernel(SampleUVT) + AdvectionRK4 + pset.Kernel(periodicBC),
#              runtime=delta(days=(len(dr)-2)),
#              dt=-delta(hours=1),
#              output_file=out,
#              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

## Circumpolar: release during first two years
# pset.execute(AdvectionRK4 + pset.Kernel(periodicBC),
#              runtime=delta(days=(365*2)),
#              dt=delta(hours=1),
#              output_file=out,
#              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
# pset.repeatdt = None
# pset.execute(AdvectionRK4 + pset.Kernel(periodicBC),
#              runtime=delta(days=(len(dr)-2-(365*2))),
#              dt=delta(hours=1),
#              output_file=out,
#              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

out.close()