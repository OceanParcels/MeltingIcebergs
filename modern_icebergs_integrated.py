'''
Code for modern Antarctic iceberg simulations
---------------------------------------------
#(!)#  marks changes required every run (filename, files)
'''

#================================== IMPORT ===================================#
from parcels import (FieldSet, ParticleSet, Variable, JITParticle, ErrorCode,
                     ParticleFile, UnitConverter, Field, VectorField)
from datetime import datetime, date
from datetime import timedelta as delta
from glob import glob

import math
import numpy as np
import xarray as xr
import pandas as pd

import sys
sys.path.append('../active/')

import icebergs_kernels as kernels
import icebergs_particleclass as particleclass



#==================== LOAD/SELECT DATA & DEFINE FIELDSET ====================#
coord_path = '/storage/shared/oceanparcels/input_data/'
data_path = coord_path
data_pathm = 'fields/'

### Ocean fields #(!)#
# ufiles = sorted(glob(data_path+'MOi/psy4v3r1/psy4v3r1-daily_U_2021-*.nc')) # daily
ufiles = sorted(glob(data_path+'MOi/psy4v3r1-monthly_U_2021-*.nc'))        # monthly
vfiles = [f.replace('_U_', '_V_') for f in ufiles]
tfiles = [f.replace('_U_', '_T_') for f in ufiles]
wfile = coord_path+'MOi/psy4v3r1/psy4v3r1-daily_W_2021-01-01.nc'

mesh_mask = coord_path + 'MOi/domain_ORCA0083-N006/coordinates.nc'

filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfile, 'data': ufiles},
             'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfile, 'data': vfiles},
             'T': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfile, 'data': tfiles}}

variables = {'U': 'vozocrtx', 'V': 'vomecrty', 'T': 'votemper'}

dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
              'T': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}}

### Bathymetry
bathymetry = coord_path + 'MOi/domain_ORCA0083-N006/bathymetry_ORCA12_V3.3.nc'

filenamesb = {'B': {'lon': bathymetry, 'lat': bathymetry, 'depth': bathymetry, 'data': bathymetry}}
variablesb = {'B': 'Bathymetry'}
dimensionsb = {'B': {'lon': 'nav_lon', 'lat': 'nav_lat'}}

### Wind #(!)#
# filesw = sorted(glob(data_path+'ERA5/reanalysis-era5-single-level_wind10m_202101*.nc')) # daily
filesw = sorted(glob(data_pathm+'ERA5/reanalysis-era5-single-level_wind10m_2021*.nc'))  # monthly

filenamesw = {'U10': {'lon': filesw[0], 'lat': filesw[0], 'data': filesw},
                'V10': {'lon': filesw[0], 'lat': filesw[0], 'data': filesw}}
variablesw = {'U10': 'u10', 'V10': 'v10'}
dimensionsw = {'U10': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'},
                 'V10': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}}

### Antarctica indices
indices = {'depth': range(0, 33),  # Note at least two depth layers need to be loaded
           'lat': np.arange(0, 800).tolist()} # -77, -50

indicesw = {'lat': np.arange(550, 721).tolist()} # -50, -90

### Time
dr = pd.date_range(start='01/01/2021', end='12/31/2021', freq='d')
dates = dr[(dr.day != 29) | (dr.month != 2)]
timestampsf = np.expand_dims(np.array(dates), axis=1)
timestamps = timestampsf
drf = dr
#(!)# add for monthly
dr = pd.date_range(start='01/01/2021', end='12/31/2021', freq='MS')
dates = dr[(dr.day != 29) | (dr.month != 2)]
timestampsf = np.expand_dims(np.array(dates), axis=1)

### Load fieldset
fieldset = FieldSet.from_nemo(filenames, variables, dimensions, indices=indices, timestamps=timestamps)
## Add bathymetry
fieldsetb = FieldSet.from_netcdf(filenamesb, variablesb, dimensionsb, indices=indices)
fieldset.add_field(fieldsetb.B)
## Add wind
fieldsetw = FieldSet.from_netcdf(filenamesw, variablesw, dimensionsw, indices=indicesw, timestamps=timestampsf, allow_time_extrapolation=True) # #(!)# w/ timestamps
fieldset.add_field(fieldsetw.U10)
fieldset.add_field(fieldsetw.V10)
UV10 = VectorField('UV10', fieldset.U10, fieldset.V10)
fieldset.add_vector_field(UV10)



#==================== DEFINE PARTICLESET ====================#
### Define particle class
pclass = particleclass.Iceberg_FW


### Define start locations
## Forwards
lons = np.array([-55., -50., -45., -40., -35., -30., -25., -20.,
                 -15., -10.,  -5.,   0.,   5.,  10.,  15.,  20.,
                  25.,  30.,  35.,  40.,  45.,  50.,  55.,  60.,
                 -57.82330253, -59.25, -60.08662877, -60.01036191, -58.82106809])
lats = np.array([-77.99537227, -78.45787684, -79.17977465, -79.75459498, -79.99027112, -80.07689735, -79.9997914 , -80.12097168,
                 -79.24559138, -78.39033646, -77.33178711, -76.13251966, -75.17986379, -74.68456009, -73.6896764 , -73.53134648,
                 -72.9724492 , -72.3117251 , -72.29006817, -72.16268516, -70.86573029, -68.58359528, -68.58359528, -66.59729004,
                 -77.        , -75.        , -73.        , -71.        , -69.        ]) + 5. # shift northwards to compensate for geography

lonsar = np.tile(lons, len(timestamps))
latsar = np.tile(lats, len(timestamps))
timear = np.repeat(timestamps, len(lons))

### Define global parameters
fieldset.add_constant('ODP_lat', -61.849083)  # Latitude ODP 696 present-day [°]
fieldset.add_constant('ODP_lon', -42.933067)  # Longitude ODP 696 present-day [°]

fieldset.add_constant('rho_i', 850.0)         # Ice density [kg/m3]
fieldset.add_constant('rho_o', 1027.5)        # Ocean density [kg/m3]
fieldset.add_constant('rho_a', 1.293)         # Air density [kg/m3]
fieldset.add_constant('Ti', -1.92)            # Freezing temperature [°C] --> average determined from salinity field 00380101 in Weddell Sea
fieldset.add_constant('sec_to_day', 86400.)   # Seconds per day [s/d]

fieldset.add_constant('halo_west', min(fieldset.U.grid.lon))
fieldset.add_constant('halo_east', max(fieldset.U.grid.lon))
fieldset.add_periodic_halo(zonal=True)

### Load particleset
pset = ParticleSet.from_list(fieldset=fieldset,
                             pclass=pclass,
                             lon=lonsar,
                             lat=latsar,
                             depth=np.full(len(latsar), [0]),
                             time=timear)



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


### Kernel groups
properties = pset.Kernel(kernels.NewProperties) + \
    pset.Kernel(kernels.SetInitialPropertiesC4)

melt_terms = pset.Kernel(kernels.BuoyantConvection) + \
    pset.Kernel(kernels.BasalMelt) + pset.Kernel(kernels.WaveErosion)

dynamic =  pset.Kernel(kernels.AdvectionEE) + pset.Kernel(periodicBC) + \
    pset.Kernel(kernels.MinimumDistance)


#==================== EXECUTE ====================#
### Define output file #(!)#
out = pset.ParticleFile(name='modern_1y_int_30d_1hdtar_C4_monthly_ate.zarr', outputdt=delta(days=30)) # monthly
# out = pset.ParticleFile(name='modern_1y_int_30d_1hdtar_C4_hourly.zarr', outputdt=delta(days=30))      # hourly

### Execute #(!)#
## Daily
# pset.execute(properties + pset.Kernel(kernels.SampleFieldsModern) + melt_terms + pset.Kernel(kernels.Grounding) + dynamic,
#              runtime=delta(days=(len(dr)-2)),
#              dt=delta(hours=1),
#              output_file=out,
#              recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
## Monthly
pset.execute(properties + pset.Kernel(kernels.SampleFieldsModern) + melt_terms + pset.Kernel(kernels.Grounding) + dynamic,
             runtime=delta(days=(len(drf)-2)),
             dt=delta(hours=1),
             output_file=out,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

out.close()