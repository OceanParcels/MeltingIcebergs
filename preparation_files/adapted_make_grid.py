"""
Created on Thu Apr 18 16:06:46 2019

@author: nooteboom
Adapted by MV Elbertsen
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pylab as plt

# Read the grid file and bathymetry file
gridfile = Dataset('grid_coordinates_pop_tx0.1_38ma.nc') # Follow steps form Nooteboom et al. (2022)

#%% Write new netcdf file with new grid
dataset = Dataset('new_grid_coordinates_pop_tx0.1_38ma.nc', 'w')

i_indexs = dataset.createDimension('i_index', 3600)
j_indexs = dataset.createDimension('j_index', 2550)
depth_ts = dataset.createDimension('depth_t', 42)
w_deps   = dataset.createDimension('w_dep', 43)

ins         = dataset.createVariable('i_index', np.float32,('i_index',))
jns         = dataset.createVariable('j_index', np.float32,('j_index',))
depth_tns   = dataset.createVariable('depth_t', np.float32,('depth_t',))
w_depns     = dataset.createVariable('w_dep', np.float32,('w_dep',))
latitudes   = dataset.createVariable('T_LAT_2D', np.float32,('j_index','i_index',))
longitudes  = dataset.createVariable('T_LON_2D', np.float32,('j_index','i_index',))
latitudes2  = dataset.createVariable('U_LAT_2D', np.float32,('j_index','i_index',))
longitudes2 = dataset.createVariable('U_LON_2D', np.float32,('j_index','i_index',))

htns   = dataset.createVariable('HTN', np.float32,('j_index','i_index',))
htes   = dataset.createVariable('HTE', np.float32,('j_index','i_index',))
huss   = dataset.createVariable('HUS', np.float32,('j_index','i_index',))
huws   = dataset.createVariable('HUW', np.float32,('j_index','i_index',))
angles = dataset.createVariable('ANGLE', np.float32,('j_index','i_index',))
tareas = dataset.createVariable('TAREA', np.float32,('j_index','i_index',))

# Write data
ins[:]         = gridfile['i_index'][:]
jns[:]         = gridfile['j_index'][:]
depth_tns[:]   = gridfile['depth_t'][:]
w_depns[:]     = gridfile['w_dep'][:]
latitudes[:]   = gridfile['T_LAT_2D'][:]
longitudes[:]  = gridfile['T_LON_2D'][:]
latitudes2[:]  = np.concatenate((gridfile['U_LAT_2D'][:,250:], gridfile['U_LAT_2D'][:,:250]), axis=1)
longitudes2[:] = np.concatenate((gridfile['U_LON_2D'][:,250:], gridfile['U_LON_2D'][:,:250]), axis=1)

htns[:]   = gridfile['HTN'][:]
htes[:]   = gridfile['HTE'][:]
huss[:]   = gridfile['HUS'][:]
huws[:]   = gridfile['HUW'][:]
angles[:] = gridfile['ANGLE'][:]
tareas[:] = gridfile['TAREA'][:]

#Attributes:
latitudes.long_name   = 'latitude on t-grid'
latitudes.units       = 'degrees N'
longitudes.long_name  = 'longitude on t-grid'
longitudes.units      = 'degrees N'
latitudes2.long_name  = 'latitude on u-grid'
latitudes2.units      = 'degrees N'
longitudes2.long_name = 'longitude on u-grid'
longitudes2.units     = 'degrees N'
w_depns.units         = 'meters'
w_depns.long_name     = 'T-grid depth'
depth_tns.units       = 'meters'
depth_tns.long_name   = 'T-grid depth'
ins.long_name         = 'i-coordinate index'
jns.long_name         = 'j-coordinate index'

dataset.close()