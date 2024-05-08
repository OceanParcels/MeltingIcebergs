"""
Created on Wed Apr 24 14:58:55 2019

@author: nooteboom
Adapted by MV Elbertsen
"""

from netCDF4 import Dataset
from scipy.interpolate import griddata
import numpy as np

#%% Read nc bathymetry file
# Load the coarse topobathymetry file: Follow steps in supplement from Baatsen et al. (2016) or
# download prepared file form Nooteboom et al. (2022)
topbathfile = Dataset('TopoBathy38.nc')
Z   = topbathfile['Z'][:]
lon = topbathfile['longitude'][:]
lat = topbathfile['latitude'][:]

# Convert to a meshgrid
lons, lats = np.meshgrid(lon, lat)

# Add a '90 degrees North' latitude band, to assure that all interpolated values are within the domain of the starting grid
Z    = np.concatenate((Z, Z[-1][np.newaxis,:]), axis = 0)
lats = np.concatenate((lats,np.full(3600, 90)[np.newaxis,:]), axis = 0)
lons = np.concatenate((lons,lons[-1][np.newaxis,:]), axis = 0)

# Use periodic boundary conditions in the longitude direction for the starting grid
Z    = np.concatenate((Z,Z,Z), axis = 1)
lats = np.concatenate((lats,lats,lats), axis = 1)
lons = np.concatenate((lons-360,lons,lons+360), axis = 1)

# Load the grid we like to interpolate to
gridfile = Dataset('new_grid_coordinates_pop_tx0.1_38ma.nc')
gridlon  = gridfile['T_LON_2D'][:]
gridlat  = gridfile['T_LAT_2D'][:]

# Concatenate the gridlon and gridlat to one array:
points = np.concatenate((lons.flatten()[::1][:,np.newaxis], lats.flatten()[::1][:,np.newaxis]), axis=1)



#%%
# Perform the interpolation 
grid_z = griddata(points, Z.flatten()[::1], (gridlon, gridlat), method='linear')

# Write bathymetry to netcdf file:
dataset = Dataset('adapted_bathymetry.nc', 'w')

i_indexs = dataset.createDimension('i_index', 3600)
j_indexs = dataset.createDimension('j_index', 2400+150)
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

baths = dataset.createVariable('Bathymetry', np.float32,('j_index','i_index',))

# Write data
ins[:]         = gridfile['i_index'][:]
jns[:]         = gridfile['j_index'][:]
depth_tns[:]   = gridfile['depth_t'][:]
w_depns[:]     = gridfile['w_dep'][:]
latitudes[:]   = gridfile['T_LAT_2D'][:]
longitudes[:]  = gridfile['T_LON_2D'][:]
latitudes2[:]  = gridfile['U_LAT_2D'][:]
longitudes2[:] = gridfile['U_LON_2D'][:]

baths[:] = grid_z

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