## MeltingIcebergs > preparation_files

#### Main Code
## :page_with_curl: Description
This folder contains the files used to construct the model grid and bathymetry and define the iceberg release locations and coastal regions.

## Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* netCDF4 (1.5.7)
* numpy (1.23.3)
* scipy (1.9.1)
* shapely (1.8.4)
* xarray (0.20.1)


## :file_folder: Folder Structure
This folder contains the model code used to prepare some basic settings for the simulations.

```bash
|   adapted_make_bathymetry.py
|   adapted_make_grid.py
|   make_regions.ipynb
|   ReleaseLocations.ipynb

```

## :computer: Code Usage
1) adapted_make_bathymetry.py: Interpolate the low-resolution bathymetry (TopoBathy38.nc) from Baatsen et al. (2020) to the high-resolution grid (new_grid_coordinates_pop_x0.1_38ma.nc) and create a high-resolution bathymetry file 'adapted_bathymetry.nc'. Adapted from a script of Nooteboom et al. (2022).
2) adapted_make_grid.py*: Align the high-resolution grid (grid_coordinates_pop_tx0.1_38ma.nc) and create the file 'new_grid_coordinates_pop_tx0.1_38ma.nc'. Adapted from Nooteboom et al. (2022).
3) make_regions.ipynb: Define a field containing the coastal regions from Carter et al. (2017) using the defined forward release locations and high-resolution bathymetry to create a file 'bathymetry_regions.py'.
4) ReleaseLocations.ipynb: Define the position of ODP Site 696 and the forward release locations along the 500m bathymetry line based on the coastal regions defined by Carter et al. (2017). Define backward release locations within the grid cell of ODP Site 696.

	
*Depending on the used depth range of the forcing model, the depth range of the model grid should be adapted to match. In the case of the late Eocene simulations, this was done using: ncks -d depth_t,,730. -d w_dep,,815. new_grid_coordinates_pop_tx0.1_38ma.nc edited_grid_coordinates_pop_tx0.1_38ma.nc
