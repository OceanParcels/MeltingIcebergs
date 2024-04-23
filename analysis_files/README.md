
## MeltingIcebergs > analysis_files

#### Main Code
## :page_with_curl: Description
During (some of) the simulations, values for the ocean salinity and Monin-Obukhov length are required. The files in this folder can be used to estimate those values based on the forcing model fields or previous studies.

## Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* numpy (1.23.3)
* shapely (1.8.4)
* xarray (0.20.1)


## :file_folder: Folder Structure
This folder contains the code to analyse specific model components (Monin-Obukhov length and salinity) in preparation for the simulations.

```bash
|   BasalMeltSensitivity.ipynb
|   Salinity.ipynb

```

## :computer: Code Usage
1) BasalMeltSensitivity.ipynb: Reproduce the figures from FitzMaurice & Stern (2017) to estimate the order of size value of the Monin-Obukhov length ($L_o$) to use in the simulations using the adapted basal melt equation.
2) Salinity.ipynb: Determine an average salinity in the Weddell Sea region from a single forcing model salinity field (model day 380101) to use in the simulations for freezing temperature ($T_f$) calculations.
