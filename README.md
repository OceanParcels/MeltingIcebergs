# MeltingIcebergs
❗ Note: This is the official code repository for Mark Elbertsen's MSc project on melting icebergs in the Eocene.

#### Main Code
## :page_with_curl: Description
This repository contains the code for Mark Elbertsen's MSc project used to perform depth-integrated Lagrangian iceberg tracing around Antarctica during the late Eocene. Using the OceanParcels framework, iceberg melting (or growth) is simulated using several kernels, including for the dominant iceberg melt terms: basal melt, buoyant convection and wave erosion. By defining kernels for five iceberg size classes, the model can be used to determine the minimum iceberg size required for icebergs to survive the late Eocene warmth.

## Dependencies
The model was run using Python 3.8.13.

The simulations (*.py files 'modern_icebergs' and those stored in 'eocene_icebergs' folder) have the following dependencies:
* parcels (2.4)
* numpy (1.23.3)
* pandas (1.4.4)
* xarray (0.20.1)

The files used in preparation for some model files (*.py and *.ipynb files in 'preparation_files' folder) have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* netCDF4 (1.5.7)
* numpy (1.23.3)
* scipy (1.9.1)
* shapely (1.8.4)
* xarray (0.20.1)

The files used to analyse some inital model files (*.ipynb files in 'analysis_files' folder) have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* numpy (1.23.3)
* shapely (1.8.4)
* xarray (0.20.1)

Finally, the files used for plotting (*.ipynb files in 'plotting_files' folder) have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* numpy (1.23.3)
* pandas (1.4.4)
* shapely (1.8.4)
* xarray (0.20.1)


## :file_folder: Folder Structure
1) preparation_files: Model files used to prepare the simulations.
2) analysis_files: Analyse specific model components.
3) plotting_files: Notebooks used to generate the plots and animations used in the thesis.
4) eocene_icebergs: Model code to simulate iceberg trajectories.
5) icebergs: Model code (kernels) to initialise and adapt the icebergs throughout the simulation.

```bash
│   LICENSE
│   modern_icebergs_integrated.py
│   output.doc
│   outputfolder.txt
│   passive_particles.py
│   README.md
│
├───analysis_files
│       BasalMeltSensitivity.ipynb
│       README.txt
│       Salinity.ipynb
│
├───eocene_icebergs
│       eocene_icebergs_integrated.py
│       eocene_icebergs_surface.py
│
├───icebergs
│       icebergs_kernels.py
│       icebergs_particleclass.py
│
├───plotting_files
│       AnimationPreparation.ipynb
│       README.txt
│       ThesisFigures.ipynb
│
└───preparation_files
        adapted_make_bathymetry.py
        adapted_make_grid.py
        make_regions.ipynb
        README.txt
        ReleaseLocations.ipynb
```

## :computer: Code Usage (Please be a little bit more specific how the code can be used from others + if you like to follow this description approach it would be better to do the same for all the code including the jupyter notebooks) 
All simulations require the files in the 'icebergs' folder to define the kernels and iceberg particleclass used in the model. In addition, the files in the 'preparation_files' folder were used to construct the Eocene model grid and bathymetry used for the forcing model of Nooteboom et al. (2022), and define iceberg release locations and coastal regions used in the simulations.

To run the Eocene simulations, one can select one or download both (depth-integrated and/or surface-only) of the files in the 'eocene_icebergs' folder.

For a modern (test) simulation, the 'modern_icebergs_integrated.py' file is required. uses MOi and ERA5 data.


1) passive_particles.py: Simulate passive particle trajectories around Antarctica.
2) icebergs_kernels.py: Parcels kernels for the iceberg simulations.
3) icebergs_particleclass.py: Parcels particle classes for the iceberg simulations.
4) eocene_icebergs_surface.py: Simulate Eocene icebergs using only surface fields.
5) eocene_icebergs_integrated.py: Simulate Eocene icebergs using depth-integrated fields.
6) modern_icebergs_integrated.py: Simulate icebergs in the present-day using depth-integrated fields.

#### Folders
1) preparation_files: Model files used to prepare the simulations.
2) analysis_files: Analyse specific model components.
3) plotting_files: Notebooks used to generate the plots and animations used in the thesis.

## :envelope_with_arrow: Contact and contribution
For questions about this repository, please contact the author(s) **{fill the corresponding person}**, or open an Issue or Pull request in this repository.

## :balance_scale: License
This repository is licensed under a MIT License. You can view the [LICENSE here](https://github.com/AristotleKandylas/MeltingIcebergs_rev/blob/main/LICENSE)

## :bookmark: Cite this repository
...To be added after Publication on Zenodo...

## Funding
This research is supported by the ERC Starting Grant 802835 (OceaNice) to Peter K. Bijl.

## Acknowledgements
We thank Michael Kliphuis for assisting with and management of the output data. We thank Anna von der Heydt and Peter Nooteboom for providing the forcing model data, and Peter Nooteboom for assisting with the initial model set-up. We also want to thank Michael Baatsen for providing climate index data used during the analysis of the results ('PaperFiguresFinal-Other.ipynb').

