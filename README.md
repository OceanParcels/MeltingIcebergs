# MeltingIcebergs
❗ Note: This is the official code repository for Mark Elbertsen's MSc project on melting icebergs in the Eocene.

#### Main Code
## :page_with_curl: Description
This repository contains the code for Mark Elbertsen's MSc project used to perform depth-integrated Lagrangian iceberg tracing around Antarctica during the late Eocene. Using the OceanParcels framework, iceberg melting (or growth) is simulated using several kernels, including for the dominant iceberg melt terms: basal melt, buoyant convection and wave erosion. By defining kernels for five different order-of-magnitude iceberg size classes, the model can be used to determine the minimum iceberg size required for icebergs to survive the late Eocene warmth.

## Dependencies
The model was run using Python 3.8.13.

The main model fields describing the icebergs (*.py files in 'icebergs' folder) have the following dependencies:
* parcels (2.4)
* numpy (1.23.3)

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

The files used to analyse some initial model files (*.ipynb files in 'analysis_files' folder) have the following dependencies:
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
│       README.md
│       Salinity.ipynb
│
├───eocene_icebergs
│       eocene_icebergs_integrated.py
│       eocene_icebergs_surface.py
│       README.md
│
├───icebergs
│       icebergs_kernels.py
│       icebergs_particleclass.py
│       README.md
│
├───plotting_files
│       AnimationPreparation.ipynb
│       PaperFigures-Backwards.ipynb
│       PaperFigures-Forwards.ipynb
│       PaperFigures-Modern.ipynb
│       PaperFigures-Other.ipynb
│       README.md
│
└───preparation_files
        adapted_make_bathymetry.py
        adapted_make_grid.py
        make_regions.ipynb
        README.md
        ReleaseLocations.ipynb
```

## :computer: Code Usage
All simulations require the files in the 'icebergs' folder to define the kernels and iceberg particle class used in the model. By adapting these files, one could, for example, define new iceberg-size classes, implement different melt parameterisations, or change the stored model variables. Some of the variables used here (salinity, Monin-Obukof length) were determined outside of the model in the files from the 'analysis_files' folder.

Using the files in the 'preparation_files' folder, the Eocene model grid and bathymetry used in the forcing model of Nooteboom et al. (2022) can be constructed. In addition, this folder contains the files to define the iceberg release locations and the coastal regions based on Carter et al. (2017) as used in the simulations.

To run the Eocene simulations, one can select one or both (depth-integrated and/or surface-only) of the files in the 'eocene_icebergs' folder. These simulations use model data from the eddy-resolving late Eocene model by Nooteboom et al. (2022) as forcing (daily fields of ocean temperature and velocity components for model years 38 to 42, and monthly fields for the surface wind stress components).

For a short modern simulation, the 'modern_icebergs_integrated.py' file is needed. These simulations use the Mercator Ocean International (MOi) hydrodynamics dataset (Gasparin et al., 2018) for ocean temperature and velocity components. For the surface wind components, ERA5 reanalysis data is used. These simulations cover only the year 2021.

To recreate the figures and animations shown in the paper and supplements, one can use the notebooks in the 'plotting_files' folder.


## :envelope_with_arrow: Contact and contribution
For questions about this repository, please contact the authors, Mark Elbertsen (m.v.elbertsen@uu.nl) or Erik van Sebille (e.vansebille@uu.nl), or open an Issue or Pull request in this repository.

## :balance_scale: License
This repository is licensed under an MIT License. You can view the [LICENSE here](https://github.com/AristotleKandylas/MeltingIcebergs_rev/blob/main/LICENSE)

## :bookmark: Cite this repository
...To be added after Publication on Zenodo...

## Funding
This research is supported by the ERC Starting Grant 802835 (OceaNice) to Peter K. Bijl.

## Acknowledgements
We thank Michael Kliphuis for assisting with and management of the output data. We thank Anna von der Heydt and Peter Nooteboom for providing the forcing model data, and Peter Nooteboom for assisting with the initial model set-up. We also want to thank Michael Baatsen for providing climate index data used during the analysis of the results.

