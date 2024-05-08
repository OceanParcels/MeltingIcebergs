# MeltingIcebergs
❗ Note: This is the official code repository for Mark Elbertsen's MSc project on melting icebergs in the Eocene.

## :page_with_curl: Description
This repository contains the code for Mark Elbertsen's MSc project used to perform depth-integrated Lagrangian iceberg tracing around Antarctica during the late Eocene. Using the OceanParcels framework, iceberg melting (or growth) is simulated using several kernels, including for the dominant iceberg melt terms: basal melt, buoyant convection and wave erosion. By defining kernels for five different order-of-magnitude iceberg size classes, the model can be used to determine the minimum iceberg size required for icebergs to survive the late Eocene warmth.

## ⚙️ Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* pandas (1.4.4)
* parcels (2.4)
* shapely (1.8.4)

### Environment file
To set up an environment with the dependencies directly, one could use the 'meltingicebergs.yml' file using:
```
conda env create -f meltingicebergs.yml
conda activate meltingicebergs
```


## :file_folder: Folder Structure
1) preparation_files: Model files used to prepare the simulations.
2) analysis_files: Analyse specific model components.
3) plotting_files: Notebooks used to generate the plots and animations used in the thesis.
4) iceberg_simulations: Model code to simulate iceberg trajectories.
5) iceberg_model: Model code (kernels) to initialise and adapt the icebergs throughout the simulation.

```bash
│   LICENSE
│   meltingicebergs.yml
│   output.doc
│   outputfolder.txt
│   README.md
│
├───analysis_files
│       BasalMeltSensitivity.ipynb
│       README.md
│       Salinity.ipynb
│
├───iceberg_simulations
│       eocene_icebergs_integrated.py
│       eocene_icebergs_surface.py
│       modern_icebergs_integrated.py
│       passive_particles.py
│       README.md
│
├───iceberg_model
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
All simulations require the files in the 'iceberg_model' folder to define the kernels and iceberg particle class used in the model. By adapting these files, one could, for example, define new iceberg-size classes, implement different melt parameterisations, or change the stored model variables. Some of the variables used here (salinity, Monin-Obukof length) were determined outside of the model in the files from the 'analysis_files' folder.

Using the files in the 'preparation_files' folder, the Eocene model grid and bathymetry used in the forcing model of Nooteboom et al. (2022) can be constructed. In addition, this folder contains the files to define the iceberg release locations and the coastal regions based on Carter et al. (2017) as used in the simulations.

To run the Eocene simulations, one can select one or both (depth-integrated and/or surface-only) of the files in the 'iceberg_simulations' folder. These simulations use model data from the eddy-resolving late Eocene model by Nooteboom et al. (2022) as forcing (daily fields of ocean temperature and velocity components for model years 38 to 42, and monthly fields for the surface wind stress components). For a short modern simulation, the 'modern_icebergs_integrated.py' file in the 'iceberg_simulations' folder is needed. These simulations use the Mercator Ocean International (MOi) hydrodynamics dataset (Gasparin et al., 2018) for ocean temperature and velocity components. For the surface wind components, ERA5 reanalysis data is used. These simulations cover only the year 2021. For simulations with passive particles in the late Eocene, the file 'passive_particles.py' from the folder 'iceberg_simulations' is needed.

To recreate the figures and animations shown in the paper and supplements, one can use the notebooks in the 'plotting_files' folder.


## :envelope_with_arrow: Contact and contribution
For questions about this repository, please contact the authors, [Mark Elbertsen](https://github.com/mvelbertsen) (m.v.elbertsen@uu.nl) or [Erik van Sebille](https://github.com/erikvansebille) (e.vansebille@uu.nl), or open an Issue or Pull request in this repository.

## :balance_scale: License
This repository is licensed under an MIT License. You can view the [LICENSE here](https://github.com/AristotleKandylas/MeltingIcebergs_rev/blob/main/LICENSE)

## :bookmark: Cite this repository
...To be added after Publication on Zenodo...

## Funding
This research is supported by the ERC Starting Grant 802835 (OceaNice) to [Peter K. Bijl](https://www.uu.nl/staff/PKBijl).

## Acknowledgements
We thank [Michael Kliphuis](https://github.com/michaelkliphuis) for assisting with and management of the output data. We thank [Anna von der Heydt](https://www.uu.nl/staff/ASvonderHeydt) and [Peter Nooteboom](https://github.com/pdnooteboom) for providing the forcing model data, and [Peter Nooteboom](https://github.com/pdnooteboom) for assisting with the initial model set-up. We also want to thank [Michael Baatsen](https://github.com/MichielBaatsen) for providing climate index data used during the analysis of the results.

