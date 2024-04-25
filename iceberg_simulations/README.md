## MeltingIcebergs > iceberg_simulations

## :page_with_curl: Description
This folder contains the files to run the Eocene, modern and passive simulations. For the Eocene, one can select one or both (depth-integrated and/or surface-only) of the files depending on the intended simulation(s) to perform. Both simulations use model data from the eddy-resolving late Eocene model by Nooteboom et al. (2022) as forcing (daily fields of ocean temperature and velocity components for model years 38 to 42, and monthly fields for the surface wind stress components).

For a short modern simulation, the 'modern_icebergs_integrated.py' file is needed. These simulations use the Mercator Ocean International (MOi) hydrodynamics dataset (Gasparin et al., 2018) for ocean temperature and velocity components. For the surface wind components, ERA5 reanalysis data is used. These simulations cover only the year 2021 and are always depth-integrated.

Simulations with passive surface particles ('passive_particles.py') were performed for the same timespan and model fields as the Eocene simulations.

## ⚙️ Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* parcels (2.4)
* numpy (1.23.3)
* pandas (1.4.4)
* xarray (0.20.1)


## :file_folder: Folder Structure
This folder contains the model code to simulate late Eocene iceberg trajectories.

```bash
│   eocene_icebergs_integrated.py
│   eocene_icebergs_surface.py
│   modern_icebergs_integrated.py
│   passive_particles.py
```

## :computer: Code Usage
1) eocene_icebergs_integrated.py: Model code for simulating forward or backward iceberg trajectories during the late Eocene using depth-integrated model forcing and iceberg melt calculations.
2) eocene_icebergs_surface.py: Model code for simulating forward or backward iceberg trajectories during the late Eocene using only surface fields for forcing and iceberg melt calculations.
3) modern_icebergs_integrated.py: Model code for simulating short, forward iceberg trajectories during the year 2021 using depth-integrated model forcing and iceberg melt calculations at different temporal resolutions of wind forcing.
4) passive_particles.py: Model code for simulating forward or backward passive particle trajectories during the late Eocene using only surface fields for forcing.
