## MeltingIcebergs > eocene_icebergs

#### Main Code
## :page_with_curl: Description
This folder contains the files to run the Eocene simulations. One can select one or both (depth-integrated and/or surface-only) of the files depending on the intended simulation(s) to perform. Both simulations use model data from the eddy-resolving late Eocene model by Nooteboom et al. (2022) as forcing (daily fields of ocean temperature and velocity components for model years 38 to 42, and monthly fields for the surface wind stress components).

## Dependencies
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

```

## :computer: Code Usage
1) eocene_icebergs_integrated.py: Model code for simulating forward or backward iceberg trajectories during the late Eocene using depth-integrated model forcing and iceberg melt calculations.
2) eocene_icebergs_surface.py: Model code for simulating forward or backward iceberg trajectories during the late Eocene using only surface field for forcing and iceberg melt calculations.
