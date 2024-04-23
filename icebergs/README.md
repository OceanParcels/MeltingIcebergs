## MeltingIcebergs > icebergs

#### Main Code
## :page_with_curl: Description
This folder contains the files to define the model kernels and iceberg particle classes used. By adapting these files, one could, for example, define new iceberg-size classes, implement different melt parameterisations, or change the stored model variables.

## Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* parcels (2.4)
* numpy (1.23.3)


## :file_folder: Folder Structure
This folder contains the model code (kernels) to initialise and adapt the icebergs throughout the simulations.

```bash
|   icebergs_kernels.py
|   icebergs_particleclass.py

```

## :computer: Code Usage
1) icebergs_kernels.py: Definition of Parcels kernels (iceberg sizes, melt terms, field sampling) required for the iceberg simulations.
2) icebergs_particleclass.py: Definition of Parcels particle classes for the iceberg simulations.
