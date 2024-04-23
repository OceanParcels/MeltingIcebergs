## MeltingIcebergs > plotting_files

#### Main Code
## :page_with_curl: Description
This folder contains the files to analyse the data from the simulations and generate the figures and animations used in the paper and supplements.

## Dependencies
The model was run using Python 3.8.13.

The files have the following dependencies:
* cartopy (0.18.0)
* cmocean (3.0.3)
* matplotlib-inline (0.1.6)
* numpy (1.23.3)
* pandas (1.4.4)
* shapely (1.8.4)
* xarray (0.20.1)


## :file_folder: Folder Structure
This folder contains the notebooks used to generate the figures and animations used in the paper.

```bash
|   AnimationPaper.ipynb
│   PaperFigures-Backwards.ipynb
│   PaperFigures-Forwards.ipynb
│   PaperFigures-Modern.ipynb
│   PaperFigures-Other.ipynb

```

## :computer: Code Usage
1) AnimationPaper.ipynb: Prepare and compile the iceberg animation.
2) PaperFigures-Backwards.ipynb: Analyse and generate the figures from backward simulations.
2) PaperFigures-Forwards.ipynb: Analyse and generate the figures from forward simulations.
2) PaperFigures-Modern.ipynb: Analyse and generate the figures from modern simulations.
2) PaperFigures-Other.ipynb: Generate additional figures (maps).
