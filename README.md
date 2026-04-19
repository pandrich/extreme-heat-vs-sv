# Extreme heat and sexual violence

This repository contains the code base used for the analysis in `Extreme heat and sexual violence victimisation amongst children, adolescents and youth in Southern Africa`.

## Data availability

Temperature and Universal Thermal Climate Index (UTCI) data collected from the [Climate Data Store](https://cds.climate.copernicus.eu/) is stored for ease of use at [this]() data repository. To run the code, unzip the file into the `data/raw` folder.

Administrative areas shapefiles were collected from the [DHS Program Spatial Repository] (https://spatialdata.dhsprogram.com/boundaries).

Access to the Violence Against Children Survey data can be requested through the Together for Girls [portal](https://www.togetherforgirls.org/en/analyzing-public-vacs-data). Note that we specifically asked for the information on the rural/urban status of the respondent location to be included in the datasets.

## Run the code

The analysis is structured in a series of notebooks covering data processing and modeling steps (as described below).
The `src` folder contains the definition of helper functions used for data manipulation and plotting.

The code was created using R version 4.5.2 and was tested using RStudio version 2025.09.2+418. renv version 1.1.5 was used to manage package dependencies, documented in the `renv.lock` file.

### Setup the coding environment

Before running the code locally create the coding environments.

1. Open the project file `extreme-heat-vs-sv.Rproj` in RStudio.
2. In the console activate the project environment
	```
	renv::init()
	```
3. Select 1 to restore the environment from the existing lock file.

CmdStan and cmdstanr need to be separately installed following the instructions [here](https://mc-stan.org/cmdstanr/)


### Notebook codebook

Running the notebooks in the order described below ensures that the processed datasets required for the modelling steps are appropriately created.

#### Data processing

| Name | Description | Notes |
|--------|---------------|---------|
| 1_preprocess_vacs_data.qmd |  | |
| 2_create_temp_datasets.qmd | |  |
| 3_create_utci_datasets.qmd | | |
| 4_preprocess_temp_utci_data.qmd |  |  |
| 5_combine_heat_viol.qmd|  |  |
| |  |   |