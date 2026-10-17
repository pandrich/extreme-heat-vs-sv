# Extreme heat and sexual violence

This repository contains the code base used for the analysis reported in the manuscript `Extreme heat and sexual violence victimisation amongst children, adolescents and youth in Southern Africa`.

## Data availability

Open access temperature and Universal Thermal Climate Index (UTCI) data collected from the [Climate Data Store](https://cds.climate.copernicus.eu/) is stored for ease of use at [this]() data repository. To run the code, unzip the file into the `data/raw` folder.

Administrative areas shapefiles were collected from the [DHS Program Spatial Repository] (https://spatialdata.dhsprogram.com/boundaries).

Access to the Violence Against Children Survey data can be requested through the Together for Girls [portal](https://www.togetherforgirls.org/en/analyzing-public-vacs-data). Note that we specifically asked for the information on the rural/urban status of the respondent location to be included in the datasets.

## Run the code

The analysis is structured in a series of notebooks covering data processing and modeling steps (as described below).
The `src` folder contains files where we define helper functions used for data manipulation and plotting.

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

Below we provide a brief description of the notebooks used to develop the analysis. These are divided into:
- Code we used to pre-process the data (`data_processing` subfolder).
- Modelling code (`modelling` subfolder).
- Code we used to explore the results (`result_analysis` subfolder).

#### Data processing

These notebooks should be run in the order described below to create the relevant intermediary datasets.

| Name | Description | 
|--------|---------------|
| 1_preprocess_vacs_data.qmd | Coding of VACS variables.  Requires VACS survey data with urban/rural classification. |
| 2_create_temp_datasets.qmd | Builds the longitudinal daily temperature dataset. Covers the 2000-2020 time period.  |
| 3_create_utci_datasets.qmd | Builds the longitudinal daily UTCI dataset. | 
| 4_preprocess_temp_utci_data.qmd | Create dataset with the 12-month look-back heat variables. Some visualizations are added for the heat variables exploration.  |
| 5_combine_heat_viol.qmd| Combine the violence survey and heat datasets. From the heat datasets the dates that match the survey dates are selected.  |

#### Modelling 

The notebooks below can be used to replicate and explore all models described in the manuscript.

| Name | Description | 
|--------|---------------|
| temp_above_1sd_np_viol.qmd | Runs the primary exposure models with non-partner sexual violence as the target variable.  | 
| temp_above_1sd_ip_viol.qmd | Runs the primary exposure models with intimate partner sexual violence as the target variable. |
| temp_above_1sd_ip_viol_only_live_with_part.qmd | Runs the primary exposure models with intimate partner sexual violence as the target variable, for the subset of respondents who currently live with a partner. |
| temp_above_1sd_nrel_viol.qmd | Runs the primary exposure models with sexual violence by a non-relative perpetrator as the target variable. | 
| temp_above_1sd_rel_viol.qmd| Runs the primary exposure models with sexual violence by a relative perpetrator as the target variable.|
| temp_above_1sd_rel_viol_only_live_with_rel.qmd| Runs the primary exposure models with sexual violence by a relative perpetrator as the target variable, for the subset of respondents who currently live with a relative.|
| temp_above_1sd_in_viol.qmd | Runs the primary exposure models with sexual violence that occurred in a family setting as the target variable. | 
| temp_above_1sd_out_viol.qmd | Runs the primary exposure models with sexual violence that occurred outside the family setting as the target variable. | 
| max_roll_temp_np_viol.qmd | Runs the models with non-partner sexual violence as the target variable and an alternative definition of high temperature exposure based on the mean maximum daily temperature over 3-day rolling windows. See the manuscript and the function `calculate_temp_12m_stats` in `src/data_processing.R` for more details. | 
| utci_above_1sd_np_viol.qmd | Runs the models with non-partner sexual violence as the target variable and an alternative definition of high heat exposure based on the Universal Thermal Comfort Index. See the manuscript and the function `calculate_utci_12m_stats` in `src/data_processing.R` for more details. | 

#### Result Analysis

We use these notebooks to compare the results of models using different exposure and target variables and to derive estimates for the expected increased victimisation resulting from high heat exposure.

| Name | Description | 
|--------|---------------|
| compare_results.qmd | Compare results across models with different exposure and target variables. |
| calculate_victimisation_statistics.qmd | Calculate estimates for table Ⅰ in the manuscript. |


