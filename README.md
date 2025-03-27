
# Project Title

Global Groundwater Recharge Model

[![DOI](https://zenodo.org/badge/797320852.svg)](https://zenodo.org/doi/10.5281/zenodo.13222685)

# Project Description

The Global Groundwater Recharge (GGR) model is a Python-based grid model designed to simulate daily rain-fed groundwater recharge using satellite imagery and environmental parameters. Comprising three layers - topsoil (root zone), subsoil, and aquifer (refer to Fig. 1) - the GGR model computes various hydrological processes, including water exchange between topsoil and atmosphere, surface runoff, topsoil recharge, water volume in soil layers, subsoil recharge from topsoil, capillary rise from subsoil to topsoil, and groundwater recharge. These calculations are performed on a daily time step and grid-based values.

The default spatial and temporal setup of the model is as follows: spatial resolution of 0.1°×0.1° and daily temporal resolution. The model encompasses a spatial extent ranging from 180.0°W to 180.0°E longitudes and 60.0°N to 60.0°S latitudes, with a temporal range spanning from January 2001 to December 2020. 

![Figure 1](https://github.com/Global-Groundwater-Model/Global_Groundwater_Recharge_Model/blob/main/figures/Figure1.png)

Fig. 1. The conceptual model of Global Groundwater Recharge (GGR). 


## Prerequisite
- [numpy](https://numpy.org/install/)
- [xarray](https://docs.xarray.dev/en/latest/getting-started-guide/installing.html)
- [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
- [geopandas](https://geopandas.org/en/stable/getting_started/install.html)
- [shapely](https://shapely.readthedocs.io/en/stable/installation.html)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)


## Getting Started
To run the Global Groundwater Recharge (GGR) model, you'll need a set of input data. Below is a table listing all required input data along with their respective sources.


Table. 1. Overview of the Input Parameters Used in the GGR Model with the Abbreviation, Spatial and Temporal Resolutions, Unit, and Source.
![Table 1](https://github.com/Global-Groundwater-Model/Global_Groundwater_Recharge_Model/blob/main/figures/Table1.png)
*These inputs must be rescaled to 0.1 degree resolution. 
**These inputs are available in the Input folder. 

To execute the GGR model, you'll need to utilize the Snakemake library in Python. The model workflow is orchestrated by the Snakefile, which automates the execution of various tasks.

Before running the Snakefile, ensure you have prepared the Input_Info.csv file. This file contains the necessary input data sources. Refer to Table 1 for a comprehensive list of input data sources and their abbreviations.

Additionally, you can visualize the workflow of the GGR model by examining the Snakemake rule graph, as illustrated in Figure. 2. This graph provides insights into the sequence of code execution within the GGR Model.

![Figure 2](https://github.com/Global-Groundwater-Model/Global_Groundwater_Recharge_Model/blob/main/figures/Figure2.png)

## Documentation

For more detailed information on the GGR model and to cite the model, please refer to the following article:

Nazari, S., Kurse, I., Moosdorf, N., "Spatiotemporal dynamics of global rain-fed groundwater recharge from 2001 to 2020. Journal of Hydrology. 650:132490.

Please consider citing this article if you use the GGR model in your research or projects.

## References to Input Data

- [Amante, C., & Eakins, B. (2009). ETOPO1 1 arc-minute global relief model: procedures, data sources and analysis. NOAA technical memorandum NESDIS NGDC-24. National Geophysical Data Center, NOAA, 10(2009), V5C8276M.](https://doi.org/10.7289/V5C8276M) 

- [Hagemann, S., & Gates, L. D. (2003). Improving a subgrid runoff parameterization scheme for climate models by the use of high resolution data derived from satellite observations. Climate Dynamics, 21(3), 349-359.](https://doi.org/10.1007/s00382-003-0349-x)
- [Hengl, T., Mendes de Jesus, J., Heuvelink, G. B., Ruiperez Gonzalez, M., Kilibarda, M., Blagotić, A., Shangguan, W., Wright, M. N., Geng, X., & Bauer-Marschallinger, B. (2017). SoilGrids250m: Global gridded soil information based on machine learning. PloS one, 12(2), e0169748.](https://doi.org/10.1371/journal.pone.0169748) 

- [Huffman, G. J., E.F. Stocker, D.T. Bolvin, E.J. Nelkin, Jackson Tan (2019). GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree V06, Edited by Andrey Savtchenko, Greenbelt, MD, Goddard Earth Sciences Data and Information Services Center (GES DISC).](https://doi.org/10.5067/GPM/IMERGDF/DAY/06)

- [Muñoz Sabater, J. (2019). ERA5-Land hourly data from 1950 to present.](https://doi.org/10.24381/cds.e2161bac)

- [Simons, G., Koster, R., & Droogers, P. (2020). HiHydroSoil v2. 0-High Resolution Soil Maps of Global Hydraulic Properties FutureWater report 134, Wageningen, The Netherlands.](https://www.futurewater.eu/projects/hihydrosoil/) 

- [Stacke, T., & Hagemann, S. (2021). HydroPy (v1. 0): A new global hydrology model written in Python. Geoscientific Model Development Discussions, 1-28.](https://zenodo.org/records/4541239) 


