# omz-bda

Code and figures for oxygen minimum zone (OMZ) analysis of World Ocean Atlas 2018 (WOA18) Dissolved Oxygen (DO) data documented in:

Yamane MT, Vierra W, Barton W, Miller M  
*Seasonal trends in Pacific oxygen minimum zones suggest expansion with warming waters*

for MS342 - Chemical and Physical Oceanography in Spring 2021

**data** stores output from OMZ-BDA runs; [WOA18 data](https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/bin/woa18oxnu.pl?parameter=o) as NetCDF files should be stored here

**figures** stores figures used in the study

**OMZBDA.py** holds the OMZ boundary detection algorithm (OMZ-BDA) and helper functions

**woa_analysis.ipynb** provides an example OMZ-BDA run (with status updates)

**woa_eda.ipynb** exploratory data analysis of the WOA18 DO dataset

**woa_figures.ipynb** creates figures used in the study

**woa_omzbda.ipynb** provides validation for major components of the OMZ-BDA

## Dependencies
* numpy
* netCDF4
* matplotlib (for figures)
* IPython (for animations)
