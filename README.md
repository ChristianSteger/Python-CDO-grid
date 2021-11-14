# Python-CDO-grid

# General Information
The module *remap_utilities_cdo.py* contains two function to facilitate CDO remapping from Python. The module *grid_utilities.py*  contains some auxiliary functions and tools to work with grid/polygon data. Application examples of the module's functions can be found in *remap_utilities_cdo_examples.py* and *grid_utilities_examples.py*.

# Getting started
The following example data is required for the scripts *remap_utilities_cdo_examples.py* and *grid_utilities_examples.py*:

- [ERA5-Land Geopotential](https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation#ERA5Land:datadocumentation-parameterlistingParameterlistings): Download *Geopotential (netCDF4)* &rarr; rename to *ERA5-Land_geopotential.nc*
- [EURO-CORDEX orographies](https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cordex-domains-single-levels?tab=form): Europe &rarr; Evaluation &rarr; 0.11 degree x 0.11 degree &rarr; Fixed &rarr; Orography &rarr; ERA-Interim &rarr; CNRM-ALADIN63, ICTP-RegCM4-6, KNMI-RACMO22E &rarr; r1i1p1 &rarr; Compressed zip file
- [Natural Earth country borders](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/): Admin 0 – Countries &rarr; *Download countries*

Next, all example data has to be moved to a working directory *work_dir*. Finally, the paths for the working directory and the functions must be adapted in *remap_utilities_cdo_examples.py* and *grid_utilities_examples.py*.

# Dependencies

The following Python packages are required:
- numpy
- xarray
- matplotlib
- cartopy
- shapely
- fiona
- descartes
- pyproj

# Support 
In case of issues or questions, please contact Christian Steger (christian.steger@env.ethz.ch).