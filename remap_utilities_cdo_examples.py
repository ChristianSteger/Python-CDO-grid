# Load modules
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import subprocess
import cartopy.crs as ccrs
import glob

mpl.style.use("classic")

# Path to folder
path_work = os.getenv("HOME") + "/Desktop/work_dir/"

# Load required functions
sys.path.append(os.getenv("HOME") + "/Desktop/Python_CDO_grid/Functions/")
import remap_utilities_cdo
from grid_utilities import gridcoord, gridframe

###############################################################################
# Source grid: regular latitude/longitude grid
###############################################################################

# Get source grid information
file_in = path_work + "ERA5-Land_geopotential.nc"
ds = xr.open_dataset(file_in)
ds = ds.sel(longitude=slice(0.0, 40.0), latitude=slice(65.0, 30.0))
lon_cent_in = ds["longitude"].values
lat_cent_in = ds["latitude"].values
ds.close()
lon_edge_in, lat_edge_in = gridcoord(lon_cent_in, lat_cent_in)

# -----------------------------------------------------------------------------
# Target grid: regular latitude/longitude grid
# -----------------------------------------------------------------------------

# Define output grid
lon_cent = np.linspace(5.0, 20.0, 301)
lat_cent = np.linspace(42.0, 52.0, 201)
lon_edge, lat_edge = gridcoord(lon_cent, lat_cent)

# Check domain coverage
geo_crs = ccrs.PlateCarree()
plt.figure()
ax = plt.axes(projection=geo_crs)
lon_frame, lat_frame = gridframe(lon_edge_in, lat_edge_in, offset=0)
poly = plt.Polygon(list(zip(lon_frame, lat_frame)), facecolor="red",
                   edgecolor="red", alpha=0.2, linewidth=1.0)
ax.add_patch(poly)
lon_frame, lat_frame = gridframe(lon_edge, lat_edge, offset=0)
poly = plt.Polygon(list(zip(lon_frame, lat_frame)), facecolor="none",
                   edgecolor="blue", alpha=1.0, linewidth=2.5)
ax.add_patch(poly)
ax.coastlines(resolution="50m", color="black", linewidth=1.0)
ax.axis([lon_cent_in[0] - 1.0, lon_cent_in[-1] + 1.0,
         lat_cent_in.min() - 1.0, lat_cent_in.max() + 1.0])

# Remap with compact CDO grid description file
file_txt = path_work + "grid_target_lonlat.txt"
remap_utilities_cdo.griddesc(file_txt, gridtype="lonlat",
                             xsize=len(lon_cent), ysize=len(lat_cent),
                             xfirst=lon_cent[0], yfirst=lat_cent[0],
                             xinc=np.diff(lon_cent).mean(),
                             yinc=np.diff(lat_cent).mean())
for i in ("remapbil", "remapcon"):
    cmd = "cdo " + i + "," + file_txt
    sf = file_in
    tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_" + i + ".nc"
    subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Alternative method: Remap with NetCDF CDO grid description file
# -> this method is useful if (I) CDO does not understand the selected grid
#    description of if (II) user-specific grid specification are required
lon_cent_2d, lat_cent_2d = np.meshgrid(lon_cent, lat_cent)
lon_edge_2d, lat_edge_2d = np.meshgrid(lon_edge, lat_edge)
file_netcdf = path_work + "grid_target_lonlat.nc"
remap_utilities_cdo.griddesc_netcdf(file_netcdf,
                                    lon_cent=lon_cent_2d, lat_cent=lat_cent_2d,
                                    lon_edge=lon_edge_2d, lat_edge=lat_edge_2d)
for i in ("remapbil", "remapcon"):
    cmd = "cdo " + i + "," + file_netcdf
    sf = file_in
    tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_" + i + "_netcdf.nc"
    subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Remove files
files_rm = glob.glob(path_work + "*geopotential_lonlat*.nc")
for i in files_rm:
    os.remove(i)

# -----------------------------------------------------------------------------
# Target grid: rotated latitude/longitude grid
# -----------------------------------------------------------------------------

# Define output grid
pole_latitude = 42.5
pole_longitude = -160.0
rlon_cent = np.linspace(-7.5, 7.5, 301)
rlat_cent = np.linspace(-5.0, 5.0, 201)
rlon_edge, rlat_edge = gridcoord(rlon_cent, rlat_cent)
rot_pole_crs = ccrs.RotatedPole(pole_latitude=pole_latitude,
                                pole_longitude=pole_longitude)

# Check domain coverage
plt.figure()
ax = plt.axes(projection=geo_crs)
lon_frame, lat_frame = gridframe(lon_edge_in, lat_edge_in, offset=0)
poly = plt.Polygon(list(zip(lon_frame, lat_frame)), facecolor="red",
                   edgecolor="red", alpha=0.2, linewidth=1.0)
ax.add_patch(poly)
rlon_frame, rlat_frame = gridframe(rlon_edge, rlat_edge, offset=0)
poly = plt.Polygon(list(zip(rlon_frame, rlat_frame)), facecolor="none",
                   edgecolor="blue", alpha=1.0, linewidth=2.5,
                   transform=rot_pole_crs)
ax.add_patch(poly)
ax.coastlines(resolution="50m", color="black", linewidth=1.0)
ax.axis([lon_cent_in[0] - 1.0, lon_cent_in[-1] + 1.0,
         lat_cent_in.min() - 1.0, lat_cent_in.max() + 1.0])

# Remap with compact CDO grid description file
file_txt = path_work + "grid_target_rot.txt"
remap_utilities_cdo.griddesc(file_txt, gridtype="projection",
                             xsize=len(rlon_cent), ysize=len(rlat_cent),
                             xfirst=rlon_cent[0], yfirst=rlat_cent[0],
                             xinc=np.diff(rlon_cent).mean(),
                             yinc=np.diff(rlat_cent).mean(),
                             grid_mapping_name="rotated_latitude_longitude",
                             grid_north_pole_longitude=pole_longitude,
                             grid_north_pole_latitude=pole_latitude)
for i in ("remapbil", "remapcon"):
    cmd = "cdo " + i + "," + file_txt
    sf = file_in
    tf = path_work + sf.split("/")[-1][:-3] + "_rot_" + i + ".nc"
    subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Remove files
files_rm = glob.glob(path_work + "*geopotential_rot*.nc")
for i in files_rm:
    os.remove(i)

# -----------------------------------------------------------------------------
# Target grid: map projection
# -----------------------------------------------------------------------------

# Define output grid
lc_crs = ccrs.LambertConformal(central_longitude=20.0, central_latitude=47.5,
                               standard_parallels=(41.5, 53.5))
x_cent = np.arange(-800000.0, 804000.0, 4000.0)
y_cent = np.arange(-400000.0, 404000.0, 4000.0)
x_edge, y_edge = gridcoord(x_cent, y_cent)

# Check domain coverage
plt.figure()
ax = plt.axes(projection=geo_crs)
lon_frame, lat_frame = gridframe(lon_edge_in, lat_edge_in, offset=0)
poly = plt.Polygon(list(zip(lon_frame, lat_frame)), facecolor="red",
                   edgecolor="red", alpha=0.2, linewidth=1.0)
ax.add_patch(poly)
rlon_frame, rlat_frame = gridframe(x_edge, y_edge, offset=0)
poly = plt.Polygon(list(zip(rlon_frame, rlat_frame)), facecolor="none",
                   edgecolor="blue", alpha=1.0, linewidth=2.5,
                   transform=lc_crs)
ax.add_patch(poly)
ax.coastlines(resolution="50m", color="black", linewidth=1.0)
ax.axis([lon_cent_in[0] - 1.0, lon_cent_in[-1] + 1.0,
         lat_cent_in.min() - 1.0, lat_cent_in.max() + 1.0])

# Remap with compact CDO grid description file
file_txt = path_work + "grid_target_proj.txt"
remap_utilities_cdo.griddesc(file_txt, gridtype="projection",
                             xsize=len(x_cent), ysize=len(y_cent),
                             xfirst=x_cent[0], yfirst=y_cent[0],
                             xinc=np.diff(x_cent).mean(),
                             yinc=np.diff(y_cent).mean(),
                             grid_mapping_name="lambert_conformal_conic",
                             proj_params=lc_crs.proj4_init)
for i in ("remapbil", "remapcon"):
    cmd = "cdo " + i + "," + file_txt
    sf = file_in
    tf = path_work + sf.split("/")[-1][:-3] + "_proj_" + i + ".nc"
    subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Remove files
files_rm = glob.glob(path_work + "*geopotential_proj*.nc")
for i in files_rm:
    os.remove(i)

###############################################################################
# Other source grids
###############################################################################

# -----------------------------------------------------------------------------
# Source grid: rotated latitude/longitude grid
# -----------------------------------------------------------------------------

# Get source grid information
file_in = path_work + "orog_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_"\
          + "KNMI-RACMO22E_v1_fx.nc"
ds = xr.open_dataset(file_in)
rlon_cent_in = ds["rlon"].values
rlat_cent_in = ds["rlat"].values
grid_mapping_name = ds["rotated_pole"].grid_mapping_name
pole_longitude = ds["rotated_pole"].grid_north_pole_longitude
pole_latitude = ds["rotated_pole"].grid_north_pole_latitude
ds.close()

# Write grid description file for source grid
file_txt = path_work + "grid_source_rot.txt"
remap_utilities_cdo.griddesc(file_txt, gridtype="projection",
                             xsize=len(rlon_cent_in), ysize=len(rlat_cent_in),
                             xfirst=rlon_cent_in[0], yfirst=rlat_cent_in[0],
                             xinc=np.diff(rlon_cent_in).mean(),
                             yinc=np.diff(rlat_cent_in).mean(),
                             grid_mapping_name="rotated_latitude_longitude",
                             grid_north_pole_longitude=pole_longitude,
                             grid_north_pole_latitude=pole_latitude)

# Bilinear interpolation
cmd = "cdo remapbil," + path_work + "grid_target_lonlat.txt"
sf = file_in
tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_remapbil.nc"
subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Conservative interpolation
# -> coordinates of grid cell edges are not provided in NetCDF file of source
#    grid. The above created CDO grid description file is thus used to add
#    this information with 'setgrid'
cmd = "cdo remapcon," + path_work + "grid_target_lonlat.txt"
sf = "-setgrid," + file_txt + " " + file_in
tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_remapcon.nc"
subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Remove files
files_rm = glob.glob(path_work + "*KNMI-RACMO22E_v1_fx_lonlat*.nc")
for i in files_rm:
    os.remove(i)

# -----------------------------------------------------------------------------
# Source grid: map projection
# -----------------------------------------------------------------------------

# Conservative interpolation of CNRM-ALADIN63 data
file_in = path_work + "orog_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_"\
          + "CNRM-ALADIN63_v1_fx.nc"
cmd = "cdo remapcon," + path_work + "grid_target_lonlat.txt"
sf = file_in
tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_remapcon.nc"
subprocess.call(cmd + " " + sf + " " + tf, shell=True)
# -> works because grid cell edges are defined

# Conservative interpolation of ICTP-RegCM4-6 data
file_in = path_work + "orog_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_"\
          + "ICTP-RegCM4-6_v1_fx.nc"
ds = xr.open_dataset(file_in)
x_cent_in = ds["x"].values
y_cent_in = ds["y"].values
proj4_params = ds["crs"].proj4_params
ds.close()
# -> requires CDO setgrid command because coordinates of grid cell edges are
#    missing

file_txt = path_work + "grid_source_proj.txt"
remap_utilities_cdo.griddesc(file_txt, gridtype="projection",
                             xsize=len(x_cent_in), ysize=len(y_cent_in),
                             xfirst=x_cent_in[0], yfirst=y_cent_in[0],
                             xinc=np.diff(x_cent_in).mean(),
                             yinc=np.diff(y_cent_in).mean(),
                             grid_mapping_name="lambert_conformal_conic",
                             proj_params=proj4_params)

cmd = "cdo remapcon," + path_work + "grid_target_lonlat.txt"
sf = "-setgrid," + file_txt + " " + file_in
tf = path_work + sf.split("/")[-1][:-3] + "_lonlat_remapcon.nc"
subprocess.call(cmd + " " + sf + " " + tf, shell=True)

# Remove files
files_rm = glob.glob(path_work + "*lonlat_remapcon.nc")
for i in files_rm:
    os.remove(i)

# Notes
# NetCDF file with remapping weights can be generated if a lot of data
# on the same grid has to be remapped:
# cdo gencon,grid_target.txt input.nc weights.nc
# cdo remap,grid_target.txt,weights.nc input.nc output.nc
