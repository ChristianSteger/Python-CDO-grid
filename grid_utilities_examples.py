# Load modules
import sys
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from shapely.geometry import Polygon
from shapely.geometry import shape
from shapely.ops import unary_union
import shapely.vectorized
import cartopy.crs as ccrs
import fiona
from descartes import PolygonPatch
from pyproj import CRS, Transformer

mpl.style.use("classic")

# Path to folders
path_work = os.getenv("HOME") + "/Desktop/work_dir/"

# Load required functions
sys.path.append(os.getenv("HOME") + "/Desktop/Python_CDO_grid/Functions/")
from grid_utilities import gridcoord, gridframe, gridpolygon, areagridcell

###############################################################################
# Load test data
###############################################################################

# Load topography, coordinates and grid information
ds = xr.open_dataset(path_work + "orog_EUR-11_ECMWF-ERAINT_evaluation_"
                     + "r1i1p1_KNMI-RACMO22E_v1_fx.nc")
rlat_cent = ds["rlat"].values  # (412)
rlon_cent = ds["rlon"].values  # (424)
pole_longitude = ds["rotated_pole"].grid_north_pole_longitude
pole_latitude = ds["rotated_pole"].grid_north_pole_latitude
orog = ds["orog"].values  # [m]

###############################################################################
# Plot grid frame
###############################################################################

# Create grid coordinates
rlon_edge, rlat_edge = gridcoord(rlon_cent, rlat_cent)

# Plot
plt.figure()
ax = plt.axes()
plt.pcolormesh(rlon_edge, rlat_edge, orog)
for i in [0, 50, 100]:
    rlon_frame, rlat_frame = gridframe(rlon_edge, rlat_edge, offset=i)
    poly = plt.Polygon(list(zip(rlon_frame, rlat_frame)), facecolor="none",
                       edgecolor="black", alpha=1.0, linewidth=2.5)
    ax.add_patch(poly)
plt.axis([-30.0, 20.0, -25.0, 25.0])

###############################################################################
# Compute intersection of grid and polygon
###############################################################################
# -> Note: computation only exactly correct for Euclidean geometry

# Get borders of certain countries
countries = ("Switzerland", "Czechia", "Austria", "Liechtenstein")
ds = fiona.open(path_work + "ne_10m_admin_0_countries/"
                + "ne_10m_admin_0_countries.shp")
geom_names = [i["properties"]["NAME"] for i in ds]
shp_geom = [shape(ds[geom_names.index(i)]["geometry"]) for i in countries]
ds.close()
shp_geom_unar = unary_union(shp_geom)  # merge all polygons
# shp_geom_unar = Polygon(shp_geom_unar
#                         .exterior.simplify(0.1, preserve_topology=True))
# optionally: simplify polygon

# Plot country borders
plt.figure()
ax = plt.axes()
for i in shp_geom:
    poly_plot = PolygonPatch(i, facecolor="blue", edgecolor="black",
                             alpha=0.5)
    ax.add_patch(poly_plot)
poly_plot = PolygonPatch(shp_geom_unar, facecolor="none",
                         edgecolor="red", alpha=1.0, lw=2.5)
ax.add_patch(poly_plot)
ax.autoscale_view()

# Transform polygon boundaries
crs_rot_pole = ccrs.RotatedPole(pole_longitude=pole_longitude,
                                pole_latitude=pole_latitude)
lon_vert, lat_vert = shp_geom_unar.exterior.coords.xy
crs_geo = ccrs.PlateCarree()
coord_rot = crs_rot_pole.transform_points(crs_geo, np.array(lon_vert),
                                          np.array(lat_vert))
rlon_vert = coord_rot[:, 0]
rlat_vert = coord_rot[:, 1]

# Compute areas of grid cells
rlon_edge_2d, rlat_edge_2d = np.meshgrid(rlon_edge, rlat_edge)
area_frac = gridpolygon(rlon_edge_2d, rlat_edge_2d, rlon_vert, rlat_vert,
                        agg_gc=np.array([50, 5]))

# Alternative method (-> checking grid cell centres)
rlon_cent_2d, rlat_cent_2d = np.meshgrid(rlon_cent, rlat_cent)
mask = shapely.vectorized.contains(Polygon(zip(rlon_vert, rlat_vert)),
                                   rlon_cent_2d, rlat_cent_2d) \
    .astype(np.float32)

# Colormap
levels = np.arange(0, 1.05, 0.05)
cmap = plt.get_cmap("YlOrRd")
norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N)

# Test plot
plt.figure()
ax = plt.axes()
plt.pcolormesh(rlon_edge, rlat_edge, area_frac, cmap=cmap, norm=norm)
# plt.pcolormesh(rlon_edge, rlat_edge, mask, cmap=cmap, norm=norm)
plt.colorbar()
poly_sh = Polygon(zip(rlon_vert, rlat_vert))
poly_plot = PolygonPatch(poly_sh, facecolor="none", edgecolor="black", lw=2.5)
ax.add_patch(poly_plot)
plt.axis([poly_sh.bounds[0] - 0.2, poly_sh.bounds[2] + 0.2,
          poly_sh.bounds[1] - 0.2, poly_sh.bounds[3] + 0.2])

###############################################################################
# Compute area of grid cells
###############################################################################

# Map-projection (Equal Area Cylindrical projection)
crs_latlong = CRS.from_proj4("+proj=latlong +ellps=sphere")
crs_eac = CRS.from_proj4("+proj=cea +ellps=sphere")
# -> equal-area map projection required: cylindrical equal-area projection
transform = Transformer.from_crs(crs_latlong, crs_eac)
x_edge, y_edge = transform.transform(rlon_edge_2d, rlat_edge_2d)

# Calculate area of grid cells and multiply with 'area fraction'
area_gc = areagridcell(x_edge, y_edge) / (1000. ** 2)  # [km2]
area_gc *= area_frac

# Check area
area_count = {"Switzerland": 41285.0, "Austria": 83879.0,
              "Czechia": 78871.0, "Liechtenstein": 160.0}  # [km2]
# Source: https://www.wikipedia.org
dev = (area_gc.sum() / sum([area_count[i] for i in countries]) - 1.0) * 100.
print("Deviation in area: " + "%.2f" % dev + " %")
