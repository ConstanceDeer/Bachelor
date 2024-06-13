# FINAL PLOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import os
from numpy.ma import masked_where
from scipy.ndimage import label
import re

directory_path = 'C:/Users/const/Documents/Bachelor/all tracks'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]


fig = plt.figure()
 
ax = fig.add_subplot( projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
ax.add_feature(cartopy.feature.LAND)
gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,alpha=0.5,linestyle="--")
gl.top_labels = False
gl.right_labels = False
ax.add_feature(cartopy.feature.BORDERS)

for i in range(len(file_paths)):
#for i in range(5):
    fn = file_paths[i]

    # Extracting the '011' substring using regular expression
    name = os.path.splitext(os.path.basename(fn))[0].split('_')[6]
    name = str(int(name))
    # Load NC datafile
    load_file = nc.Dataset(fn, "r")

    # Assuming your NetCDF file has variables 'latitude' and 'longitude'
    latitudes = load_file.variables["latitude"][:]
    longitudes = load_file.variables["longitude"][:]
    longitudes[longitudes > 180] -= 360

    # Close the NetCDF file
    load_file.close()

    # Plot data with label
    ax.plot(longitudes.flatten(), latitudes.flatten(), '*', markersize=0.15)

    # Find the index where latitude is closest to zero
    idx = np.abs(latitudes.flatten()).argmin()

    # Annotate the plot at latitude 0 with legend names
    ax.text(longitudes.flatten()[idx], latitudes.flatten()[idx], name, horizontalalignment='center', verticalalignment='center', fontsize=7, bbox=dict(facecolor='white', alpha=1,pad=0.4, edgecolor='none'))


ax.set_extent([-180, 180,80, -80], crs=ccrs.PlateCarree())
fig = plt.gcf()
fig.set_dpi(800)
plt.show()   


#%%
import re
directory_path = 'C:/Users/const/Documents/Bachelor/all tracks'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]


fig = plt.figure()
 
ax = fig.add_subplot( projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
ax.add_feature(cartopy.feature.LAND)
gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False,alpha=0.5,linestyle="--")
gl.top_labels = False
gl.right_labels = False
ax.add_feature(cartopy.feature.BORDERS)


fn = file_paths[10]

# Extracting the '011' substring using regular expression
name = os.path.splitext(os.path.basename(fn))[0].split('_')[6]

# Load NC datafile
load_file = nc.Dataset(fn, "r")
    
# Extract data variables
datakeys = load_file.variables.keys()
    
# Collect data in a python dictionary
swot = {}
for k in datakeys:
         swot[k] = load_file.variables[k][:]
    
            # Close the NetCDF file
load_file.close()
    
# Assuming your NetCDF file has variables 'latitudes', 'longitudes', and 'ssha_data'
latitudes = swot["latitude"][:]
longitudes = swot["longitude"][:]
longitudes[longitudes > 180] -= 360            
    
# Show data global
   
ax.plot(longitudes,latitudes,'*r', markersize=2, label=name)
ax.set_extent([-180, 180,80, -80], crs=ccrs.PlateCarree())
fig = plt.gcf()

ax.legend()
fig.set_dpi(800)
plt.show()   


#%%
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import netCDF4 as nc

directory_path = 'C:/Users/const/Documents/Bachelor/all tracks'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]

fig = plt.figure()
ax = fig.add_subplot(projection=ccrs.PlateCarree())
ax.coastlines(resolution='10m')
ax.add_feature(cartopy.feature.LAND)
gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False, alpha=0.5, linestyle="--")
gl.top_labels = False
gl.right_labels = False
ax.add_feature(cartopy.feature.BORDERS)

fn = file_paths[10]

# Extracting the '011' substring using regular expression
name = os.path.splitext(os.path.basename(fn))[0].split('_')[6]

# Load NC datafile
load_file = nc.Dataset(fn, "r")

# Assuming your NetCDF file has variables 'latitude' and 'longitude'
latitudes = load_file.variables["latitude"][:]
longitudes = load_file.variables["longitude"][:]

# Close the NetCDF file
load_file.close()

# Flatten latitude and longitude arrays
latitudes = latitudes.flatten()
longitudes = longitudes.flatten()
longitudes[longitudes > 180] -= 360

# Plot data with label
line, = ax.plot(longitudes, latitudes, '*r', markersize=0.5)

# Set plot extent
ax.set_extent([-180, 180, 80, -80], crs=ccrs.PlateCarree())

# Find the index where latitude is closest to zero
idx = np.abs(latitudes).argmin()

# Annotate the plot at latitude 0 with legend names
ax.text(longitudes[idx], latitudes[idx], name, horizontalalignment='center', verticalalignment='center', fontsize=6, bbox=dict(facecolor='white', alpha=1))


# Set DPI
fig.set_dpi(800)

# Show plot
plt.show()
