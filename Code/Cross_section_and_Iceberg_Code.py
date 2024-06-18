
#%% Wave pattern with cross-section 

import swot_ssh_utils as swot
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import os
import math
from numpy.ma import masked_where
from scipy.ndimage import label
import cartopy.feature as cfeature
from scipy.optimize import curve_fit
# skud ud: https://github.com/SWOT-community/SWOT-OpenToolkit/blob/main/src/swot_ssh_utils.py
fn = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/2_D_plot/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

# initialize the class SSH_L2
dd=swot.SSH_L2()

# print the variables in the file
# load data into the class
# the lat_bounds keyword is used to subset the data
# the return is dd.left and dd.right, which correspond to the left and right swath
dd.load_data(fn,lat_bounds=[-53.5,-51.5],)


# Create a map using cartopy
fig, ax = plt.subplots(figsize=(3, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))

# add land, rivers, lakes, and coastlines
ax.add_feature(cartopy.feature.LAND)

# add states and provinces
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

long, latt, ssha = [], [], []

# Loads data from the two tracks and compine into 3 arrays.
for tmp in [dd.left,dd.right]:
    #find data points in the open ocean (flag==0)
    #This flag is not always reliable, but it is the best we have for now
    surf=tmp.ancillary_surface_classification_flag==0

    # There is a flag ssh_karin_2_qual that can be used to filter out bad data
    # But it is not always reliable yet in the current products (Dec 1, 2023)
    # This flag is not yet used in the code below, but you can uncomment the following line to check it. 
    #qc_flag=tmp.ssh_karin_2_qual==0 # flag==0 means good data
    
    #mask out non-ocean and bad data
    dtm=np.where(surf,tmp.ssh_karin_2.data,np.nan) 
    
    dtm=dtm-np.nanmin(dtm.flatten()) # this can remove the cross-swath bias 
    
    lon,lat=tmp.longitude.data,tmp.latitude.data
    
    m=~np.isnan(lon+lat).flatten() # mask out the nan values
    
    mss=swot.get_mss22(lon,lat) # get the mean sea surface
    dtm-=mss # remove the mean sea surface 
    
    # the following line find the latitude limit
    # to bound the area for along-track mean removal
    # you can change the bounds to get a better fit for your region. 
    # The 10th column of the lat array is used as an approximate. 
    msk=(lat[:,10]>-53)&(lat[:,10]<-50)
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:] 
    
    # plot the data using scatter 

    # append values to lists
    long.append(lon.flatten()[m])
    latt.append(lat.flatten()[m])
    ssha.append(dtm.flatten()[m])
    
# convert lists to NumPy arrays
long = np.concatenate(long)
latt = np.concatenate(latt)
ssha = np.concatenate(ssha)

#plotting
cax=ax.scatter(long, latt, 
               c=ssha, s=0.05, vmin=-0.1,vmax=0.1, 
               cmap='Spectral' ,transform=ccrs.PlateCarree())
   
# add colorbar and gridlines
gl = ax.gridlines(draw_labels=True, color='gray', 
                  alpha=0.5, linestyle='--')
cbar = plt.colorbar(cax, ax=ax, label='(m)')
# remove the right side ylabels
gl.right_labels = False
#cbar.set_label('ssh_karin_2')
ax.set_title(fn.split('/')[-1].split('_')[7] + ' - SSHA 250m resolution')

ax.plot([-40.8, -42.54], [-52.8, -52.3], marker='o', linestyle='-', color='red', transform=ccrs.PlateCarree())
# Adjust figure size and DPI
fig = plt.gcf()
fig.set_size_inches(10, 6)  # Set the figure size in inches (width, height)
fig.set_dpi(300)  # Set the DPI (dots per inch)

# set the extent of the map
lonmin,lonmax,latmin,latmax=np.nanmin(lon-1),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat)
ax.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())




#%% 2-dimensional plots (Cross section)

a = -0.2873563218
b = 38.92413793
L=lambda x : a*x+b

mask=np.abs(latt-L(long)) < 0.004

lats=latt[mask]
longs=long[mask]
ssha1=ssha[mask]

# Maske for sortering basert på longs
sorted_indices = np.argsort(longs)


lats = lats[sorted_indices]
longs = longs[sorted_indices]
ssha1 = ssha1[sorted_indices]

#plotting

# Create a map using cartopy
fig, ax = plt.subplots(figsize=(3, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))

# add land, rivers, lakes, and coastlines
ax.add_feature(cartopy.feature.LAND)

# add states and provinces
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')


cax=ax.scatter(longs, lats, 
               c=ssha1, s=0.05, vmin=-0.1,vmax=0.1, 
               cmap='Spectral' ,transform=ccrs.PlateCarree())
   
# add colorbar and gridlines
gl = ax.gridlines(draw_labels=True, color='gray', 
                  alpha=0.5, linestyle='--')
cbar = plt.colorbar(cax, ax=ax)
# remove the right side ylabels
gl.right_labels = False
#cbar.set_label('ssh_karin_2')
ax.set_title(fn.split('/')[-1].split('_')[7] + ' - SSHA 250m resolution')

#ax.plot([-40.8, -42.54], [-52.8, -52.3], marker='o', linestyle='-', color='red', transform=ccrs.PlateCarree())
# Adjust figure size and DPI
fig = plt.gcf()
fig.set_size_inches(10, 6)  # Set the figure size in inches (width, height)
fig.set_dpi(300)  # Set the DPI (dots per inch)

# set the extent of the map
lonmin,lonmax,latmin,latmax=np.nanmin(lon-1),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat)
ax.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())



import numpy as np
import matplotlib.pyplot as plt
long1 = longs - 360

# Calculate median values
neighbor_count = 10
med_distances = []
med_ssha = []

for i in range(len(longs) - neighbor_count):
    med_dist = np.median(long1[i:i+neighbor_count])
    med_ssha_value = np.median(ssha1[i:i+neighbor_count])
    
    med_distances.append(med_dist)
    med_ssha.append(med_ssha_value)

# Create the plot
fig, ax = plt.subplots(figsize=(20, 6))
long1 = longs - 360
plt.plot(long1, ssha1, marker='.', linestyle='None', color='blue', alpha=1, label='Data points')
plt.plot(med_distances, med_ssha, color='red', label='Median')  # Add the median line

plt.xlim([-42.6, -40.75])  # Adjust the x-axis limits
plt.xlabel('Longitude (Degrees_east)')
plt.ylabel('SSH (m)')
plt.legend()
plt.grid()

plt.show()




#%% Movement of two icebergs

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

# Load elevation data
fn_elevation = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Starting coordinates
start_lon = -42.7
start_lat = -52.85

# Starting coordinates
start_lon2 = -41.35
start_lat2 = -53.355

# Create subplot
fig, ax = plt.subplots(figsize=(10, 6))

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='viridis', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Set labels and title
plt.xlabel('Longitude (degrees)', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=14)

# Set plot limits
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
#plt.ylim([-57, -48])
#plt.xlim([-52, -35])

# Add grid and colorbar
plt.grid(True)

# Plot iceberg path
iceberg_lon = np.array([-42.7, -42.52, -42.38, -42.17, -41.78])
iceberg_lat = np.array([-52.85, -52.72, -52.56, -52.32, -52.06])


iceberg_lon2 = np.array([-41.35, -41.62,  0 ,0, 0, -42.2, -42.27, -42.3, -42.27, -42.22, -42.02, 0,  -41.3, -40.8, 0])
iceberg_lat2 = np.array([-53.355, -53.27, 0 ,0, 0, -53.05, -53.0, -52.92, -52.84, -52.7, -52.52, 0, -52.17, -52.03, 0])

indexes = np.array(range(len(iceberg_lat2)))

plt.plot(iceberg_lon, iceberg_lat, marker='o', markersize=4, linestyle='-', color='red',  linewidth=4, label='Observed iceberg 1 by SWOT')

# Filter out 0 points from iceberg_lon and iceberg_lat
iceberg_lon_filtered = iceberg_lon2[iceberg_lon2 != 0]
iceberg_lat_filtered = iceberg_lat2[iceberg_lon2 != 0]
indexes_filtered = indexes[iceberg_lat2 != 0]

# Plot iceberg path with filtered points
plt.plot(iceberg_lon_filtered, iceberg_lat_filtered, marker='o', markersize=4, linestyle='-', linewidth=4, color='yellow', label="Observed iceberg 2 by SWOT")


#plt.plot(iceberg_lon2, iceberg_lat2, marker='o', markersize=3, linestyle='-', color='lime', label='Iceberg')


# Add legend
plt.legend()

# Add text annotations
for idx, (lon, lat) in enumerate(zip(iceberg_lon, iceberg_lat), 1):
    plt.text(lon, lat, str(idx-1), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))
    
    

# Plot iceberg path with filtered points and numbering starting from 0
for idx, (index, lon, lat) in enumerate(zip(indexes_filtered, iceberg_lon_filtered, iceberg_lat_filtered)):
    plt.plot(lon, lat, marker='o', markersize=4, color='black')
    plt.text(lon, lat, str(index), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))


elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar to the side
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)  # Set label font size
colorbar.ax.tick_params(labelsize=12)

# Reverse the direction of the colorbar
colorbar.invert_yaxis()


plt.show()




#%% The movement of Iceberg 1
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc

# Load elevation data
fn_elevation =  'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Define function to find the nearest non-NaN point
def find_nearest_non_nan(lon, lat, u, v, start_lon, start_lat):
    # Compute distances from the starting point to all other points
    distances = np.sqrt((lon - start_lon)**2 + (lat - start_lat)**2)
    # Find the index of the minimum distance
    min_index = np.unravel_index(np.nanargmin(distances), distances.shape)
    # Check if the minimum distance point is NaN
    while np.isnan(u[min_index]) or np.isnan(v[min_index]): 
        distances[min_index] = np.inf  # Set the distance to infinity so it's not selected again
        min_index = np.unravel_index(np.nanargmin(distances), distances.shape)
    return min_index


directory_path = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/F'

# Create a list of file paths
#file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]

file_paths  = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]
file_paths = file_paths[:18]


# Load wind data
directory_path_wind = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data2/Filer4'
file_paths_wind = [os.path.join(directory_path_wind, f"{i}.npy") for i in range(18)]  # Data for vector arrows


# Load current data from directory_path_current2
directory_path_current2 = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data/filer'
file_paths_current2 = [os.path.join(directory_path_current2, f"{i}.npy") for i in range(18)]


# Starting coordinates


#Isbjerg 1
start_lon = -42.7
start_lat = -52.85


start_lon2 = -42.7
start_lat2 = -52.85


# Assume dt is the time step (in hours)
dt = 1

# Lists to store coordinates for the path
path_lon = [start_lon]
path_lat = [start_lat]
path_lon_current2 = [start_lon2]
path_lat_current2 = [start_lat2]

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop

path_current_stopped = False
path_current2_stopped = False

for i,  (file_path_wind, file_path, file_path_current2) in enumerate(zip(file_paths_wind, file_paths, file_paths_current2)):
    #data = np.load(file_path)
    data = np.load(file_path, allow_pickle=True)
    
    # Print some information about the loaded data
    lon = data[0] 
    lat = data[1]
    u = data[2] * 3.6 
    v = data[3] * 3.6 
    
    data_wind = np.load(file_path_wind)
    
    # Print some information about the loaded wind data
    lon_wind = data_wind[0, :, :]
    lat_wind = data_wind[1, :, :]
    u_wind = data_wind[2, :, :] * 3.6  # Convert to km/h
    v_wind = data_wind[3, :, :] * 3.6  # Convert to km/h
        
    # Load current data from directory_path_current2
    data_current2 = np.load(file_path_current2)
    
    # Print some information about the loaded current data from directory_path_current2
    lon_current2 = data_current2[0, :, :]
    lat_current2 = data_current2[1, :, :]
    u_current2 = data_current2[2, :, :] * 3.6  # Convert to km/h
    v_current2 = data_current2[3, :, :] * 3.6  # Convert to km/h
    
    for t in range(25):
        if not path_current_stopped:
            # Find the index of the nearest non-NaN point
            nearest_index = find_nearest_non_nan(lon, lat, u, v, start_lon, start_lat)
            # Interpolate wind data at the nearest point
            u_ny = u[nearest_index]
            v_ny = v[nearest_index]
            
            u_wind_ny = griddata((lat_wind.flatten(), lon_wind.flatten()), u_wind.flatten(), (start_lat, start_lon), method='nearest')
            v_wind_ny = griddata((lat_wind.flatten(), lon_wind.flatten()), v_wind.flatten(), (start_lat, start_lon), method='nearest')
            
            u_combined = u_wind_ny * 0.00 + u_ny * 1
            v_combined = v_wind_ny * 0.00 + v_ny * 1
            
            # Update the coordinates based on the displacement
            delta_lon_km = u_combined * dt
            delta_lat_km = v_combined * dt
            
            delta_lon_deg = delta_lon_km / (111.32 * np.cos(np.radians(start_lat)))
            km_per_degree_lat = 111.32 * np.cos(np.radians(start_lat))
            delta_lat_deg = delta_lat_km / km_per_degree_lat
            
            start_lon = start_lon + delta_lon_deg
            start_lat = start_lat + delta_lat_deg
            
            # Check if out of bounds
            if start_lon < -43 or start_lon > -40 or start_lat < -53.5 or start_lat > -51.5:
                print(f"Stopping SWOT path as point is out of bounds: {start_lon}, {start_lat}")
                path_current_stopped = True
            
            path_lon.append(start_lon)
            path_lat.append(start_lat)
        
        if not path_current2_stopped:
            u_current2_ny = griddata((lat_current2.flatten(), lon_current2.flatten()), u_current2.flatten(), (start_lat2, start_lon2), method='nearest')
            v_current2_ny = griddata((lat_current2.flatten(), lon_current2.flatten()), v_current2.flatten(), (start_lat2, start_lon2), method='nearest')
            
            u_wind_ny2 = griddata((lat_wind.flatten(), lon_wind.flatten()), u_wind.flatten(), (start_lat2, start_lon2), method='nearest')
            v_wind_ny2 = griddata((lat_wind.flatten(), lon_wind.flatten()), v_wind.flatten(), (start_lat2, start_lon2), method='nearest')
            
            u_combined2 = u_wind_ny2 * 0.05 + u_current2_ny * 0.95
            v_combined2 = v_wind_ny2 * 0.05 + v_current2_ny * 0.95
           
            
            delta_lon_km2 = u_combined2 * dt
            delta_lat_km2 = v_combined2 * dt
            
            delta_lon_deg2 = delta_lon_km2 / (111.32 * np.cos(np.radians(start_lat2)))
            km_per_degree_lat2 = 111.32 * np.cos(np.radians(start_lat2))
            delta_lat_deg2 = delta_lat_km2 / km_per_degree_lat2
            
            start_lon2 += delta_lon_deg2
            start_lat2 += delta_lat_deg2
            
            # Check if out of bounds
            if start_lon2 < -43 or start_lon2 > -40 or start_lat2 < -53.5 or start_lat2 > -51.5:
                print(f"Stopping OSCAR path as point is out of bounds: {start_lon2}, {start_lat2}")
                path_current2_stopped = True
            
            path_lon_current2.append(start_lon2)
            path_lat_current2.append(start_lat2)    
                

# Plot the current path
plt.plot(path_lon, path_lat, marker='o', markersize=4, linestyle='-', color='lime', label='Modelled path based on SWOT')

# Plot the current path from directory_path_current2
plt.plot(path_lon_current2, path_lat_current2, marker='o', markersize=4, linestyle='-', color='red', label='Modelled path based on OSCAR')

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='viridis', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Set labels and title
plt.xlabel('Longitude (degrees)', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=14)

# Plot starting point
plt.plot(path_lon[0], path_lat[0], marker='o', markersize=2, color='red')

# Iceberg path coordinates
iceberg_lon = np.array([-42.7, -42.52, -42.38, -42.17, -41.78, 0, 0, 0])
iceberg_lat = np.array([-52.85, -52.72, -52.56, -52.32, -52.06, 0, 0, 0])
indexes = np.array(range(len(iceberg_lat)))

# Filter out 0 points from iceberg_lon and iceberg_lat
iceberg_lon_filtered = iceberg_lon[iceberg_lon != 0]
iceberg_lat_filtered = iceberg_lat[iceberg_lon != 0]
indexes_filtered = indexes[iceberg_lat != 0]

# Plot iceberg path with filtered points
plt.plot(iceberg_lon_filtered, iceberg_lat_filtered, marker='o', markersize=4, linestyle='-', color='yellow', linewidth=4, label='Observed iceberg 1 by SWOT')
#plt.plot(iceberg_lon_filtered, iceberg_lat_filtered, marker='o', markersize=4, linestyle='-', color='yellow', label='Observed Iceberg 2')

plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
# Add legend
plt.legend()

# Plot time steps starting from 0
for t in range(0, len(path_lon), 25):
    plt.plot(path_lon[t], path_lat[t], marker='o', markersize=4, color='black')
    plt.text(path_lon[t], path_lat[t], str(t // 24), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

for t in range(0, len(path_lon_current2), 25):
    plt.plot(path_lon_current2[t], path_lat_current2[t], marker='o', markersize=4, color='black')
    plt.text(path_lon_current2[t], path_lat_current2[t], str(t // 24), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

# Plot iceberg path with filtered points and numbering starting from 0
for idx, (index, lon, lat) in enumerate(zip(indexes_filtered, iceberg_lon_filtered, iceberg_lat_filtered)):
    plt.plot(lon, lat, marker='o', markersize=4, color='black')
    plt.text(lon, lat, str(index), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

# Re-plot the elevation data
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar to the side
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)
colorbar.ax.tick_params(labelsize=12)

# Reverse the direction of the colorbar
colorbar.invert_yaxis()
plt.grid(True)
plt.savefig('plot.eps', format='eps')
plt.show()





#%% The movement of Iceberg 2
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc

# Load elevation data
fn_elevation =  'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Define function to find the nearest non-NaN point
def find_nearest_non_nan(lon, lat, u, v, start_lon, start_lat):
    # Compute distances from the starting point to all other points
    distances = np.sqrt((lon - start_lon)**2 + (lat - start_lat)**2)
    # Find the index of the minimum distance
    min_index = np.unravel_index(np.nanargmin(distances), distances.shape)
    # Check if the minimum distance point is NaN
    while np.isnan(u[min_index]) or np.isnan(v[min_index]): 
        distances[min_index] = np.inf  # Set the distance to infinity so it's not selected again
        min_index = np.unravel_index(np.nanargmin(distances), distances.shape)
    return min_index


directory_path = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/F'

# Create a list of file paths
#file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]

file_paths  = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]
file_paths = file_paths[:11]


# Load wind data
directory_path_wind = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data2/Filer4'
file_paths_wind = [os.path.join(directory_path_wind, f"{i}.npy") for i in range(18)]  # Data for vector arrows


# Load current data from directory_path_current2
directory_path_current2 = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data/filer'
file_paths_current2 = [os.path.join(directory_path_current2, f"{i}.npy") for i in range(18)]


# Starting coordinates
#Isbjerg 1
#start_lon = -42.7
#start_lat = -52.85

#Isbjerg 2
start_lon = -41.35
start_lat = -53.355

start_lon2 = -41.35
start_lat2 = -53.355


# Assume dt is the time step (in hours)
dt = 1

# Lists to store coordinates for the path
path_lon = [start_lon]
path_lat = [start_lat]
path_lon_current2 = [start_lon2]
path_lat_current2 = [start_lat2]

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop

path_current_stopped = False
path_current2_stopped = False

for i,  (file_path_wind, file_path, file_path_current2) in enumerate(zip(file_paths_wind, file_paths, file_paths_current2)):
    #data = np.load(file_path)
    data = np.load(file_path, allow_pickle=True)
    
    # Print some information about the loaded data
    lon = data[0] 
    lat = data[1]
    u = data[2] * 3.6 
    v = data[3] * 3.6 
    
    data_wind = np.load(file_path_wind)
    
    # Print some information about the loaded wind data
    lon_wind = data_wind[0, :, :]
    lat_wind = data_wind[1, :, :]
    u_wind = data_wind[2, :, :] * 3.6  # Convert to km/h
    v_wind = data_wind[3, :, :] * 3.6  # Convert to km/h
        
    # Load current data from directory_path_current2
    data_current2 = np.load(file_path_current2)
    
    # Print some information about the loaded current data from directory_path_current2
    lon_current2 = data_current2[0, :, :]
    lat_current2 = data_current2[1, :, :]
    u_current2 = data_current2[2, :, :] * 3.6  # Convert to km/h
    v_current2 = data_current2[3, :, :] * 3.6  # Convert to km/h
    
    for t in range(25):
        if not path_current_stopped:
            # Find the index of the nearest non-NaN point
            nearest_index = find_nearest_non_nan(lon, lat, u, v, start_lon, start_lat)
            # Interpolate wind data at the nearest point
            u_ny = u[nearest_index]
            v_ny = v[nearest_index]
            
            u_wind_ny = griddata((lat_wind.flatten(), lon_wind.flatten()), u_wind.flatten(), (start_lat, start_lon), method='nearest')
            v_wind_ny = griddata((lat_wind.flatten(), lon_wind.flatten()), v_wind.flatten(), (start_lat, start_lon), method='nearest')
            
            u_combined = u_wind_ny * 0.01 + u_ny * 0.99
            v_combined = v_wind_ny * 0.01 + v_ny * 0.99
            
            # Update the coordinates based on the displacement
            delta_lon_km = u_combined * dt
            delta_lat_km = v_combined * dt
            
            delta_lon_deg = delta_lon_km / (111.32 * np.cos(np.radians(start_lat)))
            km_per_degree_lat = 111.32 * np.cos(np.radians(start_lat))
            delta_lat_deg = delta_lat_km / km_per_degree_lat
            
            start_lon = start_lon + delta_lon_deg
            start_lat = start_lat + delta_lat_deg
            
            # Check if out of bounds
            if start_lon < -43 or start_lon > -40 or start_lat < -53.5 or start_lat > -51.5:
                print(f"Stopping SWOT path as point is out of bounds: {start_lon}, {start_lat}")
                path_current_stopped = True
            
            path_lon.append(start_lon)
            path_lat.append(start_lat)
        
        if not path_current2_stopped:
            u_current2_ny = griddata((lat_current2.flatten(), lon_current2.flatten()), u_current2.flatten(), (start_lat2, start_lon2), method='nearest')
            v_current2_ny = griddata((lat_current2.flatten(), lon_current2.flatten()), v_current2.flatten(), (start_lat2, start_lon2), method='nearest')
            
            u_wind_ny2 = griddata((lat_wind.flatten(), lon_wind.flatten()), u_wind.flatten(), (start_lat2, start_lon2), method='nearest')
            v_wind_ny2 = griddata((lat_wind.flatten(), lon_wind.flatten()), v_wind.flatten(), (start_lat2, start_lon2), method='nearest')
            
            u_combined2 = u_wind_ny2 * 0.01 + u_current2_ny * 0.99
            v_combined2 = v_wind_ny2 * 0.01 + v_current2_ny * 0.99
            
            delta_lon_km2 = u_combined2 * dt
            delta_lat_km2 = v_combined2 * dt
            
            delta_lon_deg2 = delta_lon_km2 / (111.32 * np.cos(np.radians(start_lat2)))
            km_per_degree_lat2 = 111.32 * np.cos(np.radians(start_lat2))
            delta_lat_deg2 = delta_lat_km2 / km_per_degree_lat2
            
            start_lon2 += delta_lon_deg2
            start_lat2 += delta_lat_deg2
            
            # Check if out of bounds
            if start_lon2 < -43 or start_lon2 > -40 or start_lat2 < -53.5 or start_lat2 > -51.5:
                print(f"Stopping OSCAR path as point is out of bounds: {start_lon2}, {start_lat2}")
                path_current2_stopped = True
            
            path_lon_current2.append(start_lon2)
            path_lat_current2.append(start_lat2)    
                

# Plot the current path
plt.plot(path_lon, path_lat, marker='o', markersize=4, linestyle='-', color='lime', label='Modelled path based on SWOT')

# Plot the current path from directory_path_current2
plt.plot(path_lon_current2, path_lat_current2, marker='o', markersize=4, linestyle='-', color='red', label='Modelled path based on OSCAR')

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='viridis', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Set labels and title
plt.xlabel('Longitude (degrees)', fontsize=14)
plt.ylabel('Latitude (degrees)', fontsize=14)

# Plot starting point
plt.plot(path_lon[0], path_lat[0], marker='o', markersize=2, color='red')

# Iceberg path coordinates
iceberg_lon = np.array([-41.35, -41.62, 0, 0, 0, -42.2, -42.27, -42.3, -42.27, -42.22, -42.02, 0, -41.3, -40.8, 0])
iceberg_lat = np.array([-53.355, -53.27, 0, 0, 0, -53.05, -53.0, -52.92, -52.84, -52.7, -52.52, 0, -52.17, -52.03, 0])
indexes = np.array(range(len(iceberg_lat)))

# Filter out 0 points from iceberg_lon and iceberg_lat
iceberg_lon_filtered = iceberg_lon[iceberg_lon != 0]
iceberg_lat_filtered = iceberg_lat[iceberg_lon != 0]
indexes_filtered = indexes[iceberg_lat != 0]

# Plot iceberg path with filtered points
plt.plot(iceberg_lon_filtered, iceberg_lat_filtered, marker='o', markersize=4, linestyle='-', color='yellow', linewidth=4, label='Observed iceberg 2 by SWOT')
#plt.plot(iceberg_lon_filtered, iceberg_lat_filtered, marker='o', markersize=4, linestyle='-', color='yellow', label='Observed Iceberg 2')

plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
# Add legend
plt.legend()

# Plot time steps starting from 0
for t in range(0, len(path_lon), 25):
    plt.plot(path_lon[t], path_lat[t], marker='o', markersize=4, color='black')
    plt.text(path_lon[t], path_lat[t], str(t // 24), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

for t in range(0, len(path_lon_current2), 25):
    plt.plot(path_lon_current2[t], path_lat_current2[t], marker='o', markersize=4, color='black')
    plt.text(path_lon_current2[t], path_lat_current2[t], str(t // 24), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

# Plot iceberg path with filtered points and numbering starting from 0
for idx, (index, lon, lat) in enumerate(zip(indexes_filtered, iceberg_lon_filtered, iceberg_lat_filtered)):
    plt.plot(lon, lat, marker='o', markersize=4, color='black')
    plt.text(lon, lat, str(index), fontsize=12, color='black', ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.1'))

# Re-plot the elevation data
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar to the side
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)
colorbar.ax.tick_params(labelsize=12)

# Reverse the direction of the colorbar
colorbar.invert_yaxis()
plt.grid(True)
plt.savefig('plot.eps', format='eps')
plt.show()


#%% Estimate area iceberg 1, sigma_0 plot

h = swot.distance_between_points(-42.56, -42.56, -52.78, -52.65)/1000 
L1 = swot.distance_between_points(-42.56, -42.48, -52.78, -52.78)/1000  


print(f"{h} km")
print(f"{L1} km")


A1 = (h*L1)/2
print(f" Areal {A1} km^2")

import swot_ssh_utils as swot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.patches import Polygon

fn = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/GIF/SWOT_L2_LR_SSH_Unsmoothed_540_005_20230602T150806_20230602T155911_PIB0_01.nc'

def plot_unsmoothed_sea_ice_sig0(fn, cc, 
                                 s3sys=None,
                                 test_plot=False,
                                 for_google_earth=False):
    dd = swot.SSH_L2()
    dd.load_data(fn, s3sys=s3sys, lat_bounds=[])
    pn=fn.split('/')[-1].split('_')[6]
    cn=fn.split('/')[-1].split('_')[5]
    fnn = fn.split('/')[-1]

    lat = dd.left.latitude.data[:, 120]
    lon = dd.left.longitude.data[:, 120]
    if len(cc['lat_bounds']) == 2:
        lat0, lat1 = cc['lat_bounds']
        msk = (lat > lat0) & (lat < lat1)
    else:
        lon0, lon1 = cc['lon_bounds']
        msk = (lon > lon0) & (lon < lon1)
    # Create a figure and axis with Orthographic projection
    fig, ax = plt.subplots(1,1, figsize=(6,5),
                    subplot_kw={'projection': 
                            ccrs.LambertAzimuthalEqualArea(central_latitude=np.nanmean(dd.left.latitude.data),
                            central_longitude=np.nanmean(dd.left.longitude.data))})
    
    # Loop through the left and right data
    lons, lats = [], []
    if test_plot:
        skip=5
    else:
        skip=1
    for i, tmp in enumerate([dd.left, dd.right]):
               
        lon, lat = tmp.longitude.data[msk, :], tmp.latitude.data[msk, :]
        lons.append(lon)
        lats.append(lat)
        dtm = tmp.sig0_karin_2[msk, :].data
        dtm -= np.nanmin(dtm+1e-10)
        mm=np.isfinite(dtm)
        dtm[mm] = np.log2(dtm[mm])
        cax2 = ax.pcolor(lon[::skip,::skip], lat[::skip,::skip], 
                        dtm[::skip,::skip],
                        cmap='bone',
                        transform=ccrs.PlateCarree())    
    
    min_lon=np.nanmin(np.array(lons))
    max_lon=np.nanmax(np.array(lons))
    min_lat=np.nanmin(np.array(lats))
    max_lat=np.nanmax(np.array(lats))
    ax.set_extent([min_lon, max_lon, min_lat, max_lat ] ,crs=ccrs.PlateCarree())
    if for_google_earth: #not working yet
        # Remove the border frame
        ax.set_frame_on(False)
        # Turn off the axis
        ax.axis('off')
        ax.set_facecolor("none")
        # Ensure that the aspect ratio is equal
        ax.set_aspect('equal', 'box')
        plt.savefig('../media/figures/Unsmoothed_sig0_images/'+fn.split('/')[-1].split('.')[0]+'_noborder.png', 
                    dpi=300,bbox_inches='tight', pad_inches=0)
    else:
        gl0 = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        ax.text(0.02, 0.02, fnn,
                horizontalalignment='left',
                fontsize=6,
                transform=ax.transAxes,
                color='white',)
        cbar = plt.colorbar(cax2, ax=ax,
                        shrink=0.9, orientation='horizontal',
                        pad=0.1, aspect=50,
                        label='log2(sig0)')
        
        plt.axis('equal')
        fig = plt.gcf()
        #fig.set_size_inches(10, 10)  # Set the figure size in inches (width, height)
        fig.set_dpi(300)
    
    # Adding the polygon here
    polygon_points = [[-42.56, -42.48, -42.56], [-52.78, -52.78, -52.65]]
    #polygon_points = [[-42.56, -42.505, -42.505, -42.56], [-52.78, -52.78, -52.65, -52.65]]
    polygon = Polygon(np.array(polygon_points).T, closed=True, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree())
    ax.add_patch(polygon)

    # Return the ax object for further modification outside the function if needed
    return ax

cc={'pn':'021','lat_bounds':[-53.5,-51.5],'lon_bounds':[-43, -40],'day':'0313','region_name':'Antarctica'}

#[-53.35,-53.15], [-42, -41.6]

ax = plot_unsmoothed_sea_ice_sig0(fn, cc,
                                  for_google_earth=False,
                                  test_plot=False,
                                  s3sys=None)

plt.show()



#%% Estimate area iceberg 2, sigma_0 plot

h = swot.distance_between_points(-41.65, -41.66, -53.23, -53.28)/1000 
L1 = swot.distance_between_points(-41.65, -41.3, -53.23, -53.29)/1000 
L2 = swot.distance_between_points(-41.85, -41.5, -53.25, -53.31)/1000 

print(f"{h} km")
print(f"{L1} km")
print(f"{L2} km")

A2 = h*L2
print(f" Areal {A2} km^2")


import swot_ssh_utils as swot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.patches import Polygon

fn = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/GIF/SWOT_L2_LR_SSH_Unsmoothed_540_005_20230602T150806_20230602T155911_PIB0_01.nc'

def plot_unsmoothed_sea_ice_sig0(fn, cc, 
                                 s3sys=None,
                                 test_plot=False,
                                 for_google_earth=False):
    dd = swot.SSH_L2()
    dd.load_data(fn, s3sys=s3sys, lat_bounds=[])
    pn=fn.split('/')[-1].split('_')[6]
    cn=fn.split('/')[-1].split('_')[5]
    fnn = fn.split('/')[-1]

    lat = dd.left.latitude.data[:, 120]
    lon = dd.left.longitude.data[:, 120]
    if len(cc['lat_bounds']) == 2:
        lat0, lat1 = cc['lat_bounds']
        msk = (lat > lat0) & (lat < lat1)
    else:
        lon0, lon1 = cc['lon_bounds']
        msk = (lon > lon0) & (lon < lon1)
    # Create a figure and axis with Orthographic projection
    fig, ax = plt.subplots(1,1, figsize=(6,5),
                    subplot_kw={'projection': 
                            ccrs.LambertAzimuthalEqualArea(central_latitude=np.nanmean(dd.left.latitude.data),
                            central_longitude=np.nanmean(dd.left.longitude.data))})
    
    # Loop through the left and right data
    lons, lats = [], []
    if test_plot:
        skip=5
    else:
        skip=1
    for i, tmp in enumerate([dd.left, dd.right]):
               
        lon, lat = tmp.longitude.data[msk, :], tmp.latitude.data[msk, :]
        lons.append(lon)
        lats.append(lat)
        dtm = tmp.sig0_karin_2[msk, :].data
        dtm -= np.nanmin(dtm+1e-10)
        mm=np.isfinite(dtm)
        dtm[mm] = np.log2(dtm[mm])
        cax2 = ax.pcolor(lon[::skip,::skip], lat[::skip,::skip], 
                        dtm[::skip,::skip],
                        cmap='bone',
                        transform=ccrs.PlateCarree())    
    
    min_lon=np.nanmin(np.array(lons))
    max_lon=np.nanmax(np.array(lons))
    min_lat=np.nanmin(np.array(lats))
    max_lat=np.nanmax(np.array(lats))
    ax.set_extent([min_lon, max_lon, min_lat, max_lat ] ,crs=ccrs.PlateCarree())
    if for_google_earth: #not working yet
        # Remove the border frame
        ax.set_frame_on(False)
        # Turn off the axis
        ax.axis('off')
        ax.set_facecolor("none")
        # Ensure that the aspect ratio is equal
        ax.set_aspect('equal', 'box')
        plt.savefig('../media/figures/Unsmoothed_sig0_images/'+fn.split('/')[-1].split('.')[0]+'_noborder.png', 
                    dpi=300,bbox_inches='tight', pad_inches=0)
    else:
        gl0 = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        ax.text(0.02, 0.02, fnn,
                horizontalalignment='left',
                fontsize=6,
                transform=ax.transAxes,
                color='white',)
        cbar = plt.colorbar(cax2, ax=ax,
                        shrink=0.9, orientation='horizontal',
                        pad=0.1, aspect=50,
                        label='log2(sig0)')
        plt.axis('equal')
        fig = plt.gcf()
        #fig.set_size_inches(10, 10)  # Set the figure size in inches (width, height)
        fig.set_dpi(300)
    
    # Adding the polygon here
    polygon_points = [[-41.65, -41.85, -41.5, -41.3], [-53.23, -53.25, -53.31, -53.29]]
    #ax.plot([-41.66], [-53.28], marker='o', linestyle='-', color='black', transform=ccrs.PlateCarree())
    #polygon_points = [[-41.65, -41.81, -41.5, -41.3], [-53.19, -53.25, -53.31, -53.25]]
    #polygon_points = [[-41.81, -41.81, -41.4, -41.4], [-53.25, -53.31, -53.31, -53.25]]
    polygon = Polygon(np.array(polygon_points).T, closed=True, edgecolor='black', facecolor='none', transform=ccrs.PlateCarree())
    ax.add_patch(polygon)

    # Return the ax object for further modification outside the function if needed
    return ax

cc={'pn':'021','lat_bounds':[-53.5,-51.5],'lon_bounds':[-43, -40],'day':'0313','region_name':'Antarctica'}

ax = plot_unsmoothed_sea_ice_sig0(fn, cc,
                                  for_google_earth=False,
                                  test_plot=False,
                                  s3sys=None)
plt.show()




#%% Distance to coast iceberg 1
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc

import swot_ssh_utils as swot
L1 =swot.distance_between_points(-42.7, -42.6, -52.85, -53.25)/1000 



print(L1)

# Load elevation data
fn_elevation = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Create figure and axis with Cartopy
fig, ax = plt.subplots(figsize=(12, 6))

# Set plot limits
ax.set_ylim([-53.5, -51.5])
ax.set_xlim([-43, -40])

# Set labels
ax.set_xlabel('Longitude (degrees)', fontsize=14)
ax.set_ylabel('Latitude (degrees)', fontsize=14)

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)
colorbar.ax.tick_params(labelsize=12)



point2 = (-42.7, -52.85)
point1 = (-42.6, -53.25)

ax.plot(point1[0], point1[1], color='black', marker='o',  markersize=8)
ax.plot(point2[0], point2[1], color='red', marker='o',  markersize=8)



# Reverse the direction of the colorbar
colorbar.invert_yaxis()

# Add grid
plt.grid(True)

# Save and show plot
plt.savefig('plot.eps', format='eps')
plt.show()





#%% Distance to coast iceberg 2

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc

import swot_ssh_utils as swot
L1 =swot.distance_between_points(-41.62, -41.85, -53.27, -53.35)/1000 

print(L1)

# Load elevation data
fn_elevation = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Create figure and axis with Cartopy
fig, ax = plt.subplots(figsize=(12, 6))

# Set plot limits
ax.set_ylim([-53.5, -51.5])
ax.set_xlim([-43, -40])

# Set labels
ax.set_xlabel('Longitude (degrees)', fontsize=14)
ax.set_ylabel('Latitude (degrees)', fontsize=14)

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)
colorbar.ax.tick_params(labelsize=12)



point2 = (-41.62, -53.27)
point1 = (-41.85, -53.35)

# Plot the points
ax.plot(point1[0], point1[1], color='black', marker='o',  markersize=8)
ax.plot(point2[0], point2[1], color='red', marker='o',  markersize=8)



# Reverse the direction of the colorbar
colorbar.invert_yaxis()

# Add grid
plt.grid(True)

# Save and show plot
plt.savefig('plot.eps', format='eps')
plt.show()



#%% Wind plots 
import numpy as np
import matplotlib.pyplot as plt
import os

# Load wind data
#directory_path_wind = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data2/Filer4'
#file_paths_wind = [os.path.join(directory_path_wind, f"{i}.npy") for i in range(5)]  # Data for vector arrows

directory_path_wind = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/Data2/Filer4'
#file_paths_wind = [os.path.join(directory_path_wind, f"{i}.npy") for i in range(6, 8)]  # Data for vector arrows
file_paths_wind = [os.path.join(directory_path_wind, f"{i}.npy") for i in range(9, 13)] 

# Load elevation data
fn_elevation = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Vind/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
with nc.Dataset(fn_elevation, 'r') as f:
    elevation_var = f.variables['elevation'][:]
    lon_elevation = f.variables['lon'][:]
    lat_elevation = f.variables['lat'][:]

# Initialize lists to store wind data
lon_wind_list = []
lat_wind_list = []
u_wind_list = []
v_wind_list = []

# Load wind data from the first six files and store in lists
for file_path_wind in file_paths_wind:
    data_wind = np.load(file_path_wind)
    lon_wind_list.append(data_wind[0, :, :])
    lat_wind_list.append(data_wind[1, :, :])
    u_wind_list.append(data_wind[2, :, :] *3.6 )  # Convert to km/h
    v_wind_list.append(data_wind[3, :, :] *3.6)  # Convert to km/h

# Calculate the mean of the wind data
lon_wind_mean = np.mean(lon_wind_list, axis=0)
lat_wind_mean = np.mean(lat_wind_list, axis=0)
u_wind_mean = np.mean(u_wind_list, axis=0)
v_wind_mean = np.mean(v_wind_list, axis=0)


fig, ax = plt.subplots(figsize=(12, 6)) 

# Set plot limits
ax.set_ylim([-53.5, -51.5])
ax.set_xlim([-43, -40])

ax.set_xlabel('Longitude (degrees)', fontsize=14)
ax.set_ylabel('Latitude (degrees)', fontsize=14)

# Plot elevation
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='Blues_r', origin='lower', aspect='auto', vmin=-4500, vmax=1)

# Add colorbar
colorbar = plt.colorbar(elevation_plot, ax=ax, orientation='vertical', pad=0.05)
colorbar.set_label('Elevation (m)', fontsize=14)
colorbar.ax.tick_params(labelsize=12)


plt.title('      ')
quiver = plt.quiver(lon_wind_mean, lat_wind_mean, u_wind_mean, v_wind_mean, scale=200)

plt.quiverkey(quiver, 0.8, 0.91, 10, "10 km/h", labelpos="E", coordinates='figure')
plt.show()




#%% Calculations of geostrophic currents from SWOT:
    
#%%
import swot_ssh_utils as swot
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import netCDF4 as nc
import cartopy
import cartopy.crs as ccrs
import os
from numpy.ma import masked_where
from scipy.ndimage import label, gaussian_filter
import cartopy.feature as cfeature
from PIL import Image
import imageio
import sys
# skud ud: https://github.com/SWOT-community/SWOT-OpenToolkit/blob/main/src/swot_ssh_utils.py
from scipy import ndimage

def interpolate_data(file,value,longitude, latitude):
    from scipy.interpolate import RectBivariateSpline as rbs
    from scipy.interpolate import griddata
    from xarray import open_dataset
    import numpy as np
    import requests
    import xarray as xr 
    
    #LOAD DATA
    # Create a mask for missing data in the input longitude and latitude grids
    m = np.isnan(longitude * latitude).flatten()
    
    # Calculate the min and max values for longitude and latitude
    lon_min = np.nanmin(longitude)
    lon_max = np.nanmax(longitude)
    lat_min = np.nanmin(latitude)
    lat_max = np.nanmax(latitude)
    
    # Load geoid data from normal file
    load_file = nc.Dataset(file, "r")
     
    # Extract data variables
    datakeys = load_file.variables.keys()
     
    # Collect data in a python dictionary
    swot = {}
    for k in datakeys:
         swot[k] = load_file.variables[k][:]
     
    # Close the NetCDF file
    load_file.close()
     
     
     # Assuming your NetCDF file has variables 'latitude', 'longitude', and 'ssha'
    lat = swot["latitude"][:]
    lon = swot["longitude"][:]
    #lon[lon > 180] -= 360
    geoid = swot[value][:]

    lat_indices = (lat >= lat_min) & (lat <= lat_max)

    lon=lon[lat_indices]
    lat=lat[lat_indices]
    geoid=geoid[lat_indices]

    lon_indices = (lon >= lon_min) & (lon <= lon_max)
    
    lonm=lon[lon_indices]
    latm=lat[lon_indices]
    geoid=geoid[lon_indices]
    
    lonm = lonm.flatten()
    latm = latm.flatten()

    # Mask for missing values in the selected MSS22 data
    mm = ~np.isnan(geoid).flatten()
    
    # iNTERPOLATE, NOTE ITS FLATTENS TO 1D arrays
    # Interpolate the MSS22 data onto the input grid using griddata
    aa = griddata((latm[mm], lonm[mm]), geoid.flatten()[mm], (latitude.flatten()[~m], longitude.flatten()[~m]))
   
    
    # Initialize an array for interpolated MSS22 data
    geoid = np.ones((longitude.size))

    # Assign missing values and interpolated data to the MSS22 array
    geoid[m] = np.nan
    geoid[~m] = aa

    return geoid.reshape(longitude.shape)    


#%%
import numpy as np
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

fn = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/GIF/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'
basic = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'


#fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

#basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

#fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

#basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

#fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

# initialize the class SSH_L2
dd=swot.SSH_L2()

# print the variables in the file
# load data into the class
# the lat_bounds keyword is used to subset the data
# the return is dd.left and dd.right, which correspond to the left and right swath


dd.load_data(fn,lat_bounds=[-53.5,-51.5],)


# Create a map using cartopy
fig, ax = plt.subplots(figsize=(12, 10),subplot_kw=dict(projection=ccrs.PlateCarree()))

# add land, rivers, lakes, and coastlines
ax.add_feature(cartopy.feature.LAND)

# add states and provinces
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')



# Initialize the class SSH_L2
dd = swot.SSH_L2()

# Load data into the class
dd.load_data(fn, lat_bounds=[-53.5, -51.5])

# Create a map using cartopy
fig, ax = plt.subplots(figsize=(12, 10), subplot_kw=dict(projection=ccrs.PlateCarree()))

# Add land, rivers, lakes, and coastlines
ax.add_feature(cartopy.feature.LAND)

# Add states and provinces
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

side = ['left', 'right']
n = 0

# Visualize the data separately for left and right swaths
for tmp, side_name in zip([dd.left, dd.right], side):
    if tmp is None:
        print(f"No data available for {side_name} swath")
        continue

    surf=tmp.ancillary_surface_classification_flag==0

    # There is a flag ssh_karin_2_qual that can be used to filter out bad data
    # But it is not always reliable yet in the current products (Dec 1, 2023)
    # This flag is not yet used in the code below, but you can uncomment the following line to check it. 
    #qc_flag=tmp.ssh_karin_2_qual==0 # flag==0 means good data
    
    #mask out non-ocean and bad data
    dtm=np.where(surf,tmp.ssh_karin_2.data,np.nan) 
    
    #dtm=dtm-np.nanmin(dtm.flatten()) # this can remove the cross-swath bias 
    
    lon,lat=tmp.longitude.data,tmp.latitude.data
    
    m=~np.isnan(lon+lat).flatten() # mask out the nan values
    
 
    
    # REMOVE ICEBERGS 
    # Smooth the data to reduce noise
    dtm_smoothed = gaussian_filter(dtm, sigma=4)
    
    # Compute the gradient
    gradient_x, gradient_y = np.gradient(dtm_smoothed)
    
    # Calculate the magnitude of the gradient
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    
    # Threshold the magnitude to identify significant changes
    threshold = 0.08 # Adjust this threshold based on your data
    significant_changes_mask = gradient_magnitude > threshold
    
    # Dilate the detected regions
    from scipy.ndimage import binary_dilation
    significant_changes_mask = binary_dilation(significant_changes_mask)
    
    # Now, significant_changes_mask contains True where significant changes are detected
    # You can use this mask to locate the cluster in your data
    cluster_indices = np.where(significant_changes_mask)
    dtm[cluster_indices]=np.nan
    
    # the following line find the latitude limit
    # to bound the area for along-track mean removal
    # you can change the bounds to get a better fit for your region. 
    # The 10th column of the lat array is used as an approximate. 
    msk=(lat[:,10]>-53)&(lat[:,10]<-50)
    
    cossover=interpolate_data(basic,'height_cor_xover',lon, lat) # cross-over bias
    dtm+=cossover # remove the mean sea surface 

    #mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    #dtm-=mss # remove the mean sea surface 
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:]  # trend ssha
   

    
    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm+=mss # remove the mean sea surface 
    
    #dtm+=mss # mss add
    
    mss=interpolate_data(basic,'geoid',lon, lat) # get the geoid at the same points
    dtm-=mss # remove  geoid

    combined_array = np.stack((tmp.longitude.data, tmp.latitude.data, dtm), axis=0)

     # Save the combined array
    name = f"Slut/geoid_19-06-2023-{side_name}.npy"
    np.save(name, combined_array)

    print(f"Saved data for {side_name} swath to {name}")

    lon, lat = tmp.longitude.data, tmp.latitude.data
    m = ~np.isnan(lon + lat).flatten()  # Mask out the nan values


 # plot the data using scatter , vmin=-0.1, vmax=0.1
    cax=ax.scatter(lon.flatten()[m], lat.flatten()[m], vmin=-5, vmax=5,
                   c=dtm.flatten()[m], s=0.05, 
                   cmap='Spectral_r' ,transform=ccrs.PlateCarree()) 
    
    
# add colorbar and gridlines
gl = ax.gridlines(draw_labels=True, color='gray', 
                  alpha=0.5, linestyle='--')


cbar = plt.colorbar(cax, ax=ax,shrink=0.7)
cbar.ax.tick_params(labelsize=12, )
cbar.set_label('Ocean Surface Topography, OST [m]', rotation=270, fontsize=16,  labelpad=-55)

# remove the right side ylabels
gl.right_labels = False
gl.top_labels = False  # hide the top axis labels
#cbar.set_label('ssh_karin_2')
gl.xlabel_style = {'size': '18'}
gl.ylabel_style = {'size': '18'}



date_part = fn.split('/')[-1].split('_')[7].split('T')[0]
#ax.set_title(date_part + ' - Track 05')
ax.text(0.88, 0.02, date_part,
        horizontalalignment='left',
        fontsize=12,
        transform=ax.transAxes,
        color='grey',)



#ax.plot([-40.8, -42.54], [-52.8, -52.3], marker='o', linestyle='-', color='red', transform=ccrs.PlateCarree())
# Adjust figure size and DPI
fig = plt.gcf()
fig.set_dpi(800)  # Set the DPI (dots per inch)

# set the extent of the map
lonmin,lonmax,latmin,latmax=np.nanmin(lon-1),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat)
ax.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())

plt.show()

#%%  
import numpy as np

left = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Slut/geoid_19-06-2023-left.npy'
data_left = np.load(left)

right = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Slut/geoid_19-06-2023-right.npy'
data_right = np.load(right)

right_ssha = data_right[2, :, :]

left_ssha = data_left[2, :, :]


# Determine the maximum number of rows between the two arrays
max_rows = max(left_ssha.shape[0], right_ssha.shape[0])

# Determine the number of columns (should be the same for both arrays)
num_columns = left_ssha.shape[1]

# Pad the smaller array with zeros to match the dimensions of the larger array
left_ssha_padded = np.pad(left_ssha, ((0, max_rows - left_ssha.shape[0]), (0, 0)), mode='constant')
right_ssha_padded = np.pad(right_ssha, ((0, max_rows - right_ssha.shape[0]), (0, 0)), mode='constant')

#Flipper
left_ssha_padded = np.flipud(left_ssha_padded)
right_ssha_padded = np.flipud(right_ssha_padded )

# Concatenate the padded arrays along axis 1
combined_ssha = np.hstack((left_ssha_padded, right_ssha_padded))


right_lon = data_right[0, :, :]
left_lon = data_left[0, :, :]



# Determine the maximum number of rows between the two arrays
max_rows = max(left_lon.shape[0], right_lon.shape[0])

# Determine the number of columns (should be the same for both arrays)
num_columns = left_lon.shape[1]

# Pad the smaller array with zeros to match the dimensions of the larger array
left_lon_padded = np.pad(left_lon, ((0, max_rows - left_lon.shape[0]), (0, 0)), mode='constant')
right_lon_padded = np.pad(right_lon, ((0, max_rows - right_lon.shape[0]), (0, 0)), mode='constant')

#Flipper
#left_lon_padded = np.flipud(left_lon_padded ) - 360
left_lon_padded = np.fliplr(np.flipud(left_lon_padded)) - 360

right_lon_padded = np.flipud(right_lon_padded) - 360


# Concatenate the padded arrays along axis 1
combined_lon = np.hstack((left_lon_padded, right_lon_padded))
combined_lon = combined_lon

right_lat = data_right[1, :, :]
left_lat = data_left[1, :, :]



# Determine the maximum number of rows between the two arrays
max_rows = max(left_lat.shape[0], right_lat.shape[0])

# Determine the number of columns (should be the same for both arrays)
num_columns = left_lat.shape[1]

# Pad the smaller array with zeros to match the dimensions of the larger array
left_lat_padded = np.pad(left_lat, ((0, max_rows - left_lat.shape[0]), (0, 0)), mode='constant')
right_lat_padded = np.pad(right_lat, ((0, max_rows - right_lat.shape[0]), (0, 0)), mode='constant')


#left_lat_padded = np.flipud(left_lat_padded)
left_lat_padded = np.fliplr(np.flipud(left_lat_padded))
right_lat_padded = np.flipud(right_lat_padded) 


# Concatenate the padded arrays along axis 1
combined_lat = np.hstack((left_lat_padded, right_lat_padded)) 

#%%
step = 10

from netCDF4 import Dataset

Expert = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/SWOT_L2_LR_SSH_Expert_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
dataset = Dataset(Expert)

heading = dataset.variables['velocity_heading']
heading=np.array(heading)
#velocity_heading = dataset.variables['velocity_heading'][:]

# Keep every nth point in both directions
filtered_lat = combined_lat[::step, ::step]
filtered_lon = combined_lon[::step, ::step]
filtered_ssha = combined_ssha[::step, ::step]


# Print the new shapes
print(f"Original shape: {combined_lat.shape}")
print(f"Downsampled shape: {filtered_lat.shape}")

gradient_dist_km = 250 * step


DOT = filtered_ssha 

N, M = DOT.shape


lol=np.array(heading)
lol[lol>100]=np.nan
heading=lol


oa = np.full(9866, 90-62)
# Convert gradients to N and E coordinates
orbit_angle = np.deg2rad(oa)

orbit_angle = np.repeat(orbit_angle[:,np.newaxis], axis=1, repeats=M)

orbit_angle = np.resize(orbit_angle, DOT.shape)

## Correct direction

heading_mean = np.nanmedian(heading)

grad_along, grad_cross = np.gradient(DOT)
grad_along, grad_cross = -grad_along, -grad_cross

grad_cross /= gradient_dist_km # Convert unit to m / m
grad_along /= gradient_dist_km # Convert unit to m / m


# Angle with respect to true north of the horizontal component of the spacecraft Earth-relative velocity vector.

# Values between 0 and 90 deg indicate that the velocity vector has a northward component,

# and values between 90 and 180 deg indicate that the velocity vector has a southward component.

 

if heading_mean > 90: # Southward

    correction_angle = np.deg2rad(180) - orbit_angle

    # Along track

    grad_along_x = grad_along*np.sin(correction_angle)

    grad_along_y = -grad_along*np.cos(correction_angle)

 

    # Cross track

    grad_cross_x = -grad_cross*np.cos(correction_angle)

    grad_cross_y = -grad_cross*np.sin(correction_angle)
    print('wrong')

 

else: # Northward

    # Along track

    grad_along_x = grad_along*np.sin(orbit_angle)

    grad_along_y = grad_along*np.cos(orbit_angle)

 

    # Cross track

    grad_cross_x = grad_cross*np.cos(orbit_angle)

    grad_cross_y = -grad_cross*np.sin(orbit_angle)

 

grad_x = grad_cross_x + grad_along_x

grad_y = grad_cross_y + grad_along_y

# Surface Geostrophic current

g = 9.82 # m/s

Omega = 7.292e-5 # s^-1

f = 2 * Omega * np.sin(np.deg2rad(filtered_lat ))
f=np.array(f)



geo_u =- g/f * (grad_y) # convert gradient to m/m

geo_v = +g/f * (grad_x) # convert gradient to m/m

geo_mag = np.sqrt(geo_u**2 + geo_v**2)


geo_mag[geo_mag==np.inf]=0
leftu=geo_u[:,2:23]
leftv=geo_v[:,2:23]
llon=filtered_lon[:,2:23]
llat=filtered_lat[:,2:23]

rightu=geo_u[:,26:48]
rightv=geo_v[:,26:48]
rlon=filtered_lon[:,26:48]
rlat=filtered_lat[:,26:48]


# Concatenate the padded arrays along axis 1
lon = np.hstack((llon, rlon))
lat = np.hstack((llat, rlat))
geo_u = np.hstack((leftu, rightu))
geo_v = np.hstack((leftv, rightv))

data = np.array([lon,lat, geo_u, geo_v])

# Angiv den ønskede sti til filen
#file_path = "C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/datafiles2/Gradient/data_18.npy"
file_path = "C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/Slut/data_18.npy"

# Gem dataene i en npy-fil
np.save(file_path, data)

#%% Distance to coast SWOT
from netCDF4 import Dataset


fn2 = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

# Open the NetCDF file
dataset = Dataset(fn2)

# Print the variables contained in the file
print("Variables in the NetCDF file:")
for var in dataset.variables:
    print(var)
    

distance_to_coast = dataset.variables['distance_to_coast'][:]
latitude= dataset.variables['latitude'][:]
longitude=dataset.variables['longitude'][:] - 360


import matplotlib.pyplot as plt

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(longitude, latitude, c=distance_to_coast, cmap='viridis', marker='.', edgecolor='none', vmin=200000, vmax=400000)
plt.colorbar(label='Distance to Coast (units)')
plt.title('Distance to Coast')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.show()

#%% rad_surface_type_flag SWOT
from netCDF4 import Dataset


fn2 = 'C:/Users/mathi/OneDrive - Danmarks Tekniske Universitet/6_Semester_DTU/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

# Open the NetCDF file
dataset = Dataset(fn2)

# Print the variables contained in the file
print("Variables in the NetCDF file:")
for var in dataset.variables:
    print(var)
    

rad_surface_type_flag = dataset.variables['rad_surface_type_flag'][:]
latitude= dataset.variables['latitude'][:]
longitude=dataset.variables['longitude'][:] - 360
