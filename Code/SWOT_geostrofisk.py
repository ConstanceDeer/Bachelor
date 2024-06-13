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
from scipy.ndimage import median_filter

fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_566_005_20230628T110430_20230628T115536_PIB0_01.nc'
#fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_539_005_20230601T151728_20230601T160833_PIB0_01.nc'

basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

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

    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm-=mss # remove the mean sea surface 
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:]  # trend ssha
   

    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm+=mss # remove the mean sea surface 
    
    #dtm+=mss # mss add
    
    mss=interpolate_data(basic,'geoid',lon, lat) # get the geoid at the same points
    dtm-=mss # remove  geoid
    
    #dtm = gaussian_filter(dtm, sigma=5)
    #dtm=median_filter(dtm, size=17)

    combined_array = np.stack((tmp.longitude.data, tmp.latitude.data, dtm), axis=0)

    # Save the combined array
    name = f"Slut/geoid_19-06-2023-{side_name}.npy"
    np.save(name, combined_array)

    print(f"Saved data for {side_name} swath to {name}")

    lon, lat = tmp.longitude.data, tmp.latitude.data
    m = ~np.isnan(lon + lat).flatten()  # Mask out the nan values


 # plot the data using scatter , vmin=-0.1, vmax=0.1
    cax=ax.scatter(lon.flatten()[m], lat.flatten()[m], vmin=-1.25, vmax=-1,
                   c=dtm.flatten()[m], s=0.05, 
                   cmap='Spectral_r' ,transform=ccrs.PlateCarree()) 
    
    
# add colorbar and gridlines
gl = ax.gridlines(draw_labels=True, color='gray', 
                  alpha=0.5, linestyle='--')


cbar = plt.colorbar(cax, ax=ax,shrink=0.7)
cbar.ax.tick_params(labelsize=12, )
cbar.set_label('Ocean Surface Topography, OST [m]', rotation=270, fontsize=16,  labelpad=-70)

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
fig = plt.gcf()
fig.set_dpi(800)  # Set the DPI (dots per inch)

# set the extent of the map
lonmin,lonmax,latmin,latmax=np.nanmin(lon-1),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat)
ax.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())


plt.show()

#%% data 
import numpy as np

left = 'C:/Users/const/Documents/Bachelor/Slut/geoid_19-06-2023-left.npy'
data_left = np.load(left)

right = 'C:/Users/const/Documents/Bachelor/Slut/geoid_19-06-2023-right.npy'
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
step = 20

from netCDF4 import Dataset

Expert ='C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Expert_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
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
file_path = "C:/Users/const/Documents/Bachelor/Slut/data.npy"

# Gem dataene i en npy-fil
np.save(file_path, data)

#%


fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop
# Plot the path
plt.quiver(lon, lat,geo_u,geo_v)
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')
plt.title('Path of the red vector')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop
# Plot the path
plt.scatter(lon, lat, c=geo_v, s=100,vmin=-0.25, vmax=0.25, cmap='Spectral_r')
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')
plt.title('Updown')
plt.colorbar( label='Vertical velocity [ Pos = up, neg=down]')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop
# Plot the path
plt.scatter(lon, lat, c=geo_u, s=100,vmin=-0.25, vmax=0.25, cmap='Spectral_r')
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')
plt.title('Updown')
plt.colorbar( label='Vertical velocity [ Pos = right, neg=left]')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.show()

#%% Kurve - nearest
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as nc

# Load elevation data
fn_elevation =  'C:/Users/const/Documents/Bachelor/gebco_2023_n-48.0_s-57.0_w-52.0_e-35.0.nc'
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


directory_path = 'C:/Users/const/Documents/Bachelor/Slut/Gradient'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]

# Starting coordinates
#Isbjerg 1
start_lon = -42.7
start_lat = -52.85

#Isbjerg 2
#start_lon = -41.35
#start_lat = -53.355


# Assume dt is the time step (in hours)
dt = 1

# Lists to store coordinates for the path
path_lon = [start_lon]
path_lat = [start_lat]

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop

for i, file_path in enumerate(file_paths):
    #data = np.load(file_path)
    data = np.load(file_path, allow_pickle=True)
    
    # Print some information about the loaded data
    lon = data[0] 
    lat = data[1]
    u = data[2] * 3.6 
    v = data[3] * 3.6 
    
        
    for t in range(25):
    
        # Find the index of the nearest non-NaN point
        nearest_index = find_nearest_non_nan(lon, lat, u, v, start_lon, start_lat)
        # Interpolate wind data at the nearest point
        u_ny = u[nearest_index]
        v_ny = v[nearest_index]
        
        print(f"u_ny: {u_ny}, v_ny: {v_ny}")
        # Update the coordinates based on the displacement
        delta_lon_km = u_ny * dt  # Longitudinal displacement in kilometers
        delta_lat_km = v_ny * dt  # Latitudinal displacement in kilometers
        
        # Convert longitudinal displacement from kilometers to degrees
        delta_lon_deg = delta_lon_km / (111.32 * np.cos(np.radians(start_lat)))
        
        # Calculate km per degree latitude at the start latitude
        km_per_degree_lat = 111.32 * np.cos(np.radians(start_lat))
        # Convert latitudinal displacement from kilometers to degrees
        delta_lat_deg = delta_lat_km / km_per_degree_lat
        
        # Update the coordinates based on the displacement
        start_lon = start_lon + delta_lon_deg
        start_lat = start_lat + delta_lat_deg
        
        print(f"Koordinater efter {t} timer: {start_lon}, {start_lat}")
        
        # Add the updated coordinates to the path
        path_lon.append(start_lon)
        path_lat.append(start_lat)
    
    # Check if the red point is within the plot limits
    if start_lon < -43 or start_lon > -40 or start_lat < -53.5 or start_lat > -51.5:
        break

# Plot the path
plt.plot(path_lon, path_lat, marker='o', markersize=5, linestyle='-', color='red')
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='viridis', origin='lower', aspect='auto', vmin=-4500,vmax=1)
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')
plt.title('Path of the red vector')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.colorbar(elevation_plot, ax=ax, orientation='vertical', label='Elevation (m)')  # Add colorbar
plt.show()

#%% Kurve interpolere

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os

# Example data (replace with actual loading logic)
directory_path = 'C:/Users/const/Documents/Bachelor/Slut/Gradient'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".npy")]

# Starting coordinates
start_lon = -41.35
start_lat = -53.355

# Lists to store coordinates for the path
path_lon = [start_lon]
path_lat = [start_lat]

fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis outside the loop

for i, file_path in enumerate(file_paths):
    data = np.load(file_path, allow_pickle=True)
    
    lon = data[0]
    lat = data[1]
    u = data[2] * 3.6
    v = data[3] * 3.6
    
    # Check for NaNs and remove them
    valid_indices = ~np.isnan(lon) & ~np.isnan(lat) & ~np.isnan(u) & ~np.isnan(v)
    lon = lon[valid_indices]
    lat = lat[valid_indices]
    u = u[valid_indices]
    v = v[valid_indices]

    # Define the points for interpolation
    points = np.column_stack((lon, lat))
    
    for t in range(25):
        # Perform the interpolation directly on lon and lat
        u_ny = griddata(points, u, (start_lon, start_lat), method='linear')
        v_ny = griddata(points, v, (start_lon, start_lat), method='linear')
        
        print(f"u_ny: {u_ny}, v_ny: {v_ny}")
        
        # Update the coordinates based on the displacement
        delta_lon_km = u_ny * dt  # Longitudinal displacement in kilometers
        delta_lat_km = v_ny * dt  # Latitudinal displacement in kilometers
        
        # Convert longitudinal displacement from kilometers to degrees
        delta_lon_deg = delta_lon_km / (111.32 * np.cos(np.radians(start_lat)))
        
        # Calculate km per degree latitude at the start latitude
        km_per_degree_lat = 111.32 * np.cos(np.radians(start_lat))
        # Convert latitudinal displacement from kilometers to degrees
        delta_lat_deg = delta_lat_km / km_per_degree_lat
        
        # Update the coordinates based on the displacement
        start_lon = start_lon + delta_lon_deg
        start_lat = start_lat + delta_lat_deg
        
        print(f"Coordinates after {t} hours: {start_lon}, {start_lat}")
        
        # Add the updated coordinates to the path
        path_lon.append(start_lon)
        path_lat.append(start_lat)
    
    # Check if the point is within the plot limits
    if start_lon < -43 or start_lon > -40 or start_lat < -53.5 or start_lat > -51.5:
        break

# Plot the path
plt.plot(path_lon, path_lat, marker='o', markersize=5, linestyle='-', color='red')
elevation_plot = ax.imshow(elevation_var, extent=(lon_elevation.min(), lon_elevation.max(), lat_elevation.min(), lat_elevation.max()), cmap='viridis', origin='lower', aspect='auto', vmin=-4500,vmax=1)
plt.xlabel('Longitude (degrees_east)')
plt.ylabel('Latitude (degrees_north)')
plt.title('Path of the red vector')
plt.ylim([-53.5, -51.5])
plt.xlim([-43, -40])
plt.grid(True)
plt.colorbar(elevation_plot, ax=ax, orientation='vertical', label='Elevation (m)')  # Add colorbar
plt.show()
