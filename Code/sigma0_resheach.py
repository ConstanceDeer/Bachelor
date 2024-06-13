# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:07:14 2024

@author: const
"""
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
#%% One
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import swot_ssh_utils as swot
from scipy.ndimage import label, gaussian_filter
import xarray as xr


fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_543_005_20230605T143959_20230605T153105_PIB0_01.nc'
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
# visualize the data seperately for left and right swaths


#%
# initialize the class SSH_L2
dd=swot.SSH_L2()

#load data
dd.load_data(fn,lat_bounds=[-55.5,-53.5],)

side=['left','right']
Lo = []
La = []
OST =[]
n=0
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
    
    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm-=mss # remove the mean sea surface 
    
    # REMOVE ICEBERGS 
    # Smooth the data to reduce noise
    dtm_smoothed = gaussian_filter(dtm, sigma=4)
    
    # Compute the gradient
    gradient_x, gradient_y = np.gradient(dtm_smoothed)
    
    # Calculate the magnitude of the gradient
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    
    # Threshold the magnitude to identify significant changes
    threshold = 0.05 # Adjust this threshold based on your data
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
    msk=(lat[:,10]>-56)&(lat[:,10]<-53)
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:] 
    
    # Combine the arrays into a single 3D array
    #combined_array = np.stack((lon, lat, dtm), axis=0)
    
    # Save the combined array to a file (e.g., using NumPy's save function)
    #name = f"datafiles/geoid_27-06-2023-{side[n]}.npy"
    
    #np.save(name, combined_array)
    Lo = np.concatenate((Lo, lon.flatten()))
    La = np.concatenate((La, lat.flatten()))
    OST = np.concatenate((OST, dtm.flatten()))
    n=n+1
#%
# Create a figure and a grid of subplots
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))

# Add land, rivers, lakes, and coastlines to both subplots
for ax in axs:
    ax.add_feature(cfeature.LAND)

# Add states and provinces to both subplots
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

for ax in axs:
    ax.add_feature(states_provinces, edgecolor='gray')

# Plot the data using scatter in the first subplot
cax = axs[0].scatter(Lo, La, vmin=-0.15, vmax=0.1, c=OST, s=8, cmap='Spectral', transform=ccrs.PlateCarree())

# Add colorbar and gridlines to the first subplot
gl = axs[0].gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
# Set the extent of the map in the first subplot
axs[0].set_extent([-42, -41.4, -54.1, -53.6], crs=ccrs.PlateCarree())

# Add text annotation to the first subplot
date_part = fn.split('/')[-1].split('_')[7].split('T')[0]
axs[0].text(0.88, 0.95, date_part, horizontalalignment='left', fontsize=9, transform=axs[0].transAxes, color='black')

# Adjust figure size and DPI
fig.set_dpi(800)

axs[0].add_feature(cfeature.LAND)
axs[0].add_feature(cfeature.OCEAN)
axs[0].add_feature(cfeature.COASTLINE)
cbar = plt.colorbar(cax, ax=axs[0],
                    shrink=0.9, orientation='horizontal',
                    pad=0.1, aspect=50,
                    label='Sea Surface Height Anomaly [m]')
gl.right_labels = False  # Remove right side labels
gl.top_labels = False
#plt.axis('equal')
fig = plt.gcf()
fig.set_dpi(800)

#% Make sigma0 part of plot
import swot_ssh_utils as swot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
# Skud ud: https://github.com/SWOT-community/SWOT-OpenToolkit/blob/main/examples/unsmoothed_sea_ice_250m.ipynb

def plot_unsmoothed_sea_ice_sig0(ax, fn, cc, s3sys=None, test_plot=False, for_google_earth=False):
    dd = swot.SSH_L2()
    dd.load_data(fn, s3sys=s3sys, lat_bounds=[])
    pn = fn.split('/')[-1].split('_')[6]
    cn = fn.split('/')[-1].split('_')[5]
    fnn = fn.split('/')[-1]

    lat = dd.left.latitude.data[:, 120]
    lon = dd.left.longitude.data[:, 120]
    if len(cc['lat_bounds']) == 2:
        lat0, lat1 = cc['lat_bounds']
        msk = (lat > lat0) & (lat < lat1+0.5)
    else:
        lon0, lon1 = cc['lon_bounds']
        msk = (lon > lon0) & (lon < lon1)

    lons, lats = [], []
    if test_plot:
        skip = 5
    else:
        skip = 1
    for i, tmp in enumerate([dd.left, dd.right]):
        lon, lat = tmp.longitude.data[msk, :], tmp.latitude.data[msk, :]
        lons.append(lon)
        lats.append(lat)
        dtm = tmp.sig0_karin_2[msk, :].data
        dtm -= np.nanmin(dtm + 1e-10)
        mm = np.isfinite(dtm)
        dtm[mm] = np.log2(dtm[mm])
        cax2 = ax.pcolor(lon[::skip, ::skip], lat[::skip, ::skip], 
                         dtm[::skip, ::skip],
                         cmap='bone',
                         transform=ccrs.PlateCarree())

    ax.set_extent([cc['lon_bounds'][0], cc['lon_bounds'][1], cc['lat_bounds'][0], cc['lat_bounds'][1]], crs=ccrs.PlateCarree())
    
    if for_google_earth: 
        ax.set_frame_on(False)
        ax.axis('off')
        ax.set_facecolor("none")
        ax.set_aspect('equal', 'box')
        plt.savefig('../media/figures/Unsmoothed_sig0_images/'+fn.split('/')[-1].split('.')[0]+'_noborder.png', 
                    dpi=300, bbox_inches='tight', pad_inches=0)
    else:
        gl0 = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        #ax.text(1.84, 0.02, date_part, horizontalalignment='left', fontsize=12, transform=axs[0].transAxes, color='black')
        cbar = plt.colorbar(cax2, ax=ax, shrink=0.9, orientation='horizontal', pad=0.1, aspect=50,label='log2(sig0)')
        #plt.axis('equal')
        fig = plt.gcf()
        fig.set_dpi(800)
        gl0.right_labels = False  # Remove right side labels
        gl0.top_labels = False

# Assuming fig and axs are defined elsewhere
# Call the function to plot in the second subplot
# first: 'lat_bounds':[-52.3, -52.1],'lon_bounds':[-42, -41.675]
# Second: [-42, -41.4, -54.1, -53.6]

cc={'pn':'021','lat_bounds':[-54.1, -53.6],'lon_bounds':[-42, -41.39],'day':'0565','region_name':'Antarctica'}

plot_unsmoothed_sea_ice_sig0(axs[1], fn, cc)  # Adjust the parameters accordingly

plt.show()

#%% Two
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import swot_ssh_utils as swot
from scipy.ndimage import label, gaussian_filter
import xarray as xr


fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_541_005_20230603T145844_20230603T154949_PIB0_01.nc'
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
# visualize the data seperately for left and right swaths

#%
# initialize the class SSH_L2
dd=swot.SSH_L2()

#load data
dd.load_data(fn,lat_bounds=[-53.5,-51.5],)

side=['left','right']
Lo = []
La = []
OST =[]
n=0
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
    
    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm-=mss # remove the mean sea surface 
    
    # REMOVE ICEBERGS 
    # Smooth the data to reduce noise
    dtm_smoothed = gaussian_filter(dtm, sigma=4)
    
    # Compute the gradient
    gradient_x, gradient_y = np.gradient(dtm_smoothed)
    
    # Calculate the magnitude of the gradient
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    
    # Threshold the magnitude to identify significant changes
    threshold = 0.05 # Adjust this threshold based on your data
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
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:] 
    
    # Combine the arrays into a single 3D array
    #combined_array = np.stack((lon, lat, dtm), axis=0)
    
    # Save the combined array to a file (e.g., using NumPy's save function)
    #name = f"datafiles/geoid_27-06-2023-{side[n]}.npy"
    
    #np.save(name, combined_array)
    Lo = np.concatenate((Lo, lon.flatten()))
    La = np.concatenate((La, lat.flatten()))
    OST = np.concatenate((OST, dtm.flatten()))
    n=n+1
#%
# Create a figure and a grid of subplots
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))

# Add land, rivers, lakes, and coastlines to both subplots
for ax in axs:
    ax.add_feature(cfeature.LAND)

# Add states and provinces to both subplots
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

for ax in axs:
    ax.add_feature(states_provinces, edgecolor='gray')

# Plot the data using scatter in the first subplot
cax = axs[0].scatter(Lo, La, vmin=-0.15, vmax=0.1, c=OST, s=8, cmap='Spectral_r', transform=ccrs.PlateCarree())

# Add colorbar and gridlines to the first subplot
gl = axs[0].gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
# Set the extent of the map in the first subplot
axs[0].set_extent([-42, -41.6, -52.4, -52.1], crs=ccrs.PlateCarree())

# Add text annotation to the first subplot
date_part = fn.split('/')[-1].split('_')[7].split('T')[0]
#axs[0].text(0.88, 0.95, date_part, horizontalalignment='left', fontsize=9, transform=axs[0].transAxes, color='black')

# Adjust figure size and DPI
fig.set_dpi(800)

axs[0].add_feature(cfeature.LAND)
axs[0].add_feature(cfeature.OCEAN)
axs[0].add_feature(cfeature.COASTLINE)
cbar = plt.colorbar(cax, ax=axs[0],
                    shrink=0.9, orientation='horizontal',
                    pad=0.1, aspect=50,
                    label='Sea Surface Height Anomaly, SSHA [m]')
gl.right_labels = False  # Remove right side labels
gl.top_labels = False
#plt.axis('equal')
fig = plt.gcf()
fig.set_dpi(800)

#% Make sigma0 part of plot
import swot_ssh_utils as swot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
# Skud ud: https://github.com/SWOT-community/SWOT-OpenToolkit/blob/main/examples/unsmoothed_sea_ice_250m.ipynb

def plot_unsmoothed_sea_ice_sig0(ax, fn, cc, s3sys=None, test_plot=False, for_google_earth=False):
    dd = swot.SSH_L2()
    dd.load_data(fn, s3sys=s3sys, lat_bounds=[])
    pn = fn.split('/')[-1].split('_')[6]
    cn = fn.split('/')[-1].split('_')[5]
    fnn = fn.split('/')[-1]

    lat = dd.left.latitude.data[:, 120]
    lon = dd.left.longitude.data[:, 120]
    if len(cc['lat_bounds']) == 2:
        lat0, lat1 = cc['lat_bounds']
        msk = (lat > lat0) & (lat < lat1+0.5)
    else:
        lon0, lon1 = cc['lon_bounds']
        msk = (lon > lon0) & (lon < lon1)

    lons, lats = [], []
    if test_plot:
        skip = 5
    else:
        skip = 1
    for i, tmp in enumerate([dd.left, dd.right]):
        lon, lat = tmp.longitude.data[msk, :], tmp.latitude.data[msk, :]
        lons.append(lon)
        lats.append(lat)
        dtm = tmp.sig0_karin_2[msk, :].data
        dtm -= np.nanmin(dtm + 1e-10)
        mm = np.isfinite(dtm)
        #dtm[mm] = np.log2(dtm[mm])
        dtm[mm] = dtm[mm]
        cax2 = ax.pcolor(lon[::skip, ::skip], lat[::skip, ::skip], 
                         dtm[::skip, ::skip],
                         cmap='bone',
                         transform=ccrs.PlateCarree())

    ax.set_extent([cc['lon_bounds'][0], cc['lon_bounds'][1], cc['lat_bounds'][0], cc['lat_bounds'][1]], crs=ccrs.PlateCarree())
    
    if for_google_earth: 
        ax.set_frame_on(False)
        ax.axis('off')
        ax.set_facecolor("none")
        ax.set_aspect('equal', 'box')
        plt.savefig('../media/figures/Unsmoothed_sig0_images/'+fn.split('/')[-1].split('.')[0]+'_noborder.png', 
                    dpi=300, bbox_inches='tight', pad_inches=0)
    else:
        gl0 = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE)
        #ax.text(1.84, 0.02, date_part, horizontalalignment='left', fontsize=12, transform=axs[0].transAxes, color='black')
        cbar = plt.colorbar(cax2, ax=ax, shrink=0.9, orientation='horizontal', pad=0.1, aspect=50,label='Normalized Radar Cross Section (NRCS)')
        #plt.axis('equal')
        fig = plt.gcf()
        fig.set_dpi(800)
        gl0.right_labels = False  # Remove right side labels
        gl0.top_labels = False

# Assuming fig and axs are defined elsewhere
# Call the function to plot in the second subplot
# first: 'lat_bounds':[-52.3, -52.1],'lon_bounds':[-42, -41.675]
# Second: [-42, -41.4, -54.1, -53.6]

cc={'pn':'021','lat_bounds':[-52.4, -52.1],'lon_bounds':[-42, -41.6],'day':'0565','region_name':'Antarctica'}

plot_unsmoothed_sea_ice_sig0(axs[1], fn, cc)  # Adjust the parameters accordingly

plt.show()

#%% Front page
import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import swot_ssh_utils as swot
from scipy.ndimage import label, gaussian_filter
import xarray as xr


fn = 'C:/Users/const/Documents/Bachelor/250m/SWOT_L2_LR_SSH_Unsmoothed_541_005_20230603T145844_20230603T154949_PIB0_01.nc'
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
# visualize the data seperately for left and right swaths

#%
# initialize the class SSH_L2
dd=swot.SSH_L2()

#load data
dd.load_data(fn,lat_bounds=[-53.5,-51.5],)

side=['left','right']
Lo = []
La = []
OST =[]
n=0
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
    
    mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
    dtm-=mss # remove the mean sea surface 
    
    # REMOVE ICEBERGS 
    # Smooth the data to reduce noise
    dtm_smoothed = gaussian_filter(dtm, sigma=4)
    
    # Compute the gradient
    gradient_x, gradient_y = np.gradient(dtm_smoothed)
    
    # Calculate the magnitude of the gradient
    gradient_magnitude = np.sqrt(gradient_x**2 + gradient_y**2)
    
    # Threshold the magnitude to identify significant changes
    threshold = 0.05 # Adjust this threshold based on your data
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
    
    # remove the along-track mean, the ad hoc way of removing cross-swath bias
    dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:] 
    
    # Combine the arrays into a single 3D array
    #combined_array = np.stack((lon, lat, dtm), axis=0)
    
    # Save the combined array to a file (e.g., using NumPy's save function)
    #name = f"datafiles/geoid_27-06-2023-{side[n]}.npy"
    
    #np.save(name, combined_array)
    Lo = np.concatenate((Lo, lon.flatten()))
    La = np.concatenate((La, lat.flatten()))
    OST = np.concatenate((OST, dtm.flatten()))
    n=n+1
#%
# Create a figure and a grid of subplots
fig, axs = plt.subplots( subplot_kw=dict(projection=ccrs.PlateCarree()))

# Add land, rivers, lakes, and coastlines to both subplots

axs.add_feature(cfeature.LAND)

# Add states and provinces to both subplots
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='10m',
    facecolor='none')

axs.add_feature(states_provinces, edgecolor='gray')

# Plot the data using scatter in the first subplot
cax = axs.scatter(Lo, La, vmin=-0.15, vmax=0.1, c=OST, s=8, cmap='Spectral_r', transform=ccrs.PlateCarree())

# Add colorbar and gridlines to the first subplot
#gl = axs.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')
# Set the extent of the map in the first subplot
axs.set_extent([-42.01, -41.67, -52.34, -52.06], crs=ccrs.PlateCarree())

# Add text annotation to the first subplot
#date_part = fn.split('/')[-1].split('_')[7].split('T')[0]
#axs[0].text(0.88, 0.95, date_part, horizontalalignment='left', fontsize=9, transform=axs[0].transAxes, color='black')

# Adjust figure size and DPI
fig.set_dpi(800)

#axs[0].add_feature(cfeature.LAND)
#axs[0].add_feature(cfeature.OCEAN)
#axs[0].add_feature(cfeature.COASTLINE)
#cbar = plt.colorbar(cax, ax=axs[0], shrink=0.9, orientation='horizontal', pad=0.1, aspect=50, label='Sea Surface Height Anomaly, SSHA [m]')
#gl.right_labels = False  # Remove right side labels
#gl.top_labels = False
#plt.axis('equal')
fig = plt.gcf()
fig.set_dpi(800)