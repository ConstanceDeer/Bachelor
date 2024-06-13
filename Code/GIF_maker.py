#%%  MSSH

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

directory_path = 'C:/Users/const/Documents/Bachelor/250m'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]


#fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

# initialize the class SSH_L2
dd=swot.SSH_L2()

# print the variables in the file
# load data into the class
# the lat_bounds keyword is used to subset the data
# the return is dd.left and dd.right, which correspond to the left and right swath

image_paths = []

for i in range(len(file_paths)):
    chosen_file_index = i # Change this to the index of the file you want to load
    fn = file_paths[chosen_file_index]
    dd.load_data(fn,lat_bounds=[-55.5,-53.5],)
    #fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

    
    # Create a map using cartopy
    fig, ax = plt.subplots(figsize=(6, 12), subplot_kw=dict(projection=ccrs.PlateCarree()))
    
    # add land, rivers, lakes, and coastlines
    ax.add_feature(cartopy.feature.LAND)
    
    # add states and provinces
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')
    
    #min_lat = -62.95
    #max_lat = -60.3
    
    # visualize the data seperately for left and right swaths
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
        #sys.exit()
        # the following line find the latitude limit
        # to bound the area for along-track mean removal
        # you can change the bounds to get a better fit for your region. 
        # The 10th column of the lat array is used as an approximate. 
      
        #msk=(lat[:,10]>-53)&(lat[:,10]<-51)
       
        #dtm = dtm -  np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:] 
        
        cossover=interpolate_data(basic,'height_cor_xover',lon, lat) # cross-over bias
        dtm+=cossover # remove the mean sea surface 
        
        # plot the data using scatter , vmin=-0.1, vmax=0.1
        cax=ax.scatter(lon.flatten()[m], lat.flatten()[m], vmin=-0.1, vmax=0.1, 
                       c=dtm.flatten()[m], s=0.05, 
                       cmap='Spectral' ,transform=ccrs.PlateCarree())
        
       
    # add colorbar and gridlines
    gl = ax.gridlines(draw_labels=True, color='gray', 
                      alpha=0.5, linestyle='--')
    cbar = plt.colorbar(cax, ax=ax)
    # remove the right side ylabels
    gl.right_labels = False
    #cbar.set_label('ssh_karin_2')
    
    date_part = fn.split('/')[-1].split('_')[7].split('T')[0]
    ax.set_title(date_part + ' - Track 05')
    
    #ax.plot([-40.8, -42.54], [-52.8, -52.3], marker='o', linestyle='-', color='red', transform=ccrs.PlateCarree())
    # Adjust figure size and DPI
    fig = plt.gcf()
    fig.set_size_inches(20, 12)  # Set the figure size in inches (width, height)
    fig.set_dpi(400)  # Set the DPI (dots per inch)
    
    # set the extent of the map
    lonmin,lonmax,latmin,latmax=np.nanmin(lon-1),np.nanmax(lon),np.nanmin(lat),np.nanmax(lat)
    ax.set_extent([lonmin,lonmax,latmin,latmax], crs=ccrs.PlateCarree())
  
    
    # Save the figure as an image
    image_path = f"testtest/plot_{i}.png"
    plt.savefig(image_path)
    plt.show()
    plt.close()  # Close the figure to release memory
    image_paths.append(image_path)
    
#os.remove('output.gif')
# Create a GIF from the saved images
#%%
with imageio.get_writer('June_250m_test.gif', mode='I', duration=0.5) as writer:
    for image_path in image_paths:
        image = imageio.imread(image_path)
        writer.append_data(image)
#%%
# Optionally, remove the individual image files
for image_path in image_paths:
    os.remove(image_path)
    os.remove

#%% Inter_polate data
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

#%% Geoid

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
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
directory_path = 'C:/Users/const/Documents/Bachelor/250m'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

#fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

# initialize the class SSH_L2
dd=swot.SSH_L2()

# print the variables in the file
# load data into the class
# the lat_bounds keyword is used to subset the data
# the return is dd.left and dd.right, which correspond to the left and right swath

image_paths = []

for i in range(len(file_paths)):
    chosen_file_index = i # Change this to the index of the file you want to load
    fn = file_paths[chosen_file_index]
    dd.load_data(fn,lat_bounds=[-53.5,-51.5],)
    #fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

    
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
    
    #min_lat = -62.95
    #max_lat = -60.3
    
    # visualize the data seperately for left and right swaths
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
        
        #mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
        #dtm-=mss # remove the mean sea surface
        
        cossover=interpolate_data(basic,'height_cor_xover',lon, lat) # cross-over bias
        dtm+=cossover # remove the mean sea surface 
        
        msk=(lat[:,10]>-53)&(lat[:,10]<-50)
        
        # remove the along-track mean, the ad hoc way of removing cross-swath bias
        dtm=dtm - np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:]  # trend ssha
       
        mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
        dtm+=mss # remove the mean sea surface 
        
        
        geoid=interpolate_data(basic,'geoid',lon, lat) # get the geoid at the same points
        dtm-=geoid # remove  geoid
        
        
        # plot the data using scatter , vmin=-0.1, vmax=0.1, vmin=-1.25, vmax=-0.75
        cax=ax.scatter(lon.flatten()[m], lat.flatten()[m],  vmin=-2, vmax=5,
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
  
    
    # Save the figure as an image
    image_path = f"testtest/plot_{i}.png"
    plt.savefig(image_path)
    plt.show()
    plt.close()  # Close the figure to release memory
    image_paths.append(image_path)

#%%
# Create a GIF from the saved images
#%
with imageio.get_writer('June_250m_Geoid_5351.gif', mode='I', duration=1) as writer:
    for image_path in image_paths:
        image = imageio.imread(image_path)
        writer.append_data(image)
#%%
# Optionally, remove the individual image files
for image_path in image_paths:
    os.remove(image_path)
    os.remove
    
#%% MSS

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
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
directory_path = 'C:/Users/const/Documents/Bachelor/250m'

# Create a list of file paths
file_paths = [os.path.join(directory_path, file) for file in os.listdir(directory_path) if file.endswith(".nc")]
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'

#fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

# initialize the class SSH_L2
dd=swot.SSH_L2()

# print the variables in the file
# load data into the class
# the lat_bounds keyword is used to subset the data
# the return is dd.left and dd.right, which correspond to the left and right swath

image_paths = []

for i in range(len(file_paths)):
    chosen_file_index = i # Change this to the index of the file you want to load
    fn = file_paths[chosen_file_index]
    dd.load_data(fn,lat_bounds=[-53.5,-51.5],)
    #fn = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Unsmoothed_565_005_20230627T111353_20230627T120458_PIB0_01.nc'

    
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
    
    #min_lat = -62.95
    #max_lat = -60.3
    
    # visualize the data seperately for left and right swaths
    # Loads data from the two tracks and compine into 3 arrays.
    for tmp in [dd.left,dd.right]:
        #find data points in the open ocean (flag==0)
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
       
       mss=interpolate_data(basic,'mean_sea_surface_cnescls',lon, lat) # get the geoid at the same points
       dtm-=mss # remove the mean sea surface 
       
       cossover=interpolate_data(basic,'height_cor_xover',lon, lat) # cross-over bias
       dtm+=cossover # remove the mean sea surface

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
      
       cluster_indices = np.where(significant_changes_mask)
       dtm[cluster_indices]=np.nan
      
       # remove the along-track mean, the ad hoc way of removing cross-swath bias
       msk=(lat[:,10]>-53)&(lat[:,10]<-51)
       dtm = dtm -  np.nanmean(dtm[msk,:],axis=0)[np.newaxis,:]  
       
       # plot the data using scatter , vmin=-0.1, vmax=0.1
       cax=ax.scatter(lon.flatten()[m], lat.flatten()[m], vmin=-0.15, vmax=0.15,
                      c=dtm.flatten()[m], s=0.05, 
                      cmap='Spectral_r' ,transform=ccrs.PlateCarree())
        
        
    # add colorbar and gridlines
    gl = ax.gridlines(draw_labels=True, color='gray', 
                      alpha=0.5, linestyle='--')
    
  
    cbar = plt.colorbar(cax, ax=ax,shrink=0.7)
    cbar.ax.tick_params(labelsize=12, )
    cbar.set_label('Sea Surface Height Anomaly, SSHA [m]', rotation=270, fontsize=16,  labelpad=-75)
   
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
  
    
    # Save the figure as an image
    image_path = f"testtest/plot_{i}.png"
    plt.savefig(image_path)
    plt.show()
    plt.close()  # Close the figure to release memory
    image_paths.append(image_path)

#%%
# Create a GIF from the saved images
#%
with imageio.get_writer('June_250m_MSS_5351_0.15.gif', mode='I', duration=1) as writer:
    for image_path in image_paths:
        image = imageio.imread(image_path)
        writer.append_data(image)
#%%
# Optionally, remove the individual image files
for image_path in image_paths:
    os.remove(image_path)
    os.remove
    
    