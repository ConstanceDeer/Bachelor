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
import swot_ssh_utils as swot
import scipy
from scipy.fft import fft, ifft, fftfreq
import os
from scipy.interpolate import RectBivariateSpline as rbs
from scipy.interpolate import griddata
from xarray import open_dataset
import numpy as np
import xarray as xr
# skud ud: https://github.com/SWOT-community/SWOT-OpenToolkit/blob/main/src/swot_ssh_utils.py

#%PDS - Bjarkes functioner

def powerSpectralDensity(ssh, sampling, unit="m", tapering_f = "boxcar", tapering=8):
    # Setup
    N = ssh.size

    # Interpolate gaps
    mask = np.isnan(ssh)
    ssh[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), ssh[~mask])

    # Detrend
    ssh = scipy.signal.detrend(ssh, axis=0, type='linear')
    
    match tapering_f:
        case "boxcar":
            kernel = scipy.signal.windows.boxcar(N)
        case "tukey":
            kernel = scipy.signal.windows.tukey(N, alpha=0.1)
        case "hann":
            kernel = scipy.signal.windows.hann(N)
    ssh *= kernel

    # Perform FFT
    SSH = fft(ssh)[:N//2]
    k = fftfreq(N, sampling)[:N//2]
                
    # Calculate power spectral density
    kernel_norm = kernel.sum()/N
    psd = 2 * (sampling)/(N*kernel_norm)* ((abs(SSH))**2)
    
    # Perform unit correction
    match unit:
        case "m":
            psd *= 1
        case "cm":
            psd *= 1e4

    # Convert units
    psd /= 1000 # Convert from [unit]^2/cpm to [unit]^2/cpkm
    k *= 1000 # Convert from cpm to cpkm

    # Remove first element (is inf)
    k = k[1:]
    psd = psd[1:]
    
    return psd, k
# ========================================================================
def periodogram_list(X, conf_lvl=68, sampling=2000, unit="m", tapering_f = "boxcar", tapering=8, skip_n=1):

    # Make sure largest element comes first
    X = sorted(X, key=len, reverse=True)

    N = np.array(X[0]).size
    psd_stack = np.zeros(N//2-1)
    k = fftfreq(N, sampling)[:N//2]
    k *= 1000 # Convert from cpm to cpkm
    k = k[1:]

    for i in range(0,len(X),skip_n):
        x_profile = np.array(X[i])

        # Skip if empty segment
        if (~np.isnan(x_profile)).sum() == 0:
            continue

        # Compute PSD
        psd_, k_ = powerSpectralDensity(x_profile, sampling=sampling, unit=unit, 
                                        tapering_f=tapering_f, tapering=tapering)
        
        # Fit up to N size
        psd_temp = np.zeros(N//2-1)
        psd_temp[:] = np.nan
        psd_temp[(N//2-1-psd_.size):] = psd_

        psd_stack = np.c_[psd_stack, psd_temp]
    psd_stack = psd_stack[:,1:]

    psd = np.nanmedian(psd_stack,axis=1)
    conf_lvl = 68
    cd_interval = np.nanpercentile(psd_stack,[100-conf_lvl,conf_lvl],axis=1)

    return psd, k, cd_interval
# ========================================================================

#%% Load data

right= np.load('datafiles/geoid_27-06-2023-right.npy')
left=np.load('datafiles/geoid_27-06-2023-left.npy')

ssha_left=left[2,:,:]
ssha_right=right[2,:,:]

longitude_left=left[0,:,:]
latitude_left=left[1,:,:]

longitude_right=right[0,:,:]
latitude_right=right[1,:,:]

#%% Only part of each swath 

psdl, kl, sdl = periodogram_list(ssha_left[515:565,:], conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)
psdr, kr, sdr = periodogram_list(ssha_right[380:430,:], conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)

# Plot
fig, ax = plt.subplots(figsize=(10,4))
ax.loglog(kl, psdl, 'r',label="Left Swath")
ax.fill_between(kl, sdl[0,:], sdl[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.loglog(kr, psdr, 'b',label="Right Swath")
ax.fill_between(kr, sdr[0,:], sdr[1,:], color='b',edgecolor=None, alpha=.2,label="_nolegend_")

ax.legend()

# Setup plot settings
ax.set_yscale("log")
ax.set_xlim([50**-1,2])
ax.set_xscale("log")
ax.set_xlabel("Wavenumber (cpkm)")
ax.set_ylabel("PSD (cm$^2$/cpkm)")
ax.set_yticks(np.logspace(-6,6,13))
ax.grid(which="major", alpha=0.8)
ax.grid(which="minor", alpha=0.4, linestyle=":")
fig = plt.gcf()

fig.set_dpi(800)
# Create second axis

km_lims = [20, 0.5]
ax2 = ax.twiny()
ax2.set_xlabel("Wavelength (km)")
ax2.set_xscale("log")

# Synchronize the limits of ax2 with ax
ax2.set_xlim(ax.get_xlim())

# Define the wavelengths for the ticks
ticks_wavelengths = np.array([0.6,1,2,3, 5,10,50])  # Add any desired wavelengths here

# Calculate the positions of the ticks on the second x-axis
ticks_positions = 1 / ticks_wavelengths

# Set the ticks and labels on the second x-axis
ax2.set_xticks(ticks_positions)
ax2.set_xticklabels([str(t) + " km" for t in ticks_wavelengths], fontsize=8)

plt.plot()

#%% For entier swath 

#Load data
psdl, kl, sdl = periodogram_list(ssha_left, conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=10)
psdr, kr, sdr = periodogram_list(ssha_right, conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=10)

#Start figure 
fig, ax = plt.subplots(figsize=(10,7))

ax.loglog(kl, psdl, 'r',label="Left Swath")
ax.fill_between(kl, sdl[0,:], sdl[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.loglog(kr, psdr, 'b',label="Right Swath")
ax.fill_between(kr, sdr[0,:], sdr[1,:], color='b',edgecolor=None, alpha=.2,label="_nolegend_")

ax.legend()

# Setup plot settings
ax.set_yscale("log")
ax.set_xlim([50**-1,2])
ax.set_xscale("log")
ax.set_xlabel("Wavenumber (cpkm)")
ax.set_ylabel("PSD (cm$^2$/cpkm)")
ax.set_yticks(np.logspace(-6,6,13))
ax.grid(which="major", alpha=0.8)
ax.grid(which="minor", alpha=0.4, linestyle=":")
fig = plt.gcf()

fig.set_dpi(300)
# Create second axis

km_lims = [20, 0.5]
ax2 = ax.twiny()
ax2.set_xlabel("Wavelength (km) - Whole")
ax2.set_xscale("log")

# Synchronize the limits of ax2 with ax
ax2.set_xlim(ax.get_xlim())

# Define the wavelengths for the ticks
ticks_wavelengths = np.array([0.6,1,2,3, 5,10,50])  # Add any desired wavelengths here

# Calculate the positions of the ticks on the second x-axis
ticks_positions = 1 / ticks_wavelengths

# Set the ticks and labels on the second x-axis
ax2.set_xticks(ticks_positions)
ax2.set_xticklabels([str(t) + " km" for t in ticks_wavelengths], fontsize=8)

plt.plot()
plt.show()

#%%
def interpolate_data(file,value,longitude, latitude):
    from scipy.interpolate import RectBivariateSpline as rbs
    from scipy.interpolate import griddata
    from xarray import open_dataset
    import numpy as np
    import requests
    import xarray as xr 

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

    # Interpolate the MSS22 data onto the input grid using griddata
    aa = griddata((latm[mm], lonm[mm]), geoid.flatten()[mm], (latitude.flatten()[~m], longitude.flatten()[~m]))
   
    
    # Initialize an array for interpolated MSS22 data
    geoid = np.ones((longitude.size))

    # Assign missing values and interpolated data to the MSS22 array
    geoid[m] = np.nan
    geoid[~m] = aa

    return geoid.reshape(longitude.shape)


#%% For the SSM for baseline 
basic = 'C:/Users/const/Documents/Bachelor/SWOT_L2_LR_SSH_Basic_565_005_20230627T111642_20230627T120459_PIB0_01.nc'
lon=left[0,:,:]
lat=left[1,:,:]

geoid=interpolate_data(basic,'geoid',lon, lat) # get the mean sea surface

choosen_area=geoid[515:565,:]

psdn, kn, sdn = periodogram_list(choosen_area, conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)
psdl, kl, sdl = periodogram_list(ssha_left[515:565,:], conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)
# Plot
fig, ax = plt.subplots(figsize=(10,7))
ax.loglog(kn, psdn, 'r',label="Geoid")
ax.fill_between(kn, sdn[0,:], sdn[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.loglog(kl, psdl, 'b',label="Left Swath")
ax.fill_between(kl, sdl[0,:], sdl[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.legend()

# Setup plot settings
ax.set_yscale("log")
ax.set_xlim([50**-1,2])
ax.set_xscale("log")
ax.set_xlabel("Wavenumber (cpkm)")
ax.set_ylabel("PSD (cm$^2$/cpkm)")
ax.set_yticks(np.logspace(-6,6,13))
ax.grid(which="major", alpha=0.8)
ax.grid(which="minor", alpha=0.4, linestyle=":")
fig = plt.gcf()

fig.set_dpi(300)
# Create second axis

km_lims = [20, 0.5]
ax2 = ax.twiny()
ax2.set_xlabel("Geoid PDS comparison - Left")
ax2.set_xscale("log")

# Synchronize the limits of ax2 with ax
ax2.set_xlim(ax.get_xlim())

# Define the wavelengths for the ticks
ticks_wavelengths = np.array([0.6,1,2,3, 5,10,50])  # Add any desired wavelengths here

# Calculate the positions of the ticks on the second x-axis
ticks_positions = 1 / ticks_wavelengths

# Set the ticks and labels on the second x-axis
ax2.set_xticks(ticks_positions)
ax2.set_xticklabels([str(t) + " km" for t in ticks_wavelengths], fontsize=8)

plt.plot()

#%%
lon=right[0,:,:]
lat=right[1,:,:]
geoid=interpolate_data(basic,'geoid',lon, lat) # get the mean sea surface

choosen_area=geoid[380:430,:]

psdn, kn, sdn = periodogram_list(choosen_area, conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)
psdl, kl, sdl = periodogram_list(ssha_right[380:430,:], conf_lvl=68, tapering_f="tukey", sampling=250, unit="cm",skip_n=1)
# Plot
fig, ax = plt.subplots(figsize=(10,7))
ax.loglog(kn, psdn, 'r',label="Geoid")
ax.fill_between(kn, sdn[0,:], sdn[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.loglog(kl, psdl, 'b',label="Right Swath")
ax.fill_between(kl, sdl[0,:], sdl[1,:], color='r',edgecolor=None, alpha=.2,label="_nolegend_")

ax.legend()

# Setup plot settings
ax.set_yscale("log")
ax.set_xlim([50**-1,2])
ax.set_xscale("log")
ax.set_xlabel("Wavenumber (cpkm)")
ax.set_ylabel("PSD (cm$^2$/cpkm)")
ax.set_yticks(np.logspace(-6,6,13))
ax.grid(which="major", alpha=0.8)
ax.grid(which="minor", alpha=0.4, linestyle=":")
fig = plt.gcf()

fig.set_dpi(300)
# Create second axis

km_lims = [20, 0.5]
ax2 = ax.twiny()
ax2.set_xlabel("Geoid PDS comparison - Right")
ax2.set_xscale("log")

# Synchronize the limits of ax2 with ax
ax2.set_xlim(ax.get_xlim())

# Define the wavelengths for the ticks
ticks_wavelengths = np.array([0.6,1,2,3, 5,10,50])  # Add any desired wavelengths here

# Calculate the positions of the ticks on the second x-axis
ticks_positions = 1 / ticks_wavelengths

# Set the ticks and labels on the second x-axis
ax2.set_xticks(ticks_positions)
ax2.set_xticklabels([str(t) + " km" for t in ticks_wavelengths], fontsize=8)

plt.plot()
