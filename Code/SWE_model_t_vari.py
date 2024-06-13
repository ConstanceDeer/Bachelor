
#%% initilization for plots

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def centeredDistanceMatrix(n):
    # Make sure n is odd
    x, y = np.meshgrid(range(n), range(n))
    return np.sqrt((x - (n / 2) + 1)**2 + (y - (n / 2) + 1)**2)

def arbitraryfunction(d, y, n):
    # Interpolate the values of y for given distances
    x = np.linspace(0, np.max(d), len(y))
    f = interp1d(x, y, fill_value="extrapolate")
    return f(d.flatten()).reshape(d.shape)

#%% Zeta model Vesion 1 - STOR
import numpy as np
from scipy.special import j0
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
rho_w = 1024  # define your density of water
rho_i = 916   # define your density of ice
g = 9.8       # gravity acceleration
h0 = 5        # define your h0_d value here
d = 3000      # water depth
a = 10
b = 5 # define your b_a value here
r = 100
alpha = 1.5
vi = np.sqrt(2 * g * (h0 - b / 2))
va = vi / (1 + (2 / 3) * (rho_w / rho_i) * (a / b))
m=np.pi*rho_i*(a**2)*b

volume=2*np.pi*(b/2)**2 * a*2
mass=volume*rho_i

print(volume,mass)
#%
#%
def B_integrant(x, y):
    return np.sqrt(1 - y**2) * j0((a/d) * x * y) * y

def B(x):
    return quad(lambda y: B_integrant(x, y), 0, 1)[0]

def I_zeta_integrant(x, t):
    return j0((r/d)*x)*B(x)*np.sqrt(x*np.tanh(x))*np.sin(np.sqrt(((g*(t**2))/d)*x*np.tanh(x)))*x

def I_zeta(t):
    return quad(lambda x: I_zeta_integrant(x, t), 0, 1500)[0] 

def zeta(t):
    return (-2*(a/d)**3 * (np.sqrt(2*(h0/d)-(a/d)*(b/a))/(1+(2/3)*(rho_w/rho_i)*(a/b)) ) * I_zeta(t))*d
#%
T = np.linspace(0, 100, 200)
LL=np.zeros_like(T)
for i in range(len(T)):
    LL[i]=zeta(T[i])
#%
fig, axs = plt.subplots(figsize=(10, 5))
plt.plot(T,LL)
axs.set_xlim([0,100])
axs.set_ylabel('Sea Surface height, η [m]',fontsize=14)
axs.set_xlabel('Time from impact, t [s]',fontsize=14)

# Increase tick label size
axs.tick_params(axis='both', which='major', labelsize=14)
axs.tick_params(axis='both', which='minor', labelsize=14)  # If you have minor ticks
axs.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
fig.set_dpi(300)

# Create an interpolation function
interp_func_large = interp1d(T, LL, kind='linear')

#%% Zeta model Vesion 2 - LILLE

import numpy as np
from scipy.special import j0
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
rho_w = 1024  # define your density of water kg/m^3
rho_i = 916   # define your density of ice   kg/m^3
g = 9.8       # gravity acceleration
h0 = 1        # define your h0_d value here
d = 3000
a = 2
b = 1 # define your b_a value here
r = 25
alpha = 1.5
vi = np.sqrt(2 * g * (h0 - b / 2))
va = vi / (1 + (2 / 3) * (rho_w / rho_i) * (a / b))
m=np.pi*rho_i*(a**2)*b

volume=2*np.pi*(b/2)**2 * a*2
mass=volume*rho_i

print(volume,mass)
#%
#%
def B_integrant(x, y):
    return np.sqrt(1 - y**2) * j0((a/d) * x * y) * y

def B(x):
    return quad(lambda y: B_integrant(x, y), 0, 1)[0]

def I_zeta_integrant(x, t):
    return j0((r/d)*x)*B(x)*np.sqrt(x*np.tanh(x))*np.sin(np.sqrt(((g*(t**2))/d)*x*np.tanh(x)))*x

def I_zeta(t):
    return quad(lambda x: I_zeta_integrant(x, t), 0, 3500)[0] #3000

def zeta(t):
    return (-2*(a/d)**3 * (np.sqrt(2*(h0/d)-(a/d)*(b/a))/(1+(2/3)*(rho_w/rho_i)*(a/b)) ) * I_zeta(t))*d
#%
T = np.linspace(0, 105, 400)
LL=np.zeros_like(T)
for i in range(len(T)):
    LL[i]=zeta(T[i])
    
fig, axs = plt.subplots(figsize=(10, 5))
plt.plot(T,LL)
axs.set_xlim([0,100])
axs.set_ylabel('Sea Surface height, η [m]',fontsize=14)
axs.set_xlabel('Time from impact, t [s]',fontsize=14)

# Increase tick label size
axs.tick_params(axis='both', which='major', labelsize=14)
axs.tick_params(axis='both', which='minor', labelsize=14)  # If you have minor ticks
axs.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
fig.set_dpi(300)

# Create an interpolation function
interp_func_small = interp1d(T, LL, kind='linear')


#%% 1D model, med dynamisk model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# The following code is build on the code: https://github.com/hongkangcarl/PDE-solver-shallow-water-equations/blob/main/swe.py

# finite difference treatment of grad h = [h(x+âx) - h(x-âx)]/(2âx)
def grad(h):
    grad_h=np.gradient(h, dx)
    return  grad_h

# finite difference treatment of div v
def div(vec):
    return np.gradient(vec, dx)

# parameter values
# 20.45 kilometers , amplitude around 0.05m = 5cm
L_x = 22500/2           # Length of domain in x-direction
g = 9.8              # Acceleration of gravity [m/s^2]
H = 3000              # Depth of fluid [m]
N_x = 1000                         # Number of grid points in x-direction
dx = L_x/(N_x - 1)                 # Grid spacing in x-direction
dt = 0.1*dx/np.sqrt(g*H)       # Time step (defined from the CFL condition)
BB = 1/12         # friction coefficient [ts^-1]
xs = np.linspace(0, L_x, N_x)  # Array with x-points
N=75
n=round(N/dt)
ts = np.linspace(0, N, n)


#% Køre selve model 
x = xs

# prepare the list to store the numerical solutions
h_list = []
v_list = []

# initial condition: no initial surface distrotion, no initial velocities
h_init = np.zeros_like(x)
v_init = np.zeros_like(h_init)

# Update initial condition at x=0 with given SSH
#h_init[0] = time_step[i]

# store the initial condition into the lists
h_list.append(h_init)
v_list.append(v_init)

# start the simulation, for loop over the discretized time
for i in range(1, len(ts)):
    # calculate the velocities based on the current h(x, t)
    v = (1 - dt * BB) * v_list[-1] - dt * g * grad(h_list[-1])
    
    # calculate the next-step h(x, t+ât) based on the velocities
    h = h_list[-1] - dt * H * div(v)
    
    # Define the length of the damping zone (sponge layer)
    damping_length = 10  # Adjust as needed

    # Define the damping factor
    damping_factor = 0.1  # Adjust as needed

    if ts[i] < 100:
        h[0] = interp_func_large(ts[i])
        for j in range(damping_length):
            h[-j-1] *= (1 - damping_factor * (damping_length - j) / damping_length)
        
    else:
        h[0] = 0
        for j in range(damping_length):
            h[j] *= (1 - damping_factor * (damping_length - j) / damping_length)
            h[-j-1] *= (1 - damping_factor * (damping_length - j) / damping_length)

   
    # store the numerical results into the lists
    v_list.append(v)
    h_list.append(h)

#%% Gif - plot generator 
import imageio
image_paths = []
step = 25
for i in range(1, len(h_list), step): 
    data = h_list[i]

    # Create a circular distance matrix
    n = len(data)
    d = centeredDistanceMatrix(n)

    # Interpolate the data for circular plot
    f = arbitraryfunction(d, data, n)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(16, 6))

    # Plot the first plot
    axs[0].plot(x+r, data, label=f"t={ts[i]}")
    axs[0].legend()
    axs[0].set_xlabel('Distrance from impact, [m]')
    axs[0].set_ylabel('Surface elevation, [m]')
    #axs[0].set_title('Spacial Wave model in 1-Dimension')
    axs[0].set_ylim([-1, 1])
    axs[0].set_xlim([0, 20500/2])
    axs[0].grid(True)

    # Plot the second plot
    img = axs[1].imshow(f.T, origin='lower', interpolation='nearest', cmap='Spectral_r', extent=[-20500/2, 20500/2, -20500/2, 20500/2])
    cbar = plt.colorbar(img, ax=axs[1], label='Surface elevation, [m]')
    axs[1].set_xlabel('Distrance from impact, [m]')
    img.set_clim(-1, 1)  # Set colorbar limits
    #axs[1].set_title('Axially Symmetric Data')
    axs[1].set_aspect('equal', adjustable='box')  # Ensure aspect ratio is equal for circular plot
    fig.set_dpi(300)
    # Save the figure as an image
    image_path = f"animation_big/plot_{i}.png"
    plt.savefig(image_path)
    #plt.show()
    plt.close() 
    image_paths.append(image_path)

#%% gif maker
with imageio.get_writer('Animation_large.gif', mode='I', duration=dt*step) as writer:
    for image_path in image_paths:
        image = imageio.imread(image_path)
        writer.append_data(image)
#%% 20'th loop 
step = 10
for i in range(0+50, len(h_list), step):
    plt.figure(figsize=(8, 6))
    plt.plot(x+50, h_list[i], label=f"t={ts[i]}")
    plt.xlabel('x')
    plt.ylabel('h')
    plt.title('Shallow Water Equation Simulation (1D)')
    plt.ylim([-3, 3])
    plt.xlim([0, 130000])
    plt.legend()
    plt.grid(True)
    plt.show()
    
#%% other plot
step = 20
for i in range(0, len(h_list), step):
    # Assuming h_list[10022] contains your data
    data = h_list[i]
    
    # Create a circular distance matrix
    n = len(data)
    d = centeredDistanceMatrix(n)
    
    # Interpolate the data for circular plot
    f = arbitraryfunction(d, data, n)
    
    # Plot the data
    plt.figure(figsize=(8, 8))
    plt.imshow(f.T, origin='lower', interpolation='nearest', cmap='seismic', extent=[-6000, 6000, -6000, 6000])
    plt.colorbar(label='h')
    plt.clim(-1, 1)  # Set colorbar limits
    plt.title('Axially Symmetric Data')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is equal for circular plot
    plt.show()

#%% Plot for specific point 
#% 1D specific points
index=10022-500
data = h_list[index]
plt.figure(figsize=(10, 3))
plt.plot(x+100, data)
plt.xlabel('Distrance from impact, [m]',fontsize=12)
plt.ylabel('Surface elevation, [m]',fontsize=12)
plt.text(8000, 0.85, f"t = {ts[index]}s", horizontalalignment='left', fontsize=9, color='black')
plt.ylim([-1, 1])
plt.xlim([0, 20500/2])
fig = plt.gcf()
fig.set_dpi(600)
plt.grid(True,alpha=0.7)
plt.show()


#% 2D For one of the plot 
# Assuming h_list[10022] contains your data
data = h_list[index]

# Create a circular distance matrix
n = len(data)
d = centeredDistanceMatrix(n)

# Interpolate the data for circular plot
f = arbitraryfunction(d, data, n)

# Plot the data
plt.figure(figsize=(9, 7))
plt.imshow(f.T, origin='lower', interpolation='nearest', cmap='seismic', extent=[-22500/2, 22500/2, -22500/2, 22500/2])
plt.colorbar(label='Sea surface height, η [m]')
plt.clim(-1, 1)  # Set colorbar limits
plt.xlabel('Distance From Impact Site, r [m]',fontsize=12)
# Increase tick label size
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='minor', labelsize=12)  # If you have minor ticks
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
plt.text(4000, 11500, f"t = {ts[index]}s", horizontalalignment='left', fontsize=9, color='black')
fig = plt.gcf()
fig.set_dpi(600)
#plt.title('Axially Symmetric Data')
plt.gca().set_aspect('equal', adjustable='box')  # Ensure aspect ratio is equal for circular plot
plt.show()
#%% Subplot
import matplotlib.pyplot as plt

# Assuming h_list[10022] contains your data
data = h_list[10022]

# Create a circular distance matrix
n = len(data)
d = centeredDistanceMatrix(n)

# Interpolate the data for circular plot
f = arbitraryfunction(d, data, n)

# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(16, 6))

# Plot the first plot
axs[0].plot(x+50, data, label=f"t={ts[i]}")
axs[0].legend()
axs[0].set_xlabel('Distrance from impact, [m]')
axs[0].set_ylabel('Surface elevation, [m]')
axs[0].set_title('Spacial Wave model in 1-Dimension')
axs[0].set_ylim([-1.5, 1.5])
axs[0].set_xlim([0, 10000])
axs[0].grid(True)

# Plot the second plot
img = axs[1].imshow(f.T, origin='lower', interpolation='nearest', cmap='seismic', extent=[-10000, 10000, -10000, 10000])
cbar = plt.colorbar(img, ax=axs[1], label='Surface elevation, [m]')
img.set_clim(-0.3, 0.3)  # Set colorbar limits
axs[1].set_title('Axially Symmetric Data')
axs[1].set_aspect('equal', adjustable='box')  # Ensure aspect ratio is equal for circular plot

plt.show()
#%% Round plot 
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Sample data
X = x
y =  h_list[10022]


# Smoothing the data
f = interp1d(X, y, kind='cubic')
x_smooth = np.linspace(X.min(), X.max(), 1000)
y_smooth = f(x_smooth)

# Creating circular data
theta = np.linspace(0, 2*np.pi, 100)  # Angle
r = np.linspace(X.min(), X.max(), 100)  # Radius adjusted to match the range of x

# Create meshgrid
theta, r = np.meshgrid(theta, r)

# Convert polar coordinates to Cartesian coordinates
x_circular = r * np.cos(theta)
y_circular = r * np.sin(theta)

# Calculate elevation values for each point on the grid
z_circular = f(r.flatten())

# Reshape z to match the shape of x and y
z_circular = z_circular.reshape(x_circular.shape)

# Plot the circular disk with colormap representing elevation using pcolormesh
plt.figure(figsize=(10, 8))
ax = plt.gca()  # Get the current Axes object
ax.set_facecolor('white')  # Set background color to light blue
plt.pcolormesh(x_circular, y_circular, z_circular, cmap='seismic', shading='auto')
plt.colorbar(label='Elevation')
plt.ylim([-6050,6050])
plt.xlim([-6050,6050])
plt.title('Circular Disk Representation of Elevation')
#plt.axis('equal')  # Equal aspect ratio
plt.show()
