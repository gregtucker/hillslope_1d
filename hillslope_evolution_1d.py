# -*- coding: utf-8 -*-
"""
1D model of regolith production and hillslope evolution

This program implements a 1D model of hillslope evolution. The hillslope has
a layer of mobile regolith, of thickness H, on top of a bedrock surface. There
are two fixed boundaries, one on each side, representing the elevation of a
pair of streams cutting down at a steady rate `baselevel_rate`. The governing
equations are:

    Regolith: dH/dt = (rho_r / rho_b) P0 exp( -H / Hstar ) - dq/dx
    
    Bedrock surface: dzb/dt = baselevel_rate - P0 exp( -H / Hstar )
    
    Land surface elevation: z = zb + H
    
    Regollith flux: q = -k dz/dx (1 - exp( -H / Hq ))
    
The exponential term in the regolith-flux equation reduces flux when the 
regolith is thin, going to zero flux when the regolith thickness is zero.

Written as demo for modeling class, Feb 2016

Created on Wed Feb 10 12:59:24 2016

@author: gtucker
"""

import numpy as np
import matplotlib.pyplot as plt


# INITIALIZE

# User-defined variables
H0 = 1.0             # Starting regolith thickness, m
z0 = 0.0             # Starting elevation relative to baselevel, m
num_nodes = 101      # Number of nodes
run_duration = 5.0e6 # Run duration in years
kappa = 0.01         # Transport coefficient, m2/yr
P0 = 0.0001          # Bare-bedrock regolith production rate, m/yr
Hstar = 0.6          # Soil-production decay coefficient, m
dx = 5.0             # Node spacing, m
rhos = 1600.0        # Bulk density of soil, kg/m3
rhob = 2700.0        # Bulk density of rock, kg/m3
Hq = 0.05            # Decay coefficient for soil transport, m
baselevel_rate = 0.00005  # Rate of baselevel fall, m/yr

# Run control
plot_interval = 5.0e4  # Interval for plotting, yr

# Set derived parameters
densrat = rhob / rhos   # Density ratio
next_plot = plot_interval

# Create arrays
H = np.zeros(num_nodes)  # Mobile regolith thickness, m
H += H0
z = np.zeros(num_nodes)  # Land surface elevation, m
zb = z - H               # Bedrock surface elevation, m
x = np.arange(0.0, num_nodes * dx, dx)  # x coordinates of nodes

# Set initial conditions
# (for now, we'll keep it as zero elevation)

# Set time step size and current time
dt = 0.1 * dx * dx / kappa   # Time-step size, yr
current_time = 0.0     # Current time, yr

# Plotting parameters
yaxis_max = 0.6 * (num_nodes - 1) * dx

# Plot starting condition
plt.plot(x, z, 'b')
plt.hold(True)
plt.plot(x, zb, 'r')
plt.hold(False)
plt.ylim([0, yaxis_max])
plt.xlabel('Distance (m)')
plt.ylabel('Elevation (m)')
plt.draw()


# RUN

while current_time < run_duration:
    
    # Make sure not to overshoot
    if current_time + dt > run_duration:
        dt = run_duration - current_time

    # Calculate soil production rate
    P = P0 * np.exp( -H / Hstar )
    
    # Find slope angle
    dzdx = np.diff( z ) / dx
    
    # Find soil thickness at cell edges (use average)
    Hedge = (H[:num_nodes-1] + H[1:]) / 2.0
    
    # Calculate soil flux
    q = -kappa * dzdx * (1.0 - np.exp( -Hedge / Hq ))
    
    # Calculate flux divergence
    dqdx = np.diff(q) / dx
    
    # Update soil thickness
    H[1:num_nodes-1] += (densrat * P[1:num_nodes-1] - dqdx) * dt
    
    # Update bedrock elevation
    zb[1:num_nodes-1] += (baselevel_rate - P[1:num_nodes-1]) * dt
    
    # Update elevation
    z = zb + H
    
    # Control the boundaries
    z[0] = 0.0
    z[num_nodes-1] = 0.0
    H[0] = H[1]
    H[num_nodes-1] = H[num_nodes-2]
    zb[0] = -H[0]
    zb[num_nodes-1] = -H[num_nodes-1]
    
    # Update current time
    current_time += dt
    
    # Plot
    if current_time >= next_plot:
        plt.plot(x, z, 'k')
        plt.hold(True)
        plt.plot(x, zb, 'r')
        plt.hold(False)
        plt.ylim([0, yaxis_max])
        plt.xlabel('Distance (m)')
        plt.ylabel('Elevation (m)')
        plt.draw()
        next_plot += plot_interval
    

# FINALIZE

