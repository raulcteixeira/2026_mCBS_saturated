# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 08:58:34 2025
This file calculates the contrast taking into account the density modulation of the atoms 
due to the standing wave, but without absorption. The density modulation is calculated 
from a numerical integration of the atomic movement in the standing wave, and then used 
to calculate the contrast using the equation of I as a function of s, Omega_0 and Omega_M, 
with a Gaussian beam.
It also takes into account the finite size of the beam, but not the fact that the beam is tilted and 
not centered on the atoms.
It also takes into account the finite optical depth of the cloud and the absorption of the beam.

"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.lib.scimath as smath
import os

# Constants
c = 2.997e8
kB = 1.38e-23
m = 87.9 * 1.66e-27        # Atomic mass of 88Sr
lambda0 = 461e-9           # Blue line
k = 2 * np.pi / lambda0    # Wavevector
G_Sr = 2 * np.pi * 30.5e6  # Linewidth in Hz
I_sat = 40.5               # mW/cm^2

# Experimental parameters
h = 0.0042                 # Distance atom cloud <-> virtual mirror (m)
theta_0 = 4.3 * np.pi / 180  # Incidence angle in radians
theta_f = np.pi / (k * h * theta_0)  # Fringe period
sw_lambda = lambda0/np.cos(theta_0)  # Standing wave period
h = (h//sw_lambda)*sw_lambda         # Exact number of standing wave periods

N_point_theta = 300
theta = np.linspace(theta_0 - 3 * theta_f, theta_0 + 3 * theta_f, N_point_theta)

L = 0.6                  # Atom ↔ mirror distance (m)
tau_c = 2 * L / c         # Round trip travel time

# Atomic cloud size
sigma_x = 0.5e-3
sigma_y = 0.5e-3
sigma_z = 0.5e-3

# central optical depth
b0 = 0.6

# Laser beam waist
w0 = 2.1e-3
z0 = np.pi * w0 ** 2 / lambda0
d_l = 0; #1.08e-3 / 2         # Distance between center of cloud and center of beam

s_in = np.array([0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])  # saturation parameter of incoming beam at the center of the cloud
pulse_time = [30, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10] # us

#s_in = np.array([5.0, 9.0])
#pulse_time = [10, 10]

N_points_Omega = len(s_in)
Omega_0 = np.sqrt(s_in * G_Sr**2 / 2)  # Rabi frequency of the incoming beam, for each power
s_0 = 4*s_in         # saturation parameter at maximum of standing wave
Power = I_sat * s_in * np.pi * (w0 * 100)**2 / 2   #power of incident beam. mW


# Exp data
Delta_t_exp = np.array([30,15,8,10,10,15,15,10,12,15,15])  # Delta t used for each power, in us
P_exp = np.array([1,2,5,10,15,25,3,15,20,25,3])   # mW
P_err = 0.07*P_exp
s_exp = 2*P_exp/(np.pi*(w0*100)**2)/I_sat
s_err = 0.07*s_exp
Contrast_exp = np.array([0.74002,0.86048,0.61458,0.33236,0.1703,0.17106,0.58747,0.30815,0.25758,0.14104,0.69217])
Contrast_err = np.array([0.03602,0.04185,0.03549,0.02789,0.03394,0.03286,0.04734,0.02736,0.03638,0.02033,0.05181])

#%% Opening numerical results on density modulation


path = r"C:\Users\Raul\Documents\Repositories\2026_mCBS_saturated\mCBS_Atomic_Movement\data"

with open(os.path.join(path,f"atomic_movement_s={s_in[0]:.1f}.npy"), 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,2]
N_points_zsw = len(density_m_avg)

r_max = 3 * min(w0, max([sigma_x, sigma_y]))
z_max = 3 * min(w0, sigma_z)
N_points_r = 51
r = np.linspace(0, r_max, N_points_r)
N_points_z = 201
z = np.linspace(-z_max - h, z_max - h, N_points_z)
dz = z[2]-z[1]
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

Z, R, Zsw = np.meshgrid(z, r, zsw, indexing='ij')
rho = np.exp(-(R**2 + (Z + h)**2) / (2 * sigma_x**2))
R_sum = np.sum(R, axis=2)
rho_sum = np.sum(rho, axis=2)

# put the atoms at the same initial position of the standing wave.
Zsw = Zsw + np.pi/(2*k*np.cos(theta_0)) - Z % (np.pi/k/np.cos(theta_0))

rho_sum_BL = np.sum(rho, axis=0)
#OD = np.exp(-R**2 / (2 * sigma_x**2))*b0;     #Normalized for Beer-Lambert law
rho_center = rho/rho_sum_BL/dz*b0; #Density at the center of the cloud, normalized for Beer-Lambert law and for the integration over z to give b0

Intensities = np.zeros((N_points_Omega,N_point_theta))
Contrast = np.zeros(N_points_Omega)

Intensities_no_dephasing = np.zeros((N_points_Omega,N_point_theta))
Contrast_no_dephasing = np.zeros(N_points_Omega)

final_OD = np.zeros_like(s_in)

for j in range(N_points_Omega):
    
    s_0_Z = np.zeros_like(Z)
    s_0_Z[0,:,:] = s_in[j]*np.exp(-R[0,:,:]**2/w0**2)

    #Saturation variation for the incoming beam
    for i in range(len(z)-1):
        ds_0 = -rho_center[i,:,:]/(1+s_0_Z[i,:,:])*s_0_Z[i,:,:]*dz
        s_0_Z[i+1,:,:] = s_0_Z[i,:,:]+ds_0
    
    plt.figure(2)
    plt.plot(z, s_0_Z[:,int((N_points_r-1)/2),int((N_points_zsw-1)/2)])
    plt.xlabel("z (m)")
    plt.ylabel("Saturation parameter s_0")
    #plt.title(f"Saturation parameter for s={s_in[j]:.1f}")
    #plt.show()

    #Saturation variation for the reflected beam
    Z_reflected = np.flip(Z, 1)
    s_0_Z_reflected = np.zeros_like(Z)
    s_0_Z_reflected[0,:,:] = s_0_Z[-1,:,:]

    rho_center_reflected = np.flip(rho_center,0)

    for i in range(0, len(z)-1):
        ds_0_reflected = -rho_center_reflected[i,:,:]/(1+s_0_Z_reflected[i,:,:])*s_0_Z_reflected[i,:,:]*dz
        s_0_Z_reflected[i+1,:,:] = s_0_Z_reflected[i,:,:]+ds_0_reflected

    s_0_Z_reflected = np.flip(s_0_Z_reflected, 0)

    plt.figure(2)
    plt.plot(z, s_0_Z_reflected[:,int((N_points_r-1)/2),int((N_points_zsw-1)/2)])
    plt.xlabel("z (m)")
    plt.ylabel("Saturation parameter s_0")
    #plt.title(f"Saturation parameter for s={s_in[j]:.1f}")

    #Interference of incoming and reflected beam
    Omega_Z = np.sqrt(s_0_Z/2)*G_Sr
    Omega_Z_reflected = np.sqrt(s_0_Z_reflected/2)*G_Sr

    final_OD[j] = -np.log(s_0_Z_reflected[0,int((N_points_r-1)/2),int((N_points_zsw-1)/2)]/s_0_Z[0,int((N_points_r-1)/2),int((N_points_zsw-1)/2)])


final_OD = np.array(final_OD)
print(final_OD)
print(s_in)

plt.figure(5)
plt.plot(s_in, final_OD)
plt.xlabel("s_in")
plt.ylabel("Final OD")
plt.xlim(0,10)
plt.ylim(0,1.2)
plt.title("Final OD as a function of s_in")
plt.show()


data_to_save = np.column_stack((s_in, final_OD))
np.savetxt(os.path.join(path,"final_OD.csv"), data_to_save, delimiter=",")