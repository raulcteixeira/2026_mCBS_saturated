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
c = 3e8
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

s_in_movement = np.array([0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])  # saturation parameter of incoming beam at the center of the cloud
pulse_time = [30, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10] # us
N_points_Omega_movement = len(s_in_movement)
Omega_in_movement = np.sqrt(s_in_movement / 2)*G_Sr  # Rabi frequency of the incoming beam, for each power

#s_in = np.array([0.4, 1.2, 3.0, 7.0, 9.0]) 
#pulse_time = [30, 15, 10, 10, 10] # us

# s_in = np.array([5.0, 9.0])
# pulse_time = [10, 10]

s_in = np.linspace(0.2, 10.0, 50)  # saturation parameter of incoming beam at the center of the cloud
N_points_Omega = len(s_in)
Omega_in = np.sqrt(s_in / 2)*G_Sr  # Rabi frequency of the incoming beam, for each power
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
with open(os.path.join(path,f"atomic_movement_abs_s={s_in_movement[0]:.1f}.npy"), 'rb') as f:
    density_m = np.load(f)

N_points_zsw = np.shape(density_m)[0]
#N_points_zsw = 20

r_max = 3 * min(w0, max([sigma_x, sigma_y]))
z_max = 3 * min(w0, sigma_z)
N_points_r = 51
r = np.linspace(0, r_max, N_points_r)
N_points_z = 101
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
rho_center = rho/rho_sum_BL[0,0]/dz*b0; #Density at the center of the cloud, normalized for Beer-Lambert law and for the integration over z to give b0

Intensities = np.zeros((N_points_Omega,N_point_theta))

Contrast = np.zeros(N_points_Omega)
Contrast_0 = np.zeros(N_points_Omega)
Contrast_no_dephasing = np.zeros(N_points_Omega)
Contrast_abs = np.zeros(N_points_Omega)
Contrast_abs_density = np.zeros(N_points_Omega_movement)
Contrast_elastic = np.zeros(N_points_Omega_movement)

def sat_param(Om):
    return 2 * (np.abs(Om))**2 / G_Sr**2

def g1(Om):
    Omega_M = np.sqrt(Om**2 - G_Sr**2 / 16 + 0j)
    sat = sat_param(Om)
    P_el = 1/(1+sat)
    P_in = sat/(1+sat)
    g1_el = 1
    g1_in = (1+sat)/(2*sat)*(np.exp(-G_Sr*tau_c/2)
                        + (sat-1)/(sat+1)*np.cos(Omega_M*tau_c)*np.exp(-3*G_Sr*tau_c/4)
                        + G_Sr/(4*Omega_M)*(5*sat-1)/(sat+1)*np.sin(Omega_M*tau_c)*np.exp(-3*G_Sr*tau_c/4))
    return P_el*g1_el + P_in*g1_in

def ss(s, g1, density,theta_i):
    return np.real(np.sum(s/(1+s)*(1 + np.real(g1*np.exp(-2*1j*k*(Z+Zsw)*np.cos(theta_i))))*density,axis=2))


## Calculation of the theoretical curve with atomic movement
for j in range(N_points_Omega_movement):

    # with open(os.path.join(path,f"atomic_movement_s={s_in[j]:.1f}.npy"), 'rb') as f:
    with open(os.path.join(path,f"atomic_movement_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
        density_m = np.load(f)
    with open(os.path.join(path,f"atomic_movement_abs_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
        density_m_abs = np.load(f)
    with open(os.path.join(path,f"atomic_movement_abs_2_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
        density_m_abs_2 = np.load(f)

    density_m = density_m / np.mean(density_m, axis=0)  # Normalize each column
    density_m_abs = density_m_abs / np.mean(density_m_abs, axis=0)  # Normalize each column
    density_m_abs_2 = density_m_abs_2 / np.mean(density_m_abs_2, axis=0)  # Normalize each column

    # For including all times, we verage over the first 20 time steps, which correspond to 10 us.
    density_m_t = density_m[:,0:20] + density_m_abs[:,0:20] + density_m_abs_2[:,0:20]
    density_m_avg = np.mean(density_m_t,axis=1)

    # If we consider the highest impact of movement, we take fixed t = 1 us.
    # density_m_avg = density_m[:,2] + density_m_abs[:,2] + density_m_abs_2[:,2]
    

    plt.figure(1)
    plt.plot(zsw, density_m_avg)
    plt.xlabel("zsw (m)")
    plt.ylabel("Density modulation")
    plt.title(f"Density modulation for s={s_in[j]:.1f}")
    #plt.show()

    #Omega versus R - without absortion, the Rabi frequency only depends on the local intensity of the beam, which has a Gaussian profile.
    Omega = Omega_in_movement[j]*np.exp(-R**2/w0**2)

    #Omega versus R and Z - to account for absorption, we need to calculate the saturation parameter
    # at each position in the cloud, which depends on the local intensity of the beam.
    # The local intensity depends on the local Rabi frequency, which depends on the local saturation parameter. We can calculate the local saturation parameter by integrating the absorption of the beam as it propagates through the cloud, taking into account the density modulation of the atoms.
    s_0_Z = np.zeros_like(Z)
    s_0_Z[0,:,:] = s_in_movement[j]*np.exp(-2*R[0,:,:]**2/w0**2)

    #Saturation variation for the incoming beam
    for i in range(len(z)-1):
        ds_0 = -rho_center[i,:,:]/(1+s_0_Z[i,:,:])*s_0_Z[i,:,:]*dz
        s_0_Z[i+1,:,:] = s_0_Z[i,:,:]+ds_0

    #Saturation variation for the reflected beam
    Z_reflected = np.flip(Z, 0)
    s_0_Z_reflected = np.zeros_like(Z)
    s_0_Z_reflected[0,:,:] = s_0_Z[-1,:,:]

    rho_center_reflected = np.flip(rho_center,0)

    for i in range(0, len(z)-1):
        ds_0_reflected = -rho_center_reflected[i,:,:]/(1+s_0_Z_reflected[i,:,:])*s_0_Z_reflected[i,:,:]*dz
        s_0_Z_reflected[i+1,:,:] = s_0_Z_reflected[i,:,:]+ds_0_reflected

    #Interference of incoming and reflected beam
    Omega_Z = np.sqrt(s_0_Z/2)*G_Sr
    Omega_Z_reflected = np.flip(np.sqrt(s_0_Z_reflected/2)*G_Sr, 0)


    Omega_0_total = 2*Omega_in_movement[j]*np.cos(k*(Z+Zsw)*np.cos(theta_0))
    Omega_total = 2*Omega*np.cos(k*(Z+Zsw)*np.cos(theta_0))
    Omega_total_abs = np.abs(Omega_Z*np.exp(1j*k*(Z+Zsw)*np.cos(theta_0)) + Omega_Z_reflected*np.exp(-1j*k*(Z+Zsw)*np.cos(theta_0)))

    s = sat_param(Omega_total)
    g1_total = g1(Omega_total)

    s_abs = sat_param(Omega_total_abs)
    g1_total_abs = g1(Omega_total_abs)
    g1_total_elastic = 1

    print(s_in_movement[j])

    not1, not2, density_mesh = np.meshgrid(np.ones(len(z)), np.ones(len(r)), density_m_avg, indexing='ij')

    I_abs_density = np.zeros(N_point_theta)
    I_elastic = np.zeros(N_point_theta)

    Omega_M = np.sqrt(Omega_total**2 - G_Sr**2 / 16 + 0j)
    for i in range(N_point_theta):
        ss1 = s/(4.*(1.+s))*(2./(1.+s) + np.exp(-G_Sr*tau_c/2.) + (s-1.)/(s+1.)*np.cos(Omega_M*tau_c)*np.exp(-3.*G_Sr*tau_c/4.) + G_Sr/(4.*Omega_M)*(5.*s-1.)/(s+1.)*np.sin(Omega_M*tau_c)*np.exp(-3.*G_Sr*tau_c/4.))
        ss_avg_abs_density = ss(s_abs, g1_total_abs, density_mesh, theta[i])
        ss_avg_elastic = ss(s_abs, g1_total_elastic, density_mesh, theta[i])
        I_abs_density[i] = np.sum(np.sum(R_sum * rho_sum * ss_avg_abs_density))
        I_elastic[i] = np.sum(np.sum(R_sum * rho_sum * ss_avg_elastic))

    Contrast_abs_density[j] = (np.max(I_abs_density) - np.min(I_abs_density)) / (np.max(I_abs_density) + np.min(I_abs_density)) * 2
    Contrast_elastic[j] = (np.max(I_elastic) - np.min(I_elastic)) / (np.max(I_elastic) + np.min(I_elastic)) * 2


with open("Contrast_movement.txt", "w") as fid:
    fid.write(f"s_in_movement, Contrast_gaussian_beam_abs_density_mod,  Contrast_gaussian_beam_abs_density_mod_g1_1\n")
    for n in range(N_points_Omega_movement):
        fid.write(f"{s_in_movement[n]:.3f},  {Contrast_abs_density[n]:.3f},  {Contrast_elastic[n]:.3f}\n")



## Calculation of the theoretical curve for other cases
for j in range(N_points_Omega):


    #Omega versus R - without absortion, the Rabi frequency only depends on the local intensity of the beam, which has a Gaussian profile.
    Omega = Omega_in[j]*np.exp(-R**2/w0**2)

    #Omega versus R and Z - to account for absorption, we need to calculate the saturation parameter
    # at each position in the cloud, which depends on the local intensity of the beam.
    # The local intensity depends on the local Rabi frequency, which depends on the local saturation parameter. We can calculate the local saturation parameter by integrating the absorption of the beam as it propagates through the cloud, taking into account the density modulation of the atoms.
    s_0_Z = np.zeros_like(Z)
    s_0_Z[0,:,:] = s_in[j]*np.exp(-2*R[0,:,:]**2/w0**2)

    #Saturation variation for the incoming beam
    for i in range(len(z)-1):
        ds_0 = -rho_center[i,:,:]/(1+s_0_Z[i,:,:])*s_0_Z[i,:,:]*dz
        s_0_Z[i+1,:,:] = s_0_Z[i,:,:]+ds_0

    plt.figure(2)
    plt.plot(z, s_0_Z[:,0,int((N_points_zsw-1)/2)])
    plt.xlabel("z (m)")
    plt.ylabel("Saturation parameter s_0")
    plt.title(f"Saturation parameter versus z")
    #plt.show()

    #Saturation variation for the reflected beam
    Z_reflected = np.flip(Z, 0)
    s_0_Z_reflected = np.zeros_like(Z)
    s_0_Z_reflected[0,:,:] = s_0_Z[-1,:,:]

    rho_center_reflected = np.flip(rho_center,0)

    for i in range(0, len(z)-1):
        ds_0_reflected = -rho_center_reflected[i,:,:]/(1+s_0_Z_reflected[i,:,:])*s_0_Z_reflected[i,:,:]*dz
        s_0_Z_reflected[i+1,:,:] = s_0_Z_reflected[i,:,:]+ds_0_reflected

    #Interference of incoming and reflected beam
    Omega_Z = np.sqrt(s_0_Z/2)*G_Sr
    Omega_Z_reflected = np.flip(np.sqrt(s_0_Z_reflected/2)*G_Sr, 0)


    Omega_0_total = 2*Omega_in[j]*np.cos(k*(Z+Zsw)*np.cos(theta_0))
    Omega_total = 2*Omega*np.cos(k*(Z+Zsw)*np.cos(theta_0))
    Omega_total_abs = np.abs(Omega_Z*np.exp(1j*k*(Z+Zsw)*np.cos(theta_0)) + Omega_Z_reflected*np.exp(-1j*k*(Z+Zsw)*np.cos(theta_0)))

    s_0 = sat_param(Omega_0_total)
    g1_total_0 = g1(Omega_0_total)

    s = sat_param(Omega_total)
    g1_total = g1(Omega_total)

    s_abs = sat_param(Omega_total_abs)
    g1_total_abs = g1(Omega_total_abs)
    g1_total_elastic = 1

    print(s_in[j])

    plt.figure(10)
    plt.plot(Omega_0_total[int((N_points_r-1)/2),int((N_points_z-1)/2),:])
    plt.plot(Omega_total[int((N_points_r-1)/2),int((N_points_z-1)/2),:])
    plt.plot(Omega_total_abs[int((N_points_r-1)/2),int((N_points_z-1)/2),:])
    #plt.show()

    I_0 = np.zeros(N_point_theta)
    I = np.zeros(N_point_theta)
    I_abs = np.zeros(N_point_theta)

    Omega_M = np.sqrt(Omega_total**2 - G_Sr**2 / 16 + 0j)
    for i in range(N_point_theta):
        ss_avg_0 = ss(s_0, g1_total_0, np.ones_like(Z), theta[i])
        ss1 = s/(4.*(1.+s))*(2./(1.+s) + np.exp(-G_Sr*tau_c/2.) + (s-1.)/(s+1.)*np.cos(Omega_M*tau_c)*np.exp(-3.*G_Sr*tau_c/4.) + G_Sr/(4.*Omega_M)*(5.*s-1.)/(s+1.)*np.sin(Omega_M*tau_c)*np.exp(-3.*G_Sr*tau_c/4.))
        ss_avg = ss(s, g1_total, np.ones_like(Z), theta[i])
        ss_avg = np.sum(s/(1.+s) + 2.*np.real(np.exp(-1j*2.*k*(Z+Zsw)*np.cos(theta[i]))*ss1),axis=2)
        ss_avg_abs = ss(s_abs, g1_total_abs, np.ones_like(Z), theta[i])
        I_0[i] = np.sum(np.sum(R_sum * rho_sum * ss_avg_0))  # Intensity as a function of theta (fringes);
        I[i] = np.sum(np.sum(R_sum * rho_sum * ss_avg))
        I_abs[i] = np.sum(np.sum(R_sum * rho_sum * ss_avg_abs))


    plt.figure(20)
    plt.plot(theta - theta_0, I / np.mean(I))
    Intensities[j, :] = I
    Contrast_0[j] = (np.max(I_0) - np.min(I_0)) / (np.max(I_0) + np.min(I_0)) * 2
    Contrast[j] = (np.max(I) - np.min(I)) / (np.max(I) + np.min(I)) * 2
    Contrast_abs[j] = (np.max(I_abs) - np.min(I_abs)) / (np.max(I_abs) + np.min(I_abs)) * 2


with open("Fringes_centered_beam.txt", "w") as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.3f}  ")
        fid.write("\n")

with open("Contrast.txt", "w") as fid:
    fid.write(f"s_in,  Contrast_plane_wave,  Contrast_gaussian_beam,  Contrast_gaussian_beam_abs \n")
    for n in range(N_points_Omega):
        fid.write(f"{s_in[n]:.3f},  {Contrast_0[n]:.3f},  {Contrast[n]:.3f},  {Contrast_abs[n]:.3f}\n")

plt.figure(30)
plt.plot(s_in, Contrast_0, label='Plane wave')
plt.plot(s_in, Contrast, label='Gaussian beam')
plt.plot(s_in, Contrast_abs, label='Gaussian beam with absorption')
plt.plot(s_in_movement, Contrast_abs_density, label='Gaussian beam with absorption and density modulation')
plt.plot(s_in_movement, Contrast_elastic, label='Gaussian beam with absorption and density modulation, g1 = 1')
plt.errorbar(s_exp, Contrast_exp, xerr=s_err, yerr=Contrast_err, fmt='o', label='Experimental')
plt.xlabel("s")
plt.ylabel("Contrast")
plt.title('Contrast as a function of saturation parameter s')
plt.legend()
plt.show()