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

from Contraste_waist_Complet_Gaussian_4 import N_per_avg

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
print(f'Rabi frequency for Ω₀ τ_c = 1: {1 / tau_c / 2 / np.pi / 1e6:.2f} MHz')

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

s_in = [0.4, 0.8, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]  # saturation parameter of incoming beam at the center of the cloud
N_points_Omega = len(s_in)
Omega_0 = np.sqrt(s_in * G_Sr**2 / 2)  # Rabi frequency of the incoming beam, for each power
s_0 = 4*s_in         # saturation parameter at maximum of standing wave
Power = I_sat * s_in * np.pi * (w0 * 100)**2 / 2   #power of incident beam. mW
pulse_time = [30, 15, 15, 10, 10, 10, 10, 10, 10, 10, 10] # us

# Exp data
Delta_t_exp = np.array([30,15,8,10,10,15,15,10,12,15,15])  # Delta t used for each power, in us
P_exp = np.array([1,2,5,10,15,25,3,15,20,25,3])   # mW
P_err = 0.07*P_exp
s_exp = 2*P_exp/(np.pi*(w0*100)**2)/I_sat
s_err = 0.07*s_exp
Contrast_exp = np.array([0.74002,0.86048,0.61458,0.33236,0.1703,0.17106,0.58747,0.30815,0.25758,0.14104,0.69217])
Contrast_err = np.array([0.03602,0.04185,0.03549,0.02789,0.03394,0.03286,0.04734,0.02736,0.03638,0.02033,0.05181])

#%% Opening numerical results on density modulation


path = "G:\\Meu Drive\\Documentos\\Academics\\Sr 1\\projetos\\2020 - saturated mCBS\\Scripts\\Results_density_modulation\\"

with open(path+f"atomic_movement_s_={s_in[0]:.1f}.npy", 'rb') as f:
    density_m = np.load(f)

density_m_avg = np.mean(density_m,axis=1)
N_points_zsw = len(density_m_avg)

r_max = 3 * min(w0, max([sigma_x, sigma_y]))
z_max = 3 * min(w0, sigma_z)
N_points_r = 51
r = np.linspace(0, r_max, N_points_r)
N_points_z = 101
z = np.linspace(-z_max - h, z_max - h, N_points_z)
dz = z[2]-z[1]
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

R, Z, Zsw = np.meshgrid(r, z, zsw, indexing='ij')
rho = np.exp(-(R**2 + (Z + h)**2) / (2 * sigma_x**2))
R_sum = np.sum(R, axis=2)
rho_sum = np.sum(rho, axis=2)

OD = np.exp(-R**2 / (2 * sigma_x**2))*b0;     #Normalized for Beer-Lambert law

Intensities = np.zeros (N_points_Omega,N_point_theta)
Contrast = np.zeros(N_points_Omega)

Intensities_no_dephasing = np.zeros (N_points_Omega,N_point_theta)
Contrast_no_dephasing = np.zeros(N_points_Omega)

for j in range(N_points_Omega):
    
    with open(path+f"atomic_movement_s_={s_in[j]:.1f}.npy", 'rb') as f:
        density_m = np.load(f)

    density_m_avg = np.mean(density_m,axis=1)
    #density_m_avg = density_m[:,instant_mod]
    N_points_zsw = len(density_m_avg)

    #Omega versus Z - to account for absorption, we need to calculate the saturation parameter 
    # at each position in the cloud, which depends on the local intensity of the beam. 
    # The local intensity depends on the local Rabi frequency, which depends on the local saturation parameter. We can calculate the local saturation parameter by integrating the absorption of the beam as it propagates through the cloud, taking into account the density modulation of the atoms.
    s_0_Z = np.zeros(np.size(Z));
    s_0_Z[:,0,:] = 2*(Omega_0[j]/G_Sr)**2.*np.exp(-R[:,0,:]**2/w0**2)

    #Saturation variation for the incoming beam
    for i in range(1, len(z)-1):
        ds_0 = -OD[:,i,:]/(1+s_0_Z[:,i,:])*s_0_Z[:,i,:]*dz
        s_0_Z[:,i+1,:] = s_0_Z[:,i,:]+ds_0
    
    #Saturation variation for the reflected beam
    Z_reflected = np.flip(Z, 1)
    s_0_Z_reflected = np.zeros_like(Z)
    s_0_Z_reflected[:,0,:] = s_0_Z[:,-1,:]

    OD_reflected = np.flip(OD, 1)

    for i in range(1, len(z)-1):
        ds_0_reflected = -OD_reflected[:,i,:]/(1+s_0_Z_reflected[:,i,:])*s_0_Z_reflected[:,i,:]*dz
        s_0_Z_reflected[:,i+1,:] = s_0_Z_reflected[:,i,:]+ds_0_reflected

    #Interference of incoming and reflected beam
    Omega_Z = np.sqrt(s_0_Z/2)*G_Sr
    Omega_Z_reflected = np.fliplip(np.sqrt(s_0_Z_reflected/2)*G_Sr, 1)

    Omega_total = np.abs(Omega_Z*np.exp(np.sqrt(-1)*k*(Z+Zsw)*np.cos(theta_0)) + Omega_Z_reflected*np.exp(-np.sqrt(-1)*k*(Z+Zsw)*np.cos(theta_0)))
    

    Omega = Omega_total
    Omega_M = np.sqrt(Omega**2 - G_Sr**2 / 16 + 0.j)
    s = 2 * (np.abs(Omega))**2 / G_Sr**2

    P_el = 1/(1+s)
    P_in = s/(1+s)
    g1_el = 1
    g1_in = 1
    g1_total = P_el*g1_el + P_in*g1_in
    g1_no_dephasing = 1


    plt.figure(10)
    plt.plot(Omega[int((N_points_r-1)/2),int((N_points_z-1)/2),:])
    

    not1, not2, density_mesh = np.meshgrid(np.ones(len(r)), np.ones(len(z)), density_m_avg, indexing='ij')
    
    I = np.zeros(N_point_theta)
    I_no_dephasing = np.zeros(N_point_theta)

    for i in range(N_point_theta):
        ss_avg = np.sum(s/(1+s)*(1 + g1_total*np.cos(2*k*(Z+Zsw)*np.cos(theta[i])))*density_mesh,axis=2);
        I[i] = np.sum(R_sum * rho_sum * ss_avg);  # Intensity as a function of theta (fringes);   
        ss_avg_no_dephasing = np.sum(s/(1+s)*(1 + g1_no_dephasing*np.cos(2*k*(Z+Zsw)*np.cos(theta[i])))*density_mesh,axis=2);
        I_no_dephasing[i] = np.sum(R_sum * rho_sum * ss_avg_no_dephasing);  # Intensity as a function of theta (fringes);


    plt.figure(20)
    plt.plot(theta - theta_0, I / np.mean(I))
    Intensities[j, :] = I
    Contrast[j] = (np.max(I) - np.min(I)) / (np.max(I) + np.min(I)) * 2
    Intensities_no_dephasing[j, :] = I_no_dephasing
    Contrast_no_dephasing[j] = (np.max(I_no_dephasing) -    np.min(I_no_dephasing)) / (np.max(I_no_dephasing) + np.min(I_no_dephasing)) * 2

with open("Fringes_centered_beam.txt", "w") as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.1f}  ")
        fid.write("\n")



#plt.figure(1)
#plt.plot(Omega_0 * tau_c, Contrast, 'g')
#plt.plot(Omega_0 * tau_c, Contrast_no_dephasing, 'm')
#plt.xlabel("Ω₀ τ_c")
#plt.ylabel("Contrast")
#plt.legend(["Gaussian beam, centered on the atoms"])

#plt.figure(2)
#plt.plot(Power, Contrast, 'g')
#plt.errorbar(P_exp, Contrast_exp, xerr=P_err,yerr=Contrast_err,fmt='o')
#plt.xlabel("Incident laser power (mW)")
#plt.ylabel("Contrast")
#plt.title('density modulation at 2us, same power beams')

plt.figure(3)
plt.plot(s_in, Contrast, 'g', label='With dephasing')
plt.plot(s_in, Contrast_no_dephasing, 'm', label='Without dephasing')
plt.errorbar(s_exp, Contrast_exp, xerr=s_err, yerr=Contrast_err, fmt='o', label='Experimental')
plt.xlabel("s")
plt.ylabel("C")
plt.title('density modulation at 2us, same power beams')
plt.legend()