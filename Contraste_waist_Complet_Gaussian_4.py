# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 08:58:34 2025

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.lib.scimath as smath

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

# Laser beam waist
w0 = 2.1e-3
z0 = np.pi * w0 ** 2 / lambda0
d_l = 1.08e-3 / 2         # Distance between center of cloud and center of beam

N_points_Omega = 50
Omega_0 = np.linspace(0.2, 2.0, N_points_Omega) / tau_c
s_0 = 8 * Omega_0**2 / G_Sr**2          # saturation parameter at maximum of standing wave
Power = I_sat * s_0 / 4 * np.pi * (w0 * 100)**2 / 2   #power of incident beam. mW

# Exp data
Delta_t_exp = np.array([30,15,8,10,10,15,15,10,12,15,15])  # Delta t used for each power, in us
P_exp = np.array([1,2,5,10,15,25,3,15,20,25,3])   # mW
P_err = 0.07*P_exp
s_exp = 2*P_exp/(np.pi*(w0*100)**2)/I_sat
s_err = 0.07*s_exp
Contrast_exp = np.array([0.74002,0.86048,0.61458,0.33236,0.1703,0.17106,0.58747,0.30815,0.25758,0.14104,0.69217])
Contrast_err = np.array([0.03602,0.04185,0.03549,0.02789,0.03394,0.03286,0.04734,0.02736,0.03638,0.02033,0.05181])

N_points_zsw = 50
N_per_avg = int(100/N_points_zsw)
print(N_per_avg)

s_dens_vec = np.array([0.2,1.0,3.0,5.0,9.0])
instant_mod = 5
dens_mod = np.zeros([N_points_zsw,len(s_dens_vec)])

#%% Opening numerical results on density modulation



path = "G:\\Meu Drive\\Documentos\\Academics\\Sr 1\\projetos\\2020 - saturated mCBS\\Scripts\\Results_density_modulation\\"

with open(path+"atomic_movement_s_=_0_2.npy", 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,instant_mod]
for ii in range(N_points_zsw):
    dens_mod[ii,0] = np.mean(density_m_avg[N_per_avg*ii:N_per_avg*(ii+1)-1])


with open(path+"atomic_movement_s = 1.npy", 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,instant_mod]
for ii in range(N_points_zsw):
    dens_mod[ii,1] = np.mean(density_m_avg[N_per_avg*ii:N_per_avg*(ii+1)-1])


with open(path+"atomic_movement_s = 3.npy", 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,instant_mod]
for ii in range(N_points_zsw):
    dens_mod[ii,2] = np.mean(density_m_avg[N_per_avg*ii:N_per_avg*(ii+1)-1])


with open(path+"atomic_movement_s = 5.npy", 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,instant_mod]
for ii in range(N_points_zsw):
    dens_mod[ii,3] = np.mean(density_m_avg[N_per_avg*ii:N_per_avg*(ii+1)-1])


with open(path+"atomic_movement_s = 9.npy", 'rb') as f:
    density_m = np.load(f)

#density_m_avg = np.mean(density_m,axis=1)
density_m_avg = density_m[:,instant_mod]
for ii in range(N_points_zsw):
    dens_mod[ii,4] = np.mean(density_m_avg[N_per_avg*ii:N_per_avg*(ii+1)-1])
    
    
plt.figure(10)

for jj in range(5):
    plt.plot(dens_mod[:,jj],label='s = '+str(s_dens_vec[jj]))
plt.legend()

#%% Section 1 - Contrast calculated using C = abs(J2())exp()



r_max = 3 * min(w0, max([sigma_x, sigma_y]))
z_max = 3 * min(w0, sigma_z)
N_points_r = 50
r = np.linspace(0, r_max, N_points_r)
N_points_z = 100
z = np.linspace(-z_max - h, z_max - h, N_points_z)
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

R, Z, Zsw = np.meshgrid(r, z, zsw, indexing='ij')
rho = np.exp(-(R**2 + (Z + h)**2) / (2 * sigma_x**2))
R_sum = np.sum(R, axis=2)
rho_sum = np.sum(rho, axis=2)

Intensities = np.zeros((N_points_Omega, N_point_theta))
Contraste_1 = np.zeros(N_points_Omega)

for j in range(N_points_Omega):
    Omega = 2 * Omega_0[j] * np.cos(k * (Z + Zsw) * np.cos(theta_0))
    Omega_M = np.sqrt(Omega**2 - G_Sr**2 / 16 + 0.j)
    s = 2 * (np.abs(Omega))**2 / G_Sr**2

    ss = s / (4 * (1 + s)) * (
         2 / (1 + s)
         + np.exp(-G_Sr * tau_c / 2)
         + (s - 1) / (s + 1) * np.cos(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
         + G_Sr / (4 * Omega_M) * (5 * s - 1) / (s + 1) * np.sin(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
    )
    
    # here we make the hypothesis of a broadened spectrum with same width as the Rabi freq
    # ss = s / (4 * (1 + s)) * (
    #     2 / (1 + s) + np.exp(-G_Sr * tau_c / 2) - np.exp(-G_Sr * np.sqrt(s + 1) * tau_c / 2) )


    I = np.zeros(N_point_theta)
    for i in range(N_point_theta):
        phase = np.exp(-1j * 2 * k * (Z + Zsw) * np.cos(theta[i]))
        ss_avg = np.sum(s / (1 + s) + 2 * np.real(phase * ss), axis=2)
        I[i] = np.sum(R_sum * rho_sum * ss_avg)

    Intensities[j, :] = I
    plt.figure(10)
    plt.plot(theta - theta_0, I / np.mean(I))
    Contraste_1[j] = (np.max(I) - np.min(I)) / (np.max(I) + np.min(I)) * 2

with open("Fringes_plane_wave.txt", "w") as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.1f}  ")
        fid.write("\n")

plt.figure(2)
plt.plot(Power, Contraste_1, 'r')
plt.xlabel("Incident laser power (mW)")
plt.ylabel("C")

plt.figure(3)
plt.plot(s_0 / 4, Contraste_1, 'r')
plt.xlabel("s")
plt.ylabel("C")

#%% Section 2 - Reduction of the contrast due to the amplitude variation of the probe beam.
# Calculation with an incident laser beam centered on the atomic cloud and
# as if the amplitude was not tilted
# Calculations done using the equation of I as a function of s, Omega_0 and Omega_M

N_points_r = 51
N_points_z = 101
r_max = 3 * min(w0, max([sigma_x, sigma_y]))
# calculate zmax such as all z are multiple of sw_lambda
z_max = ((3 * min(w0, sigma_z))//((N_points_z-1)/2*sw_lambda))*(N_points_z-1)/2*sw_lambda 

r = np.linspace(0, r_max, N_points_r)
z = np.linspace(-z_max - h, z_max - h, N_points_z)
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

R, Z, Zsw = np.meshgrid(r, z, zsw, indexing='ij')
rho = np.exp(-(R**2 + (Z + h)**2) / (2 * sigma_x**2))
R_sum = np.sum(R, axis=2)
rho_sum = np.sum(rho, axis=2)

Intensities = np.zeros((N_points_Omega, N_point_theta))
Contraste_2 = np.zeros(N_points_Omega)

Intensities0 = np.zeros((N_points_Omega, N_point_theta))
Contraste0_2 = np.zeros(N_points_Omega)

for j in range(N_points_Omega):
    Omega = 2 * Omega_0[j] * np.cos(k * (Z + Zsw) * np.cos(theta_0)) * np.exp(-R**2 / w0**2)
    Omega_M = np.sqrt(Omega**2 - G_Sr**2 / 16 + 0.j)
    s = 2 * (np.abs(Omega))**2 / G_Sr**2

    ss = s / (4 * (1 + s)) * (
        2 / (1 + s)
        + np.exp(-G_Sr * tau_c / 2)
        + (s - 1) / (s + 1) * np.cos(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
        + G_Sr / (4 * Omega_M) * (5 * s - 1) / (s + 1) * np.sin(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
    )
    
    ss0 = s / (2 * (1 + s))
    
    I = np.zeros(N_point_theta)
    I0 = I
    for i in range(N_point_theta):
        phase = np.exp(-1j * 2 * k * (Z + Zsw) * np.cos(theta[i]))
        ss_avg = np.sum((s / (1 + s) + 2 * np.real(phase * ss)), axis=2) 
        ss0_avg = np.sum((s / (1 + s) + 2 * np.real(phase * ss0)), axis=2)
        I[i] = np.sum(R_sum * rho_sum * ss_avg)
        I0[i] = np.sum(R_sum * rho_sum * ss0_avg)



    plt.figure(20)
    plt.plot(theta - theta_0, I / np.mean(I))
    Intensities[j, :] = I
    Intensities0[j, :] = I0
    Contraste0_2[j] = (np.max(I0) - np.min(I0)) / (np.max(I0) + np.min(I0)) * 2

with open("Fringes_centered_beam.txt", "w") as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.1f}  ")
        fid.write("\n")
    
    plt.figure(1)
    plt.plot(Omega_0 * tau_c, Contraste_2, 'g')
    plt.plot(Omega_0 * tau_c, Contraste0_2, 'g')
    plt.legend(["Gaussian beam, centered on the atoms"])

    plt.figure(2)
    plt.plot(Power, Contraste_2, 'g')
    plt.plot(Power, Contraste0_2, 'g')
    plt.errorbar(P_exp, Contrast_exp, xerr=P_err,yerr=Contrast_err,fmt='o')
    plt.xlabel("Incident laser power (mW)")
    plt.ylabel("Contrast")
    #plt.title('no density modulation')
    plt.title('density modulation at 2us, same power beams')
    #plt.title('BF = 0.5, 80% cut on heaviside')

    plt.figure(3)
    plt.plot(s_0 / 4, Contraste_2, 'g')
    plt.plot(s_0 / 4, Contraste0_2, 'g')
    plt.errorbar(s_exp, Contrast_exp, xerr=s_err,yerr=Contrast_err,fmt='o')
    plt.xlabel("s")
    plt.ylabel("C")
    #plt.title('no density modulation')
    plt.title('density modulation at 2us, same power beams')

# here we make the hypothesis of a broadened spectrum with same width as the Rabi freq
#    ss = s / (4 * (1 + s)) * (
#        2 / (1 + s) + np.exp(-G_Sr * tau_c / 2) - np.exp(-G_Sr * np.sqrt(s + 1) * tau_c / 2) )


#%% Section 3 - Same as before, but now we admit some density modulation

# option 1: some artificial modulation    

# we must here choose the density modulation for the atoms
#    BF = 0.5  # bunching factor for atoms
#    (1-BF* s/(1+s) ) 
#    (1 - 0.8*np.heaviside(np.cos(k*(Z+Zsw)*np.cos(theta_0))+BF,0.5) )
#    (1 - BF*np.cos(k*(Z+Zsw)*np.cos(theta_0)))
#    (1-s/(1 + s))
#    (1 - (s/(1+s))**2)

# option 2: modulation calculated from numerical integration

N_points_r = 51
N_points_z = 101
r_max = 3 * min(w0, max([sigma_x, sigma_y]))
# calculate zmax such as all z are multiple of sw_lambda
z_max = ((3 * min(w0, sigma_z))//((N_points_z-1)/2*sw_lambda))*(N_points_z-1)/2*sw_lambda 

r = np.linspace(0, r_max, N_points_r)
z = np.linspace(-z_max - h, z_max - h, N_points_z)
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

R, Z, Zsw = np.meshgrid(r, z, zsw, indexing='ij')
rho = np.exp(-(R**2 + (Z + h)**2) / (2 * sigma_x**2))
R_sum = np.sum(R, axis=2)
rho_sum = np.sum(rho, axis=2)

Intensities = np.zeros((N_points_Omega, N_point_theta))
Contraste_2 = np.zeros(N_points_Omega)

for j in range(N_points_Omega):
    Omega = 2 * Omega_0[j] * np.cos(k * (Z + Zsw) * np.cos(theta_0)) * np.exp(-R**2 / w0**2)
    Omega_M = np.sqrt(Omega**2 - G_Sr**2 / 16 + 0.j)
    s = 2 * (np.abs(Omega))**2 / G_Sr**2

    ss = s / (4 * (1 + s)) * (
        2 / (1 + s)
        + np.exp(-G_Sr * tau_c / 2)
        + (s - 1) / (s + 1) * np.cos(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
        + G_Sr / (4 * Omega_M) * (5 * s - 1) / (s + 1) * np.sin(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
    )

    plt.figure(10)
    plt.plot(Omega[int((N_points_r-1)/2),int((N_points_z-1)/2),:])
    
    if s_0[j]/4 < s_dens_vec[0]:
        density = dens_mod[:,0]
    elif s_0[j]/4 < s_dens_vec[1]:
        density = dens_mod[:,0] + (dens_mod[:,1]-dens_mod[:,0])/(s_dens_vec[1]-s_dens_vec[0])*(s_0[j]/4-s_dens_vec[0])
    elif s_0[j]/4 < s_dens_vec[2]:
        density = dens_mod[:,1] + (dens_mod[:,2]-dens_mod[:,1])/(s_dens_vec[2]-s_dens_vec[1])*(s_0[j]/4-s_dens_vec[1])
    elif s_0[j]/4 < s_dens_vec[3]:
        density = dens_mod[:,2] + (dens_mod[:,3]-dens_mod[:,2])/(s_dens_vec[3]-s_dens_vec[2])*(s_0[j]/4-s_dens_vec[2])
    elif s_0[j]/4 < s_dens_vec[4]:
        density = dens_mod[:,3] + (dens_mod[:,4]-dens_mod[:,3])/(s_dens_vec[4]-s_dens_vec[3])*(s_0[j]/4-s_dens_vec[3])
    else:
        density = dens_mod[:,4]
        
    not1, not2, density_mesh = np.meshgrid(np.ones(len(r)), np.ones(len(z)), density, indexing='ij')
    
    I = np.zeros(N_point_theta)
    for i in range(N_point_theta):
        phase = np.exp(-1j * 2 * k * (Z + Zsw) * np.cos(theta[i]))
        ss_avg = np.sum((s / (1 + s) + 2 * np.real(phase * ss)) * density_mesh, axis=2)      
        I[i] = np.sum(R_sum * rho_sum * ss_avg)



    plt.figure(20)
    plt.plot(theta - theta_0, I / np.mean(I))
    Intensities[j, :] = I
    Contraste_2[j] = (np.max(I) - np.min(I)) / (np.max(I) + np.min(I)) * 2

with open("Fringes_centered_beam.txt", "w") as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.1f}  ")
        fid.write("\n")



plt.figure(1)
plt.plot(Omega_0 * tau_c, Contraste_2, 'g')
plt.legend(["Gaussian beam, centered on the atoms"])

plt.figure(2)
plt.plot(Power, Contraste_2, 'g')
plt.errorbar(P_exp, Contrast_exp, xerr=P_err,yerr=Contrast_err,fmt='o')
plt.xlabel("Incident laser power (mW)")
plt.ylabel("Contrast")
#plt.title('no density modulation')
plt.title('density modulation at 2us, same power beams')
#plt.title('BF = 0.5, 80% cut on heaviside')

plt.figure(3)
plt.plot(s_0 / 4, Contraste_2, 'g')
plt.errorbar(s_exp, Contrast_exp, xerr=s_err,yerr=Contrast_err,fmt='o')
plt.xlabel("s")
plt.ylabel("C")
#plt.title('no density modulation')
plt.title('density modulation at 2us, same power beams')


#%%
# Section 3 - Reduction of the contrast due to the intensity variation of the probe beam.
# Calculation with an incident laser beam NOT centered on the atomic cloud
# and the real atomic density distribution

x_max = 3 * min(w0, sigma_x)    # Integration limit
y_max = 3 * min(w0, sigma_y)
z_max = 3 * min(w0, sigma_z)

N_points_x = 50
N_points_y = 50
x = np.linspace(-x_max, x_max, N_points_x)
y = np.linspace(-y_max, y_max, N_points_y)

N_points_z = 100
z = np.linspace(-z_max - h, z_max - h, N_points_z)

N_points_zsw = 20
zsw = np.linspace(0, lambda0/np.cos(theta_0) * (N_points_zsw - 1) / N_points_zsw, N_points_zsw)

X, Y, Z, Zsw = np.meshgrid(x, y, z, zsw, indexing='ij')

rho = np.exp(-X**2 / (2 * sigma_x**2) - Y**2 / (2 * sigma_y**2) - (Z + h)**2 / (2 * sigma_z**2))  # Atomic density profile

X_sum = np.sum(X, axis=3)
Y_sum = np.sum(Y, axis=3)
rho_sum = np.sum(rho, axis=3)

Intensities = np.zeros((N_points_Omega, N_point_theta))
Contraste_3 = np.zeros(N_points_Omega)

for j in range(N_points_Omega):
    # Laser beam path calculated with theta_0, beams not centered on the atoms
    Omega = Omega_0[j] * (
        np.exp(-(X**2 + (Y * np.cos(theta_0) + (Z + Zsw) * np.sin(theta_0))**2) / w0**2) *
        np.exp(1j * k * (-Y * np.sin(theta_0) + (Z + Zsw) * np.cos(theta_0)))
        +
        np.exp(-(X**2 + (-Y * np.cos(theta_0) + (Z + Zsw) * np.sin(theta_0))**2) / w0**2) *
        np.exp(1j * k * (-Y * np.sin(theta_0) - (Z + Zsw) * np.cos(theta_0)))
    )

    Omega_M = np.sqrt(np.abs(Omega)**2 - G_Sr**2 / 16 +0.j)
    s = 2 * np.abs(Omega)**2 / G_Sr**2

    ss = s / (4 * (1 + s)) * (
        2 / (1 + s)
        + np.exp(-G_Sr * tau_c / 2)
        + (s - 1) / (s + 1) * np.cos(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
        + G_Sr / (4 * Omega_M) * (5 * s - 1) / (s + 1) * np.sin(Omega_M * tau_c) * np.exp(-3 * G_Sr * tau_c / 4)
    )

    I = np.zeros(N_point_theta)
    for i in range(N_point_theta):
        phase = np.exp(-1j * 2 * k * (Z + Zsw) * np.cos(theta[i]))
        ss_avg = np.sum(s / (1 + s) + 2 * np.real(phase * ss), axis=3)
        I[i] = np.sum(rho_sum * ss_avg)

    plt.figure(30)
    plt.plot(theta - theta_0, I / np.mean(I))
    Intensities[j, :] = I
    Contraste_3[j] = (np.max(I) - np.min(I)) / (np.max(I) + np.min(I)) * 2  # Contrast calculation

# Write to file
with open('Fringes_separated_beams.txt', 'w') as fid:
    for n in range(N_points_Omega):
        for m in range(N_point_theta):
            fid.write(f"{Intensities[n, m]:.1f}  ")
        fid.write("\n")

# Plot contrast vs power
plt.figure(2)
plt.plot(Power, Contraste_3, 'm')
plt.xlabel("Incident laser power (mW)", fontsize=16)
plt.ylabel("C", fontsize=16)
plt.legend(['Plane wave', 'Gaussian beam', 'Gaussian beam, tilted, not centered on the atoms, θ₀'])

# Plot contrast vs s
plt.figure(3)
plt.plot(s_0 / 4, Contraste_3, 'm')
plt.xlabel("s", fontsize=16)
plt.ylabel("C", fontsize=16)
plt.legend(['Plane wave', 'Gaussian beam', 'Gaussian beam, tilted, not centered on the atoms, θ₀'])