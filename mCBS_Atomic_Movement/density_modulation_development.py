# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:32:19 2025

@author: Raul
"""

def local_s(s,z,k_vec):
    
    s_l = 2*s*(1+np.cos(2*k_vec*z))
    
    return s_l

def photon_drag(s,z,k_vec,sqr_alph):
    
    ph_drag = s*((1-sqr_alph)**2 + 2*sqr_alph*(1-sqr_alph)*(np.cos(k_vec*z))**2)
    
    return ph_drag


import numpy as np
import matplotlib.pyplot as plt
import random

lambda0 = 461e-9;
k = 2*np.pi/lambda0

s_in = 1.0        # saturation parameter of incoming beam
alpha = 1.0         # power of reflected beam relative to incoming
sqr_alpha = np.sqrt(alpha)

m = 87.9 * 1.66e-27 
hbar = 6.63e-34/2/np.pi
kB = 1.38e-23
T0 = 2e-6
vT0 = np.sqrt(kB*T0/m)
vr = hbar*k/m
print(vT0,vr)
G_Sr = 2 * np.pi * 30.5e6  # Linewidth in Hz

Nat = 1000000

N_x = 100
N_steps_T = 21
T_max = 10e-6
dt = 1/G_Sr
T_step = T_max/(N_steps_T-1)
delta_t_max = T_max/(N_steps_T-1)/10

final_x = np.zeros((N_x,N_steps_T+1))

for ii in range(Nat):
    print(ii)
    x = random.random()*lambda0
    v = random.gauss(mu=0.,sigma=vT0)
    #print(x/lambda0)
    time = 0
    x_index = int((x/lambda0*N_x)%N_x)
    final_x[x_index,0] = final_x[x_index,0] + 1
    for jj in range(N_steps_T):
        T_limit= (jj+1)*T_step
        while time < T_limit:
            s_now = local_s(s_in,k,x)
            # sp_em_rate: spontaneous emission rate
            sp_em_rate = 1/2 * G_Sr * s_now/(1 + s_now) 
            # photon_abs_rate: average force from beam imbalance
            photon_abs_rate = 1/2 * G_Sr * photon_drag(s_in,k,x,sqr_alpha)/(1 + s_now)
            time = time + dt
            x = x + v*dt
            sort = random.random()
            if sort < 1/6*sp_em_rate*dt:
                v = v + vr
            elif sort < 1/3*sp_em_rate*dt:
                v = v - vr
            elif sort > (1 - photon_abs_rate*dt):   # average force from beam imbalance
                v = v + vr
        x_index = int((x/lambda0*N_x)%N_x)
        final_x[x_index,jj+1] = final_x[x_index,jj+1] + 1


x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)
pop_ex = local_s(s_in,k,x_vector)/(local_s(s_in,k,x_vector)+1)


time_interval = 1

plt.figure(1)
for jj in range(int((N_steps_T-1)/time_interval)+1):
    time_index = jj*time_interval
    plt.plot(x_vector/lambda0,final_x[:,time_index],label=str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,Nat/N_x*pop_ex,label = 'population')
plt.title('s = 1.0, alpha = 0.9')
plt.legend()
plt.xlabel('Position (lambda)')
plt.ylabel('Frequency')

#%%

time_interval = 2

plt.figure(2)
for jj in range(int((N_steps_T-1)/time_interval)+1):
    time_index = jj*time_interval
    plt.plot(x_vector/lambda0,final_x[:,time_index],label=str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,Nat/N_x*pop_ex,label = 'population')
plt.title('s = 1.0, alpha = 0.9')
plt.legend()
plt.xlabel('Position (lambda)')
plt.ylabel('Frequency')

with open('test.npy', 'wb') as f:

    np.save(f, final_x)


# plt.hist(final_x[:,0],bins = 20,histtype='step',label = 't = 0')
# plt.hist(final_x[:,1],bins = 20,histtype='step',label = 't = 0.1us')
# plt.hist(final_x[:,2],bins = 20,histtype='step',label = 't = 0.2us')
# plt.hist(final_x[:,3],bins = 20,histtype='step',label = 't = 0.3us')
# plt.hist(final_x[:,5],bins = 20,histtype='step',label = 't = 0.5us')
# plt.hist(final_x[:,7],bins = 20,histtype='step',label = 't = 0.7us')
# plt.hist(final_x[:,10],bins = 20,histtype='step',label = 't = 1us')
# plt.hist(final_x[:,15],bins = 20,histtype='step',label = 't = 1.5us')
# plt.hist(final_x[:,20],bins = 20,histtype='step',label = 't = 2us')
# plt.title('s = 5')
# plt.legend()