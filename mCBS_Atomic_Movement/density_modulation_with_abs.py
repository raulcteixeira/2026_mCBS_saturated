# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 15:32:19 2025

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import os


def local_s(s,z,k_vec,alph):
    
    s_l = s*(1 + alph + 2*np.sqrt(alph)*np.cos(2*k_vec*z))
    
    return s_l

lambda0 = 461e-9
k = 2*np.pi/lambda0

m = 87.9 * 1.66e-27 
hbar = 6.63e-34/2/np.pi
kB = 1.38e-23
T0 = 2e-6
vT0 = np.sqrt(kB*T0/m)
vr = hbar*k/m
print(vT0,vr)
G_Sr = 2 * np.pi * 30.5e6  # Linewidth in Hz

path = r"C:\Users\Raul\Documents\Repositories\2026_mCBS_saturated\mCBS_Atomic_Movement\data"
data_b0 = np.loadtxt(os.path.join(path,"Attenuation_Raul.txt"))
s_in = data_b0[:,0]
b0 = (data_b0[:,1] + data_b0[:,2])/2
alpha = np.exp(-b0)

plt.figure()
plt.plot(s_in,b0,label='b0')
plt.xlabel('saturation parameter at the center of the cloud')
plt.ylabel('optical thickness b0')
plt.legend()
plt.show()



Nat = 100000

N_x = 100
N_steps_T = 21
T_max = 10e-6
dt = 1/G_Sr
T_step = T_max/(N_steps_T-1)
delta_t_max = T_max/(N_steps_T-1)/10

for kk in range(len(s_in)):
    is_excited = False
    final_x = np.zeros((N_x,N_steps_T+1))
    for ii in range(Nat):
        if ii%100 ==0: 
            print(f"s_in = {s_in[kk]}: Processing atom {ii}")
        x = random.random()*lambda0
        v = random.gauss(mu=0.,sigma=vT0)
        #print(x/lambda0)
        time = 0
        x_index = int((x/lambda0*N_x)%N_x)
        final_x[x_index,0] = final_x[x_index,0] + 1
        for jj in range(N_steps_T):
            T_limit= (jj+1)*T_step
            while time < T_limit:
                time = time + dt
                x = x + v*dt
                sort = random.random()
                if not is_excited:
                    s_now = local_s(s_in[kk],k,x,alpha[kk])
                    photon_rate = 1/2 * G_Sr * s_now/(s_now+1)
                    if sort < photon_rate*dt:
                        is_excited = True
                        sort2 = random.random()
                        if sort2 < 1/(1+alpha[kk]):
                            v = v + vr
                        else:
                            v = v - vr
                else:
                    if sort < G_Sr*dt:
                        is_excited = False
                        if sort < 1/6*G_Sr*dt:
                            v = v + vr
                        elif sort < 1/3*G_Sr*dt:
                            v = v - vr
            x_index = int((x/lambda0*N_x)%N_x)
            final_x[x_index,jj+1] = final_x[x_index,jj+1] + 1
            
    with open(f'atomic_movement_abs_2_s={s_in[kk]:.1f}.npy', 'wb') as f:
        np.save(f, final_x)



x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)
pop_ex = local_s(s_in[kk],k,x_vector,alpha[kk])/(local_s(s_in[kk],k,x_vector,alpha[kk])+1)


time_interval = 1

for jj in range(int((N_steps_T-1)/time_interval)+1):
    time_index = jj*time_interval
    plt.plot(x_vector/lambda0,final_x[:,time_index],label=str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,Nat/N_x*pop_ex,label = 'population')
plt.title(f's = {s_in[kk]:.1f}')
plt.legend()
plt.xlabel('Position (lambda)')
plt.ylabel('Frequency')
plt.show()



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

#%%

#time_interval = 2

#for jj in range(int((N_steps_T-1)/time_interval)+1):
#    time_index = jj*time_interval
#    plt.plot(x_vector/lambda0,final_x[:,time_index],label=str(round(time_index*T_step*1e6,1))+' us')
#plt.plot(x_vector/lambda0,Nat/N_x*pop_ex,label = 'population')
#plt.title('s = 5')
#plt.legend()
#plt.xlabel('Position (lambda)')
#plt.ylabel('Frequency')

#with open('test.npy', 'wb') as f:

#    np.save(f, final_x)