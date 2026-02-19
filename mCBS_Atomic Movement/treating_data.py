# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:52 2026

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random

def local_s(s,z,k_vec):
    
    s_l = 2*s*(1+np.cos(2*k_vec*z))
    
    return s_l

T_max = 10e-6

## s = 0.2
s_in = 0.2
data = np.load('atomic_movement_s = 0,2.npy',encoding = 'ASCII')
N_x = data.shape[0]
Nat = int(np.sum(data,0)[0])
N_steps_T = data.shape[1]-1
x_vector = data[:,0]
density_dist = data[:,1:]/Nat*N_x
T_step = T_max/(N_steps_T-1)

lambda0 = 461e-9
k = 2*np.pi/lambda0
x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)

pop_ex = local_s(s_in,k,x_vector)/(local_s(s_in,k,x_vector)+1)/2

plt.figure()
for jj in [0,2,4,8,20]:
    time_index = jj
    plt.plot(x_vector/lambda0,density_dist[:,time_index],label="Δt = "+str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,pop_ex,label = 'exc. population',linestyle ='--')
plt.title('(a) s = '+str(s_in),fontsize = 14)
plt.legend(fontsize = 12)
plt.xlabel('z/$\lambda$',fontsize = 14)
plt.ylabel('Normalized density',fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim([-0.05,2.0])
plt.xlim([-0.05,1.1])
plt.savefig('movement_s = '+str(s_in)+'.svg',format = 'svg')

##%% s = 1
s_in = 1
data = np.load('atomic_movement_s = 1.npy',encoding = 'ASCII')
N_x = data.shape[0]
Nat = int(np.sum(data,0)[0])
N_steps_T = data.shape[1]-1
x_vector = data[:,0]
density_dist = data[:,1:]/Nat*N_x
T_step = T_max/(N_steps_T-1)

lambda0 = 461e-9
k = 2*np.pi/lambda0
x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)

pop_ex = local_s(s_in,k,x_vector)/(local_s(s_in,k,x_vector)+1)/2

plt.figure()
for jj in [0,2,4,8,20]:
    time_index = jj
    plt.plot(x_vector/lambda0,density_dist[:,time_index],label="Δt = "+str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,pop_ex,label = 'exc. population',linestyle = '--')
plt.title('(b) s = '+str(s_in),fontsize = 14)
plt.legend(fontsize = 12)
plt.xlabel('z/$\lambda$',fontsize = 14)
plt.ylabel('Normalized density',fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim([-0.05,2.0])
plt.xlim([-0.05,1.1])
plt.savefig('movement_s = '+str(s_in)+'.svg',format = 'svg')

##%% s = 5
s_in = 5
data = np.load('atomic_movement_s = 5.npy',encoding = 'ASCII')
N_x = data.shape[0]
Nat = int(np.sum(data,0)[0])
N_steps_T = data.shape[1]-1
x_vector = data[:,0]
density_dist = data[:,1:]/Nat*N_x
T_step = T_max/(N_steps_T-1)

lambda0 = 461e-9
k = 2*np.pi/lambda0
x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)

pop_ex = local_s(s_in,k,x_vector)/(local_s(s_in,k,x_vector)+1)/2

plt.figure()
for jj in [0,2,4,8,20]:
    time_index = jj
    plt.plot(x_vector/lambda0,density_dist[:,time_index],label="Δt = "+str(round(time_index*T_step*1e6,1))+' us')
plt.plot(x_vector/lambda0,pop_ex,label = 'exc. population',linestyle = '--')
plt.title('(c) s = '+str(s_in),fontsize = 14)
plt.legend(fontsize = 12)
plt.xlabel('z/$\lambda$',fontsize = 14)
plt.ylabel('Normalized density',fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim([-0.05,2.0])
plt.xlim([-0.05,1.1])
plt.savefig('movement_s = '+str(s_in)+'.svg',format = 'svg')

