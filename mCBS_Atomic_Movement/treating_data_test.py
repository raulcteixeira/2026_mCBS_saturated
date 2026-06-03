# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:52 2026

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random



T_max = 15e-6

s_in = 1.2
data = np.load(f'mCBS_Atomic_Movement/data/atomic_movement_s={s_in:.1f}.npy', encoding='ASCII')
#data = np.load(f'mCBS_Atomic_Movement/old/atomic_movement_s = 0,2.npy', encoding='ASCII')
print(data.shape)
N_x = data.shape[0]
Nat = int(np.sum(data,0)[0])
N_steps_T = data.shape[1]-2
#x_vector = data[:,0]
#density_dist = data[:,1:]/Nat*N_x
print(N_steps_T)
density_dist = data/Nat*N_x

T_step = T_max/N_steps_T

lambda0 = 461e-9
k = 2*np.pi/lambda0
x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)


plt.figure()
for jj in [0,2,4,8,20]:
    time_index = jj
    plt.plot(x_vector/lambda0,density_dist[:,time_index],label="$\Delta t = $"+str(round(time_index*T_step*1e6,1))+' $\mu$s')
plt.title('(b) $s_0 = $'+str(s_in),fontsize = 18)
plt.legend(fontsize = 16, loc = 'upper right')
plt.xlabel('$z/\lambda$',fontsize = 18)
#plt.ylabel('Normalized density',fontsize = 18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylim([0.6,1.9])
plt.locator_params(axis='y', nbins=6)
plt.xlim([-0.05,1.1])
plt.tight_layout()
#plt.savefig('movement_s = '+str(s_in)+'.svg',format = 'svg')
plt.show()