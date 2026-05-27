# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:52 2026

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random



T_max = 15e-6

s_in = 0.8
data = np.load('test.npy',encoding = 'ASCII')
print(data.shape)
N_x = data.shape[0]
Nat = int(np.sum(data,0)[0])
N_steps_T = data.shape[1]-1
x_vector = data[:,0]
density_dist = data[:,1:]/Nat*N_x
T_step = T_max/(N_steps_T-1)

lambda0 = 461e-9
k = 2*np.pi/lambda0
x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)


plt.figure()
for jj in [0,2,4,8,9]:
    time_index = jj
    plt.plot(x_vector/lambda0,density_dist[:,time_index],label="Δt = "+str(round(time_index*T_step*1e6,1))+' us')
plt.title('(a) $s_0 = $'+str(s_in),fontsize = 14)
plt.legend(fontsize = 12)
plt.xlabel('z/$\lambda$',fontsize = 14)
plt.ylabel('Normalized density',fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim([-0.05,2.0])
plt.xlim([-0.05,1.1])
#plt.savefig('movement_s = '+str(s_in)+'.svg',format = 'svg')
plt.show()