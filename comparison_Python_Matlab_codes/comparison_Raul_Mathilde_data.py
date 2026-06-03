# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:52 2026

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import os

path = r"C:\Users\Raul\Documents\Repositories\2026_mCBS_saturated"
data_Raul = np.loadtxt(os.path.join(path,"Contrast.txt"),skiprows = 1,delimiter=',')
data_Mathilde = np.loadtxt(os.path.join(path,"Donnees_Mathilde.txt"))
data_Mathilde_2 = np.loadtxt(os.path.join(path,"Contrast_Mathilde.txt"))
data_Mathilde_2_Dirac = np.loadtxt(os.path.join(path,"Contrast_Mathilde_Dirac.txt"))
attenuation_Raul = np.loadtxt(os.path.join(path,"Attenuation_Raul.txt"))
attenuation_Mathilde = np.loadtxt(os.path.join(path,"Attenuation_Mathilde.txt"))

s_in_Raul = data_Raul[:,0]
C_Raul = data_Raul[:,1]
C_Raul_el = data_Raul[:,2]

s_in_Mathilde = data_Mathilde[:,0]
C_Mathilde = data_Mathilde[:,1]
C_Mathilde_el = data_Mathilde[:,2]

s_in_Mathilde_2 = data_Mathilde_2[:,0]
C_Mathilde_2 = data_Mathilde_2[:,1]
C_Mathilde_2_el = data_Mathilde_2_Dirac[:,1]

plt.figure()
plt.plot(s_in_Raul, C_Raul, label="Raul")
plt.plot(s_in_Mathilde, C_Mathilde, label="Mathilde")
plt.plot(s_in_Raul, C_Raul_el, label="Raul elastic", linestyle="dashed")
plt.plot(s_in_Mathilde, C_Mathilde_el, label="Mathilde elastic", linestyle="dashed")
plt.plot(s_in_Mathilde_2, C_Mathilde_2, label="Mathilde 2")
plt.plot(s_in_Mathilde_2, C_Mathilde_2_el, label="Mathilde 2 elastic", linestyle="dashed")
plt.xlabel("saturation parameter s")
plt.ylabel("Contrast")
plt.legend()
plt.title("Comparison of contrast vs saturation parameter")
plt.show()

plt.figure()
plt.plot(attenuation_Raul[:,0], attenuation_Raul[:,1], label="Raul")
plt.plot(attenuation_Raul[:,0], attenuation_Raul[:,2], label="Raul")
plt.plot(attenuation_Mathilde[:,0], attenuation_Mathilde[:,1], label="Mathilde")
plt.plot(attenuation_Mathilde[:,0], attenuation_Mathilde[:,2], label="Mathilde")
plt.xlabel("saturation parameter s")
plt.ylabel("Attenuation (OD)")
plt.legend()
plt.title("Comparison of attenuation vs saturation parameter")
plt.show()
