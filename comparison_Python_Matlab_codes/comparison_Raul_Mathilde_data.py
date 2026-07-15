# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 16:34:52 2026

@author: Raul
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import os

path = r"C:\Users\Raul\Documents\Repositories\2026_mCBS_saturated\comparison_Python_Matlab_codes"
data_Raul = np.loadtxt(os.path.join(path,"Contrast.txt"),skiprows = 1,delimiter=',')
data_Raul_100 = np.loadtxt(os.path.join(path,"Contrast_Raul_100.txt"),skiprows = 1,delimiter=',')
data_Mathilde = np.loadtxt(os.path.join(path,"Donnees_Mathilde.txt"))
data_Mathilde_2 = np.loadtxt(os.path.join(path,"Contrast_Mathilde.txt"))
data_Mathilde_2_Dirac = np.loadtxt(os.path.join(path,"Contrast_Mathilde_Dirac.txt"))
data_Mathilde_3 = np.loadtxt(os.path.join(path,"Contraste_Gaussian_Attenuation.txt"))
attenuation_Raul = np.loadtxt(os.path.join(path,"Attenuation_Raul.txt"))
attenuation_Mathilde = np.loadtxt(os.path.join(path,"Attenuation_Mathilde.txt"))
Omega_Raul = np.loadtxt(os.path.join(path,"Omega_Raul.txt"))
Omega_Mathilde = np.loadtxt(os.path.join(path,"Omega_Mathilde.txt"))
ss_Raul = np.loadtxt(os.path.join(path,"ss_Raul.txt"))
ss_Raul1 = np.loadtxt(os.path.join(path,"ss_Raul1.txt"))
ss_Raul2 = np.loadtxt(os.path.join(path,"ss_Raul2.txt"))
ss_Mathilde = np.loadtxt(os.path.join(path,"ss_Mathilde.txt"))
ss_Mathilde_elastic = np.loadtxt(os.path.join(path,"ss_Mathilde_elastic.txt"))

s_in_Raul = data_Raul[:,0]
C_Raul = data_Raul[:,1]
C_Raul_el = data_Raul[:,2]

s_in_Mathilde = data_Mathilde[:,0]
C_Mathilde = data_Mathilde[:,1]
C_Mathilde_el = data_Mathilde[:,2]

s_in_Mathilde_2 = data_Mathilde_2[:,0]
C_Mathilde_2 = data_Mathilde_2[:,1]
C_Mathilde_2_el = data_Mathilde_2_Dirac[:,1]

s_in_Mathilde_3 = data_Mathilde_3[:,0]
C_Mathilde_3 = data_Mathilde_3[:,1]

plt.figure()
plt.plot(ss_Raul[:,0], ss_Raul[:,1], label="Raul")
plt.plot(ss_Raul[:,0], ss_Raul[:,2], label="Raul elastic")
plt.plot(ss_Mathilde[:,0], ss_Mathilde[:,1], label="Mathilde", linestyle="dashed")
plt.plot(ss_Mathilde_elastic[:,0], ss_Mathilde_elastic[:,1], label="Mathilde elastic", linestyle="dashed")
plt.xlabel("theta (degrees)")
plt.ylabel("ss")
plt.legend()
plt.title("Comparison of ss vs theta")
plt.show()


# plt.figure()
# plt.plot(ss_Raul[0,1:], label="Raul")
# plt.plot(ss_Mathilde_elastic[0,1:], label="Mathilde", linestyle="dashed")
# plt.plot(ss_Raul[13,1:], label="Raul 1")
# plt.plot(ss_Raul1[13,1:], label="Raul 1a", linestyle="dashed")
# plt.plot(ss_Raul2[13,1:], label="Raul 1b", linestyle="dashed")
# plt.plot(ss_Mathilde_elastic[13,1:], label="Mathilde 1", linestyle="dashed")
# plt.xlabel("z (degrees)")
# plt.ylabel("ss")
# plt.legend()
# plt.title("Comparison of ss vs theta")
# plt.show()

plt.figure()
plt.plot(Omega_Raul[0,1:], label="Raul")
plt.plot(Omega_Mathilde[0,1:], label="Mathilde", linestyle="dashed")
plt.xlabel("Position")
plt.ylabel("Omega")
plt.legend()
plt.title("Comparison of Omega vs Position")
plt.show()

plt.figure()
plt.plot(s_in_Raul, C_Raul, label="Raul")
plt.plot(s_in_Mathilde, C_Mathilde, label="Mathilde",linestyle="dashed")
plt.plot(s_in_Raul, C_Raul_el, label="Raul elastic")
plt.plot(s_in_Mathilde, C_Mathilde_el, label="Mathilde elastic", linestyle="dashed")
plt.plot(s_in_Mathilde_2, C_Mathilde_2, label="Mathilde 2", linestyle="dashed")
plt.plot(s_in_Mathilde_2, C_Mathilde_2_el, label="Mathilde 2 elastic", linestyle="dashed")
plt.plot(s_in_Mathilde_3, C_Mathilde_3, label="Mathilde 3", linestyle="dashdot")
plt.plot(s_in_Raul, data_Raul_100[:,1], label="Raul with 100 points")
plt.plot(s_in_Raul, data_Raul_100[:,2], label="Raul elastic with 100 points")
plt.xlabel("saturation parameter s")
plt.ylabel("Contrast")
plt.legend()
plt.title("Comparison of contrast vs saturation parameter")
plt.show()

plt.figure()
plt.plot(attenuation_Raul[:,0], attenuation_Raul[:,1], label="Raul in")
plt.plot(attenuation_Raul[:,0], attenuation_Raul[:,2], label="Raul refl")
plt.plot(attenuation_Mathilde[:,0], attenuation_Mathilde[:,1], label="Mathilde in",linestyle="dashed")
plt.plot(attenuation_Mathilde[:,0], attenuation_Mathilde[:,2], label="Mathilde refl",linestyle="dashed")
plt.xlabel("saturation parameter s")
plt.ylabel("Attenuation (OD)")
plt.legend()
plt.title("Comparison of attenuation vs saturation parameter")
plt.show()
