from operator import lt

import numpy as np
import matplotlib.pyplot as plt
import os

def local_s(s,z,k_vec):
    
    s_l = 2*s*(1+np.cos(2*k_vec*z))
    
    return s_l

def plot_columns(path, delimiter=','):

    # List of line styles to cycle through for each y column
    linestyles = ['-', '--', '-.', ':', (0, (1, 1)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    # '-' solid, '--' dashed, '-.' dash-dot, ':' dotted,
    # (0, (1, 1)) densely dotted, (0, (5, 1)) long dash, (0, (3, 1, 1, 1)) dash-dot-dot

    lambda0 = 461e-9
    k = 2*np.pi/lambda0

    time_step = 0.5e-6 # Time step in seconds

    s_in_movement = np.array([0.4, 1.2, 5.0])

    f, ax = plt.subplots(1, 3, figsize=(16, 5))

    for j in range(len(s_in_movement)):
        with open(os.path.join(path,f"atomic_movement_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
            density_m = np.load(f)
        with open(os.path.join(path,f"atomic_movement_abs_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
            density_m_abs = np.load(f)
        with open(os.path.join(path,f"atomic_movement_abs_2_s={s_in_movement[j]:.1f}.npy"), 'rb') as f:
            density_m_abs_2 = np.load(f)

        N_x = density_m.shape[0]

        density_m = density_m / np.mean(density_m, axis=0)  # Normalize each column
        density_m_abs = density_m_abs / np.mean(density_m_abs, axis=0)  # Normalize each column
        density_m_abs_2 = density_m_abs_2 / np.mean(density_m_abs_2, axis=0)  # Normalize each column

        density_dist = (density_m[:,0:21] + density_m_abs[:,0:21] + density_m_abs_2[:,0:21])/3

        x_vector = np.linspace(0,lambda0*(N_x-1)/N_x,N_x)

        pop_ex = local_s(s_in_movement[j],k,x_vector)/(local_s(s_in_movement[j],k,x_vector)+1)/2

        for jj in [0,2,4,8,20]:
            time = jj*time_step
            ax[j].plot(x_vector/lambda0,density_dist[:,jj],label="$\Delta t = $"+str(round(time*1e6,1))+' $\mu$s',linewidth=2)
        
        ax[j].plot(x_vector/lambda0,pop_ex,label = 'exc. population',linestyle = '--',linewidth=2)

        ax[j].set_xlabel("z/$\lambda$",fontsize=20)
        ax[j].tick_params(labelsize=16)
        ax[j].set_title(f"$s_0 = {s_in_movement[j]:.1f}$",fontsize=20)
        ax[j].set_ylim([-0.05,2.05])
        
    ax[0].set_ylabel("Relative density",fontsize=20)
    ax[2].legend(bbox_to_anchor=(1.6, 1),loc='upper right', borderaxespad=0.,fontsize = 16)
    plt.tight_layout()


if __name__ == "__main__":
    path = 'C:\\Users\\Raul\\Documents\\Repositories\\2026_mCBS_saturated\\mCBS_Atomic_Movement\\data'

    plot_columns(path)
    plt.show()
