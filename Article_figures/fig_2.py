import numpy as np
import matplotlib.pyplot as plt
import os

def plot_columns(path, delimiter=','):

    # List of line styles to cycle through for each y column
    linestyles = ['-', '--', '-.', ':', (0, (1, 1)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    # '-' solid, '--' dashed, '-.' dash-dot, ':' dotted,
    # (0, (1, 1)) densely dotted, (0, (5, 1)) long dash, (0, (3, 1, 1, 1)) dash-dot-dot


    filename = os.path.join(path,"Contrast.txt")
    # Read header (column names)
    with open(filename, 'r') as f:
        header = f.readline().strip().split(delimiter)

    # Read numeric data (skip header row)
    data = np.loadtxt(filename, skiprows=1,delimiter=delimiter)

    # Handle case of only 2 columns (1D array issue)
    if data.ndim == 1:
        data = data.reshape(-1, len(header))

    size_data = data.shape[1]
    x = data[:, 0]
    x_label = header[0]

    header[1] = "Plane wave"
    header[2] = "Finite size of beam"
    header[3] = "Finite size of beam, absorption"

    plt.figure(1,figsize=(8, 5))
    for i in range(1, data.shape[1]):

        plt.plot(x, data[:, i], label=header[i], linestyle=linestyles[i-1], linewidth=2)

    filename = os.path.join(path,"Contrast_movement.txt")
    # Read header (column names)
    with open(filename, 'r') as f:
        header = f.readline().strip().split(delimiter)

    # Read numeric data (skip header row)
    data = np.loadtxt(filename, skiprows=1,delimiter=delimiter)

    # Handle case of only 2 columns (1D array issue)
    if data.ndim == 1:
        data = data.reshape(-1, len(header))

    x = data[:, 0]
    x_label = header[0]

    header[1] = "Finite size of beam, absorption, motion"
    header[2] = "Only elastic scattering"

    for i in range(1, data.shape[1]):
        style = linestyles[(i - 1 + size_data) % len(linestyles)] 
        plt.plot(x, data[:, i], label=header[i], linestyle=style, linewidth=2)

    plt.xlabel("$s_0$",fontsize=20)
    plt.ylabel("Contrast C",fontsize=20)
    plt.tick_params(labelsize=16)
    plt.xlim([0, 10])
    plt.ylim([0, 1])
    #plt.title("Data plot")
    plt.tight_layout()
    #plt.grid(True)

if __name__ == "__main__":
    path = 'C:\\Users\\Raul\\Documents\\Repositories\\2026_mCBS_saturated'

    # Laser beam waist
    w0 = 2.1e-3
    I_sat = 40.5 # mW/cm^2

    # Exp data
    Delta_t_exp = np.array([30,15,8,10,10,15,15,10,12,15,15])  # Delta t used for each power, in us
    P_exp = np.array([1,2,5,10,15,25,3,15,20,25,3])   # mW
    P_err = 0.07*P_exp
    s_exp = 2*P_exp/(np.pi*(w0*100)**2)/I_sat
    s_err = 0.07*s_exp
    Contrast_exp = np.array([0.74002,0.86048,0.61458,0.33236,0.1703,0.17106,0.58747,0.30815,0.25758,0.14104,0.69217])
    Contrast_err = np.array([0.03602,0.04185,0.03549,0.02789,0.03394,0.03286,0.04734,0.02736,0.03638,0.02033,0.05181])

    plot_columns(path)
    plt.errorbar(s_exp, Contrast_exp, xerr=s_err, yerr=Contrast_err, fmt='o', label='Experimental data')
    plt.legend(fontsize=14)
    plt.show()
