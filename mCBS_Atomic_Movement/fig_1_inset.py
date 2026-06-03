import numpy
import matplotlib.pyplot as plt
import os
import scipy.ndimage

path = 'G:\\Meu Drive\\Documentos\\Academics\\Sr 1\\projetos\\2020 - saturated mCBS'

x = []
y = []

with open(os.path.join(path,"fringes_1_mW_exp.txt"), "r") as file:
    for line in file:
        # Skip empty lines
        if line.strip() == "":
            continue

        # Split the line into two columns
        col1, col2 = line.split()

        # Convert to float (or int if needed)
        x.append(float(col1))
        y.append(float(col2))

y = scipy.ndimage.gaussian_filter1d(y, sigma=1)

plt.plot(x, y)
#plt.xlabel("X")
#plt.ylabel("Y")
plt.xlim([2770,3420])
plt.ylim([0.1,0.6])
plt.xticks([])
plt.yticks([])
#plt.title("Fringes for 1 mW")
plt.show()
