import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import pearsonr




val,dt = np.genfromtxt('data.txt', unpack=True, usecols=(0,1), max_rows=1) #Extrae valor de dt.

x, z = np.genfromtxt('data.txt', unpack=True, usecols=(0,1),skip_header=1)


fig, axes = plt.subplots(figsize=(6, 6))
axes.plot(x, z, '.', color='black', label=r'$(x,z). dt_{adecuado} = {}$'.format(dt))

# Se ajustan demás detalles del gráfico.
axes.set_xlabel('x', fontsize=12)
axes.set_ylabel('z', fontsize=12)

axes.legend(loc='upper left')   #
axes.grid(True, linestyle='--')

axes.set_title("z vs x. Alcance balón = {}".format(x[-1]), fontsize=14)
plt.tight_layout()
fig.savefig('pPlot.png')
