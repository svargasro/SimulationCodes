import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import pearsonr







t, ybola, ypared = np.genfromtxt('data.txt', unpack=True, usecols=(0,1,2))

initialData = 9*int(len(t)/10)
t= t[initialData:]
ybola = ybola[initialData:]
ypared = ypared[initialData:]

fig, axes = plt.subplots(figsize=(6, 6))

axes.plot(t, ybola, '-', color='black', label=r'$(x,y)_{b}$')
axes.plot(t, ypared, '-', color='blue', label=r'$(x,y)_{p}$')


# Se ajustan demás detalles del gráfico.
axes.set_xlabel('t', fontsize=12)
axes.set_ylabel('y', fontsize=12)

axes.legend(loc='upper left')   #
axes.grid(True, linestyle='--')

axes.set_title("y vs x.", fontsize=14)
plt.tight_layout()
fig.savefig('initialPlot.png')
