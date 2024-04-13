import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import pearsonr




val,dt = np.genfromtxt('data.txt', unpack=True, usecols=(0,1), max_rows=1)

print(val)
print(dt)
t, y = np.genfromtxt('data.txt', unpack=True, usecols=(0,1),skip_header=1)


peaks, _ = find_peaks(y) #Se encuentran los máximos locales.

tpeak = t[peaks] #Se hallan los tiempos para dichos picos.
tpeaksD = np.roll(t[peaks], 1) #Se corren los valores del arreglo en una posición.
periodos = tpeak - tpeaksD #Se hallan el tiempo que hay entre cada pico.
periodos = periodos[1:] #Se descarta el primer valor puesto que no aporta información.
periodo = np.mean(periodos)


fig, axes = plt.subplots(figsize=(6, 6))

axes.plot(t, y, '.', color='black', label=r'$(x,y). dt = {}$'.format(dt))


# Se ajustan demás detalles del gráfico.
axes.set_xlabel('x', fontsize=12)
axes.set_ylabel('y', fontsize=12)

axes.legend(loc='upper left')   #
axes.grid(True, linestyle='--')

axes.set_title("y vs x. Periodo = {}".format(periodo), fontsize=14)
plt.tight_layout()
fig.savefig('pPlot.png')
