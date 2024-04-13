import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import pearsonr



def lineal_model(x, a, b):
    return a*x + b

frecuencies = np.array((0.1,0.2,0.4,0.8,1.6,2.0,3.2))

numberFrec = len(frecuencies)
x = np.zeros(numberFrec)
y = np.zeros(numberFrec)

for i in range(numberFrec):
    x[i],y[i] = np.genfromtxt(f'ydesvF{i}.txt', unpack=True, usecols=(0,1))


fig, axes = plt.subplots(figsize=(6, 6))



ydesv = np.log(-y)
frecuencies = np.log(frecuencies)

parameters, covarian_matrix = curve_fit(lineal_model, frecuencies, ydesv)
a, b = parameters

yajustado = lineal_model(frecuencies, a, b)

r2, _ = pearsonr(ydesv, yajustado)

axes.plot(frecuencies, ydesv, '.', color='black', label=r'$(log(f),log(y_{desv}))$')

axes.plot(frecuencies, yajustado, '-', color='blue', label=f'Ajuste lineal R^2 ={round(r2,4)}')

# Se ajustan demás detalles del gráfico.
axes.set_xlabel(r'$log(f)$', fontsize=12)
axes.set_ylabel(r'$log(y_{desv})$', fontsize=12)
axes.legend(loc='upper left')   #
axes.grid(True, linestyle='--')

axes.set_title(r"$log(y_{{desv}})$ vs $log(frecuencia)$. Ley potencias: $y_{{desv}} = {}* f^{{{}}} $".format(round(np.exp(b),3),round(a,3)), fontsize=10)
plt.tight_layout()
fig.savefig('pPlot.png')
