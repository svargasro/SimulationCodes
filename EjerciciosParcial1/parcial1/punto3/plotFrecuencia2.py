import numpy as np
import matplotlib.pyplot as plt


x, y = np.genfromtxt('dataXYF2.txt', unpack=True, usecols=(0,1))

frec, xdesv, ydesv = np.genfromtxt('TablaydesvVSF_SinTorque.txt', unpack=True, usecols=(0,1,2),skip_header=1)

ydesv2 = ydesv[5] #Siempre es el quinto valor.

fig, axes = plt.subplots(figsize=(6, 6))
axes.plot(x, y, '.', color='black', label=r'$(x,y). y_{{desv}} = {}$'.format(ydesv2))

# Se ajustan demás detalles del gráfico.
axes.set_xlabel('x', fontsize=12)
axes.set_ylabel('y', fontsize=12)

axes.legend(loc='upper left')   #
axes.grid(True, linestyle='--')

axes.set_title(r"y vs x. Fuerza de Magnus. $f=2 H_z$", fontsize=12)
plt.tight_layout()
fig.savefig('yVSxYFrecuencia2.png')
