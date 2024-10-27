import numpy as np
import matplotlib.pyplot as plt
import os

# Ruta del directorio que quieres verificar
directorio = './output/'
# Obtener la lista de archivos en el directorio
archivos = os.listdir(directorio)

cantidadArchivos = len(archivos)

#print(archivos)

w = np.zeros(cantidadArchivos)

Fy = np.zeros(cantidadArchivos)

Fm = np.zeros(cantidadArchivos)

for i, archivo in enumerate(archivos):

    wArray, fyArray, fmArray = np.genfromtxt(f'./output/{archivo}',delimiter=' ', usecols=(0,1,2),unpack=True)
    w[i] = wArray[0]
    Fm[i] = fmArray[0]
    Fy[i] = np.mean(fyArray[100:])


plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')

fig, axes = plt.subplots(1, 1, figsize=(7, 6))


#Gráfico.

axes.plot(w, Fm, "bo", label=r"$F_m$ vs $\omega$")
axes.plot(w, Fy, "ro", label=r"$F_y$ vs $\omega$")
#Se ajustan demás detalles del gráfico.
axes.set_xlabel(r'$\omega$', fontsize=12)
axes.set_ylabel(r'$F_y$,$F_m$',fontsize=12)
axes.legend()
axes.grid(True, linestyle='--')
axes.set_title(r"Fuerza en y vs $\omega$", fontsize=14)


plt.tight_layout()
plt.savefig(f'FyvsFm.png')


#plt.show()
