import numpy as np
import matplotlib.pyplot as plt
import os

# Ruta del directorio que quieres verificar
directorio = './output2/'
# Obtener la lista de archivos en el directorio
archivos = os.listdir(directorio)


for archivo in archivos:
    tArray, fyArray = np.genfromtxt(f'./output2/{archivo}',delimiter=' ', usecols=(0,1),unpack=True)


plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')

fig, axes = plt.subplots(1, 1, figsize=(7, 6))


#Gráfico.

axes.plot(tArray, fyArray ,"bo", label=r"$F_y$ vs $t$. $U_{fan}=0.1$")
#Se ajustan demás detalles del gráfico.
axes.set_xlabel(r'$t$', fontsize=12)
axes.set_ylabel(r'$F_y$',fontsize=12)
axes.legend()
axes.grid(True, linestyle='--')
axes.set_title("Fuerza en y vs. tiempo", fontsize=14)


plt.tight_layout()
plt.savefig(f'Fyvst.png')
#plt.show()
