import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.stats import pearsonr




def lineal_model(x, a, b):
    return a*x + b

# Ruta del directorio que quieres verificar
directorio = './output/'
# Obtener la lista de archivos en el directorio
archivos = os.listdir(directorio)
# Contar la cantidad de archivos
cantidadArchivos = len(archivos)

#print(archivos)

Re = np.zeros(cantidadArchivos)

cA = np.zeros(cantidadArchivos)


for i, archivo in enumerate(archivos):

    ReArray, cAArray = np.genfromtxt(f'./output/{archivo}',delimiter=' ', usecols=(0,1),unpack=True)

    Re[i] = ReArray[0]
    cA[i] = np.mean(cAArray[100:]) #Se toma desde el dato 100 porque es el dato aproximado en el que se estabiliza.


plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')

fig, axes = plt.subplots(2, 1, figsize=(7, 6))

# print(Re)
# print(cA)

ReSinLog = Re
cASinLog = cA

Re = np.log(Re)
cA = np.log(cA)

#Ajuste lineal:
parameters, covarian_matrix = curve_fit(lineal_model, Re, cA)
a, b = parameters

cA_Adjusted = lineal_model(Re, a, b)

r2, _ = pearsonr(cA, cA_Adjusted)


#Gráficos.
#
axes[0].plot(Re, cA,"bo", label="Log(cA) vs. Log(Re)")
axes[0].plot(Re, cA_Adjusted, label=r"Ajuste lineal ax + b. $R^2=${}, $a={}$ , $b={}$".format(round(r2,2),round(a,2),round(b,2)))
#Se ajustan demás detalles del gráfico.
axes[0].set_xlabel('Log(Re)', fontsize=12)
axes[0].set_ylabel('Log(cA)',fontsize=12)
axes[0].legend()
axes[0].grid(True, linestyle='--')
axes[0].set_title("Gráfico loglog Coeficiente de arrastre vs. Número de Reynolds", fontsize=14)

axes[1].plot(ReSinLog, cASinLog,"bo", label="cA vs. Re")
#Se ajustan demás detalles del gráfico.
axes[1].set_xlabel('Re', fontsize=12)
axes[1].set_ylabel('cA',fontsize=12)
axes[1].legend()
axes[1].grid(True, linestyle='--')
axes[1].set_title("Coeficiente de arrastre vs. Número de Reynolds", fontsize=14)


plt.tight_layout()
plt.savefig(f'cAvsRe.png')
#plt.show()
