import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

# Parámetros de la elipse
icenter = 64
Ly = 128
W = 32

# Leer los datos del archivo Waves2D.dat
data = np.loadtxt("Waves2D.dat")

# Extraer las columnas de datos
x_data = data[:, 0]
y_data = data[:, 1]
magnitude = data[:, 2]

# Crear una figura y un eje
fig, ax = plt.subplots()

# Graficar los datos con un mapa de colores
sc = ax.scatter(x_data, y_data, c=magnitude, s=10, cmap='viridis', edgecolor='none')
fig.colorbar(sc, ax=ax, label='Magnitud')

# Crear la elipse
theta = np.linspace(0, 2 * np.pi, 100)
x_ellipse = icenter + W * np.cos(theta)
y_ellipse = (Ly / 2) + (Ly / 2) * np.sin(theta)

# Graficar la elipse
ax.plot(x_ellipse, y_ellipse, color='red', label='Elipse')

# Configurar etiquetas y título
ax.set_xlabel('X', fontsize=14)
ax.set_ylabel('Y', fontsize=14)
ax.legend()
ax.set_aspect('equal')  # Mantener la proporción de la elipse

ax.axvline(x=161, color='blue', linestyle='--', label='Rayo Vertical (x=160)')

# Guardar la gráfica
plt.savefig("lenteBiesferica2.jpg")

# Mostrar la gráfica
#plt.show()
