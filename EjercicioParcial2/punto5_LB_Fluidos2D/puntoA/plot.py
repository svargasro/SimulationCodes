import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del archivo .dat
data = np.loadtxt('wind.dat')

# Asignar los datos a las variables correspondientes
x = data[:, 0]
y = data[:, 1]
u = data[:, 2]
v = data[:, 3]

# Filtrar los puntos donde el vector es cero (u == 0 y v == 0)
mask = (u != 0) | (v != 0)
x_filtered = x[mask]
y_filtered = y[mask]
u_filtered = u[mask]
v_filtered = v[mask]

# Crear el gráfico del campo vectorial
plt.figure(figsize=(8, 6))
plt.quiver(x_filtered, y_filtered, u_filtered, v_filtered, angles='xy', scale_units='xy', scale=1.5, color='purple', linewidth=1.5)
plt.xlim([x.min()-1, x.max()+1])
plt.ylim([y.min()-1, y.max()+1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Campo de velocidades punto A')
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('plotA.pdf')


