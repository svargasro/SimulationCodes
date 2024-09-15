import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Leer los datos del archivo .dat
data = np.loadtxt('FxFy.dat')

# Asignar los datos a las variables correspondientes
b = data[:, 0]
Fx = data[:, 1]

# Realizar ajuste lineal usando linregress
slope, intercept, r_value, p_value, std_err = linregress(b, Fx)

# Calcular R^2
r_squared = r_value**2

# Crear una figura
plt.figure(figsize=(10, 5))

# Graficar los datos originales
plt.scatter(b, Fx, color='b', label='Fx(b)')

# Graficar la línea de ajuste
plt.plot(b, slope * b + intercept, color='r', label=f'Ajuste lineal {slope:.2f}b+{intercept:.2f} (R² = {r_squared:.4f})')

# Etiquetas y título
plt.xlabel('b')
plt.ylabel('Fx(b)')
plt.title('Fx vs b con ajuste lineal')
plt.grid(True)
plt.legend()

# Ajustar el espacio para asegurar que todo se ajuste bien
plt.tight_layout()

# Guardar la figura en un archivo
plt.savefig('puntoE.png')

# Imprimir R²
