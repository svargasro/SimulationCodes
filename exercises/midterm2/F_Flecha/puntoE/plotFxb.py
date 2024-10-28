import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Leer los datos del archivo .dat
data = np.loadtxt('Fxb.dat')

# Asignar los datos a las variables correspondientes
b = data[:, 0]
Fx = data[:, 1]

# Crear una figura con dos subplots: uno para Fx vs b y otro para Log(Fx) vs Log(b)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5))

# -----------------------
# Gráfico 1: Fx vs b
# -----------------------

# Realizar ajuste lineal en Fx vs b

# Graficar los datos originales
ax1.scatter(b, Fx, color='b', label='Fx(b)')

# Etiquetas y título para Fx vs b
ax1.set_xlabel('b')
ax1.set_ylabel('Fx(b)')
ax1.set_title('Fx vs b')
ax1.grid(True)
ax1.legend()

# -----------------------
# Gráfico 2: Log(Fx) vs Log(b)
# -----------------------

# Aplicar logaritmos a los datos
log_b = b
log_Fx = np.log(Fx)

# Realizar ajuste lineal en Log(Fx) vs Log(b)
slope2, intercept2, r_value2, p_value2, std_err2 = linregress(log_b, log_Fx)
r_squared2 = r_value2**2

# Graficar los datos originales en el espacio logarítmico
ax2.scatter(log_b, log_Fx, color='g', label='Log(Fx) vs b')
ax2.plot(log_b, slope2 * log_b + intercept2, color='r', label=f'Log(Fx) = {slope2:.2f}b + ({intercept2:.2f}) -> Fx = {np.exp(intercept2):.2f}e^({slope2:.2f}b)  (R² = {r_squared2:.4f})')

# Etiquetas y título para Log(Fx) vs Log(b)
ax2.set_xlabel('b')
ax2.set_ylabel('Log(Fx)')
ax2.set_title('Log(Fx) vs b con ajuste lineal')
ax2.grid(True)
ax2.legend()

# Ajustar el espacio para asegurar que todo se ajuste bien
plt.tight_layout()

# Guardar la figura en un archivo
plt.savefig('puntoE.png')
